"""apogee_daily — the APOGEE daily reduction on Airflow (card O3).

ONE DAG, observatories chained SERIALLY (LCO then APO), restoring the
production ar_main.py structure per AKS review 2026-09-03: the Airflow
metadata DB is SQLite (single-writer), so parallel per-observatory DAG runs
contend; serial is deliberate.

  stale_run_guard -> update(repo -> sdsscore -> sync_logs)
    -> setup(mjd, date_mjd)
    -> [ LCO group ] -> [ APO group ] -> completion_notification

Each observatory group:
  resolve_night -> wait_for_raw -> tunnel_check -> almanac -> runlist
    -> start_notification -> select_mode -> [ LOCAL chain | slurm_submit ]
    -> join -> metrics_append (N1) -> daily_summary (N3)   (+ chain_status)

LOCAL chain (testing): p3d2d -> quartz(runlist,traces,extract,relflux)
             -> dome(...) -> p2d1d_full -> plots -> dashboard
             -> madgics_gate -> madgics_pipeline -> spectra_workup
             -> workup_cleanup (W3 validate, then rm raw batch dirs;
                legacy workup.jl retired 2026-09-04 — single workup)

SLURM mode (DEFAULT): a single sbatchAKS-style submission of
scripts/daily/run_all.sh per observatory (production dailies as run
historically). AKS launches the Airflow instance himself, so the chain is
human-initiated. Use AR_AIRFLOW_MODE=local (or conf {"mode": "local"}) for
on-node testing.

A no-data observatory (e.g. LCO weathered out) skips its group cleanly and
does NOT block the other: each group's resolve_night runs on none_failed,
so upstream skips pass through while a real LCO failure still blocks APO
(ar_main.py's serial semantics).

Trigger a manual/backfill run with e.g.:
  airflow dags trigger apogee_daily --conf '{"mjd": 61284}'
Re-run a FAILED night with a fresh trigger (not by clearing tasks; see
airflow/INSTALL.md, "Re-running a failed night"):
  airflow dags trigger apogee_daily --logical-date <now, ISO> \\
      --run-id manual__rerun_<N>_<stamp> \\
      --conf '{"mjd": N, "almanac_clobber": true}'
Smoke-test conf (MUST force local mode — the default is slurm — and the
dev channel; summaries are the only posts in summary_only):
  {"mjd": ..., "mode": "local", "outdir": ".../o3_dag_test/",
   "slack_mode": "summary_only", "test_label": true, "run_kind": "test",
   "run_madgics": false, "slack_channel": "C07KQ7BJY5P",
   "raw_nowait": true}
"""

from __future__ import annotations

import os
import subprocess
import sys
import time
from datetime import datetime, timedelta, timezone

# make sibling ar_common importable regardless of how the bundle is parsed
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from airflow import DAG
try:  # Airflow 3 Task SDK location; airflow.exceptions is deprecated there
    from airflow.sdk.exceptions import (AirflowFailException,
                                        AirflowSkipException)
except ImportError:  # pragma: no cover - older Airflow
    from airflow.exceptions import AirflowFailException, AirflowSkipException
from airflow.providers.standard.operators.bash import BashOperator
from airflow.providers.standard.operators.empty import EmptyOperator
from airflow.providers.standard.operators.python import (
    BranchPythonOperator,
    PythonOperator,
    ShortCircuitOperator,
)
from airflow.providers.standard.sensors.python import PythonSensor
from airflow.utils.task_group import TaskGroup
from airflow.utils.trigger_rule import TriggerRule

import ar_common as C


# ---------------------------------------------------------------------------
# setup group — night selection, VERBATIM ar_main.py semantics
# ---------------------------------------------------------------------------

def compute_mjd(**context):
    """ar_main.py setup.mjd, verbatim semantics:
        int(Time(data_interval_start).mjd) - 1
    (stdlib mjd_of replaces astropy). The `mjd` param (>= 0) overrides for
    manual/past-day runs (data already on disk) and flows through the whole
    chain via this task's xcom.

    WHY THIS TIMING WORKS (AKS): the Globus Utah->CCA raw sync is a fixed
    daily task at 4:00-4:30 am, so by the 7 am ET trigger the night's data
    has already landed — the run processes the data that arrived that
    morning. wait_for_raw downstream is cheap insurance only (instant pass
    normally).

    OPERATIONAL SUBTLETY (validated against history, 2026-09-03): on this
    Airflow 3 instance a cron-string schedule is a CronTriggerTimetable —
    data_interval_start == the trigger time itself (zero-width interval).
    So this expression selects int(mjd(7am-ET-of-run-day)) - 1, and the 300
    historical daily almanac files confirm production behaved exactly so
    (203/205 scheduled runs at delta -1 from the run morning; the only
    outliers are the two very first runs on 2025-09-02, instance stand-up).
    Under Airflow 2's CronDataIntervalTimetable the SAME expression would
    have named one day earlier — do not port this line to an Airflow 2
    instance unchanged."""
    p = context["params"]
    override = int(p.get("mjd", -1))
    if override >= 0:
        print(f"mjd from param override: {override}")
        return override
    # A manual trigger without --logical-date has no data interval in
    # Airflow 3 (2026-11-17: KeyError 'data_interval_start'); fall back to
    # logical_date, then run_after, then now (C.run_reference_time).
    ref, src = C.run_reference_time(context)
    mjd = int(C.mjd_of(ref)) - 1
    print(f"mjd from {src}={ref}: {mjd}")
    return mjd


def compute_date_mjd(**context):
    """ar_main.py setup.date_mjd, verbatim semantics:
        (data_interval_start astimezone America/New_York
         - timedelta(days=2)).strftime('%Y-%m-%d')
    (zoneinfo replaces pytz).

    Manual runs (2026-11-17): when conf gives `mjd`, the label is derived
    from that MJD (C.date_label_from_mjd, the same convention), not from
    the trigger time, which named the wrong night for past-day re-runs.
    Without a data interval (a trigger without --logical-date, which made
    this task fail with KeyError 'data_interval_start') the reference time
    falls back to logical_date, then run_after, then now."""
    label, src = C.compute_date_label(context)
    print(f"date_mjd from {src}: {label}")
    return label


# ---------------------------------------------------------------------------
# per-observatory callables
# ---------------------------------------------------------------------------

def resolve_night(tele, **context):
    """Derive this observatory's paths from setup.mjd. Gating on raw-data
    availability lives in the wait_for_raw sensor downstream."""
    p = context["params"]
    mjd = context["ti"].xcom_pull(task_ids="setup.mjd")
    outdir = p.get("outdir") or C.AR_OUTDIR
    paths = C.night_paths(tele, mjd, outdir)
    paths["night_date"] = (context["ti"].xcom_pull(task_ids="setup.date_mjd")
                           or C.night_date_str(mjd))
    print(f"night resolution: tele={tele} mjd={mjd} outdir={outdir} "
          f"raw_status={C.raw_night_status(tele, mjd)}")
    os.makedirs(paths["logdir"], exist_ok=True)
    # Rotate a leftover steps.csv from a previous run of this night (retrigger
    # / past-day rerun): metrics_append censuses this file, and accumulated
    # rows made a clean rerun report phantom failures (seen live on lco 61286:
    # "7 steps, 1 failed" where the failure was the PREVIOUS night's exit-127
    # almanac row). History is preserved in a timestamped sidecar.
    steps = os.path.join(paths["logdir"], "steps.csv")
    if os.path.exists(steps):
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        os.replace(steps, os.path.join(paths["logdir"], f"steps_prev_{stamp}.csv"))
        print(f"rotated stale steps.csv -> steps_prev_{stamp}.csv")
    return paths


def raw_transfer_ready(tele, **context):
    """Sensor poke: gate reduction on raw-transfer COMPLETION, not the date
    (ported semantics of ar_main.py's TransferFileSensor: an ancient night is
    assumed transferred; otherwise wait for the completion marker — there the
    Globus transfer-<mjd>.done file, here the <mjd>.md5sum in the mirror).

    Why ar_main.py worked WITHOUT wiring this sensor: the Globus Utah->CCA
    sync is a fixed daily task at 4:00-4:30 am, so the 7 am ET schedule
    always postdates both observatories' transfers — the data was simply
    already there (confirmed in the raw mirror: transfer markers land the
    UTC day matching the folder name for 365/368 recent apo and 366/366
    lco nights). The sensor is therefore a fast no-op in the normal case
    (first poke returns True) and only earns its keep on the rare late
    transfer (3/368 apo nights lagged by 1+ days).

    params['raw_nowait'] (smoke tests / dags test, where a reschedule wait
    would spin): skip immediately instead of waiting when data is absent."""
    p = context["params"]
    paths = context["ti"].xcom_pull(task_ids=f"{tele}.resolve_night")
    mjd = paths["mjd"]
    if mjd < 59148:  # ar_main.py: "ancient, assume transfer complete"
        print(f"mjd {mjd} < 59148: ancient, assuming transfer complete")
        return True
    status = C.raw_night_status(tele, mjd)
    print(f"raw status for {tele} {mjd}: {status}")
    if status == "ready":
        return True
    if p.get("raw_nowait"):
        raise AirflowSkipException(
            f"raw_nowait: {tele} {mjd} is '{status}' — skipping this "
            "observatory immediately instead of waiting")
    # No-data fast path: if the sync round that would have delivered this
    # night has already run (other observatory's same-night marker, or a
    # later night's marker for this one), a frameless dir means the
    # observatory took no APOGEE data — skip now instead of holding the
    # serial chain for the 12 h timeout. "incomplete" (frames present, no
    # md5sum yet) still waits: that is a transfer in flight.
    if status == "no_data" and C.raw_sync_has_run(tele, mjd):
        raise AirflowSkipException(
            f"{tele} {mjd}: no APOGEE frames, and the nightly sync round "
            "for this night has already completed — observatory took no "
            "data; skipping")
    return False


def select_mode(tele, **context):
    mode = context["params"].get("mode") or C.AR_MODE_DEFAULT
    if mode == "slurm":
        return f"{tele}.slurm_submit"
    return f"{tele}.local.p3d2d"


def madgics_enabled(**context):
    """RUN_MADGICS gate for the local chain. Default TRUE — production ran
    arMADGICS + its workup unconditionally every night (AKS 2026-09-03).
    One knob only: the run_madgics param (no env fallback)."""
    p = context["params"]
    if not p.get("run_madgics", True):
        print("arMADGICS disabled (run_madgics=false)")
        return False
    if not os.path.isdir(C.AR_MADGICS_DIR):
        print(f"run_madgics=true but arMADGICS not found at "
              f"{C.AR_MADGICS_DIR}; skipping")
        return False
    return True


def slurm_submit_and_wait(tele, **context):
    """SLURM execution mode: one sbatch of run_all.sh per observatory,
    restoring ar_main.py's submit_and_wait mechanics (AKS 2026-09-03):
    the sbatch_prefix (-vvv, -D outdir, --mail-type=ALL, the slurm-notify
    Slack-ingest --mail-user), SLACK_CHANNEL exported into the julia layer,
    and the 5-second `squeue -j` poll loop (his explicit call).

    Addition kept beyond ar_main.py: a final `sacct` state check so a
    FAILED/TIMEOUT job fails this task (ar_main.py had no verdict — failure
    surfaced only via the job's own mail). Because accounting lags the
    queue, the verdict (C.slurm_job_verdict) re-queries sacct until the
    state is terminal (>= 30 s apart, up to SACCT_SETTLE_TIMEOUT_S).

    NOTE the current run_all.sh takes the NEW 12-arg signature (R1 layout;
    ar_main.py drove the old 5-arg script whose 5th arg was
    update_sdsscore). Values reproduce production semantics:
    update_sdsscore=false HERE because the DAG's update.sdsscore task
    already ran this morning (double-update avoided, matching his old
    pattern of doing it in the update group); run_madgics defaults true
    (production ran arMADGICS + workup every night)."""
    p = context["params"]
    paths = context["ti"].xcom_pull(task_ids=f"{tele}.resolve_night")
    outdir = paths["outdir"]
    mjd = paths["mjd"]
    logroot = os.path.join(outdir, "slurm_logs")
    os.makedirs(logroot, exist_ok=True)

    cmd = [
        "sbatch", "-vvv", "--parsable",
        "-D", outdir,
        "--mail-type=ALL",
        "--mail-user=slurm_notify-aaaaq7zc7ou7enlnbsexygd3he@sdss5.slack.com",
        f"--job-name=ar_daily_{tele}_{mjd}",
        f"{C.AR_REPO}/scripts/daily/run_all.sh",
        # 12-arg signature of the current run_all.sh:
        tele,                                  # 1  tele
        str(mjd),                              # 2  mjd
        outdir,                                # 3  outdir
        C.AR_MADGICS_DIR,                      # 4  path2arMADGICS
        "false",                               # 5  update_sdsscore (DAG did it)
        p.get("checkpoint_mode", "commit_exists"),  # 6  checkpoint_mode
        "false",                               # 7  run_2d_only
        C.CALDIR_DARKS,                        # 8  caldir_darks
        C.CALDIR_FLATS,                        # 9  caldir_flats
        C.GAIN_READ_CAL_DIR,                   # 10 gain_read_cal_dir
        "false",                               # 11 almanac_clobber_mode
        "true" if p.get("run_madgics", True) else "false",  # 12 run_madgics
    ]
    # SLACK_CHANNEL flows into the sbatch'd julia layer (ar_main.py parity).
    env = os.environ.copy()
    env["SLACK_CHANNEL"] = p.get("slack_channel") or C.AR_SLACK_CHANNEL
    # Exposure-type classifier artifact, passed like the other calibration
    # inputs but by env rather than a 13th positional: run_all.sh's positional
    # signature is consumed by other callers, and SLACK_CHANNEL already
    # establishes the env route into the sbatch'd layer. Setting it (even to "")
    # is what run_all.sh keys on, so production never depends on a default
    # buried in the Julia arg table.
    env["AR_EXP_CLASS_MODEL"] = C.EXP_CLASS_MODEL
    # chip-transform calibration for the sky wavelength fit ("" = off, as before)
    env["AR_CHIP_TRANSFORM_CAL"] = C.CHIP_TRANSFORM_CAL
    print("submitting:", " ".join(cmd))
    res = subprocess.run(cmd, capture_output=True, text=True, env=env)
    if res.returncode != 0:
        raise RuntimeError(f"sbatch failed: {res.stderr}")
    job_id = res.stdout.strip().split(";")[0]
    print(f"submitted Slurm job {job_id}")

    # sbatchAKS-style bookkeeping: snapshot args next to the logs
    with open(os.path.join(logroot, f"slurm-{job_id}.args"), "w") as f:
        f.write("\n".join(cmd) + "\n")

    # ar_main.py wait_for_slurm: poll until the job leaves the queue.
    while True:
        q = subprocess.run(["squeue", "-j", job_id, "-h"],
                           capture_output=True, text=True)
        if q.returncode != 0 or not q.stdout.strip():
            break
        time.sleep(5)  # ar_main.py's 5 s cadence (AKS's explicit call)
    # Accounting lags the queue: re-query sacct (>= 30 s apart, bounded)
    # until the state is terminal instead of trusting one immediate call.
    return C.slurm_job_verdict(job_id)


def start_notification(tele, **context):
    """ar_main.py's per-observatory initial_notification: start-of-reduction
    message with the public exposure-list link — built from PUBLIC_URL_SLUG
    instead of the old token-as-slug coupling."""
    p = context["params"]
    paths = context["ti"].xcom_pull(task_ids=f"{tele}.resolve_night")
    mjd, date = paths["mjd"], paths.get("night_date", "?")
    url = C.public_log_url(tele, mjd)
    link = f" Exposure list can be found <{url}|here>." if url else ""
    label = "[TEST — please ignore] " if p.get("test_label") else ""
    text = (f"{label}Starting reduction for {tele} for SJD {mjd} "
            f"(night of {date}).{link}")
    if p.get("slack_mode", "full") == "full":
        C.slack_post(text, channel=p.get("slack_channel") or None)
    else:
        print(f"[slack suppressed by slack_mode="
              f"{p.get('slack_mode')}] {text}")


def metrics_append(tele, **context):
    paths = context["ti"].xcom_pull(task_ids=f"{tele}.resolve_night")
    if not paths:
        print("no resolved night; no metrics row written")
        return None
    if not C.read_steps_csv(paths["steps_csv"]):
        # Observatory was skipped before any chain step ran (no data /
        # no exposures): recording a 0-step row would (a) pollute the N1
        # tables and (b) make both_observatories_done treat the skipped
        # observatory as cleanly done (caught in the structural test).
        print(f"no chain steps ran for {tele} {paths['mjd']} (observatory "
              "skipped); no metrics row written")
        return None
    return C.collect_night_metrics(paths, context)


def daily_summary(tele, **context):
    """N3 per-observatory summary one-liner (success/warnings census)."""
    p = context["params"]
    paths = context["ti"].xcom_pull(task_ids=f"{tele}.resolve_night")
    label = "[TEST — please ignore] " if p.get("test_label") else ""
    channel = p.get("slack_channel") or None
    if not paths:
        print("summary: night skipped — not posting to Slack")
        return
    row = context["ti"].xcom_pull(task_ids=f"{tele}.metrics_append") or {}
    if not row:
        print("summary: no metrics row (observatory skipped) — not posting")
        return
    mjd = paths["mjd"]
    nfail = row.get("n_steps_failed", "?")
    wall = row.get("total_wall_s", 0)
    # wall = sum of THIS run's step walls; a checkpoint-skip rerun sums to
    # seconds, which the old h/m format displayed as a confusing "0h00m"
    # (AKS 2026-09-04). Show seconds below an hour, and add the DAG's own
    # elapsed time so reruns read sensibly.
    wall_str = (f"{wall // 3600}h{(wall % 3600) // 60:02d}m" if wall >= 3600
                else f"{wall // 60}m{wall % 60:02d}s")
    dag_run = context.get("dag_run")
    start = getattr(dag_run, "start_date", None)
    if start is not None:
        el = int((datetime.now(timezone.utc) - start).total_seconds())
        wall_str += (f" (dag {el // 3600}h{(el % 3600) // 60:02d}m)"
                     if el >= 3600 else f" (dag {el // 60}m{el % 60:02d}s)")
    warn = (row.get("warn_no_wavecal", 0), row.get("warn_no_relflux", 0),
            row.get("warn_other", 0))
    prod = (row.get("n_ar2D", 0), row.get("n_ar1Dcal", 0),
            row.get("n_plots", 0))
    ok = ":white_check_mark:" if nfail == 0 else ":warning:"
    text = (f"{label}{ok} `apogee_daily` {tele} {mjd}: "
            f"{row.get('n_steps', '?')} steps, {nfail} failed, "
            f"wall {wall_str} | warnings "
            f"nowavecal={warn[0]} norelflux={warn[1]} other={warn[2]} | "
            f"products ar2D={prod[0]} ar1Dcal={prod[1]} plots={prod[2]} "
            f"dashboard={'yes' if row.get('dashboard_ok') else 'NO'}")
    if p.get("slack_mode", "full") in ("full", "summary_only"):
        posted = C.slack_post(text, channel=channel)
        print(f"summary posted={posted}: {text}")
    else:
        print(f"[slack off] {text}")


def completion_notification(**context):
    """Single end-of-DAG roll-up after BOTH observatories (restoring
    ar_main.py's completion_notification), implemented via the N1 metrics
    table: posts only when both observatories have clean rows for the
    night."""
    p = context["params"]
    mjd = context["ti"].xcom_pull(task_ids="setup.mjd")
    date = context["ti"].xcom_pull(task_ids="setup.date_mjd")
    label = "[TEST — please ignore] " if p.get("test_label") else ""
    if mjd is None:
        print("completion: no setup.mjd; nothing to do")
        return
    if not C.both_observatories_done(mjd, p.get("run_kind", "daily"),
                                     p.get("metrics_dir") or None):
        print(f"completion: metrics table lacks clean rows for both "
              f"observatories for {mjd}; not posting")
        return
    text = (f"{label}ApogeeReduction pipeline completed successfully for "
            f"SJD {mjd} (night of {date}). Both observatories processed.")
    if p.get("slack_mode", "full") == "full":
        C.slack_post(text, channel=p.get("slack_channel") or None)
    else:
        print(f"[slack suppressed by slack_mode={p.get('slack_mode')}] "
              f"{text}")


# ---------------------------------------------------------------------------
# robustness: stale-run guard, almanac step, tunnel check
# ---------------------------------------------------------------------------

def stale_run_guard(**context):
    """First task of every run: fail a SCHEDULED run whose logical date is
    older than params.stale_run_max_age_h (default AR_STALE_RUN_MAX_AGE_H,
    36 h) before anything else runs. On 2026-11-17 the scheduled run of
    2026-09-27, still 'running' when the service came back, resumed next to
    the day's own run. AirflowFailException: no retries for this verdict.

    Limitation: only a run that had not yet passed this task is caught; a
    run interrupted further along resumes where it stopped."""
    p = context["params"]
    dag_run = context.get("dag_run")
    run_type = getattr(dag_run, "run_type", None)
    logical = context.get("logical_date") or getattr(dag_run,
                                                     "logical_date", None)
    max_age = float(p.get("stale_run_max_age_h", C.STALE_RUN_MAX_AGE_H))
    reason = C.stale_run_reason(run_type, logical,
                                datetime.now(timezone.utc), max_age)
    if reason:
        raise AirflowFailException(reason)
    print(f"run type {run_type}, logical date {logical}: not stale "
          f"(limit {max_age} h for scheduled runs)")


_STMT_H = C.ALMANAC_STATEMENT_TIMEOUT_S / 3600
_EXEC_H = C.ALMANAC_EXECUTION_TIMEOUT_S / 3600
ALMANAC_ANNOUNCE = (
    f"starting. It can sit at 0% CPU for hours in the APOGEE-ID database "
    f"lookup when the DB cache is cold (about 4 h measured 2026-11-17, "
    f"40 s warm); that is NOT a hang - check pg_stat_activity for "
    f"DataFileRead before killing anything. Limits: statement_timeout "
    f"{_STMT_H:g} h (PGOPTIONS), execution_timeout {_EXEC_H:g} h, "
    f"retries {C.ALMANAC_RETRIES}.")


def almanac_cmd(tele: str) -> str:
    """The almanac step (run inside C.step_cmd's logging wrapper).

    * Reuse an existing output only if it really holds this night
      (check_almanac_output.py: raw/<tele>/<mjd>/exposures + fibers), not
      merely because the file exists: a killed run once left a 96-byte stub
      that made every retry skip.
    * Write to "$AF.partial" and move it into place only on rc 0, so a
      killed or failed run never leaves a stub at the real path.
    * PGOPTIONS statement_timeout: a statement stuck server-side errors
      with a clear message instead of waiting until execution_timeout.
    The command is wrapped in bash -c '...': no single quotes inside."""
    af = C.xn("almanac_file", tele)
    mjd = C.xn("mjd", tele)
    check = (f"{C.ALMANAC_PYTHON} "
             f"{C.AR_REPO}/airflow/scripts/check_almanac_output.py")
    return (
        "bash -c 'set -uo pipefail\n"
        f'AF="{af}"\n'
        'TMP="$AF.partial"\n'
        f'CHECK="{check}"\n'
        f'NIGHT="{tele} {mjd}"\n'
        'CLOBBER="{{ params.almanac_clobber }}"\n'
        'if [ "$CLOBBER" != "True" ] && [ "$CLOBBER" != "true" ] '
        '&& [ -e "$AF" ]; then\n'
        '  if $CHECK "$AF" $NIGHT --fibers; then\n'
        '    echo "almanac output already holds this night (skip; set '
        'almanac_clobber to redo)"; exit 0\n'
        '  fi\n'
        '  echo "existing almanac output does not hold this night: '
        'rebuilding it"\n'
        'fi\n'
        'mkdir -p "$(dirname "$AF")"\n'
        'cleanup() { rm -f "$TMP"; }\n'
        'trap cleanup EXIT\n'
        'trap "exit 143" TERM INT\n'
        'rm -f "$TMP"\n'
        # DB identity comes from ~/.almanac/config.yaml + ~/.pgpass;
        # provenance: log the clone's HEAD like the julia pipelines do
        f'echo "almanac clone: $(git -C {C.ALMANAC_DIR} '
        'rev-parse HEAD 2>/dev/null || echo UNKNOWN)"\n'
        f'export PGOPTIONS="${{PGOPTIONS:-}} {C.ALMANAC_PGOPTIONS}"\n'
        'echo "$(date -u +%FT%TZ) almanac start, PGOPTIONS=$PGOPTIONS; the '
        'APOGEE-ID lookup may take hours at 0% CPU on a cold DB cache"\n'
        f'{C.ALMANAC_BIN} -p 12 -v '
        f'--{tele} --mjd-start {mjd} --mjd-end {mjd} '
        '--output "$TMP" --fibers\n'
        'rc=$?\n'
        'echo "$(date -u +%FT%TZ) almanac exited rc=$rc"\n'
        'if [ "$rc" -ne 0 ]; then\n'
        '  echo "almanac failed: partial output discarded, $AF untouched"\n'
        '  exit "$rc"\n'
        'fi\n'
        '$CHECK "$TMP" $NIGHT --fibers || echo "WARNING: almanac succeeded '
        'but its output lacks this night (no exposures in the DB?); '
        'installing it anyway, runlist decides"\n'
        'mv -f "$TMP" "$AF" || exit 1\n'
        'echo "installed $AF"\'')


# tunnel_check: a real `select 1` through the forward, with timeouts. The
# master and the forward belong to the tunnel systemd unit; if the forward
# is missing while the master is up (a tunnel unit installed before it
# carried the LocalForward), add it once over the master as a stopgap and
# say so.
TUNNEL_CHECK_CMD = (
    "set -uo pipefail\n"
    f"UNIT={C.TUNNEL_UNIT}\n"
    f"HOST={C.TUNNEL_HOST}\n"
    f"PORT={C.DB_FORWARD_PORT}\n"
    "fatal() {\n"
    '  echo "FATAL: $1"\n'
    '  echo "The ssh ControlMaster and the DB port-forward are owned by the '
    'systemd user unit $UNIT on $(hostname -s). Check:"\n'
    '  echo "  systemctl --user status $UNIT"\n'
    '  echo "  journalctl --user -u $UNIT -n 50"\n'
    '  echo "If it hit its start limit (e.g. after a password change), fix '
    'the cause, then: systemctl --user reset-failed $UNIT && '
    'systemctl --user start $UNIT"\n'
    "  exit 1\n"
    "}\n"
    "db_ok() {\n"
    f"  PGCONNECT_TIMEOUT=10 PGOPTIONS='-c statement_timeout=15000' "
    f"timeout 30 {C.PSQL_BIN} -X -w -At -h localhost -p \"$PORT\" "
    f"-U {C.DB_USER} -d {C.DB_NAME} -c 'select 1' 2>&1\n"
    "}\n"
    "command -v "
    f"{C.PSQL_BIN} >/dev/null || fatal '{C.PSQL_BIN} not found on PATH "
    "(needed for the DB check; set AR_PSQL)'\n"
    'ssh -O check "$HOST" 2>&1 || fatal "the ssh ControlMaster to $HOST is '
    'down."\n'
    'out=$(db_ok); rc=$?\n'
    'if [ "$rc" -ne 0 ] || [ "$out" != "1" ]; then\n'
    '  echo "select 1 via localhost:$PORT failed (rc=$rc): $out"\n'
    "  if ! timeout 5 bash -c \"echo > /dev/tcp/localhost/$PORT\" "
    "2>/dev/null; then\n"
    '    echo "WARNING: port $PORT is not forwarded although the master is '
    'up. The installed $UNIT predates the LocalForward; reinstall it from '
    'airflow/systemd. Adding the forward over the master as a stopgap."\n'
    f'    ssh -O forward -L "$PORT:{C.DB_REMOTE}" "$HOST" || true\n'
    "    sleep 2\n"
    "    out=$(db_ok); rc=$?\n"
    "  fi\n"
    "fi\n"
    'if [ "$rc" -ne 0 ] || [ "$out" != "1" ]; then\n'
    '  fatal "select 1 through localhost:$PORT failed (rc=$rc): $out"\n'
    "fi\n"
    'echo "tunnel OK: select 1 through localhost:$PORT returned 1"\n')


# ---------------------------------------------------------------------------
# observatory group builder
# ---------------------------------------------------------------------------

def build_observatory_group(tele: str) -> TaskGroup:
    with TaskGroup(group_id=tele) as group:
        t_resolve = PythonOperator(
            task_id="resolve_night",
            python_callable=resolve_night,
            op_kwargs={"tele": tele},
            # none_failed: a SKIPPED upstream observatory passes through
            # (LCO weathered out must not block APO), while a FAILED one
            # still blocks — ar_main.py's serial semantics.
            trigger_rule=TriggerRule.NONE_FAILED,
        )

        t_wait_raw = PythonSensor(
            task_id="wait_for_raw",
            python_callable=raw_transfer_ready,
            op_kwargs={"tele": tele},
            mode="reschedule",
            poke_interval=1800,
            timeout=12 * 3600,
            soft_fail=True,  # timeout (never transferred) -> clean skip
        )

        t_tunnel = BashOperator(
            task_id="tunnel_check",
            bash_command=TUNNEL_CHECK_CMD,
        )

        t_almanac = BashOperator(
            task_id="almanac",
            bash_command=C.step_cmd("almanac", almanac_cmd(tele), tele,
                                    announce=ALMANAC_ANNOUNCE),
            execution_timeout=timedelta(seconds=C.ALMANAC_EXECUTION_TIMEOUT_S),
            retries=C.ALMANAC_RETRIES,
            retry_delay=timedelta(seconds=C.ALMANAC_RETRY_DELAY_S),
        )

        t_runlist = BashOperator(
            task_id="runlist",
            # make_runlist_all.jl exit 16 == "no exposures tonight" -> skip
            bash_command=C.step_cmd(
                "runlist",
                C.julia("scripts/bulk/make_runlist_all.jl",
                        f"--tele {tele} "
                        f"--almanac_file {C.xn('almanac_file', tele)} "
                        f"--output {C.xn('runlist', tele)}"),
                tele,
                skip_exit_codes={16: "no exposures"},
            ),
        )

        t_start_notify = PythonOperator(
            task_id="start_notification",
            python_callable=start_notification,
            op_kwargs={"tele": tele},
        )

        t_select = BranchPythonOperator(
            task_id="select_mode",
            python_callable=select_mode,
            op_kwargs={"tele": tele},
        )

        # ------------------------- LOCAL chain -------------------------
        with TaskGroup(group_id="local"):
            t_p3d2d = BashOperator(
                task_id="p3d2d",
                bash_command=C.step_cmd(
                    "p3d2d",
                    C.julia("pipeline.jl",
                            f"--tele {tele} "
                            f"--runlist {C.xn('runlist', tele)} "
                            f"--outdir {C.xn('outdir', tele)} "
                            f"--runname {C.xn('runname', tele)} "
                            f"--chips RGB "
                            f"--caldir_darks {C.CALDIR_DARKS} "
                            f"--caldir_flats {C.CALDIR_FLATS} "
                            f"--cluster cca "
                            f"--gain_read_cal_dir {C.GAIN_READ_CAL_DIR} "
                            f"--exp_class_model {C.EXP_CLASS_MODEL} "
                            "--checkpoint_mode {{ params.checkpoint_mode }} "
                            "--workers_per_node {{ params.workers }}"),
                    tele,
                ),
            )

            # Fold the post-2D classifier verdicts into the almanac before the
            # flat runlists are built, mirroring scripts/daily/run_all.sh on the
            # SLURM chain. Without it this chain would classify but never filter.
            t_decorate = BashOperator(
                task_id="decorate_exptype",
                bash_command=C.step_cmd(
                    "decorate_exptype",
                    C.julia("scripts/cal/decorate_almanac_exptype.jl",
                            f"--almanac_file {C.xn('almanac_file', tele)} "
                            f"--apred_dir {C.xn('outdir', tele)}apred"),
                    tele,
                ),
            )
            t_p3d2d >> t_decorate
            prev = t_decorate
            for flat_type in ("quartz", "dome"):
                with TaskGroup(group_id=f"{flat_type}_flats") as fg:
                    flatrunlist = (f"{C.xn('outdir', tele)}almanac/runlist_"
                                   f"{flat_type}_{C.xn('runname', tele)}.h5")
                    t_frl = BashOperator(
                        task_id="runlist",
                        bash_command=C.step_cmd(
                            f"runlist_{flat_type}",
                            C.julia("scripts/cal/make_runlist_fiber_flats.jl",
                                    "--almanac_file "
                                    f"{C.xn('almanac_file', tele)} "
                                    f"--tele {tele} --output {flatrunlist} "
                                    f"--flat_type {flat_type}"),
                            tele,
                        ),
                    )
                    t_traces = BashOperator(
                        task_id="traces",
                        bash_command=C.step_cmd(
                            f"traces_{flat_type}",
                            f"bash -c 'mkdir -p {C.xn('outdir', tele)}"
                            f"{flat_type}_flats && "
                            + C.julia(
                                "scripts/cal/make_traces_from_flats.jl",
                                f"--tele {tele} "
                                f"--trace_dir {C.xn('outdir', tele)} "
                                f"--runlist {flatrunlist} "
                                f"--flat_type {flat_type} --slack_quiet true "
                                # FPI guide fibers for the dome-flat trace plots
                                # are read from the night's configurations here,
                                # not from a hardcoded per-telescope pair.
                                f"--almanac_file {C.xn('almanac_file', tele)} "
                                "--checkpoint_mode "
                                "{{ params.checkpoint_mode }}")
                            + "'",
                            tele,
                        ),
                    )
                    t_extract = BashOperator(
                        task_id="extract",
                        bash_command=C.step_cmd(
                            f"extract_{flat_type}",
                            f"bash -c 'mkdir -p {C.xn('outdir', tele)}"
                            "apredrelflux && "
                            + C.julia(
                                "pipeline_2d_1d.jl",
                                f"--tele {tele} --runlist {flatrunlist} "
                                f"--outdir {C.xn('outdir', tele)} "
                                f"--runname {C.xn('runname', tele)} "
                                "--relFlux false --waveSoln false "
                                "--checkpoint_mode "
                                "{{ params.checkpoint_mode }} "
                                "--workers_per_node {{ params.workers }}")
                            + "'",
                            tele,
                        ),
                    )
                    t_relflux = BashOperator(
                        task_id="relflux",
                        bash_command=C.step_cmd(
                            f"relflux_{flat_type}",
                            C.julia("scripts/cal/make_relFlux.jl",
                                    f"--trace_dir {C.xn('outdir', tele)} "
                                    f"--runlist {flatrunlist} "
                                    f"--runname {C.xn('runname', tele)} "
                                    f"--tele {tele}"),
                            tele,
                        ),
                    )
                    t_frl >> t_traces >> t_extract >> t_relflux
                prev >> fg
                prev = fg

            t_p2d1d = BashOperator(
                task_id="p2d1d_full",
                bash_command=C.step_cmd(
                    "p2d1d_full",
                    C.julia("pipeline_2d_1d.jl",
                            f"--tele {tele} "
                            f"--runlist {C.xn('runlist', tele)} "
                            f"--outdir {C.xn('outdir', tele)} "
                            f"--runname {C.xn('runname', tele)} "
                            f"--chip_transform_cal '{C.CHIP_TRANSFORM_CAL}' "
                            "--checkpoint_mode {{ params.checkpoint_mode }} "
                            "--workers_per_node {{ params.workers }}"),
                    tele,
                ),
            )
            t_plots = BashOperator(
                task_id="plots",
                bash_command=C.step_cmd(
                    "plots",
                    C.julia("scripts/daily/plot_all.jl",
                            f"--tele {tele} "
                            f"--runlist {C.xn('runlist', tele)} "
                            f"--outdir {C.xn('outdir', tele)} "
                            f"--runname {C.xn('runname', tele)} "
                            "--chips RGB"),
                    tele,
                ),
            )
            t_dashboard = BashOperator(
                task_id="dashboard",
                bash_command=C.step_cmd(
                    "dashboard",
                    C.julia("scripts/daily/generate_dashboard.jl",
                            f"--mjd {C.xn('mjd', tele)} "
                            f"--outdir {C.xn('outdir', tele)}"),
                    tele,
                ),
            )

            t_madgics_gate = ShortCircuitOperator(
                task_id="madgics_gate",
                python_callable=madgics_enabled,
                ignore_downstream_trigger_rules=False,
            )
            t_madgics = BashOperator(
                task_id="madgics_pipeline",
                bash_command=C.step_cmd(
                    "madgics_pipeline",
                    f"julia +{C.JULIA_VERSION} "
                    f"--project={C.AR_MADGICS_DIR} "
                    f"{C.AR_MADGICS_DIR}pipeline.jl "
                    f"--redux_base {C.xn('outdir', tele)} "
                    f"--almanac_file {C.xn('almanac_file', tele)} "
                    f"--outdir {C.xn('outdir', tele)}arMADGICS/"
                    f"raw_{tele}_{C.xn('mjd', tele)}/",
                    tele,
                ),
            )
            # Spectra workup — the SINGLE workup (AKS 2026-09-04: legacy
            # workup.jl retired from the chain — no row contract, and it
            # deleted the raw batches before this step could read them).
            t_spectra_workup = BashOperator(
                task_id="spectra_workup",
                bash_command=C.step_cmd(
                    "spectra_workup",
                    "bash -c '"
                    f"RW={C.AR_MADGICS_DIR}workup/run_workup.sh; "
                    'if [ ! -f "$RW" ]; then echo "WARNING: spectra workup '
                    "entrypoint $RW not found (arMADGICS PR #26 not in this "
                    'checkout?); skipping"; exit 99; fi; '
                    f"bash $RW {C.xn('outdir', tele)}arMADGICS/"
                    f"raw_{tele}_{C.xn('mjd', tele)}/ "
                    f"{C.xn('outdir', tele)} "
                    f"{C.xn('outdir', tele)}arMADGICS/"
                    f"workup_{tele}_{C.xn('mjd', tele)}/'",
                    tele,
                    skip_exit_codes={99: "entrypoint not present yet"},
                ),
            )
            # Validated cleanup — replaces workup.jl's rm role: the W3
            # validator must PASS (nonzero exit fails this task, batches
            # stay intact) before the raw batch fiber dirs are deleted.
            # AR_KEEP_MADGICS_BATCHES=true keeps them for debugging. A
            # 0-batch night (cal-only) skips both validator and rm.
            t_workup_cleanup = BashOperator(
                task_id="workup_cleanup",
                bash_command=C.step_cmd(
                    "workup_cleanup",
                    "bash -c '"
                    f"RAW={C.xn('outdir', tele)}arMADGICS/"
                    f"raw_{tele}_{C.xn('mjd', tele)}/; "
                    'if [ ! -f "${RAW}batch_info.txt" ] || '
                    'grep -qE "^# Total batches: 0$" "${RAW}batch_info.txt"; '
                    'then echo "no batches (cal-only/empty night) - nothing '
                    'to validate or clean"; exit 0; fi; '
                    f"julia +{C.JULIA_VERSION} "
                    f"--project={C.AR_MADGICS_DIR}workup "
                    f"{C.AR_MADGICS_DIR}workup/validate_workup.jl "
                    f'--rawdir "$RAW" --redux {C.xn("outdir", tele)} --K 500 '
                    f"--out {C.xn('outdir', tele)}arMADGICS/"
                    f"workup_{tele}_{C.xn('mjd', tele)}/W3_report.md; "
                    'if [ "${AR_KEEP_MADGICS_BATCHES:-false}" != "true" ]; '
                    'then echo "W3 passed - cleaning raw batch fiber dirs"; '
                    'rm -rf "$RAW"[0-9][0-9][0-9]; '
                    'else echo "AR_KEEP_MADGICS_BATCHES=true: keeping raw '
                    'batch files"; fi'
                    "'",
                    tele,
                ),
            )
            (prev >> t_p2d1d >> t_plots >> t_dashboard
             >> t_madgics_gate >> t_madgics >> t_spectra_workup
             >> t_workup_cleanup)

        # ------------------------- SLURM mode --------------------------
        t_slurm = PythonOperator(
            task_id="slurm_submit",
            python_callable=slurm_submit_and_wait,
            op_kwargs={"tele": tele},
        )

        # ------------------------- tail --------------------------------
        t_join = EmptyOperator(
            task_id="join",
            trigger_rule=TriggerRule.NONE_FAILED_MIN_ONE_SUCCESS,
        )
        # Gate mirroring the chain outcome: metrics/summary run all_done, so
        # without this the DagRun would be successful even when the chain
        # failed; as a group leaf it also carries the serial edge to the
        # next observatory group (upstream_failed here blocks it).
        t_chain_status = EmptyOperator(
            task_id="chain_status",
            trigger_rule=TriggerRule.NONE_FAILED_MIN_ONE_SUCCESS,
        )
        t_metrics = PythonOperator(
            task_id="metrics_append",
            python_callable=metrics_append,
            op_kwargs={"tele": tele},
            trigger_rule=TriggerRule.ALL_DONE,
        )
        t_summary = PythonOperator(
            task_id="daily_summary",
            python_callable=daily_summary,
            op_kwargs={"tele": tele},
            trigger_rule=TriggerRule.ALL_DONE,
        )

        (t_resolve >> t_wait_raw >> t_tunnel >> t_almanac >> t_runlist
         >> t_start_notify >> t_select)
        t_select >> t_p3d2d
        t_select >> t_slurm
        # workup_cleanup is the local madgics-chain leaf (single-workup design)
        [t_dashboard, t_workup_cleanup, t_slurm] >> t_join
        [t_dashboard, t_workup_cleanup, t_slurm] >> t_chain_status
        t_join >> t_metrics >> t_summary

    return group


# ---------------------------------------------------------------------------
# the DAG
# ---------------------------------------------------------------------------

params = {
    # -1 => auto: ar_main.py's convention, int(Time(data_interval_start).mjd)-1
    "mjd": -1,
    "outdir": C.AR_OUTDIR,
    "mode": C.AR_MODE_DEFAULT,          # "slurm" (default) | "local"
    "workers": C.AR_WORKERS_DEFAULT,
    "update_sdsscore": True,            # O1 punch list #2: default ON
    "sync_logs": True,                  # rsync .log.html mirror
    "almanac_clobber": False,
    "checkpoint_mode": "commit_exists",
    # Production behavior (AKS 2026-09-03): arMADGICS + its workup ran
    # unconditionally every night -> default TRUE, and this param is the
    # ONLY knob (no RUN_MADGICS env fallback).
    "run_madgics": True,
    # slack_mode: "full" (step posts + failure alerts + summaries),
    # "summary_only" (per-obs summaries only), "off"
    "slack_mode": "full",
    "slack_channel": C.AR_SLACK_CHANNEL,
    "test_label": False,
    "run_kind": "daily",                # metrics tag: daily|test|backfill
    "metrics_dir": C.AR_METRICS_DIR,
    # smoke tests: skip a no-data observatory immediately instead of
    # letting wait_for_raw poke/reschedule
    "raw_nowait": False,
    # stale_run_guard: a SCHEDULED run older than this (hours) fails fast;
    # manual triggers are never stale. 0 disables. Env default
    # AR_STALE_RUN_MAX_AGE_H (36).
    "stale_run_max_age_h": C.STALE_RUN_MAX_AGE_H,
}

# The first tasks of a run retry briefly: on 2026-11-17 runs launched ~20 s
# after `airflow standalone` started failed after "Pre Execute", before any
# command ran (most likely before the API server was up). They are all
# idempotent. See ar_common.STARTUP_RETRIES.
STARTUP_RETRY_ARGS = {
    "retries": C.STARTUP_RETRIES,
    "retry_delay": timedelta(seconds=C.STARTUP_RETRY_DELAY_S),
}

with DAG(
    dag_id="apogee_daily",
    description="APOGEE daily reduction, LCO then APO serially "
                "(ApogeeReduction.jl run_all.sh chain)",
    start_date=datetime(2026, 9, 1),
    # ar_main.py's schedule: 7 am ET, which empirically postdates both
    # observatories' overnight transfers (why it never needed the sensor).
    schedule="0 7 * * *",
    catchup=False,
    # 1, not ar_main.py's 2 (2026-11-17): the chain is serial anyway and
    # the metadata DB is SQLite; two active runs (a stale one + today's)
    # ran side by side after the outage and doubled the load of the slow
    # almanac DB lookup. A run still going at the next 7 am makes that
    # run QUEUE; if it waits past the stale limit, stale_run_guard fails it.
    max_active_runs=1,
    dagrun_timeout=timedelta(hours=20),  # two serial observatories
    params=params,
    default_args={
        "retries": 0,
        "on_failure_callback": C.notify_task_failure,
    },
    on_failure_callback=C.notify_dag_failure,
    tags=["apogee", "daily"],
) as dag:
    # ---- stale-run guard: the very first task, before update/setup ----
    t_stale = PythonOperator(
        task_id="stale_run_guard",
        python_callable=stale_run_guard,
        **STARTUP_RETRY_ARGS,
    )

    # ---- update group (ar_main.py template: runs first, once) ----
    with TaskGroup(group_id="update") as group_update:
        t_repo = BashOperator(
            task_id="repo",
            **STARTUP_RETRY_ARGS,
            bash_command=(
                # DELIBERATE NEW BEHAVIOR (AKS 2026-09-03): actually pull the
                # production checkouts. ar_main.py's repo task was status-only.
                #
                # AKS 2026-09-08: all three clones now sit on C.PROD_BRANCH,
                # not main, so merging to main no longer moves production the
                # next morning. Promotion is the explicit
                #     git push origin main:<PROD_BRANCH>
                # done only after validation. Assert the branch before pulling:
                # a clone on anything else means production was moved by hand,
                # and pulling would compound it rather than surface it.
                "prod_pull() {\n"
                '  dir="$1"; name="$2"\n'
                '  cd "$dir" || { echo "FATAL: missing clone $dir"; return 1; }\n'
                f'  want="{C.PROD_BRANCH}"\n'
                '  have=$(git rev-parse --abbrev-ref HEAD)\n'
                '  echo "=== $name: $dir (on $have)"\n'
                '  if [ "$have" != "$want" ]; then\n'
                '    echo "FATAL: $name clone is on \'$have\', expected production'
                ' branch \'$want\'."\n'
                '    echo "Production must not track main. Promote with:'
                ' git push origin main:$want"\n'
                '    return 1\n'
                '  fi\n'
                '  git log --oneline -3\n'
                '  git pull --ff-only origin "$want"\n'
                '  echo "$name now at $(git rev-parse --short HEAD)"\n'
                "}\n"
                f"prod_pull {C.AR_REPO} ApogeeReduction || exit 1\n"
                # arMADGICS clone too (first-night finding 2026-09-04: the
                # madgics step ran a months-stale arM main that predated the
                # raw/-layout almanac handling; AR was pulled, arM never was).
                f"prod_pull {C.AR_MADGICS_DIR} arMADGICS || exit 1\n"
                # almanac clone too (AKS 2026-09-04: run from a local clone
                # like the other pipelines; hash logged at invocation). The
                # editable install tracks code changes automatically; the
                # `uv pip install -e` re-sync only matters when almanac's
                # DEPENDENCIES change, and is a fast no-op otherwise.
                f"prod_pull {C.ALMANAC_DIR} almanac || exit 1\n"
                f"uv pip install -p {os.path.dirname(C.ALMANAC_BIN)}/python "
                f"-e {C.ALMANAC_DIR} --quiet 2>&1 | tail -1 || "
                "echo 'WARNING: almanac env re-sync failed (env unchanged)'\n"),
        )

        t_sdsscore = BashOperator(
            task_id="sdsscore",
            **STARTUP_RETRY_ARGS,
            bash_command=(
                'if [ "{{ params.update_sdsscore }}" != "True" ] && '
                '[ "{{ params.update_sdsscore }}" != "true" ]; then\n'
                '  echo "update_sdsscore disabled"; exit 0\nfi\n'
                f"cd {C.SDSSCORE_DIR}\n"
                "./update.sh\n"
                # update.sh's `submodule foreach` only maintains submodules
                # that are already initialized; a NEW config-range submodule
                # (e.g. lco 10021XXX, found 2026-09-03) never bootstraps
                # without an explicit --init. The .gitmodules URLs are SSH,
                # and this account has no GitHub SSH key, so rewrite to
                # HTTPS. Newly initialized submodules then need the same
                # main-tip checkout update.sh gives the existing ones.
                "git -c url.\"https://github.com/\".insteadOf="
                "\"git@github.com:\" submodule update --init || true\n"
                "git -c url.\"https://github.com/\".insteadOf="
                "\"git@github.com:\" submodule foreach "
                "'git checkout -q main 2>/dev/null; git pull -q || true'\n"),
        )

        t_synclogs = BashOperator(
            task_id="sync_logs",
            **STARTUP_RETRY_ARGS,
            bash_command=(
                'if [ "{{ params.sync_logs }}" != "True" ] && '
                '[ "{{ params.sync_logs }}" != "true" ]; then\n'
                '  echo "sync_logs disabled"; exit 0\nfi\n'
                f"mkdir -p {C.APOGEE_LOGS_DIR}/lco/ "
                f"{C.APOGEE_LOGS_DIR}/apo/\n"
                "time rsync -a --include='*/' --include='*.log.html' "
                f"--exclude='*' {C.RAW_ROOT}/lco/ {C.APOGEE_LOGS_DIR}/lco/\n"
                "time rsync -a --include='*/' --include='*.log.html' "
                f"--exclude='*' {C.RAW_ROOT}/apo/ "
                f"{C.APOGEE_LOGS_DIR}/apo/\n"),
        )

        t_repo >> t_sdsscore >> t_synclogs

    # ---- setup group (ar_main.py: mjd + date_mjd) ----
    with TaskGroup(group_id="setup") as group_setup:
        PythonOperator(task_id="mjd", python_callable=compute_mjd,
                       **STARTUP_RETRY_ARGS)
        PythonOperator(task_id="date_mjd", python_callable=compute_date_mjd,
                       **STARTUP_RETRY_ARGS)

    # ---- observatories, LCO first, serial (ar_main.py order; reason:
    # SQLite metadata DB is single-writer, serial avoids contention) ----
    group_lco = build_observatory_group("lco")
    group_apo = build_observatory_group("apo")

    # none_failed trigger — ar_main.py-parity leaf, and load-bearing:
    # this is the DAG's only leaf, so DagRun state follows it. On a chain
    # failure the upstream chain_status goes upstream_failed -> this leaf
    # does too -> the run is FAILED (with all_done here, a failed chain
    # still produced a "successful" run — caught live in the structural
    # test). Skipped observatories pass through (none_failed).
    t_completion = PythonOperator(
        task_id="completion_notification",
        python_callable=completion_notification,
        trigger_rule=TriggerRule.NONE_FAILED,
    )

    (t_stale >> group_update >> group_setup >> group_lco >> group_apo
     >> t_completion)
