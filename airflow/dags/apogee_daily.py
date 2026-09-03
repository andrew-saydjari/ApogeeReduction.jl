"""apogee_daily — the APOGEE daily reduction on Airflow (card O3).

ONE DAG, observatories chained SERIALLY (LCO then APO), restoring the
production ar_main.py structure per AKS review 2026-09-03: the Airflow
metadata DB is SQLite (single-writer), so parallel per-observatory DAG runs
contend; serial is deliberate.

  update(repo -> sdsscore -> sync_logs) -> setup(mjd, date_mjd)
    -> [ LCO group ] -> [ APO group ] -> completion_notification

Each observatory group:
  resolve_night -> wait_for_raw -> tunnel_check -> almanac -> runlist
    -> start_notification -> select_mode -> [ LOCAL chain | slurm_submit ]
    -> join -> metrics_append (N1) -> daily_summary (N3)   (+ chain_status)

LOCAL chain (testing): p3d2d -> quartz(runlist,traces,extract,relflux)
             -> dome(...) -> p2d1d_full -> plots -> dashboard
             -> madgics_gate -> madgics_pipeline -> madgics_workup
             -> spectra_workup (arMADGICS PR #26 entrypoint)

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
from datetime import datetime, timedelta
from zoneinfo import ZoneInfo

# make sibling ar_common importable regardless of how the bundle is parsed
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from airflow import DAG
from airflow.exceptions import AirflowSkipException
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
    dis = context["data_interval_start"]
    mjd = int(C.mjd_of(dis)) - 1
    print(f"mjd from data_interval_start={dis}: {mjd}")
    return mjd


def compute_date_mjd(**context):
    """ar_main.py setup.date_mjd, verbatim semantics:
        (data_interval_start astimezone America/New_York
         - timedelta(days=2)).strftime('%Y-%m-%d')
    (zoneinfo replaces pytz)."""
    dis = context["data_interval_start"]
    date = (dis.astimezone(ZoneInfo("America/New_York"))
            - timedelta(days=2)).strftime("%Y-%m-%d")
    print(f"date_mjd from data_interval_start={dis}: {date}")
    return date


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
    surfaced only via the job's own mail).

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
    acct = subprocess.run(
        ["sacct", "-j", job_id, "-n", "-X", "-o", "State,ExitCode"],
        capture_output=True, text=True)
    state = acct.stdout.strip().split()
    print(f"job {job_id} final: {acct.stdout.strip()}")
    if not state or not state[0].startswith("COMPLETED"):
        raise RuntimeError(f"Slurm job {job_id} ended in state "
                           f"{state[0] if state else 'UNKNOWN'}")
    return job_id


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
    warn = (row.get("warn_no_wavecal", 0), row.get("warn_no_relflux", 0),
            row.get("warn_other", 0))
    prod = (row.get("n_ar2D", 0), row.get("n_ar1Dcal", 0),
            row.get("n_plots", 0))
    ok = ":white_check_mark:" if nfail == 0 else ":warning:"
    text = (f"{label}{ok} `apogee_daily` {tele} {mjd}: "
            f"{row.get('n_steps', '?')} steps, {nfail} failed, "
            f"wall {wall // 3600}h{(wall % 3600) // 60:02d}m | warnings "
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
            bash_command=(
                "set -uo pipefail\n"
                "# The almanac DB (operations.sdss.org:5432) is reached via\n"
                "# localhost:63333 forwarded over the persistent mwm\n"
                "# ControlMaster. Only AKS can (re)open the master itself.\n"
                "if ! ssh -O check mwm; then\n"
                "  echo 'FATAL: mwm ControlMaster is DOWN. Only AKS can "
                "reopen it: sshpass -f ~/pass_file ssh -fN mwm'\n"
                "  exit 1\n"
                "fi\n"
                "if ! timeout 5 bash -c 'echo > /dev/tcp/localhost/63333' "
                "2>/dev/null; then\n"
                "  echo 'port 63333 not forwarded yet; adding forward via "
                "ControlMaster'\n"
                "  ssh -O forward -L 63333:operations.sdss.org:5432 mwm "
                "|| true\n"
                "  sleep 2\n"
                "  timeout 5 bash -c 'echo > /dev/tcp/localhost/63333' || {\n"
                "    echo 'FATAL: could not establish the 63333 forward "
                "over the mwm ControlMaster'; exit 1; }\n"
                "fi\n"
                "echo 'tunnel OK'\n"),
        )

        t_almanac = BashOperator(
            task_id="almanac",
            bash_command=C.step_cmd(
                "almanac",
                'bash -c \'set -uo pipefail\n'
                f'AF="{C.xn("almanac_file", tele)}"\n'
                'if [ -f "$AF" ] && [ "{{ params.almanac_clobber }}" != '
                '"True" ] && [ "{{ params.almanac_clobber }}" != "true" ]; '
                'then\n'
                '  echo "almanac file exists: $AF (skip; set almanac_clobber '
                'to redo)"; exit 0\nfi\n'
                'mkdir -p "$(dirname "$AF")"\n'
                # DB identity comes from ~/.almanac/config.yaml + ~/.pgpass
                f'uvx --from {C.ALMANAC_SOURCE} almanac -p 12 -v '
                f'--{tele} --mjd-start {C.xn("mjd", tele)} '
                f'--mjd-end {C.xn("mjd", tele)} '
                '--output "$AF" --fibers\'',
                tele,
            ),
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
                            "--checkpoint_mode {{ params.checkpoint_mode }} "
                            "--workers_per_node {{ params.workers }}"),
                    tele,
                ),
            )

            prev = t_p3d2d
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
            t_madgics_workup = BashOperator(
                task_id="madgics_workup",
                bash_command=C.step_cmd(
                    "madgics_workup",
                    f"julia +{C.JULIA_VERSION} "
                    f"--project={C.AR_MADGICS_DIR} "
                    f"{C.AR_MADGICS_DIR}workup.jl "
                    f"--outdir {C.xn('outdir', tele)}arMADGICS/"
                    f"raw_{tele}_{C.xn('mjd', tele)}/",
                    tele,
                ),
            )
            # Spectra workup — first-class chain step per AKS 2026-09-03.
            # Wired against the PLANNED contract `run_workup.sh <rawdir>
            # <redux> <outdir>` (arMADGICS PR #26); existence-guarded.
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
            (prev >> t_p2d1d >> t_plots >> t_dashboard
             >> t_madgics_gate >> t_madgics >> t_madgics_workup
             >> t_spectra_workup)

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
        [t_dashboard, t_spectra_workup, t_slurm] >> t_join
        [t_dashboard, t_spectra_workup, t_slurm] >> t_chain_status
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
    max_active_runs=2,                   # ar_main.py's setting
    dagrun_timeout=timedelta(hours=20),  # two serial observatories
    params=params,
    default_args={
        "retries": 0,
        "on_failure_callback": C.notify_task_failure,
    },
    on_failure_callback=C.notify_dag_failure,
    tags=["apogee", "daily"],
) as dag:
    # ---- update group (ar_main.py template: runs first, once) ----
    with TaskGroup(group_id="update") as group_update:
        t_repo = BashOperator(
            task_id="repo",
            bash_command=(
                f"cd {C.AR_REPO}\n"
                "echo '=== Git Status ==='\n"
                "git status\n"
                "echo '=== Git Log (last 3 commits) ==='\n"
                "git log --oneline -3\n"
                "echo '=== Git Remote Status ==='\n"
                "git remote -v\n"
                # DELIBERATE NEW BEHAVIOR (AKS 2026-09-03): actually pull the
                # production checkout on its branch. ar_main.py's repo task
                # was status-only (pulling was effectively off).
                "echo '=== git pull --ff-only ==='\n"
                "git pull --ff-only\n"),
        )

        t_sdsscore = BashOperator(
            task_id="sdsscore",
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
        PythonOperator(task_id="mjd", python_callable=compute_mjd)
        PythonOperator(task_id="date_mjd", python_callable=compute_date_mjd)

    # ---- observatories, LCO first, serial (ar_main.py order; reason:
    # SQLite metadata DB is single-writer, serial avoids contention) ----
    group_lco = build_observatory_group("lco")
    group_apo = build_observatory_group("apo")

    t_completion = PythonOperator(
        task_id="completion_notification",
        python_callable=completion_notification,
        trigger_rule=TriggerRule.ALL_DONE,
    )

    group_update >> group_setup >> group_lco >> group_apo >> t_completion
