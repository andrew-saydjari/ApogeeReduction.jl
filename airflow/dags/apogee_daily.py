"""apogee_daily_{apo,lco} — the APOGEE daily reduction on Airflow (card O3).

One DAG per observatory, tasks split at run_all.sh's natural boundaries
(the per-step table in 2026_08_31/daily_smoke/O1_REPORT.md):

  resolve_night -> update_sdsscore -> tunnel_check -> almanac -> runlist
    -> select_mode -> [ LOCAL chain | slurm_submit ]
    -> join -> metrics_append (N1) -> daily_summary (N3)

LOCAL chain (default; ccalin051 is the group's dedicated node and the chain
fits): p3d2d -> quartz(runlist,traces,extract,relflux)
             -> dome(...) -> p2d1d_full -> plots -> dashboard
             -> madgics_gate -> madgics_pipeline -> madgics_workup

SLURM mode: a single sbatchAKS-style submission of scripts/daily/run_all.sh
(the whole chain in one job). WRITTEN but not enabled by default —
automating sbatch is AKS's call per cluster etiquette; flip
AR_AIRFLOW_MODE=slurm (or conf {"mode": "slurm"}) deliberately.

Trigger a manual/backfill run with e.g.:
  airflow dags trigger apogee_daily_apo --conf '{"mjd": 61284}'
Smoke-test conf (at most ONE Slack message, from the summary task):
  {"mjd": ..., "outdir": ".../o3_dag_test/", "slack_mode": "summary_only",
   "test_label": true, "run_kind": "test", "run_madgics": false}
"""

from __future__ import annotations

import os
import subprocess
import sys
import time
from datetime import datetime, timedelta

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
from airflow.utils.task_group import TaskGroup
from airflow.utils.trigger_rule import TriggerRule

import ar_common as C


# ---------------------------------------------------------------------------
# Python callables
# ---------------------------------------------------------------------------

def resolve_night(tele, **context):
    """Compute the night to reduce + all derived paths; skip cleanly when the
    raw mirror has no (complete) data for it yet."""
    p = context["params"]
    mjd = int(p.get("mjd", -1))
    if mjd < 0:
        mjd = C.auto_mjd()
    outdir = p.get("outdir") or C.AR_OUTDIR
    paths = C.night_paths(tele, mjd, outdir)
    status = C.raw_night_status(tele, mjd)
    print(f"night resolution: tele={tele} mjd={mjd} outdir={outdir} "
          f"raw_status={status}")
    if status == "no_data":
        raise AirflowSkipException(
            f"No raw frames for {tele} {mjd} yet — skipping cleanly.")
    if status == "incomplete":
        raise AirflowSkipException(
            f"Raw transfer for {tele} {mjd} incomplete (no md5sum marker) — "
            "skipping; rerun manually or let O4 backfill catch it.")
    os.makedirs(paths["logdir"], exist_ok=True)
    return paths


def select_mode(**context):
    mode = context["params"].get("mode") or C.AR_MODE_DEFAULT
    if mode == "slurm":
        return "slurm_submit"
    return "local.p3d2d"


def madgics_enabled(**context):
    """RUN_MADGICS gate for the local chain (mirror of run_all.sh arg 12)."""
    p = context["params"]
    run = p.get("run_madgics")
    if run is None:
        run = os.environ.get("RUN_MADGICS", "true").lower() == "true"
    if not run:
        print("arMADGICS disabled (run_madgics=false)")
        return False
    if not os.path.isdir(C.AR_MADGICS_DIR):
        print(f"run_madgics=true but arMADGICS not found at "
              f"{C.AR_MADGICS_DIR}; skipping")
        return False
    return True


def slurm_submit_and_wait(tele, **context):
    """SLURM execution mode: one sbatchAKS-style submission of run_all.sh.

    !!! NOT executed by default and never by automated smoke tests — local
    mode is the recommended default (dedicated node, chain fits). Enabling
    this (mode="slurm") means Airflow submits Slurm jobs unattended, which is
    AKS's call per cluster etiquette. Polling is deliberately gentle:
    squeue -j <id> every 5 minutes (a single job id, read-only).
    """
    p = context["params"]
    paths = context["ti"].xcom_pull(task_ids="resolve_night")
    outdir = paths["outdir"]
    mjd = paths["mjd"]
    logroot = os.path.join(outdir, "slurm_logs")
    os.makedirs(logroot, exist_ok=True)

    # run_all.sh positional args; update_sdsscore/almanac already done by the
    # upstream tasks (checkpoint-aware: run_all.sh skips existing files).
    cmd = [
        "sbatch", "--parsable",
        "--mail-type=ALL",
        "--mail-user=slurm_notify-aaaaq7zc7ou7enlnbsexygd3he@sdss5.slack.com",
        "-D", outdir,
        f"--job-name=ar_daily_{tele}_{mjd}",
        f"{C.AR_REPO}/scripts/daily/run_all.sh",
        tele, str(mjd), outdir, C.AR_MADGICS_DIR,
        "false",                      # update_sdsscore (done upstream)
        p.get("checkpoint_mode", "commit_exists"),
        "false",                      # run_2d_only
        C.CALDIR_DARKS, C.CALDIR_FLATS, C.GAIN_READ_CAL_DIR,
        "false",                      # almanac_clobber_mode
        "true" if p.get("run_madgics", True) else "false",
    ]
    print("submitting:", " ".join(cmd))
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"sbatch failed: {res.stderr}")
    job_id = res.stdout.strip().split(";")[0]
    print(f"submitted Slurm job {job_id}")

    # mimic sbatchAKS bookkeeping: snapshot script + args next to the logs
    with open(os.path.join(logroot, f"slurm-{job_id}.args"), "w") as f:
        f.write("\n".join(cmd) + "\n")

    while True:
        time.sleep(300)  # 5 min; keep polling load on Slurm minimal
        q = subprocess.run(["squeue", "-j", job_id, "-h"],
                           capture_output=True, text=True)
        if q.returncode != 0 or not q.stdout.strip():
            break
    # final state via sacct (read-only, single job)
    acct = subprocess.run(
        ["sacct", "-j", job_id, "-n", "-X", "-o", "State,ExitCode"],
        capture_output=True, text=True)
    state = acct.stdout.strip().split()
    print(f"job {job_id} final: {acct.stdout.strip()}")
    if not state or not state[0].startswith("COMPLETED"):
        raise RuntimeError(f"Slurm job {job_id} ended in state "
                           f"{state[0] if state else 'UNKNOWN'}")
    return job_id


def metrics_append(**context):
    paths = context["ti"].xcom_pull(task_ids="resolve_night")
    if not paths:
        print("no resolved night (run skipped before resolve_night xcom); "
              "no metrics row written")
        return None
    return C.collect_night_metrics(paths, context)


def daily_summary(tele, **context):
    """N3 daily-summary one-liner at DAG end (success/warnings census)."""
    p = context["params"]
    paths = context["ti"].xcom_pull(task_ids="resolve_night")
    label = "[TEST — please ignore] " if p.get("test_label") else ""
    channel = p.get("slack_channel") or None
    if not paths:
        print("summary: night skipped (no data) — not posting to Slack")
        return
    row = context["ti"].xcom_pull(task_ids="metrics_append") or {}
    mjd = paths["mjd"]
    nfail = row.get("n_steps_failed", "?")
    wall = row.get("total_wall_s", 0)
    warn = (row.get("warn_no_wavecal", 0), row.get("warn_no_relflux", 0),
            row.get("warn_other", 0))
    prod = (row.get("n_ar2D", 0), row.get("n_ar1Dcal", 0),
            row.get("n_plots", 0))
    ok = ":white_check_mark:" if nfail == 0 else ":warning:"
    text = (f"{label}{ok} `{context['dag_run'].dag_id}` {tele} {mjd}: "
            f"{row.get('n_steps', '?')} steps, {nfail} failed, "
            f"wall {wall // 3600}h{(wall % 3600) // 60:02d}m | warnings "
            f"nowavecal={warn[0]} norelflux={warn[1]} other={warn[2]} | "
            f"products ar2D={prod[0]} ar1Dcal={prod[1]} plots={prod[2]} "
            f"dashboard={'yes' if row.get('dashboard_ok') else 'NO'}")
    if p.get("slack_mode", "full") in ("full", "summary_only"):
        slack_ok = C.slack_post(text, channel=channel)
        print(f"summary posted={slack_ok}: {text}")
    else:
        print(f"[slack off] {text}")


# ---------------------------------------------------------------------------
# DAG factory
# ---------------------------------------------------------------------------

def build_daily_dag(tele: str, schedule: str) -> DAG:
    params = {
        # -1 => auto (int(current UTC MJD) - 1, i.e. the night that just ended)
        "mjd": -1,
        "outdir": C.AR_OUTDIR,
        "mode": C.AR_MODE_DEFAULT,          # "local" | "slurm"
        "workers": C.AR_WORKERS_DEFAULT,
        "update_sdsscore": True,            # O1 punch list #2: default ON
        "almanac_clobber": False,
        "checkpoint_mode": "commit_exists",
        "run_madgics": True,                # RUN_MADGICS gate
        # slack_mode: "full" (step posts + failure alerts + summary),
        # "summary_only" (exactly one summary message), "off"
        "slack_mode": "full",
        "slack_channel": C.AR_SLACK_CHANNEL,
        "test_label": False,
        "run_kind": "daily",                # metrics tag: daily|test|backfill
        "metrics_dir": C.AR_METRICS_DIR,
    }

    dag = DAG(
        dag_id=f"apogee_daily_{tele}",
        description=f"APOGEE daily reduction for {tele} "
                    "(ApogeeReduction.jl run_all.sh chain)",
        start_date=datetime(2026, 9, 1),
        schedule=schedule,
        catchup=False,
        max_active_runs=1,
        dagrun_timeout=timedelta(hours=16),  # Airflow-3 SLA substitute
        params=params,
        default_args={
            "retries": 0,
            "on_failure_callback": C.notify_task_failure,
        },
        on_failure_callback=C.notify_dag_failure,
        tags=["apogee", "daily", tele],
    )

    with dag:
        t_resolve = PythonOperator(
            task_id="resolve_night",
            python_callable=resolve_night,
            op_kwargs={"tele": tele},
        )

        t_sdsscore = BashOperator(
            task_id="update_sdsscore",
            bash_command=(
                'if [ "{{ params.update_sdsscore }}" != "True" ] && '
                '[ "{{ params.update_sdsscore }}" != "true" ]; then\n'
                '  echo "update_sdsscore disabled"; exit 0\nfi\n'
                f"cd {C.SDSSCORE_DIR}\n"
                "./update.sh\n"),
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
                f'AF="{C.xn("almanac_file")}"\n'
                'if [ -f "$AF" ] && [ "{{ params.almanac_clobber }}" != '
                '"True" ] && [ "{{ params.almanac_clobber }}" != "true" ]; '
                'then\n'
                '  echo "almanac file exists: $AF (skip; set almanac_clobber '
                'to redo)"; exit 0\nfi\n'
                'mkdir -p "$(dirname "$AF")"\n'
                # DB identity comes from ~/.almanac/config.yaml + ~/.pgpass
                f'uvx --from {C.ALMANAC_SOURCE} almanac -p 12 -v '
                f'--{tele} --mjd-start {C.xn("mjd")} --mjd-end {C.xn("mjd")} '
                '--output "$AF" --fibers\'',
            ),
        )

        t_runlist = BashOperator(
            task_id="runlist",
            # make_runlist_all.jl exit 16 == "no exposures tonight" -> skip
            bash_command=C.step_cmd(
                "runlist",
                C.julia("scripts/bulk/make_runlist_all.jl",
                        f"--tele {tele} "
                        f"--almanac_file {C.xn('almanac_file')} "
                        f"--output {C.xn('runlist')}"),
                skip_exit_codes={16: "no exposures"},
            ),
        )

        t_select = BranchPythonOperator(
            task_id="select_mode",
            python_callable=select_mode,
        )

        # ------------------------- LOCAL chain -------------------------
        with TaskGroup(group_id="local") as local_group:
            t_p3d2d = BashOperator(
                task_id="p3d2d",
                bash_command=C.step_cmd(
                    "p3d2d",
                    C.julia("pipeline.jl",
                            f"--tele {tele} --runlist {C.xn('runlist')} "
                            f"--outdir {C.xn('outdir')} "
                            f"--runname {C.xn('runname')} --chips RGB "
                            f"--caldir_darks {C.CALDIR_DARKS} "
                            f"--caldir_flats {C.CALDIR_FLATS} "
                            f"--cluster cca "
                            f"--gain_read_cal_dir {C.GAIN_READ_CAL_DIR} "
                            "--checkpoint_mode {{ params.checkpoint_mode }} "
                            "--workers_per_node {{ params.workers }}"),
                ),
            )

            prev = t_p3d2d
            flat_groups = []
            for flat_type in ("quartz", "dome"):
                with TaskGroup(group_id=f"{flat_type}_flats") as fg:
                    flatrunlist = (f"{C.xn('outdir')}almanac/runlist_"
                                   f"{flat_type}_{C.xn('runname')}.h5")
                    t_frl = BashOperator(
                        task_id="runlist",
                        bash_command=C.step_cmd(
                            f"runlist_{flat_type}",
                            C.julia("scripts/cal/make_runlist_fiber_flats.jl",
                                    f"--almanac_file {C.xn('almanac_file')} "
                                    f"--tele {tele} --output {flatrunlist} "
                                    f"--flat_type {flat_type}"),
                        ),
                    )
                    t_traces = BashOperator(
                        task_id="traces",
                        bash_command=C.step_cmd(
                            f"traces_{flat_type}",
                            f"bash -c 'mkdir -p {C.xn('outdir')}"
                            f"{flat_type}_flats && "
                            + C.julia(
                                "scripts/cal/make_traces_from_flats.jl",
                                f"--tele {tele} --trace_dir {C.xn('outdir')} "
                                f"--runlist {flatrunlist} "
                                f"--flat_type {flat_type} --slack_quiet true "
                                "--checkpoint_mode "
                                "{{ params.checkpoint_mode }}")
                            + "'",
                        ),
                    )
                    t_extract = BashOperator(
                        task_id="extract",
                        bash_command=C.step_cmd(
                            f"extract_{flat_type}",
                            f"bash -c 'mkdir -p {C.xn('outdir')}"
                            "apredrelflux && "
                            + C.julia(
                                "pipeline_2d_1d.jl",
                                f"--tele {tele} --runlist {flatrunlist} "
                                f"--outdir {C.xn('outdir')} "
                                f"--runname {C.xn('runname')} "
                                "--relFlux false --waveSoln false "
                                "--checkpoint_mode "
                                "{{ params.checkpoint_mode }} "
                                "--workers_per_node {{ params.workers }}")
                            + "'",
                        ),
                    )
                    t_relflux = BashOperator(
                        task_id="relflux",
                        bash_command=C.step_cmd(
                            f"relflux_{flat_type}",
                            C.julia("scripts/cal/make_relFlux.jl",
                                    f"--trace_dir {C.xn('outdir')} "
                                    f"--runlist {flatrunlist} "
                                    f"--runname {C.xn('runname')} "
                                    f"--tele {tele}"),
                        ),
                    )
                    t_frl >> t_traces >> t_extract >> t_relflux
                prev >> fg
                prev = fg
                flat_groups.append(fg)

            t_p2d1d = BashOperator(
                task_id="p2d1d_full",
                bash_command=C.step_cmd(
                    "p2d1d_full",
                    C.julia("pipeline_2d_1d.jl",
                            f"--tele {tele} --runlist {C.xn('runlist')} "
                            f"--outdir {C.xn('outdir')} "
                            f"--runname {C.xn('runname')} "
                            "--checkpoint_mode {{ params.checkpoint_mode }} "
                            "--workers_per_node {{ params.workers }}"),
                ),
            )
            t_plots = BashOperator(
                task_id="plots",
                bash_command=C.step_cmd(
                    "plots",
                    C.julia("scripts/daily/plot_all.jl",
                            f"--tele {tele} --runlist {C.xn('runlist')} "
                            f"--outdir {C.xn('outdir')} "
                            f"--runname {C.xn('runname')} --chips RGB"),
                ),
            )
            t_dashboard = BashOperator(
                task_id="dashboard",
                bash_command=C.step_cmd(
                    "dashboard",
                    C.julia("scripts/daily/generate_dashboard.jl",
                            f"--mjd {C.xn('mjd')} "
                            f"--outdir {C.xn('outdir')}"),
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
                    f"--redux_base {C.xn('outdir')} "
                    f"--almanac_file {C.xn('almanac_file')} "
                    f"--outdir {C.xn('outdir')}arMADGICS/"
                    f"raw_{tele}_{C.xn('mjd')}/"),
            )
            t_madgics_workup = BashOperator(
                task_id="madgics_workup",
                bash_command=C.step_cmd(
                    "madgics_workup",
                    f"julia +{C.JULIA_VERSION} "
                    f"--project={C.AR_MADGICS_DIR} "
                    f"{C.AR_MADGICS_DIR}workup.jl "
                    f"--outdir {C.xn('outdir')}arMADGICS/"
                    f"raw_{tele}_{C.xn('mjd')}/"),
            )
            (prev >> t_p2d1d >> t_plots >> t_dashboard
             >> t_madgics_gate >> t_madgics >> t_madgics_workup)

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
        # Leaf mirroring the chain outcome: metrics/summary run all_done, so
        # without this the DagRun would be marked successful even when the
        # chain failed (DagRun state follows leaf-task states). This leaf
        # goes upstream_failed on any chain failure -> DagRun fails ->
        # notify_dag_failure fires.
        t_chain_status = EmptyOperator(
            task_id="chain_status",
            trigger_rule=TriggerRule.NONE_FAILED_MIN_ONE_SUCCESS,
        )
        t_metrics = PythonOperator(
            task_id="metrics_append",
            python_callable=metrics_append,
            trigger_rule=TriggerRule.ALL_DONE,
        )
        t_summary = PythonOperator(
            task_id="daily_summary",
            python_callable=daily_summary,
            op_kwargs={"tele": tele},
            trigger_rule=TriggerRule.ALL_DONE,
        )

        (t_resolve >> t_sdsscore >> t_tunnel >> t_almanac >> t_runlist
         >> t_select)
        t_select >> local_group
        t_select >> t_slurm
        [t_dashboard, t_madgics_workup, t_slurm] >> t_join
        [t_dashboard, t_madgics_workup, t_slurm] >> t_chain_status
        t_join >> t_metrics >> t_summary

    return dag


# APO transfer usually completes by mid-morning ET; LCO a bit later.
# Times are in the instance default_timezone (America/New_York).
dag_apo = build_daily_dag("apo", "0 9 * * *")
dag_lco = build_daily_dag("lco", "0 10 * * *")
