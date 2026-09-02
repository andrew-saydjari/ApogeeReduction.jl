"""Shared helpers for the ApogeeReduction.jl Airflow DAGs (card O3).

Everything here is import-safe at DAG-parse time: no secrets are required to
parse, and nothing talks to Slack/Slurm/the filesystem until a task actually
runs. Secrets (SLACK_TOKEN, PUBLIC_URL_SLUG) come from the environment of the
Airflow *worker* process (see airflow/INSTALL.md for how they get there).

Slack posting uses stdlib urllib against chat.postMessage so the DAGs have no
hard dependency on apache-airflow-providers-slack or an Airflow Connection.
"""

from __future__ import annotations

import csv
import glob
import json
import os
import re
import time
import urllib.request
from datetime import datetime, timezone

# ---------------------------------------------------------------------------
# Parse-time configuration (env overridable; all have production defaults)
# ---------------------------------------------------------------------------

# Repo checkout the DAG drives. Production: the sdssv service-account clone.
# For staging/tests point AR_REPO at a personal worktree.
AR_REPO = os.environ.get(
    "AR_REPO", "/mnt/home/sdssv/gitcode/ApogeeReduction.jl")

# Production daily output root (nightly products land in <outdir>/apred/<mjd>).
AR_OUTDIR = os.environ.get("AR_OUTDIR", "/mnt/ceph/users/sdssv/work/daily/")

# arMADGICS checkout (run_all.sh default is the sibling of the AR repo).
AR_MADGICS_DIR = os.environ.get(
    "AR_MADGICS_DIR",
    os.path.join(os.path.dirname(AR_REPO.rstrip("/")), "arMADGICS.jl") + "/")

# Execution mode default: "local" runs the chain directly on this node
# (ccalin051, the dedicated group node); "slurm" wraps the chain in a single
# sbatch of scripts/daily/run_all.sh. Automating sbatch is AKS's call
# (cluster etiquette): flip AR_AIRFLOW_MODE=slurm only deliberately.
AR_MODE_DEFAULT = os.environ.get("AR_AIRFLOW_MODE", "local")

# Worker cap for local mode (feeds --workers_per_node and AR_LOCAL_WORKERS).
AR_WORKERS_DEFAULT = int(os.environ.get("AR_WORKERS_DEFAULT", "16"))

# Slack channels: dev is the default everywhere; production posting is an
# explicit opt-in (set AR_SLACK_CHANNEL=C08B7FKMP16 in the Airflow env).
SLACK_CHANNEL_DEV = "C07KQ7BJY5P"   # apogee-reduction-jl-dev
SLACK_CHANNEL_PROD = "C08B7FKMP16"  # apogee-reduction-jl
AR_SLACK_CHANNEL = os.environ.get("AR_SLACK_CHANNEL", SLACK_CHANNEL_DEV)

# Metrics table location (N1 seed). Schema: airflow/README.md + metrics dir.
AR_METRICS_DIR = os.environ.get(
    "AR_METRICS_DIR", "/mnt/ceph/users/sdssv/work/daily/metrics/")

# Heartbeat file the apogee_heartbeat DAG touches; lives on ceph so the
# cross-machine watchdog (airflow/scripts/check_airflow_heartbeat.sh) can see
# it from any node.
AR_HEARTBEAT_FILE = os.environ.get(
    "AR_HEARTBEAT_FILE",
    os.path.join(AR_METRICS_DIR, "airflow_heartbeat.txt"))

RAW_ROOT = "/mnt/ceph/users/sdssv/raw/APOGEE"
SDSSCORE_DIR = os.path.join(RAW_ROOT, "sdsscore")

JULIA_VERSION = "1.11.0"  # keep in sync with run_all.sh

# almanac source pin — keep in sync with scripts/daily/run_all.sh (pin by
# COMMIT, not branch: the uvx wheel cache keys on the source string).
ALMANAC_SOURCE = ("git+https://github.com/andrew-saydjari/almanac.git"
                  "@61e0d51186e172c17c3b33cd586de1b5a61dd2cb")

CALDIR_DARKS = "/mnt/ceph/users/sdssv/work/asaydjari/2025_07_31/outdir_ref/"
CALDIR_FLATS = "/mnt/ceph/users/sdssv/work/asaydjari/2025_07_31/outdir_ref/"
GAIN_READ_CAL_DIR = "/mnt/ceph/users/sdssv/work/asaydjari/2025_07_31/pass_clean/"

# Hints appended to failure notifications for tasks whose fix needs a human.
FAILURE_HINTS = {
    "tunnel_check": ("The Utah ControlMaster (ssh alias `mwm`) is down. Only "
                     "AKS can reopen it: `sshpass -f ~/pass_file ssh -fN mwm` "
                     "on ccalin051."),
    "update_sdsscore": ("sdsscore git pull failed — check "
                        f"{SDSSCORE_DIR} (submodule auth / disk)."),
    "almanac": ("almanac failed — check the DB tunnel (port 63333), "
                "~/.almanac/config.yaml + ~/.pgpass identity, and whether "
                "confSummary files exist (sdsscore staleness)."),
}


# ---------------------------------------------------------------------------
# Night / path resolution
# ---------------------------------------------------------------------------

def current_mjd() -> float:
    """UTC MJD right now (no astropy needed)."""
    return time.time() / 86400.0 + 40587.0


def auto_mjd() -> int:
    """Default night to reduce: the night that ended this morning.

    Empirically (2026-09): when the daily fires in the ET morning of calendar
    day D, the newest complete transfer under RAW_ROOT/<tele>/ is directory
    int(MJD_utc(D)) - 1 for both observatories.
    """
    return int(current_mjd()) - 1


def night_paths(tele: str, mjd: int, outdir: str) -> dict:
    outdir = outdir if outdir.endswith("/") else outdir + "/"
    runname = f"allobs_{tele}_{mjd}"
    return {
        "tele": tele,
        "mjd": mjd,
        "outdir": outdir,
        "runname": runname,
        "raw_dir": f"{RAW_ROOT}/{tele}/{mjd}/",
        "almanac_file": f"{outdir}almanac/{runname}.h5",
        "runlist": f"{outdir}almanac/runlist_{runname}.h5",
        "logdir": f"{outdir}logs/{tele}_{mjd}/",
        "steps_csv": f"{outdir}logs/{tele}_{mjd}/steps.csv",
        "apred_dir": f"{outdir}apred/{mjd}/",
        "plots_dir": f"{outdir}plots/{mjd}/",
    }


def raw_night_status(tele: str, mjd: int) -> str:
    """'ready' | 'no_data' | 'incomplete' for the raw mirror of a night."""
    raw_dir = f"{RAW_ROOT}/{tele}/{mjd}"
    if not os.path.isdir(raw_dir):
        return "no_data"
    frames = glob.glob(os.path.join(raw_dir, "a?R-*.apz"))
    if not frames:
        return "no_data"
    # transfer-complete marker written at the end of the nightly transfer
    if not os.path.exists(os.path.join(raw_dir, f"{mjd}.md5sum")):
        return "incomplete"
    return "ready"


# ---------------------------------------------------------------------------
# Slack (N3)
# ---------------------------------------------------------------------------

def slack_post(text: str, channel: str | None = None,
               token: str | None = None) -> bool:
    """Post to Slack chat.postMessage; degrade to a log line without a token.

    Returns True iff Slack accepted the message.
    """
    token = token or os.environ.get("SLACK_TOKEN", "")
    channel = channel or AR_SLACK_CHANNEL
    if not token:
        print(f"[slack degraded to log] channel={channel} text={text}")
        return False
    payload = json.dumps({"channel": channel, "text": text}).encode()
    req = urllib.request.Request(
        "https://slack.com/api/chat.postMessage",
        data=payload,
        headers={"Authorization": f"Bearer {token}",
                 "Content-Type": "application/json; charset=utf-8"})
    try:
        with urllib.request.urlopen(req, timeout=30) as resp:
            body = json.loads(resp.read().decode())
    except Exception as e:  # network problems must not crash callbacks
        print(f"[slack post failed: {e}] text={text}")
        return False
    if not body.get("ok"):
        print(f"[slack API error: {body.get('error')}] text={text}")
        return False
    print(f"[slack posted] channel={channel} ts={body.get('ts')}")
    return True


def _ctx_params(context) -> dict:
    try:
        return dict(context.get("params") or {})
    except Exception:
        return {}


def notify_task_failure(context):
    """default_args on_failure_callback for the daily DAGs (N3).

    Respects params['slack_mode']: only 'full' posts to Slack; 'summary_only'
    and 'off' degrade to log lines (used by smoke tests so a run produces at
    most ONE Slack message, from the summary task).
    """
    p = _ctx_params(context)
    ti = context.get("task_instance") or context.get("ti")
    task_id = getattr(ti, "task_id", "?")
    dag_id = getattr(getattr(ti, "dag_id", None), "__str__", lambda: "?")()
    dag_id = getattr(ti, "dag_id", "?")
    try_number = getattr(ti, "try_number", "?")
    mjd = p.get("mjd", "?")
    hint = FAILURE_HINTS.get(task_id, "")
    label = "[TEST] " if p.get("test_label") else ""
    text = (f"{label}:rotating_light: `{dag_id}` task `{task_id}` FAILED "
            f"(mjd={mjd}, try={try_number}). {hint}")
    if p.get("slack_mode", "full") == "full":
        slack_post(text, channel=p.get("slack_channel") or None)
    else:
        print(f"[slack suppressed by slack_mode={p.get('slack_mode')}] {text}")


def notify_dag_failure(context):
    """DAG-level on_failure_callback (fires on run failure/timeout).

    NOTE (Airflow 3): sla_miss_callback was REMOVED in Airflow 3.0, so the
    SLA layer here is (a) dagrun_timeout -> this callback, and (b) the
    cross-machine watchdog script checking metrics freshness (INSTALL.md).
    """
    p = _ctx_params(context)
    dag_run = context.get("dag_run")
    dag_id = getattr(dag_run, "dag_id", "?")
    run_id = getattr(dag_run, "run_id", "?")
    label = "[TEST] " if p.get("test_label") else ""
    text = (f"{label}:rotating_light: DAG run `{dag_id}` `{run_id}` failed "
            f"or timed out (dagrun_timeout is the Airflow-3 SLA substitute).")
    if p.get("slack_mode", "full") == "full":
        slack_post(text, channel=p.get("slack_channel") or None)
    else:
        print(f"[slack suppressed by slack_mode={p.get('slack_mode')}] {text}")


# ---------------------------------------------------------------------------
# Local-mode step command builder
# ---------------------------------------------------------------------------

# Jinja fragment pulling the resolve_night xcom dict.
XN = "{{ ti.xcom_pull(task_ids='resolve_night') }}"


def xn(key: str) -> str:
    """Jinja accessor into the resolve_night xcom payload."""
    return f"{{{{ ti.xcom_pull(task_ids='resolve_night')['{key}'] }}}}"


STEP_PREAMBLE = r"""set -uo pipefail
# Scrub any inherited Slurm context so Julia pipelines take the local
# addprocs() branch, never SlurmClusterManager.
unset SLURM_NTASKS SLURM_JOB_ID SLURM_NNODES SLURM_CPUS_ON_NODE 2>/dev/null || true
# Slack from inside the Julia steps only in 'full' mode (plot_all/relFlux/
# dashboard post whenever a token is present).
if [ "{{ params.slack_mode }}" != "full" ]; then unset SLACK_TOKEN; fi
export SLACK_CHANNEL="{{ params.slack_channel }}"
export JULIA_NUM_THREADS=1
export OMP_NUM_THREADS=1
export AR_LOCAL_WORKERS="{{ params.workers }}"
LOGDIR="%(logdir)s"
mkdir -p "$LOGDIR"
"""


def step_cmd(step: str, cmd: str, skip_exit_codes: dict | None = None) -> str:
    """Wrap a chain step: log to <logdir>/<step>.log, append
    '<step>,<rc>,<wall_s>' to steps.csv (N1 raw material), tail the log into
    the Airflow log, propagate rc. skip_exit_codes maps a tool exit code to
    99 (BashOperator's skip_on_exit_code) so e.g. 'no exposures' skips
    cleanly instead of failing."""
    remap = ""
    for code, _ in (skip_exit_codes or {}).items():
        remap += f'if [ "$rc" -eq {code} ]; then rc=99; fi\n'
    preamble = STEP_PREAMBLE % {"logdir": xn("logdir")}
    return (
        preamble
        + f'START=$(date +%s)\n'
        + f'nice -n 10 {cmd} > "$LOGDIR/{step}.log" 2>&1\n'
        + 'rc=$?\n'
        + 'END=$(date +%s)\n'
        + remap
        + f'echo "{step},$rc,$((END-START))" >> "$LOGDIR/steps.csv"\n'
        + f'echo "--- last 40 lines of {step}.log ---"\n'
        + f'tail -n 40 "$LOGDIR/{step}.log"\n'
        + 'exit $rc\n')


def julia(script_rel: str, args: str) -> str:
    """A julia invocation matching run_all.sh conventions."""
    return (f"julia +{JULIA_VERSION} --project={AR_REPO} "
            f"{AR_REPO}/{script_rel} {args}")


# ---------------------------------------------------------------------------
# N1 metrics
# ---------------------------------------------------------------------------

NIGHT_METRICS_FIELDS = [
    "ts_utc", "run_kind", "dag_id", "run_id", "tele", "mjd", "mode",
    "workers", "n_steps", "n_steps_failed", "total_wall_s",
    "warn_no_wavecal", "warn_no_relflux", "warn_other",
    "n_ar2D", "n_ar2Dcal", "n_ar1Dcal", "n_ar1Dunical", "n_plots",
    "dashboard_ok",
]

STEP_METRICS_FIELDS = ["ts_utc", "run_id", "tele", "mjd", "step",
                       "exit_code", "wall_s"]


def _append_csv(path: str, fields: list[str], row: dict):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    new = not os.path.exists(path)
    with open(path, "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        if new:
            w.writeheader()
        w.writerow(row)


def read_steps_csv(steps_csv: str) -> list[dict]:
    rows = []
    if os.path.exists(steps_csv):
        with open(steps_csv) as f:
            for line in f:
                parts = line.strip().split(",")
                if len(parts) == 3:
                    rows.append({"step": parts[0],
                                 "exit_code": int(parts[1]),
                                 "wall_s": int(parts[2])})
    return rows


WARN_CLASSES = {
    "warn_no_wavecal": re.compile(r"No wavecal found", re.I),
    "warn_no_relflux": re.compile(r"Could not find any useful relfluxing",
                                  re.I),
}


def census_warnings(logdir: str) -> dict:
    counts = {"warn_no_wavecal": 0, "warn_no_relflux": 0, "warn_other": 0}
    for log in sorted(glob.glob(os.path.join(logdir, "*.log"))):
        try:
            with open(log, errors="replace") as f:
                for line in f:
                    low = line.lower()
                    if "warning" not in low and "warn:" not in low:
                        continue
                    for key, rx in WARN_CLASSES.items():
                        if rx.search(line):
                            counts[key] += 1
                            break
                    else:
                        counts["warn_other"] += 1
        except OSError:
            pass
    return counts


def product_counts(paths: dict) -> dict:
    apred = paths["apred_dir"]
    plots = paths["plots_dir"]

    def n(pat):
        return len(glob.glob(os.path.join(apred, pat)))

    return {
        "n_ar2D": n("ar2D_*"),
        "n_ar2Dcal": n("ar2Dcal_*"),
        "n_ar1Dcal": n("ar1Dcal_*"),
        "n_ar1Dunical": n("ar1Dunical_*"),
        "n_plots": len(glob.glob(os.path.join(plots, "*.png"))),
        "dashboard_ok": int(os.path.exists(
            os.path.join(plots, "dashboard.html"))),
    }


def collect_night_metrics(paths: dict, context) -> dict:
    """Assemble the per-night metrics row + step rows; append both tables."""
    p = _ctx_params(context)
    dag_run = context.get("dag_run")
    run_id = getattr(dag_run, "run_id", "?")
    dag_id = getattr(dag_run, "dag_id", "?")
    ts = datetime.now(timezone.utc).isoformat(timespec="seconds")

    steps = read_steps_csv(paths["steps_csv"])
    warn = census_warnings(paths["logdir"])
    prod = product_counts(paths)

    metrics_dir = p.get("metrics_dir") or AR_METRICS_DIR
    row = {
        "ts_utc": ts, "run_kind": p.get("run_kind", "daily"),
        "dag_id": dag_id, "run_id": run_id,
        "tele": paths["tele"], "mjd": paths["mjd"],
        "mode": p.get("mode", AR_MODE_DEFAULT),
        "workers": p.get("workers", AR_WORKERS_DEFAULT),
        "n_steps": len(steps),
        "n_steps_failed": sum(1 for s in steps
                              if s["exit_code"] not in (0, 99)),
        "total_wall_s": sum(s["wall_s"] for s in steps),
        **warn, **prod,
    }
    _append_csv(os.path.join(metrics_dir, "daily_metrics.csv"),
                NIGHT_METRICS_FIELDS, row)
    for s in steps:
        _append_csv(os.path.join(metrics_dir, "step_metrics.csv"),
                    STEP_METRICS_FIELDS,
                    {"ts_utc": ts, "run_id": run_id, "tele": paths["tele"],
                     "mjd": paths["mjd"], **s})
    print("metrics row:", json.dumps(row))
    return row
