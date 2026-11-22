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

# Execution mode default: "slurm" wraps the chain in a single sbatch of
# scripts/daily/run_all.sh (production dailies, as run historically — AKS's
# explicit call 2026-09-03: he launches the Airflow instance himself, so the
# automation chain is human-initiated and consistent with cluster
# conventions). "local" runs the chain directly on this node (ccalin051) —
# used for testing (AR_AIRFLOW_MODE=local).
AR_MODE_DEFAULT = os.environ.get("AR_AIRFLOW_MODE", "slurm")

# Worker cap for local mode (feeds --workers_per_node and AR_LOCAL_WORKERS).
AR_WORKERS_DEFAULT = int(os.environ.get("AR_WORKERS_DEFAULT", "16"))

# Slack channels: PROD is the posting default (restored 2026-09-03 per AKS —
# ar_main.py always posted to #apogee-reduction-jl), and as of the 2026-09-06
# promotion the deployed airflow_env.sh sets AR_SLACK_CHANNEL to prod
# explicitly, so both the env and the default now say prod. Export
# AR_SLACK_CHANNEL=C07KQ7BJY5P (dev), or pass conf {"slack_channel": ...},
# to demote these routine posts when testing. Note that this knob does NOT
# move the outage alerts (airflow-failure-notify.service,
# check_airflow_heartbeat.sh) — those read AR_SLACK_ALERT_CHANNEL and
# default to prod, so "the orchestrator is down" cannot land in dev.
SLACK_CHANNEL_DEV = "C07KQ7BJY5P"   # apogee-reduction-jl-dev
SLACK_CHANNEL_PROD = "C08B7FKMP16"  # apogee-reduction-jl
AR_SLACK_CHANNEL = os.environ.get("AR_SLACK_CHANNEL", SLACK_CHANNEL_PROD)

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

# Telescope-log mirror (ported from the production ar_main.py update.sync_logs
# task): nightly .log.html pages rsynced here, publicly served through the
# public_www slug dir (public_www/<PUBLIC_URL_SLUG>/apogee_logs -> this dir).
APOGEE_LOGS_DIR = "/mnt/ceph/users/sdssv/work/apogee_logs"

JULIA_VERSION = "1.11.0"  # keep in sync with run_all.sh

# almanac: run from a LOCAL git clone with an editable-install env, exactly
# like the AR/arM pipelines (AKS 2026-09-04). This replaces the uvx
# commit-pin: update.repo pulls the clone daily, invocations log
# `git rev-parse HEAD` for provenance, and the uvx wheel-cache trap (source
# string unchanged -> stale cached build served) cannot occur. One-time
# setup in INSTALL.md (clone + `uv pip install -e`).
ALMANAC_DIR = os.environ.get(
    "ALMANAC_DIR", "/mnt/home/sdssv/gitcode/almanac")
ALMANAC_BIN = os.environ.get(
    "ALMANAC_BIN", "/mnt/home/sdssv/uv_env/almanac_env/bin/almanac")

# Production deployment branch (AKS 2026-09-08). All three pipeline clones
# (AR, arMADGICS, almanac) are checked out on this branch, NOT on main.
#
# Why: `update.repo` git-pulls each clone every morning, so while the clones
# tracked main, any merge to main reached production the next morning with no
# validation step. That is how the 2026-09-08 daily broke — the LCO
# calibration guard merged to main while GAIN_READ_CAL_DIR still pointed at
# the directory the guard rejects, and the two arrived in production together.
#
# Promotion is now explicit and manual: validate, then fast-forward the branch
#     git push origin main:airflow-prod
# and the next `update.repo` picks it up. Nothing else moves production.
# update.repo asserts each clone is actually on this branch and fails loudly
# otherwise, so a clone left on main is reported rather than silently pulled.
PROD_BRANCH = os.environ.get("AR_PROD_BRANCH", "airflow-prod")

CALDIR_DARKS = "/mnt/ceph/users/sdssv/work/asaydjari/2025_07_31/outdir_ref/"
CALDIR_FLATS = "/mnt/ceph/users/sdssv/work/asaydjari/2025_07_31/outdir_ref/"
# 2026-09-06: repointed from 2025_07_31/pass_clean/, whose LCO gain and
# read-noise maps were byte-identical copies of APO's (Utah notebook cell 29
# read `outlst_apo` while writing the `lco` files). Every LCO reduction from
# 2025-06-11 to 2026-09-06 ran its 2D error model on APO detector constants.
# This directory carries APO's maps unchanged (all six bit-identical to the
# old set) plus LCO's real recovered maps; provenance in its MANIFEST.md.
# ar3D.assert_calib_map_telescope_specific hard-fails on the old directory.
GAIN_READ_CAL_DIR = "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_06/pass_clean/"

# Exposure-type classifier artifact. A calibration input like the three paths
# above, and configured in the same place for the same reason: swapping the
# model is a config change here, not a code change in ApogeeReduction.
#
# THE VERSION IS PINNED (classifier_2d v2, the 2D content forest that replaced
# v6 on AKS 2026-11-22, retrained without the read counts; sha256 in metadata/exposure_classifier_2d.sha256; its two
# stored references are pipeline.jl defaults under data/exposure_classifier_2d/).
# Never a glob and never "newest wins": the class
# list, feature layout and decision thresholds all belong to one trained forest,
# and a reduction has to be able to say which forest judged it. Retraining means
# editing this line deliberately.
#
# Set to "" to disable the post-2D exposure-type check. A path that does not
# exist is a hard error at pipeline.jl startup, never a silent skip.
#
# 2026-09-08: MOVED here from 2026_07_14/meta/, a dated analysis directory that
# was a natural scratch-cleanup target. The move was verified byte-identical
# (md5 d669241f030817ab4e68bfd32b117c0d, 175,330,224 bytes) and the OLD PATH NO
# LONGER EXISTS, so any caller still pointing at it fails loudly at startup
# rather than silently skipping the check. See mved.txt in the old directory.
# Only the v6 artifact moved; v1-v5 and the sweep outputs stay behind.
EXP_CLASS_MODEL = ("/mnt/ceph/users/sdssv/work/asaydjari/cal_ref/"
                   "exposure_classifier/exposure_classifier_2d_v2.jld2")

# Comb-measured chip-transform calibration (HDF5 built by
# scripts/cal/chip_transform/build_chip_transform.jl): the sky-line wavelength fit
# fixes the R/B chip scales to the comb's per-fiber values instead of 1. ON by
# default, like EXP_CLASS_MODEL; set to "none" to turn it off deliberately (chip scale 1;
# also the bootstrap mode for rebuilding the file from scratch, see the README below).
# A path that does not exist or cannot be read is a hard error at pipeline_2d_1d.jl
# startup, never a silent fallback. Rebuild before each release (see the README in
# scripts/cal/chip_transform/) and point this at the new versioned file.
CHIP_TRANSFORM_CAL = ("/mnt/ceph/users/sdssv/work/asaydjari/cal_ref/"
                      "chip_transform/chiptransform_pass1dr21_v2_allblock.h5")

# ---------------------------------------------------------------------------
# Robustness knobs (2026-11-17 restart after the 09-26 -> 11-17 outage)
# ---------------------------------------------------------------------------
#
# All are parse-time, env-overridable (set them in airflow_env.sh). A value
# that does not parse falls back to the default with a warning in the DAG
# processor log, rather than breaking the DAG import (a broken import would
# silently stop ALL scheduling, which is worse than a default timeout).


def parse_duration_s(text) -> int:
    """Parse a duration with a unit into whole seconds.

    Accepts PostgreSQL-style units: ``ms``, ``s``, ``min``, ``h``, ``d``
    (e.g. ``"5h"``, ``"300min"``, ``"90s"``), optionally with a space before
    the unit. A bare number is REJECTED on purpose: PostgreSQL reads a bare
    statement_timeout as milliseconds, Airflow reads seconds, and guessing
    wrong here is a 1000x error. Zero and negative values are rejected too:
    zero disables statement_timeout in PostgreSQL, which is the failure mode
    these knobs exist to prevent."""
    m = re.fullmatch(r"\s*(\d+(?:\.\d+)?)\s*(ms|s|min|h|d)\s*", str(text))
    if not m:
        raise ValueError(
            f"duration {text!r}: expected a number with a unit "
            "(ms, s, min, h, d), e.g. '5h' or '300min'")
    scale = {"ms": 0.001, "s": 1, "min": 60, "h": 3600, "d": 86400}
    seconds = float(m.group(1)) * scale[m.group(2)]
    if seconds < 1:
        raise ValueError(f"duration {text!r} must be at least 1 s")
    return int(round(seconds))


def pg_statement_timeout_option(text) -> str:
    """The PGOPTIONS fragment that sets a server-side statement_timeout.

    Normalized to integer milliseconds (PostgreSQL's native unit) so the
    server never has to parse our unit spelling."""
    return f"-c statement_timeout={parse_duration_s(text) * 1000}"


def check_almanac_timeouts(statement_timeout_s: int,
                           execution_timeout_s: int) -> None:
    """The statement timeout must fire before the task's execution_timeout:
    a statement_timeout produces a clear PostgreSQL error in almanac.log,
    whereas execution_timeout just SIGTERMs the task (the "rc -15, no
    message" failure seen on 2026-11-17)."""
    if statement_timeout_s >= execution_timeout_s:
        raise ValueError(
            f"almanac statement_timeout ({statement_timeout_s} s) must be "
            f"shorter than its execution_timeout ({execution_timeout_s} s)")


def _env_duration_s(name: str, default: str) -> int:
    raw = os.environ.get(name, default)
    try:
        return parse_duration_s(raw)
    except ValueError as e:
        print(f"WARNING: {name}={raw!r} ignored ({e}); using {default}")
        return parse_duration_s(default)


def _env_int(name: str, default: int, minimum: int = 0) -> int:
    raw = os.environ.get(name, str(default))
    try:
        val = int(raw)
        if val < minimum:
            raise ValueError(f"must be >= {minimum}")
        return val
    except ValueError as e:
        print(f"WARNING: {name}={raw!r} ignored ({e}); using {default}")
        return default


def _almanac_timeouts() -> tuple[int, int]:
    """(statement_timeout_s, execution_timeout_s) from env, defaults if the
    pair is inconsistent."""
    stmt = _env_duration_s("AR_ALMANAC_STATEMENT_TIMEOUT",
                           ALMANAC_STATEMENT_TIMEOUT_DEFAULT)
    exe = _env_duration_s("AR_ALMANAC_EXECUTION_TIMEOUT",
                          ALMANAC_EXECUTION_TIMEOUT_DEFAULT)
    try:
        check_almanac_timeouts(stmt, exe)
    except ValueError as e:
        print(f"WARNING: {e}; using the defaults "
              f"{ALMANAC_STATEMENT_TIMEOUT_DEFAULT} / "
              f"{ALMANAC_EXECUTION_TIMEOUT_DEFAULT}")
        stmt = parse_duration_s(ALMANAC_STATEMENT_TIMEOUT_DEFAULT)
        exe = parse_duration_s(ALMANAC_EXECUTION_TIMEOUT_DEFAULT)
    return stmt, exe


# almanac timeouts, chosen from the measured lookup timings (DAILY_OPS
# 2026-11-17): the slow step is ONE statement, the APOGEE-ID lookup
# (create_apogee_id_lookup), which takes ~40 s warm, 11-40 min on the
# partly-warm September nights, and an estimated 3.5-4 h fully cold (after
# the 2026-11-16 DB restart). Nights without old-catalogid targets finish
# in 6-19 s.
#   statement_timeout 5h  : above the worst cold estimate (4 h) with 25 %
#                           margin, so a legitimate cold run completes; it
#                           only fires for a statement that is truly stuck
#                           server-side, with a clear PostgreSQL error.
#   execution_timeout 6h  : the statement budget plus 1 h for the rest of
#                           almanac (confSummary parsing, other queries,
#                           HDF5 write). Revisit (e.g. 1h / 45min) once the
#                           DB-side covering index or precomputed lookup
#                           lands and the cold case is gone.
#   retries 1, 10 min     : one retry covers a transient tunnel drop; the
#                           first attempt has partly warmed the DB cache,
#                           so the retry is faster, not a repeat.
# Note the outer bound: dagrun_timeout (20 h) covers both observatories, so
# two worst-case almanac attempts on each would exceed it. In practice the
# LCO lookup warms the cache for APO (same query).
ALMANAC_STATEMENT_TIMEOUT_DEFAULT = "5h"
ALMANAC_EXECUTION_TIMEOUT_DEFAULT = "6h"
ALMANAC_STATEMENT_TIMEOUT_S, ALMANAC_EXECUTION_TIMEOUT_S = _almanac_timeouts()
ALMANAC_PGOPTIONS = pg_statement_timeout_option(
    f"{ALMANAC_STATEMENT_TIMEOUT_S}s")
ALMANAC_RETRIES = _env_int("AR_ALMANAC_RETRIES", 1)
ALMANAC_RETRY_DELAY_S = _env_duration_s("AR_ALMANAC_RETRY_DELAY", "10min")

# Python of the almanac env (has h5py) — runs the almanac output check.
ALMANAC_PYTHON = os.environ.get(
    "ALMANAC_PYTHON", os.path.join(os.path.dirname(ALMANAC_BIN), "python"))

# Startup race (2026-11-17): runs launched ~20 s after `airflow standalone`
# started failed in update.repo after "Pre Execute", before running any
# command — most likely before the API server (which Airflow 3 tasks talk
# to) was ready. The first tasks of the run are idempotent (git pull
# --ff-only, sdsscore update.sh, rsync, pure-Python setup), so they retry:
# 3 x 60 s rides out a slow API-server start; a genuine failure (clone on
# the wrong branch, ...) is reported ~3 min later than before.
STARTUP_RETRIES = _env_int("AR_STARTUP_RETRIES", 3)
STARTUP_RETRY_DELAY_S = _env_duration_s("AR_STARTUP_RETRY_DELAY", "60s")

# Stale runs (2026-11-17): the scheduled run of 2026-09-27 was still
# 'running' when the service came back 51 days later and resumed. A
# SCHEDULED run whose logical date is older than this fails fast in
# stale_run_guard, before update/setup do anything. 36 h lets a run that
# was delayed by an outage of up to a day still go (the 7 am run of day D
# is still processed by early evening of D+1), but never a previous
# night's run on top of today's. Manual triggers are never stale: they are
# how past nights are re-run on purpose. 0 disables the guard.
STALE_RUN_MAX_AGE_H = _env_int("AR_STALE_RUN_MAX_AGE_H", 36)


def stale_run_reason(run_type, logical_date, now,
                     max_age_h: float) -> str | None:
    """Why this run must not proceed, or None if it may.

    Only runs whose type is 'scheduled' are checked (``run_type`` may be a
    DagRunType enum or its string value); manual, backfill and asset runs
    pass. ``logical_date`` and ``now`` must be timezone-aware."""
    rt = str(getattr(run_type, "value", run_type) or "").lower()
    if rt != "scheduled" or not max_age_h or max_age_h <= 0:
        return None
    if logical_date is None:
        return None
    age_h = (now - logical_date).total_seconds() / 3600.0
    if age_h <= max_age_h:
        return None
    return (f"stale scheduled run: logical date {logical_date.isoformat()} "
            f"is {age_h:.1f} h old (limit {max_age_h} h, "
            "AR_STALE_RUN_MAX_AGE_H). This run outlived an outage; it is "
            "failed without running anything so it cannot run on top of "
            "today's run. To process that night, trigger a fresh run: "
            "airflow dags trigger apogee_daily --logical-date <now ISO> "
            "--run-id manual__rerun_<MJD>_<stamp> --conf "
            "'{\"mjd\": <MJD>, \"almanac_clobber\": true}' "
            "(see airflow/INSTALL.md, 'Re-running a failed night').")


# DB tunnel. The almanac DB (operations.sdss.org:5432) is reached through
# localhost:<port>, forwarded by the mwm-tunnel.service systemd user unit,
# which owns the ssh ControlMaster to `mwm` in the foreground and
# establishes the forward itself (ExitOnForwardFailure). tunnel_check runs a
# real `select 1` through it. DB user/name must match the account's
# ~/.almanac/config.yaml (password from ~/.pgpass or PGPASSFILE).
TUNNEL_HOST = os.environ.get("AR_TUNNEL_HOST", "mwm")
TUNNEL_UNIT = os.environ.get("AR_TUNNEL_UNIT", "mwm-tunnel.service")
DB_FORWARD_PORT = _env_int("AR_DB_FORWARD_PORT", 63333, minimum=1)
DB_REMOTE = os.environ.get("AR_DB_REMOTE", "operations.sdss.org:5432")
DB_USER = os.environ.get("AR_DB_USER", "sdss")
DB_NAME = os.environ.get("AR_DB_NAME", "sdss5db")
PSQL_BIN = os.environ.get("AR_PSQL", "psql")

# Hints appended to failure notifications for tasks whose fix needs a human.
FAILURE_HINTS = {
    "tunnel_check": (f"The DB tunnel is down. It is owned by the systemd "
                     f"user unit `{TUNNEL_UNIT}` on ccalin051: check "
                     f"`systemctl --user status {TUNNEL_UNIT}` and "
                     f"`journalctl --user -u {TUNNEL_UNIT}`. If it hit its "
                     "start limit (e.g. a changed password), fix the cause, "
                     f"then `systemctl --user reset-failed {TUNNEL_UNIT}` "
                     f"and `systemctl --user start {TUNNEL_UNIT}`."),
    "stale_run_guard": ("A scheduled run outlived an outage and was failed "
                        "on purpose without running anything. To process "
                        "that night, trigger a fresh run with conf "
                        "{\"mjd\": N, \"almanac_clobber\": true} "
                        "(INSTALL.md, 'Re-running a failed night')."),
    "update_sdsscore": ("sdsscore git pull failed — check "
                        f"{SDSSCORE_DIR} (submodule auth / disk)."),
    "almanac": ("almanac failed. A slow run is NOT a hang: the APOGEE-ID "
                "lookup can sit at 0% CPU for hours on a cold DB cache "
                "(check pg_stat_activity for DataFileRead before killing "
                "anything). rc -15 / SIGTERM with no error = the task's "
                "execution_timeout; 'canceling statement due to statement "
                "timeout' in almanac.log = the PGOPTIONS statement_timeout. "
                "Re-run the night with a fresh trigger, conf "
                "{\"mjd\": N, \"almanac_clobber\": true} (INSTALL.md, "
                "'Re-running a failed night'). Exit 127 means "
                "the PATH in airflow_env.sh doesn't match the account "
                "running Airflow (first-night finding 2026-09-04); "
                "otherwise check the DB tunnel (port 63333), "
                "~/.almanac/config.yaml + ~/.pgpass identity, and whether "
                "confSummary files exist: sdsscore staleness, an "
                "uninitialized new config-range submodule, or an "
                "upstream-empty range repo (seen for lco 10021XXX)."),
}


# ---------------------------------------------------------------------------
# Night / path resolution
# ---------------------------------------------------------------------------

def current_mjd() -> float:
    """UTC MJD right now (no astropy needed)."""
    return time.time() / 86400.0 + 40587.0


def mjd_of(dt) -> float:
    """UTC MJD of a tz-aware datetime — stdlib equivalent of
    astropy `Time(dt).mjd` as used by ar_main.py."""
    return dt.timestamp() / 86400.0 + 40587.0


def night_date_str(mjd: int) -> str:
    """Civil UTC date (YYYY-MM-DD) on which the night's MJD began."""
    from datetime import date, timedelta as td
    return str(date(1858, 11, 17) + td(days=int(mjd)))


def _ctx_get(context, key):
    """context[key] or None: Airflow 3 omits data_interval_* for runs
    without a data interval, and Context objects may raise KeyError."""
    try:
        return context.get(key)
    except (KeyError, AttributeError):
        try:
            return context[key]
        except (KeyError, TypeError):
            return None


def run_reference_time(context, now=None):
    """(datetime, source) the run's night arithmetic is based on:
    data_interval_start (scheduled runs: the trigger time on this
    CronTriggerTimetable instance), else logical_date, else the dag_run's
    run_after, else now (UTC). A manual `airflow dags trigger` without
    --logical-date has no data interval in Airflow 3."""
    for key in ("data_interval_start", "logical_date"):
        val = _ctx_get(context, key)
        if val is not None:
            return val, key
    dag_run = _ctx_get(context, "dag_run")
    for key in ("run_after", "logical_date"):
        val = getattr(dag_run, key, None)
        if val is not None:
            return val, f"dag_run.{key}"
    return (now or datetime.now(timezone.utc)), "now"


def date_label_from_mjd(mjd: int) -> str:
    """setup.date_mjd's label for a given night MJD.

    Same convention as the scheduled arithmetic: a run at time T processes
    mjd = int(mjd(T)) - 1 and labels it (T in America/New_York) - 2 days.
    For the 7 am ET trigger that label is the civil date on which MJD-1
    began, e.g. mjd 61360 (run of 2026-11-17) -> '2026-11-15'."""
    return night_date_str(int(mjd) - 1)


def compute_date_label(context, now=None) -> tuple[str, str]:
    """(label, source description) for setup.date_mjd: from conf `mjd`
    when given (>= 0), else ar_main.py's
    (reference time in America/New_York - 2 days)."""
    from datetime import timedelta as td
    from zoneinfo import ZoneInfo
    params = _ctx_params(context)
    try:
        override = int(params.get("mjd", -1))
    except (TypeError, ValueError):
        override = -1
    if override >= 0:
        return date_label_from_mjd(override), f"param mjd={override}"
    ref, src = run_reference_time(context, now=now)
    label = (ref.astimezone(ZoneInfo("America/New_York"))
             - td(days=2)).strftime("%Y-%m-%d")
    return label, f"{src}={ref}"


def public_log_url(tele: str, mjd: int) -> str | None:
    """Public exposure-list URL for a night, built from PUBLIC_URL_SLUG.

    Ported from ar_main.py's initial_notification, which interpolated
    SLACK_TOKEN into the URL back when the public_www dir was named after
    the token; the dir is now named after the decoupled PUBLIC_URL_SLUG,
    and public_www/<slug>/apogee_logs symlinks to APOGEE_LOGS_DIR.
    """
    slug = os.environ.get("PUBLIC_URL_SLUG", "")
    if not slug:
        return None
    return (f"https://users.flatironinstitute.org/~asaydjari/{slug}/"
            f"apogee_logs/{tele}/{mjd}/{mjd}.log.html")


def both_observatories_done(mjd: int, run_kind: str,
                            metrics_dir: str | None = None) -> bool:
    """True when daily_metrics.csv holds clean rows for BOTH observatories
    for this mjd/run_kind — the metrics-backed implementation behind the
    single completion_notification at the end of the serial DAG."""
    metrics_dir = metrics_dir or AR_METRICS_DIR
    path = os.path.join(metrics_dir, "daily_metrics.csv")
    seen = set()
    try:
        with open(path, newline="") as f:
            for row in csv.DictReader(f):
                if (row.get("mjd") == str(mjd)
                        and row.get("run_kind") == run_kind
                        and row.get("n_steps_failed") == "0"
                        # a skipped observatory must not count as done
                        and row.get("n_steps") not in ("0", "", None)):
                    seen.add(row.get("tele"))
    except OSError:
        return False
    return {"apo", "lco"} <= seen


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


def raw_sync_has_run(tele: str, mjd: int) -> bool:
    """True when the nightly Utah->CCA sync round that would have carried
    this night has demonstrably already completed, so a frameless night dir
    means the observatory took no APOGEE data (not that the transfer is
    late). Evidence, either of: the OTHER observatory's transfer-complete
    marker for the SAME mjd, or this observatory's marker for a LATER night
    (mjd+1/mjd+2 — the run processes yesterday's night, so the next night's
    data normally landed this same morning). Verified live on lco 61284/
    61285 (du Pont took no APOGEE exposures; Utah canonical + staging dirs
    both empty while apo delivered normally)."""
    other = "lco" if tele == "apo" else "apo"
    if os.path.exists(f"{RAW_ROOT}/{other}/{mjd}/{mjd}.md5sum"):
        return True
    for later in (mjd + 1, mjd + 2):
        if os.path.exists(f"{RAW_ROOT}/{tele}/{later}/{later}.md5sum"):
            return True
    return False


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
    dag_id = getattr(ti, "dag_id", "?")
    try_number = getattr(ti, "try_number", "?")
    mjd = p.get("mjd", "?")
    try:  # prefer the resolved night over the raw param (-1 = auto);
        # tasks live in per-observatory TaskGroups, so derive the group's
        # resolve_night task id from this task's prefix.
        grp = task_id.split(".")[0] if "." in str(task_id) else None
        if grp in ("apo", "lco"):
            paths = ti.xcom_pull(task_ids=f"{grp}.resolve_night")
            if paths:
                mjd = paths.get("mjd", mjd)
    except Exception:
        pass
    # hints are keyed by the bare task name (group prefixes stripped)
    hint = FAILURE_HINTS.get(str(task_id).split(".")[-1], "")
    label = "[TEST] " if p.get("test_label") else ""
    text = (f"{label}:rotating_light: `{dag_id}` task `{task_id}` FAILED "
            f"(mjd={mjd}, try={try_number}). :picard_facepalm: {hint}")
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

def xn(key: str, tele: str) -> str:
    """Jinja accessor into the observatory group's resolve_night xcom."""
    return (f"{{{{ ti.xcom_pull(task_ids='{tele}.resolve_night')"
            f"['{key}'] }}}}")


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


def step_cmd(step: str, cmd: str, tele: str,
             skip_exit_codes: dict | None = None,
             announce: str | None = None) -> str:
    """Wrap a chain step: log to <logdir>/<step>.log, append
    '<step>,<rc>,<wall_s>' to steps.csv (N1 raw material), tail the log into
    the Airflow log, propagate rc. skip_exit_codes maps a tool exit code to
    99 (BashOperator's skip_on_exit_code) so e.g. 'no exposures' skips
    cleanly instead of failing.

    announce: text echoed into the AIRFLOW task log before the step starts
    (the step's own output only reaches the Airflow log at the end, as a
    tail). Must not contain double quotes."""
    remap = ""
    for code, _ in (skip_exit_codes or {}).items():
        remap += f'if [ "$rc" -eq {code} ]; then rc=99; fi\n'
    preamble = STEP_PREAMBLE % {"logdir": xn("logdir", tele)}
    if announce:
        if '"' in announce:
            raise ValueError("announce text must not contain double quotes")
        preamble += (f'echo "$(date -u +%FT%TZ) {step}: {announce}"\n'
                     f'echo "{step}: full log (follow with tail -f): '
                     f'$LOGDIR/{step}.log"\n')
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


def census_warnings(logdir: str, since: float | None = None) -> dict:
    """Warning census over the night's step logs. `since` (epoch seconds):
    only logs modified at/after it are counted — step logs persist across
    retriggers of a night (only steps.csv is rotated), so without this a
    rerun's census would include stale logs from failed earlier attempts
    (cousin of the steps.csv phantom-failure bug, 2026-09-04)."""
    counts = {"warn_no_wavecal": 0, "warn_no_relflux": 0, "warn_other": 0}
    for log in sorted(glob.glob(os.path.join(logdir, "*.log"))):
        try:
            if since is not None and os.path.getmtime(log) < since:
                continue
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
    tele = paths["tele"]

    # apred/<mjd>/ is shared by both observatories, so count only this
    # observatory's products (filenames carry the telescope).
    def n(prefix):
        return len(glob.glob(os.path.join(apred, f"{prefix}_{tele}_*")))

    return {
        "n_ar2D": n("ar2D"),
        "n_ar2Dcal": n("ar2Dcal"),
        "n_ar1Dcal": n("ar1Dcal"),
        "n_ar1Dunical": n("ar1Dunical"),
        # plots/<mjd>/ is also shared by both observatories (found live on
        # lco 61286: its row counted apo's 1,135 plots too) — filter by the
        # telescope carried in every plot filename.
        "n_plots": len(glob.glob(os.path.join(plots, f"*_{tele}_*.png"))),
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
    start = getattr(dag_run, "start_date", None)
    warn = census_warnings(paths["logdir"],
                           since=start.timestamp() if start else None)
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
