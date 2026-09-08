#!/bin/bash
# submit_testbed.sh — run the bulk stage chain over the DR21 200-MJD testbed set.
#
# This is a THIN WRAPPER over submit_goldens.sh. It contains no pipeline logic:
# it sets the testbed-specific inputs, refuses the destructive mistakes, and
# delegates. The body stays in exactly one file on purpose — the SLURM_NTASKS
# 2x bug (fixed 2026-09-07) came from copying a launcher body between scripts
# and letting the copies drift.
#
# WHAT THIS IS NOT
# ----------------
# It is not the golden set. The goldens are 6 hand-picked (tele, mjd) pairs used
# for per-day byte diffs; this is 290 tele-nights sampled from the DR21 range
# (seed 20260903, StableRNG — see 2026_09_03/testbed_prep/testbed_manifest.toml)
# used to exercise the pipeline at scale. Same chain, different question, and a
# different output tree. submit_goldens.sh's own defaults are all golden-shaped
# — golden root, the May bulk almanac, the 6 fixed days, the job name, the
# manifest header — which is why running the testbed through it directly means
# overriding four environment variables correctly every time. That is the
# fragility this wrapper removes.
#
# Submit (from the AR checkout root, so slurm_logs/ lands beside the repo):
#   AR_OUTROOT=/mnt/ceph/users/sdssv/work/asaydjari/<date>/testbed_run/ \
#       sbatchAKS --nodes=8 test/regression/submit_testbed.sh
#
# Preview without Slurm:
#   AR_OUTROOT=/tmp/x test/regression/submit_testbed.sh --dry-run
#
# AR_OUTROOT is REQUIRED and has no default. Everything else below is a default
# that can still be overridden from the environment.
#
# Resuming: submit_goldens.sh runs with checkpoint_mode=commit_exists, so
# resubmitting into the SAME AR_OUTROOT reuses every product that already
# exists and reprocesses only what is missing. The freshness guard below
# recognises a tree this wrapper created (TESTBED_RUN stamp) and allows it;
# it refuses any OTHER non-empty tree.
#
# ------------------------------------------------------------------------------
#SBATCH --partition=cca
#SBATCH --nodes=8
#SBATCH --constraint="[genoa|icelake|rome]"
#SBATCH --mem=900G
#SBATCH --time=1-12:00
#SBATCH --job-name=ar_testbed
#SBATCH --output=slurm_logs/%x_%j.out
# ------------------------------------------------------------------------------

set -e
set -o pipefail

# ---- locate the harness ------------------------------------------------------
# Same three cases as submit_goldens.sh, and the resolved root is exported so the
# delegate skips its own derivation (whose scontrol branch would otherwise read
# THIS script's path, or an orchestrator's, and be right only by accident).
if [ -n "${AR_BASE_DIR:-}" ]; then
    harness_dir="$(cd "$AR_BASE_DIR/test/regression" && pwd)"
elif [ -n "${SLURM_JOB_ID:-}" ]; then
    script_path=$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')
    harness_dir="$(cd "$(dirname "$script_path")" && pwd)"
else
    harness_dir="$(cd "$(dirname "$(realpath "$0")")" && pwd)"
fi
AR_BASE_DIR="$(dirname "$(dirname "$harness_dir")")"
export AR_BASE_DIR

# ---- testbed inputs ----------------------------------------------------------
AR_TESTBED_PREP=${AR_TESTBED_PREP:-"/mnt/ceph/users/sdssv/work/asaydjari/2026_09_03/testbed_prep"}
AR_ALMANAC_SRC=${AR_ALMANAC_SRC:-"${AR_TESTBED_PREP}/almanac_testbed_dr21.h5"}
AR_TESTBED_DAYS_ENV=${AR_TESTBED_DAYS_ENV:-"${AR_TESTBED_PREP}/testbed_days.env"}
AR_RUN_LABEL=${AR_RUN_LABEL:-"# DR21 200-MJD testbed MANIFEST"}
export AR_ALMANAC_SRC AR_RUN_LABEL

# The 290 tele-nights. Sourced, not inlined, so the selection lives in exactly
# one place next to the manifest that records its seed and sha.
if [ -z "${AR_DAYS:-}" ]; then
    if [ ! -f "$AR_TESTBED_DAYS_ENV" ]; then
        echo "REFUSING: day list not found: $AR_TESTBED_DAYS_ENV" >&2
        echo "Set AR_TESTBED_DAYS_ENV, or AR_DAYS directly." >&2
        exit 4
    fi
    # shellcheck disable=SC1090
    . "$AR_TESTBED_DAYS_ENV"
fi
export AR_DAYS

# ---- guards: the mistakes that cost 10 hours or a baseline -------------------
if [ -z "${AR_OUTROOT:-}" ]; then
    echo "REFUSING: AR_OUTROOT is not set." >&2
    echo "submit_goldens.sh would default it to \${AR_GOLDEN_ROOT}/ApogeeReduction.jl@<sha>," >&2
    echo "i.e. it would write ~27 TB of testbed products INTO the golden baseline tree." >&2
    echo "Set AR_OUTROOT to a fresh dated directory, e.g." >&2
    echo "  AR_OUTROOT=/mnt/ceph/users/sdssv/work/asaydjari/\$(date +%Y_%m_%d)/testbed_run/" >&2
    exit 6
fi
case "$AR_OUTROOT" in */) : ;; *) AR_OUTROOT="${AR_OUTROOT}/" ;; esac
export AR_OUTROOT

_golden_root=${AR_GOLDEN_ROOT:-"/mnt/ceph/users/sdssv/work/asaydjari/2026_08_31/golden"}
_golden_abs=$(readlink -m "$_golden_root")
_out_abs=$(readlink -m "$AR_OUTROOT")
if [ "$_out_abs" = "$_golden_abs" ] || case "$_out_abs/" in "$_golden_abs"/*) true;; *) false;; esac; then
    echo "REFUSING: AR_OUTROOT ($_out_abs) is inside the golden baseline tree ($_golden_abs)." >&2
    echo "The goldens are the per-day reference the bulk diffs against; a testbed run" >&2
    echo "must not be written into them." >&2
    exit 7
fi

_stamp="${AR_OUTROOT}TESTBED_RUN"
if [ -e "$_out_abs" ]; then
    if [ ! -d "$_out_abs" ]; then
        echo "REFUSING: AR_OUTROOT exists and is not a directory: $_out_abs" >&2
        exit 8
    elif [ -f "$_stamp" ]; then
        echo "resuming into an existing testbed outroot (checkpoint_mode reuses products):"
        sed 's/^/    /' "$_stamp"
    elif [ -z "$(find "$_out_abs" -mindepth 1 -maxdepth 1 -print -quit 2>/dev/null)" ]; then
        : # exists but empty — fine
    else
        echo "REFUSING: AR_OUTROOT is a non-empty directory that this wrapper did not create:" >&2
        echo "  $_out_abs" >&2
        echo "There is no ${_stamp##*/} stamp in it, so it is somebody else's tree — and" >&2
        echo "checkpoint_mode=commit_exists would silently REUSE whatever products are" >&2
        echo "already there, including any built with superseded calibrations." >&2
        echo "Choose a fresh AR_OUTROOT, or move the existing tree aside." >&2
        exit 9
    fi
fi

# ---- stamp, then delegate ----------------------------------------------------
if [ "${1:-}" != "--dry-run" ]; then
    mkdir -p "$AR_OUTROOT"
    if [ ! -f "$_stamp" ]; then
        {
            echo "run: DR21 200-MJD testbed"
            echo "created: $(date -Is) by ${USER}"
            echo "job: ${SLURM_JOB_ID:-local}"
            echo "code: $(git -C "$AR_BASE_DIR" rev-parse --abbrev-ref HEAD 2>/dev/null) @ $(git -C "$AR_BASE_DIR" rev-parse HEAD 2>/dev/null)"
            echo "almanac: $AR_ALMANAC_SRC"
            echo "days: $AR_TESTBED_DAYS_ENV"
        } > "$_stamp"
    fi
fi

echo "submit_testbed.sh -> $harness_dir/submit_goldens.sh"
echo "  AR_BASE_DIR=$AR_BASE_DIR"
echo "  AR_OUTROOT=$AR_OUTROOT"
echo "  AR_ALMANAC_SRC=$AR_ALMANAC_SRC"
echo "  AR_RUN_LABEL=$AR_RUN_LABEL"
echo "  AR_DAYS: $(echo "$AR_DAYS" | tr ';' '\n' | grep -c .) tele-nights from $AR_TESTBED_DAYS_ENV"

exec "$harness_dir/submit_goldens.sh" "$@"
