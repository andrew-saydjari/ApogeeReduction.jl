#!/bin/bash
# run_testday.sh — golden-diff regression runner: reduce one (tele, mjd) with a
# given ApogeeReduction.jl checkout, on a workstation (no Slurm), reading raw
# data + almanac from the 2026_05_01 bulk-run inputs.
#
# Usage:
#   ./run_testday.sh <AR_checkout_dir> <tele> <mjd> <outdir>
#
# This mirrors scripts/daily/run_all.sh but:
#   * does NOT run almanac (no Utah tunnel needed): the bulk raw/-layout
#     almanac file (AR_ALMANAC_SRC, one file covering all days) is consumed
#     DIRECTLY — post-R1 (#365) read_almanac_exp_df reads the raw/ layout
#     natively; it is symlinked to <outdir>almanac/allobs_<tele>_<mjd>.h5
#     (the path every downstream stage derives from outdir+runname) and the
#     day is selected via the runlist makers' --mjd flag;
#   * does NOT update sdsscore;
#   * skips plots / dashboard / arMADGICS (science products only — those are
#     what the golden diff compares);
#   * runs either on a workstation (local Distributed workers; Slurm env
#     scrubbed) or inside an sbatch allocation (SlurmClusterManager spanning
#     the whole allocation, mirroring scripts/bulk/run_bulk.sh) — see AR_SLURM
#     below and test/regression/submit_goldens.sh.
#
# Configuration (env vars, all optional; a file named by AR_TESTDAY_CONFIG is
# sourced first if set):
#   AR_ALMANAC_SRC        raw/-layout almanac file containing the day (bulk
#                         multi-day or per-day both work; read directly)
#                         [default: 2026_05_01 bulk allobs_57600_61160.h5]
#   AR_RAW_CLUSTER        --cluster arg for pipeline.jl: "cca" resolves to the
#                         raw mirror /mnt/ceph/users/sdssv/raw/APOGEE; any
#                         other value is interpreted as a base path [cca]
#   AR_CALDIR_DARKS       dark cal dir  [2025_07_31/outdir_ref/]
#   AR_CALDIR_FLATS       flat cal dir  [2025_07_31/outdir_ref/]
#   AR_GAIN_READ_CAL_DIR  gain/readnoise cal dir [2025_07_31/pass_clean/]
#   AR_WORKERS            Distributed workers for pipeline.jl / pipeline_2d_1d
#                         in local mode [24 — headroom on 32-core ccalin051]
#   AR_SLURM              auto | true | false — see mode block below [auto]
#   AR_JULIA_VERSION      juliaup channel [1.11.0, matching run_all.sh]
#   AR_CHECKPOINT_MODE    clobber | commit_exists | commit_same [commit_exists]
#   AR_CHIPS              chips to reduce [RGB]
#   AR_EXP_CLASS_MODEL    exposure-classifier artifact; empty = skip check []
#
# Exit code: 0 on success; the step that failed is the last "--------- x ---------"
# banner in the log (${outdir}logs/run_testday_<tele>_<mjd>.log).

set -e
set -o pipefail

if [ $# -ne 4 ]; then
    echo "usage: $0 <AR_checkout_dir> <tele> <mjd> <outdir>" >&2
    exit 2
fi

harness_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ---- config ----------------------------------------------------------------
if [ -n "${AR_TESTDAY_CONFIG:-}" ]; then
    # shellcheck disable=SC1090
    source "$AR_TESTDAY_CONFIG"
fi
AR_ALMANAC_SRC=${AR_ALMANAC_SRC:-"/mnt/ceph/users/sdssv/work/asaydjari/2026_05_01/outdir/almanac/allobs_57600_61160.h5"}
AR_RAW_CLUSTER=${AR_RAW_CLUSTER:-"cca"}
AR_CALDIR_DARKS=${AR_CALDIR_DARKS:-"/mnt/ceph/users/sdssv/work/asaydjari/2025_07_31/outdir_ref/"}
AR_CALDIR_FLATS=${AR_CALDIR_FLATS:-"/mnt/ceph/users/sdssv/work/asaydjari/2025_07_31/outdir_ref/"}
AR_GAIN_READ_CAL_DIR=${AR_GAIN_READ_CAL_DIR:-"/mnt/ceph/users/sdssv/work/asaydjari/2025_07_31/pass_clean/"}
AR_WORKERS=${AR_WORKERS:-24}
AR_JULIA_VERSION=${AR_JULIA_VERSION:-"1.11.0"}
AR_CHECKPOINT_MODE=${AR_CHECKPOINT_MODE:-"commit_exists"}
AR_CHIPS=${AR_CHIPS:-"RGB"}
AR_EXP_CLASS_MODEL=${AR_EXP_CLASS_MODEL:-""}

base_dir="$(cd "$1" && pwd)"
tele=$2
mjd=$3
outdir=$4
# AR code concatenates paths as outdir * "apred/..." — the trailing slash matters.
case "$outdir" in */) : ;; *) outdir="${outdir}/" ;; esac
mkdir -p "$outdir"
outdir="$(cd "$outdir" && pwd)/"

# Two execution modes (AR_SLURM=auto|true|false, default auto = detect):
#   false — workstation: scrub the Slurm env so the pipelines never try
#           SlurmClusterManager, and pass --workers_per_node explicitly.
#   true  — inside an sbatch allocation (submit_goldens.sh): mirror
#           scripts/bulk/run_bulk.sh exactly — export SLURM_NTASKS =
#           CPUS_ON_NODE * NNODES so every julia stage spans the whole
#           allocation via SlurmClusterManager; --workers_per_node is NOT
#           passed (pipeline default -1 = keep all Slurm workers).
AR_SLURM=${AR_SLURM:-auto}
if [ "$AR_SLURM" = "auto" ]; then
    if [ -n "${SLURM_JOB_ID:-}" ]; then AR_SLURM=true; else AR_SLURM=false; fi
fi
if [ "$AR_SLURM" = "true" ]; then
    SLURM_NTASKS=$(($SLURM_CPUS_ON_NODE * $SLURM_NNODES))
    export SLURM_NTASKS
    workers_args=()
else
    unset SLURM_NTASKS SLURM_JOB_ID SLURM_NNODES SLURM_CPUS_ON_NODE SLURM_JOB_NODELIST
    workers_args=(--workers_per_node "$AR_WORKERS")
fi

runname="allobs_${tele}_${mjd}"
almanac_file=${outdir}almanac/${runname}.h5
runlist=${outdir}almanac/runlist_${runname}.h5
flat_types=("quartz" "dome")
logdir=${outdir}logs
mkdir -p "${outdir}almanac" "$logdir"
logfile=${logdir}/run_testday_${tele}_${mjd}.log

# Log everything (also to the terminal).
exec > >(tee -a "$logfile") 2>&1

echo "run_testday.sh starting $(date -Is) on $(hostname)"
echo "  AR checkout: $base_dir"
echo "  git: $(git -C "$base_dir" rev-parse HEAD 2>/dev/null || echo 'not a git checkout') ($(git -C "$base_dir" rev-parse --abbrev-ref HEAD 2>/dev/null || true))"
echo "  tele=$tele mjd=$mjd outdir=$outdir"
echo "  almanac_src=$AR_ALMANAC_SRC"
echo "  raw cluster=$AR_RAW_CLUSTER darks=$AR_CALDIR_DARKS flats=$AR_CALDIR_FLATS"
echo "  gain/read=$AR_GAIN_READ_CAL_DIR workers=$AR_WORKERS julia=$AR_JULIA_VERSION"
echo "  checkpoint_mode=$AR_CHECKPOINT_MODE chips=$AR_CHIPS exp_class_model='${AR_EXP_CLASS_MODEL}'"
if [ "$AR_SLURM" = "true" ]; then
    echo "  mode=slurm SLURM_NTASKS=$SLURM_NTASKS nodes=${SLURM_NNODES:-?} nodelist=${SLURM_JOB_NODELIST:-?}"
else
    echo "  mode=local workers=$AR_WORKERS"
fi

juliaup add "$AR_JULIA_VERSION" 2>/dev/null || true

LAST_TIME=""
print_elapsed_time() {
    current_seconds=$SECONDS
    if [ -n "$LAST_TIME" ]; then
        diff_seconds=$((current_seconds - LAST_TIME))
        printf 'Time since last step: %dd %dh:%dm:%ds\n' $((diff_seconds/86400)) $((diff_seconds%86400/3600)) $((diff_seconds%3600/60)) $((diff_seconds%60))
    fi
    printf 'Elapsed time: %dd %dh:%dm:%ds\n' $((current_seconds/86400)) $((current_seconds%86400/3600)) $((current_seconds%3600/60)) $((current_seconds%60))
    echo
    echo "--------- $1 ---------"
    echo
    LAST_TIME=$current_seconds
}

# ---- almanac staging (replaces the almanac/tunnel step) ---------------------
# The bulk raw/-layout almanac is consumed directly (all stages open it
# read-only): symlink it to the outdir+runname-derived path that pipeline.jl /
# pipeline_2d_1d.jl / make_relFlux.jl hardcode, and select the day with the
# runlist makers' --mjds flag.
print_elapsed_time "Linking almanac (raw/ layout, read directly)"
if [ ! -f "$AR_ALMANAC_SRC" ]; then
    echo "ERROR: AR_ALMANAC_SRC not found: $AR_ALMANAC_SRC" >&2
    exit 5
fi
ln -sfn "$(realpath "$AR_ALMANAC_SRC")" "$almanac_file"

# ---- runlist ----------------------------------------------------------------
print_elapsed_time "Building Runlist"
set +e
julia +"$AR_JULIA_VERSION" --project="$base_dir" "$base_dir/scripts/bulk/make_runlist_all.jl" \
    --tele "$tele" --almanac_file "$almanac_file" --output "$runlist" --mjds "$mjd"
exit_code=$?
set -e
if [ $exit_code -eq 16 ]; then
    echo "No exposures found for this night. Exiting gracefully."
    exit 0
elif [ $exit_code -ne 0 ]; then
    echo "ERROR: make_runlist_all.jl failed with exit code $exit_code"
    exit $exit_code
fi

# ---- 3D -> 2D / 2Dcal -------------------------------------------------------
print_elapsed_time "Running 3D->2D/2Dcal Pipeline for $tele"
julia +"$AR_JULIA_VERSION" --project="$base_dir" "$base_dir/pipeline.jl" \
    --tele "$tele" --runlist "$runlist" --outdir "$outdir" --runname "$runname" \
    --chips "$AR_CHIPS" --caldir_darks "$AR_CALDIR_DARKS" --caldir_flats "$AR_CALDIR_FLATS" \
    --cluster "$AR_RAW_CLUSTER" --gain_read_cal_dir "$AR_GAIN_READ_CAL_DIR" \
    --checkpoint_mode "$AR_CHECKPOINT_MODE" "${workers_args[@]}" \
    ${AR_EXP_CLASS_MODEL:+--exp_class_model "$AR_EXP_CLASS_MODEL"}

# ---- traces + relFluxing per flat type -------------------------------------
for flat_type in "${flat_types[@]}"; do
    flatrunlist=${outdir}almanac/runlist_${flat_type}_${runname}.h5
    print_elapsed_time "Making runlist for $flat_type Flats"
    julia +"$AR_JULIA_VERSION" --project="$base_dir" "$base_dir/scripts/cal/make_runlist_fiber_flats.jl" \
        --almanac_file "$almanac_file" --tele "$tele" --output "$flatrunlist" --flat_type "$flat_type" \
        --mjds "$mjd"

    print_elapsed_time "Fitting Traces from $flat_type Flats for $tele"
    mkdir -p "${outdir}${flat_type}_flats"
    julia +"$AR_JULIA_VERSION" --project="$base_dir" "$base_dir/scripts/cal/make_traces_from_flats.jl" \
        --tele "$tele" --trace_dir "$outdir" --runlist "$flatrunlist" --flat_type "$flat_type" \
        --slack_quiet true --checkpoint_mode "$AR_CHECKPOINT_MODE"

    print_elapsed_time "Running 2D->1D Pipeline without relFlux for $flat_type Flats for $tele"
    mkdir -p "${outdir}apredrelflux"
    julia +"$AR_JULIA_VERSION" --project="$base_dir" "$base_dir/pipeline_2d_1d.jl" \
        --tele "$tele" --runlist "$flatrunlist" --outdir "$outdir" --runname "$runname" \
        --relFlux false --waveSoln false --checkpoint_mode "$AR_CHECKPOINT_MODE" \
        "${workers_args[@]}"

    print_elapsed_time "Making relFlux for $flat_type Flats"
    julia +"$AR_JULIA_VERSION" --project="$base_dir" "$base_dir/scripts/cal/make_relFlux.jl" \
        --trace_dir "$outdir" --runlist "$flatrunlist" --runname "$runname" --tele "$tele"
done

# ---- 2D -> 1D ---------------------------------------------------------------
print_elapsed_time "Running 2D->1D Pipeline for $tele"
# (The historical mkdir -p ${outdir}wavecal workaround is gone: R1/#365 landed
# the proper mkpath inside skyline_medwavecal_skyline_dither, src/wavecal.jl.)
julia +"$AR_JULIA_VERSION" --project="$base_dir" "$base_dir/pipeline_2d_1d.jl" \
    --tele "$tele" --runlist "$runlist" --outdir "$outdir" --runname "$runname" \
    --checkpoint_mode "$AR_CHECKPOINT_MODE" "${workers_args[@]}"

# ---- warnings census (regression metric, cf. REFACTOR_PLAN v1 §0) ----------
print_elapsed_time "Warnings census"
for pat in \
    "No good pixels found for fiber" \
    "Non-unique or unsorted wavelengths" \
    "no useful arclamp peaks" \
    "Could not find nightly average wave soln" \
    "No fluxing file available" \
    "no useful relfluxing files" \
    "Problem with getting fiber type information" \
    "Failed to get fiber type information" \
    "Exposure-type check" \
    ; do
    n=$(grep -c "$pat" "$logfile" || true)
    echo "census: ${n}x \"$pat\""
done

print_elapsed_time "Test-day run completed"
