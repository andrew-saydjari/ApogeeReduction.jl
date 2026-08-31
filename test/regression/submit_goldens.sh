#!/bin/bash
# submit_goldens.sh — generate the T1 golden-diff baselines for all test days
# in ONE multinode Slurm allocation, mirroring scripts/bulk/run_bulk.sh
# (SLURM_NTASKS = CPUS_ON_NODE * NNODES exported once; every julia stage spans
# the allocation via SlurmClusterManager).
#
# Submit (the user does this, from the AR checkout root so slurm_logs/ lands
# beside the repo):
#   sbatchAKS test/regression/submit_goldens.sh
#
# Validate without Slurm:
#   bash -n test/regression/submit_goldens.sh
#   test/regression/submit_goldens.sh --dry-run
#
# Prerequisite (already staged): per-day almanac inputs extracted from the
# 2026_05_01 bulk almanac by extract_almanac_day.jl, under
# $AR_ALMANAC_STAGE/allobs_<tele>_<mjd>.h5 — the job has zero DB/tunnel
# dependency.
#
# Outputs: $AR_GOLDEN_ROOT/ApogeeReduction.jl@<shortsha>/<tele>_<mjd>/
# plus a MANIFEST.md at the @sha root recording branch/sha/date/config and a
# per-day exit-status table (a failed day does NOT abort the remaining days).
#
# Decisions baked in (2026-08-31): goldens are generated WITHOUT the
# exposure-classifier model (matches run_all.sh production behavior; the
# classifier enters later via R1/O5) — AR_EXP_CLASS_MODEL stays empty.

# ------------------------------------------------------------------------------
#SBATCH --partition=cca
#SBATCH --qos=cca
#SBATCH --nodes=2
#SBATCH --constraint="[genoa|icelake|rome]"
#SBATCH --mem=900G
#SBATCH --time=1-12:00
#SBATCH --job-name=ar_goldens
#SBATCH --output=slurm_logs/%x_%j.out
# ------------------------------------------------------------------------------

DRY_RUN=false
if [ "${1:-}" = "--dry-run" ]; then
    DRY_RUN=true
fi

set -e
set -o pipefail

# mirror run_bulk.sh: locate the script (sbatch copies it into spool, so use
# scontrol under Slurm), then the repo root two levels up from test/regression/
if [ -n "${SLURM_JOB_ID:-}" ]; then
    script_path=$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')
    hostname
    echo "$SLURM_JOB_NODELIST"
else
    script_path=$(realpath "$0")
fi
harness_dir="$(cd "$(dirname "$script_path")" && pwd)"
base_dir="$(dirname "$(dirname "$harness_dir")")"
echo "base_dir: $base_dir"

# ---- config (env-overridable) ----------------------------------------------
AR_GOLDEN_ROOT=${AR_GOLDEN_ROOT:-"/mnt/ceph/users/sdssv/work/asaydjari/2026_08_31/golden"}
AR_ALMANAC_STAGE=${AR_ALMANAC_STAGE:-"${AR_GOLDEN_ROOT}/almanac_inputs"}
AR_JULIA_VERSION=${AR_JULIA_VERSION:-"1.11.0"}
AR_CHECKPOINT_MODE=${AR_CHECKPOINT_MODE:-"commit_exists"}
export AR_JULIA_VERSION AR_CHECKPOINT_MODE
export AR_EXP_CLASS_MODEL=""   # decision: goldens without the classifier

# The test days (REFACTOR_PLAN v1 §4.2): 6 runs covering the 5 fixed test days
# + both telescopes at 60000. The A3 gapped-almanac day joins later (SUBMIT.md).
days=(
    "apo 60000"
    "lco 60000"
    "apo 58588"
    "apo 59337"
    "apo 58011"
    "apo 59429"
)

# sha of the code the run will execute = the golden identity
if ! git -C "$base_dir" diff --quiet || ! git -C "$base_dir" diff --cached --quiet; then
    echo "WARNING: $base_dir has uncommitted changes — goldens should come from a clean, pushed commit." >&2
    $DRY_RUN || { echo "Refusing to generate goldens from a dirty tree." >&2; exit 3; }
fi
sha=$(git -C "$base_dir" rev-parse HEAD)
branch=$(git -C "$base_dir" rev-parse --abbrev-ref HEAD)
shortsha=${sha:0:7}
outroot=${AR_GOLDEN_ROOT}/ApogeeReduction.jl@${shortsha}
manifest=${outroot}/MANIFEST.md

# ---- validation common to dry-run and real run ------------------------------
missing=0
for day in "${days[@]}"; do
    read -r tele mjd <<< "$day"
    staged=${AR_ALMANAC_STAGE}/allobs_${tele}_${mjd}.h5
    if [ ! -f "$staged" ]; then
        echo "MISSING staged almanac: $staged" >&2
        missing=1
    fi
done
if [ $missing -ne 0 ]; then
    echo "Stage them first: julia --project=$harness_dir $harness_dir/extract_almanac_day.jl <bulk_almanac> <tele> <mjd> <staged_path>" >&2
    exit 4
fi

if $DRY_RUN; then
    echo "=== DRY RUN — commands that the Slurm job will execute ==="
    echo "git identity: ${branch} @ ${sha}"
    echo "golden root:  ${outroot}"
    echo "write MANIFEST: ${manifest}"
    for day in "${days[@]}"; do
        read -r tele mjd <<< "$day"
        echo "AR_SLURM=true AR_ALMANAC_SRC=${AR_ALMANAC_STAGE}/allobs_${tele}_${mjd}.h5 \\"
        echo "  AR_CHECKPOINT_MODE=${AR_CHECKPOINT_MODE} AR_JULIA_VERSION=${AR_JULIA_VERSION} AR_EXP_CLASS_MODEL= \\"
        echo "  ${harness_dir}/run_testday.sh ${base_dir} ${tele} ${mjd} ${outroot}/${tele}_${mjd}/"
    done
    echo "=== end dry run (nothing executed) ==="
    exit 0
fi

if [ -z "${SLURM_JOB_ID:-}" ]; then
    echo "Not inside a Slurm allocation. Submit with: sbatchAKS test/regression/submit_goldens.sh" >&2
    echo "(or pass --dry-run to preview; workstation runs use run_testday.sh directly)" >&2
    exit 2
fi

# mirror run_bulk.sh: expose the whole allocation to every julia stage
SLURM_NTASKS=$(($SLURM_CPUS_ON_NODE * $SLURM_NNODES))
export SLURM_NTASKS
env | grep SLURM | while read -r line; do echo "$line"; done

juliaup add "$AR_JULIA_VERSION" 2>/dev/null || true

mkdir -p "$outroot"
{
    echo "# Golden baselines MANIFEST"
    echo
    echo "- repo: ApogeeReduction.jl"
    echo "- branch: ${branch}"
    echo "- sha: ${sha}"
    echo "- generated: $(date -Is) by ${USER} (slurm job ${SLURM_JOB_ID}, nodes ${SLURM_JOB_NODELIST})"
    echo "- julia: ${AR_JULIA_VERSION}; checkpoint_mode: ${AR_CHECKPOINT_MODE}"
    echo "- exp_class_model: (none — production run_all.sh behavior)"
    echo "- almanac inputs: ${AR_ALMANAC_STAGE} (extracted from the 2026_05_01 bulk almanac)"
    echo "- raw: cca mirror /mnt/ceph/users/sdssv/raw/APOGEE"
    echo "- cals: darks/flats $(printenv AR_CALDIR_DARKS || echo '/mnt/ceph/users/sdssv/work/asaydjari/2025_07_31/outdir_ref/ (default)'); gain/read $(printenv AR_GAIN_READ_CAL_DIR || echo '/mnt/ceph/users/sdssv/work/asaydjari/2025_07_31/pass_clean/ (default)')"
    echo
    echo "## Per-day status"
    echo
    echo "| tele | mjd | exit code | wall time |"
    echo "|---|---|---|---|"
} > "$manifest"

# Progress is fully observable from the filesystem (no Slurm queries needed):
#   stdout sentinels        "=== GOLDEN BEGIN/END <tele> <mjd> ... ===" in the
#                           sbatch log (slurm_logs/ar_goldens_<jobid>.out)
#   <tele>_<mjd>/STATUS     written when that day finishes (exit, elapsed, times)
#   $outroot/RUNSTATE       rewritten at every day boundary (what is running,
#                           how many done)
#   $outroot/ALL_DONE       written only when the whole job finishes
ndays=${#days[@]}
ndone=0
overall=0
rm -f "${outroot}/ALL_DONE"
update_runstate() {
    # $1 = free-text current-state line
    {
        echo "job: ${SLURM_JOB_ID} on ${SLURM_JOB_NODELIST}"
        echo "updated: $(date -Is)"
        echo "state: $1"
        echo "done: ${ndone}/${ndays} days"
    } > "${outroot}/RUNSTATE"
}

for day in "${days[@]}"; do
    read -r tele mjd <<< "$day"
    outdir=${outroot}/${tele}_${mjd}/
    tstart=$(date -Is)
    echo
    echo "=== GOLDEN BEGIN ${tele} ${mjd} ${tstart} ==="
    echo "================ GOLDEN RUN ${tele} ${mjd} -> ${outdir} ================"
    update_runstate "running ${tele} ${mjd} (started ${tstart})"
    t0=$SECONDS
    set +e
    AR_SLURM=true \
        AR_ALMANAC_SRC=${AR_ALMANAC_STAGE}/allobs_${tele}_${mjd}.h5 \
        "$harness_dir/run_testday.sh" "$base_dir" "$tele" "$mjd" "$outdir"
    rc=$?
    set -e
    dt=$((SECONDS - t0))
    wall=$(printf '%dh:%dm:%ds' $((dt/3600)) $((dt%3600/60)) $((dt%60)))
    tend=$(date -Is)
    ndone=$((ndone + 1))
    echo "=== GOLDEN END ${tele} ${mjd} exit=${rc} elapsed=${wall} ${tend} ==="
    mkdir -p "$outdir"
    {
        echo "tele: ${tele}"
        echo "mjd: ${mjd}"
        echo "exit: ${rc}"
        echo "elapsed: ${wall}"
        echo "started: ${tstart}"
        echo "finished: ${tend}"
    } > "${outdir}STATUS"
    update_runstate "finished ${tele} ${mjd} exit=${rc}; between days"
    echo "| ${tele} | ${mjd} | ${rc} | ${wall} |" >> "$manifest"
    if [ $rc -ne 0 ]; then
        overall=1
        echo "GOLDEN RUN FAILED (rc=${rc}) for ${tele} ${mjd} — continuing with remaining days"
    fi
done

{
    echo
    if [ $overall -eq 0 ]; then
        echo "All golden runs completed with exit 0."
    else
        echo "AT LEAST ONE GOLDEN RUN FAILED — see the table above and the per-day logs (<tele>_<mjd>/logs/)."
    fi
} >> "$manifest"

{
    echo "finished: $(date -Is)"
    echo "overall_exit: ${overall}"
    echo "days: ${ndone}/${ndays}"
    echo "manifest: ${manifest}"
} > "${outroot}/ALL_DONE"
update_runstate "ALL DONE (overall_exit=${overall})"

echo
echo "MANIFEST written to $manifest"
exit $overall
