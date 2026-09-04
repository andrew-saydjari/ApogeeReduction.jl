#!/bin/bash
# submit_goldens.sh — generate the T1 golden-diff baselines as ONE true bulk
# run over all test days, mirroring scripts/bulk/run_bulk.sh exactly (single
# combined runlist -> per-tele pipeline.jl -> flat-type cal loops -> per-tele
# pipeline_2d_1d.jl), NOT per-day harness invocations and NOT a job array.
# Outputs land in one outdir whose apred/<mjd>/ subdirs are per-day as usual.
#
# Submit (the user does this, from the AR checkout root so slurm_logs/ lands
# beside the repo):
#   sbatchAKS test/regression/submit_goldens.sh
#
# Validate without Slurm:
#   bash -n test/regression/submit_goldens.sh
#   test/regression/submit_goldens.sh --dry-run
#
# Local (workstation) execution — used for harness validation runs only:
#   AR_SLURM=false AR_WORKERS=10 AR_DAYS="apo 59429;apo 58011" \
#       AR_OUTROOT=<outdir> test/regression/submit_goldens.sh
#
# Prerequisite: a raw/-layout almanac file covering all test days
# (AR_ALMANAC_SRC; default the 2026_05_01 bulk allobs_57600_61160.h5),
# consumed directly — the job has zero DB/tunnel dependency.
#
# Telescope/runlist interaction (mirrors run_bulk.sh): ONE combined runlist
# carries both telescopes; pipeline.jl / pipeline_2d_1d.jl /
# make_traces_from_flats.jl self-filter on their --tele argument, and
# make_relFlux.jl is called once per flat type over the combined flat runlist
# (it iterates the unique teles itself and writes tele/mjd-keyed entries into
# a single runname-keyed valid-flats file — per-tele calls would clobber it).
# Because the test days are (tele, mjd) PAIRS (lco only at 60000), the
# combined runlists are built per tele (each with its own --mjds) and merged
# by concat_runlists.jl.
#
# Progress is fully observable from the filesystem (no Slurm queries needed):
#   stdout sentinels    "=== STAGE BEGIN/END/FAILED <name> ... ===" in the
#                       sbatch log (slurm_logs/ar_goldens_<jobid>.out)
#   $outroot/RUNSTATE   rewritten at every stage boundary
#   $outroot/status/STATUS_<tele>_<mjd>
#                       per-day bookkeeping written by the post-run step:
#                       product counts per class + missing-file check
#   $outroot/ALL_DONE   written only when the whole job finishes
# NOTE: bulk mode is one chain — a stage failure aborts the remaining stages
# (no per-day continue-on-failure; that was a property of the old sequential
# mode). The last "=== STAGE BEGIN" without a matching END/FAILED marks where
# it died.
#
# Decisions baked in (2026-08-31): goldens are generated WITHOUT the
# exposure-classifier model (matches run_all.sh production behavior; the
# classifier enters later via R1/O5) — AR_EXP_CLASS_MODEL stays empty.

# ------------------------------------------------------------------------------
#SBATCH --partition=cca
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
AR_ALMANAC_SRC=${AR_ALMANAC_SRC:-"/mnt/ceph/users/sdssv/work/asaydjari/2026_05_01/outdir/almanac/allobs_57600_61160.h5"}
AR_JULIA_VERSION=${AR_JULIA_VERSION:-"1.11.0"}
AR_CHECKPOINT_MODE=${AR_CHECKPOINT_MODE:-"commit_exists"}
AR_RAW_CLUSTER=${AR_RAW_CLUSTER:-"cca"}
AR_CALDIR_DARKS=${AR_CALDIR_DARKS:-"/mnt/ceph/users/sdssv/work/asaydjari/2025_07_31/outdir_ref/"}
AR_CALDIR_FLATS=${AR_CALDIR_FLATS:-"/mnt/ceph/users/sdssv/work/asaydjari/2025_07_31/outdir_ref/"}
AR_GAIN_READ_CAL_DIR=${AR_GAIN_READ_CAL_DIR:-"/mnt/ceph/users/sdssv/work/asaydjari/2025_07_31/pass_clean/"}
AR_WORKERS=${AR_WORKERS:-24}    # local (AR_SLURM=false) mode only
AR_EXP_CLASS_MODEL=""           # decision: goldens without the classifier

# The test days (REFACTOR_PLAN v1 §4.2): 6 (tele, mjd) pairs covering the 5
# fixed test days + both telescopes at 60000. The A3 gapped-almanac day joins
# later (SUBMIT.md). AR_DAYS overrides for subsets, semicolon-separated pairs:
#   AR_DAYS="apo 60000;apo 58011"
if [ -n "${AR_DAYS:-}" ]; then
    IFS=';' read -r -a days <<< "$AR_DAYS"
else
    days=(
        "apo 60000"
        "lco 60000"
        "apo 58588"
        "apo 59337"
        "apo 58011"
        "apo 59429"
    )
fi

# per-tele mjd lists (the test days are tele-specific pairs)
declare -A tele_mjds
teles=()
for day in "${days[@]}"; do
    read -r tele mjd <<< "$day"
    if [ -z "${tele_mjds[$tele]:-}" ]; then
        teles+=("$tele")
        tele_mjds[$tele]="$mjd"
    else
        tele_mjds[$tele]="${tele_mjds[$tele]},$mjd"
    fi
done
mjd_min=$(printf '%s\n' "${days[@]}" | awk '{print $2}' | sort -n | head -1)
mjd_max=$(printf '%s\n' "${days[@]}" | awk '{print $2}' | sort -n | tail -1)

# ---- git identity: sha of the code the run will execute = golden identity ---
if ! git -C "$base_dir" diff --quiet || ! git -C "$base_dir" diff --cached --quiet; then
    echo "WARNING: $base_dir has uncommitted changes — goldens should come from a clean, pushed commit." >&2
    if ! $DRY_RUN && [ "${AR_ALLOW_DIRTY:-0}" != "1" ]; then
        echo "Refusing to generate goldens from a dirty tree (AR_ALLOW_DIRTY=1 overrides — validation runs only)." >&2
        exit 3
    fi
fi
sha=$(git -C "$base_dir" rev-parse HEAD)
branch=$(git -C "$base_dir" rev-parse --abbrev-ref HEAD)
shortsha=${sha:0:7}
# AR_OUTROOT overrides the sha-derived golden root (redos into an existing
# tree when only harness files changed between shas, or validation runs).
outroot=${AR_OUTROOT:-${AR_GOLDEN_ROOT}/ApogeeReduction.jl@${shortsha}}
case "$outroot" in */) : ;; *) outroot="${outroot}/" ;; esac
manifest=${outroot}MANIFEST.md

runname="allobs_${mjd_min}_${mjd_max}"
almanac_file=${outroot}almanac/${runname}.h5
runlist=${outroot}almanac/runlist_${runname}.h5
flat_types=("quartz" "dome")

# ---- validation common to dry-run and real run ------------------------------
if [ ! -f "$AR_ALMANAC_SRC" ]; then
    echo "MISSING almanac file: $AR_ALMANAC_SRC" >&2
    echo "Point AR_ALMANAC_SRC at a raw/-layout almanac file covering the test days." >&2
    exit 4
fi

if $DRY_RUN; then
    echo "=== DRY RUN — bulk-mode stage sequence the job will execute ==="
    echo "git identity: ${branch} @ ${sha}"
    echo "golden root:  ${outroot}"
    echo "runname:      ${runname} (almanac symlink -> ${AR_ALMANAC_SRC})"
    echo "teles/days:"
    for tele in "${teles[@]}"; do
        echo "  ${tele}: --mjds ${tele_mjds[$tele]}"
    done
    echo "stages:"
    echo "  1. link almanac -> ${almanac_file}"
    echo "  2. runlist: make_runlist_all.jl per tele (--mjds as above) + concat_runlists.jl -> ${runlist}"
    for tele in "${teles[@]}"; do
        echo "  3. pipeline.jl --tele ${tele} --runlist <combined> --runname ${runname} (3D->2D/2Dcal)"
    done
    for flat_type in "${flat_types[@]}"; do
        echo "  4. [${flat_type}] make_runlist_fiber_flats.jl per tele + concat; per tele: make_traces_from_flats.jl + pipeline_2d_1d.jl (relFlux/waveSoln false); then make_relFlux.jl once (combined)"
    done
    for tele in "${teles[@]}"; do
        echo "  5. pipeline_2d_1d.jl --tele ${tele} --runlist <combined> (final 2D->1D)"
    done
    echo "  6. warnings census; per-day bookkeeping -> ${outroot}status/STATUS_<tele>_<mjd>, MANIFEST, ALL_DONE"
    echo "=== end dry run (nothing executed) ==="
    exit 0
fi

# ---- execution mode (mirrors run_testday.sh) --------------------------------
#   AR_SLURM=true  — inside sbatch: run_bulk.sh pattern, SLURM_NTASKS =
#                    CPUS_ON_NODE * NNODES, every julia stage spans the
#                    allocation via SlurmClusterManager.
#   AR_SLURM=false — workstation: scrub the Slurm env, pass --workers_per_node.
#   auto (default) — true inside an allocation; outside, refuse unless
#                    AR_SLURM=false was set explicitly (accidental login-node
#                    bulk runs are expensive).
AR_SLURM=${AR_SLURM:-auto}
if [ "$AR_SLURM" = "auto" ]; then
    if [ -n "${SLURM_JOB_ID:-}" ]; then
        AR_SLURM=true
    else
        echo "Not inside a Slurm allocation. Submit with: sbatchAKS test/regression/submit_goldens.sh" >&2
        echo "(--dry-run previews; AR_SLURM=false runs the chain locally, e.g. for harness validation)" >&2
        exit 2
    fi
fi
if [ "$AR_SLURM" = "true" ]; then
    SLURM_NTASKS=$(($SLURM_CPUS_ON_NODE * $SLURM_NNODES))
    export SLURM_NTASKS
    env | grep SLURM | while read -r line; do echo "$line"; done
    workers_args=()
else
    unset SLURM_NTASKS SLURM_JOB_ID SLURM_NNODES SLURM_CPUS_ON_NODE SLURM_JOB_NODELIST
    workers_args=(--workers_per_node "$AR_WORKERS")
fi

juliaup add "$AR_JULIA_VERSION" 2>/dev/null || true

mkdir -p "${outroot}almanac" "${outroot}status" "${outroot}logs"
logfile=${outroot}logs/submit_goldens_${runname}.log
exec > >(tee -a "$logfile") 2>&1

echo "submit_goldens.sh (bulk mode) starting $(date -Is) on $(hostname)"
echo "  git: ${branch} @ ${sha}"
echo "  outroot=$outroot runname=$runname"
echo "  almanac_src=$AR_ALMANAC_SRC"
echo "  days: ${days[*]}"
echo "  mode=$([ "$AR_SLURM" = "true" ] && echo "slurm SLURM_NTASKS=$SLURM_NTASKS" || echo "local workers=$AR_WORKERS")"

# ---- stage sentinels + RUNSTATE ---------------------------------------------
STAGE_NAME="(startup)"
STAGE_T0=$SECONDS
update_runstate() {
    {
        echo "job: ${SLURM_JOB_ID:-local} on ${SLURM_JOB_NODELIST:-$(hostname)}"
        echo "updated: $(date -Is)"
        echo "state: $1"
    } > "${outroot}RUNSTATE"
}
stage_begin() {
    STAGE_NAME="$1"
    STAGE_T0=$SECONDS
    echo
    echo "=== STAGE BEGIN ${STAGE_NAME} $(date -Is) ==="
    update_runstate "running stage: ${STAGE_NAME}"
}
stage_end() {
    local dt=$((SECONDS - STAGE_T0))
    echo "=== STAGE END ${STAGE_NAME} elapsed=$(printf '%dh:%dm:%ds' $((dt / 3600)) $((dt % 3600 / 60)) $((dt % 60))) $(date -Is) ==="
}
trap 'echo "=== STAGE FAILED ${STAGE_NAME} $(date -Is) ==="; update_runstate "FAILED in stage: ${STAGE_NAME}"' ERR

rm -f "${outroot}ALL_DONE"

# ---- stage: almanac link ----------------------------------------------------
stage_begin "link almanac"
ln -sfn "$(realpath "$AR_ALMANAC_SRC")" "$almanac_file"
stage_end

# ---- stage: combined runlist (per tele --mjds, then concat) -----------------
stage_begin "build combined runlist"
rm -f "$runlist"
tele_runlists=()
for tele in "${teles[@]}"; do
    tele_runlist=${outroot}almanac/runlist_${runname}_${tele}.h5
    rm -f "$tele_runlist"
    julia +"$AR_JULIA_VERSION" --project="$base_dir" "$base_dir/scripts/bulk/make_runlist_all.jl" \
        --tele "$tele" --almanac_file "$almanac_file" --output "$tele_runlist" \
        --mjds "${tele_mjds[$tele]}"
    tele_runlists+=("$tele_runlist")
done
julia +"$AR_JULIA_VERSION" --project="$base_dir" "$harness_dir/concat_runlists.jl" \
    "$runlist" "${tele_runlists[@]}"
stage_end

# ---- stage: 3D -> 2D / 2Dcal per tele ---------------------------------------
for tele in "${teles[@]}"; do
    stage_begin "3D->2D/2Dcal pipeline ${tele}"
    julia +"$AR_JULIA_VERSION" --project="$base_dir" "$base_dir/pipeline.jl" \
        --tele "$tele" --runlist "$runlist" --outdir "$outroot" --runname "$runname" \
        --chips "RGB" --caldir_darks "$AR_CALDIR_DARKS" --caldir_flats "$AR_CALDIR_FLATS" \
        --cluster "$AR_RAW_CLUSTER" --gain_read_cal_dir "$AR_GAIN_READ_CAL_DIR" \
        --checkpoint_mode "$AR_CHECKPOINT_MODE" "${workers_args[@]}" \
        ${AR_EXP_CLASS_MODEL:+--exp_class_model "$AR_EXP_CLASS_MODEL"}
    stage_end
done

# ---- stages: traces + relFluxing per flat type ------------------------------
for flat_type in "${flat_types[@]}"; do
    flatrunlist=${outroot}almanac/runlist_${flat_type}_${runname}.h5
    stage_begin "runlist ${flat_type} flats"
    rm -f "$flatrunlist"
    tele_flatrunlists=()
    for tele in "${teles[@]}"; do
        tele_flatrunlist=${outroot}almanac/runlist_${flat_type}_${runname}_${tele}.h5
        rm -f "$tele_flatrunlist"
        julia +"$AR_JULIA_VERSION" --project="$base_dir" "$base_dir/scripts/cal/make_runlist_fiber_flats.jl" \
            --almanac_file "$almanac_file" --tele "$tele" --output "$tele_flatrunlist" \
            --flat_type "$flat_type" --mjds "${tele_mjds[$tele]}"
        tele_flatrunlists+=("$tele_flatrunlist")
    done
    set +e
    julia +"$AR_JULIA_VERSION" --project="$base_dir" "$harness_dir/concat_runlists.jl" \
        "$flatrunlist" "${tele_flatrunlists[@]}"
    concat_rc=$?
    set -e
    stage_end
    if [ $concat_rc -eq 16 ]; then
        echo "No ${flat_type} flats on any test day — empty combined runlist written; downstream stages handle it."
    elif [ $concat_rc -ne 0 ]; then
        false   # trip the ERR trap + set -e with the stage recorded
    fi

    for tele in "${teles[@]}"; do
        stage_begin "traces from ${flat_type} flats ${tele}"
        mkdir -p "${outroot}${flat_type}_flats"
        julia +"$AR_JULIA_VERSION" --project="$base_dir" "$base_dir/scripts/cal/make_traces_from_flats.jl" \
            --tele "$tele" --trace_dir "$outroot" --runlist "$flatrunlist" \
            --flat_type "$flat_type" --slack_quiet true --checkpoint_mode "$AR_CHECKPOINT_MODE"
        stage_end

        stage_begin "2D->1D (no relFlux) ${flat_type} flats ${tele}"
        mkdir -p "${outroot}apredrelflux"
        julia +"$AR_JULIA_VERSION" --project="$base_dir" "$base_dir/pipeline_2d_1d.jl" \
            --tele "$tele" --runlist "$flatrunlist" --outdir "$outroot" --runname "$runname" \
            --relFlux false --waveSoln false --checkpoint_mode "$AR_CHECKPOINT_MODE" \
            "${workers_args[@]}"
        stage_end
    done

    stage_begin "relFlux ${flat_type} flats (all teles)"
    # once per flat type over the COMBINED runlist, like run_bulk.sh —
    # make_relFlux iterates the runlist's unique teles and writes the single
    # runname-keyed valid-flats file (per-tele calls would clobber it)
    julia +"$AR_JULIA_VERSION" --project="$base_dir" "$base_dir/scripts/cal/make_relFlux.jl" \
        --trace_dir "$outroot" --runlist "$flatrunlist" --runname "$runname"
    stage_end
done

# ---- stage: final 2D -> 1D per tele -----------------------------------------
for tele in "${teles[@]}"; do
    stage_begin "final 2D->1D pipeline ${tele}"
    julia +"$AR_JULIA_VERSION" --project="$base_dir" "$base_dir/pipeline_2d_1d.jl" \
        --tele "$tele" --runlist "$runlist" --outdir "$outroot" --runname "$runname" \
        --checkpoint_mode "$AR_CHECKPOINT_MODE" "${workers_args[@]}"
    stage_end
done

# ---- stage: warnings census (regression metric, cf. REFACTOR_PLAN v1 §0) ----
stage_begin "warnings census"
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
stage_end

# ---- stage: per-day bookkeeping ---------------------------------------------
stage_begin "per-day bookkeeping"
# per-(tele,mjd) exposure counts from the combined runlist -> expected products
nexp_table=$(julia +"$AR_JULIA_VERSION" --project="$base_dir" -e '
using JLD2
t = load(ARGS[1], "tele"); m = load(ARGS[1], "mjd")
for (tt, mm) in sort(unique(collect(zip(t, m))))
    println(tt, " ", mm, " ", count((t .== tt) .& (m .== mm)))
end' "$runlist")
echo "runlist exposure counts (tele mjd nexp):"
echo "$nexp_table"

# MANIFEST (append-aware: redo runs add a new section)
if [ -f "$manifest" ]; then
    { echo; echo "---"; echo; } >> "$manifest"
    echo "# Redo/partial run ($(date -Is))" >> "$manifest"
else
    echo "# Golden baselines MANIFEST" >> "$manifest"
fi
{
    echo
    echo "- repo: ApogeeReduction.jl"
    echo "- branch: ${branch}"
    echo "- sha: ${sha}"
    echo "- generated: $(date -Is) by ${USER} ($([ "$AR_SLURM" = "true" ] && echo "slurm job ${SLURM_JOB_ID}, nodes ${SLURM_JOB_NODELIST}" || echo "local on $(hostname), ${AR_WORKERS} workers"))"
    echo "- mode: ONE bulk run over all test days (run_bulk.sh sequence), runname ${runname}"
    echo "- julia: ${AR_JULIA_VERSION}; checkpoint_mode: ${AR_CHECKPOINT_MODE}"
    echo "- exp_class_model: (none — production run_all.sh behavior)"
    echo "- almanac: ${AR_ALMANAC_SRC} (bulk raw/-layout file, consumed directly; days selected via the runlist makers' --mjds per tele)"
    echo "- raw: cluster ${AR_RAW_CLUSTER}"
    echo "- cals: darks ${AR_CALDIR_DARKS}; flats ${AR_CALDIR_FLATS}; gain/read ${AR_GAIN_READ_CAL_DIR}"
    echo
    echo "## Per-day products (from filesystem, post-run)"
    echo
    echo "| tele | mjd | nexp | ar2D | ar2Dcal | ar1Dcal | ar1Dunical | verdict |"
    echo "|---|---|---|---|---|---|---|---|"
} >> "$manifest"

count_products() {
    # $1 = class prefix, $2 = tele, $3 = mjd. Glob count, NOT `ls | wc -l`:
    # under this script's `set -e -o pipefail`, a day with no products dir
    # (e.g. a dark-only night with zero runlist entries — apo 59136 killed
    # testbed job 6980442 here) made the failed `ls` abort the whole
    # bookkeeping loop.
    local g=("${outroot}apred/$3/$1_$2_$3_"*)
    if [ -e "${g[0]}" ]; then echo "${#g[@]}"; else echo 0; fi
}
overall=0
for day in "${days[@]}"; do
    read -r tele mjd <<< "$day"
    nexp=$(echo "$nexp_table" | awk -v t="$tele" -v m="$mjd" '$1 == t && $2 == m {print $3}')
    nexp=${nexp:-0}
    n2d=$(count_products ar2D "$tele" "$mjd")
    n2dcal=$(count_products ar2Dcal "$tele" "$mjd")
    n1dcal=$(count_products ar1Dcal "$tele" "$mjd")
    n1duni=$(count_products ar1Dunical "$tele" "$mjd")
    # missing-file check: every runlist exposure must have 3 chips of ar2D;
    # a day with exposures must have produced 1D products
    verdict=ok
    if [ "$nexp" -eq 0 ]; then
        verdict="NO-EXPOSURES"
    elif [ "$n2d" -ne $((3 * nexp)) ]; then
        verdict="MISSING-ar2D($n2d/$((3 * nexp)))"
    elif [ "$n1dcal" -eq 0 ] || [ "$n1duni" -eq 0 ]; then
        verdict="MISSING-1D"
    fi
    [ "$verdict" = "ok" ] || overall=1
    {
        echo "tele: ${tele}"
        echo "mjd: ${mjd}"
        echo "nexp_runlist: ${nexp}"
        echo "ar2D: ${n2d} (expected $((3 * nexp)))"
        echo "ar2Dcal: ${n2dcal}"
        echo "ar1Dcal: ${n1dcal}"
        echo "ar1Dunical: ${n1duni}"
        echo "verdict: ${verdict}"
        echo "finished: $(date -Is)"
        echo "diff base: <golden>/${tele}_${mjd}/apred/${mjd} vs ${outroot}apred/${mjd}"
    } > "${outroot}status/STATUS_${tele}_${mjd}"
    echo "| ${tele} | ${mjd} | ${nexp} | ${n2d} | ${n2dcal} | ${n1dcal} | ${n1duni} | ${verdict} |" >> "$manifest"
done

{
    echo
    if [ $overall -eq 0 ]; then
        echo "All per-day product checks passed."
    else
        echo "AT LEAST ONE PER-DAY PRODUCT CHECK FAILED — see the table above and the stage log (${logfile})."
    fi
} >> "$manifest"

{
    echo "finished: $(date -Is)"
    echo "overall_exit: ${overall}"
    echo "manifest: ${manifest}"
} > "${outroot}ALL_DONE"
update_runstate "ALL DONE (overall_exit=${overall})"
stage_end

echo
echo "MANIFEST written to $manifest"
exit $overall
