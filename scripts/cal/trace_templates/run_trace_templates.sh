#!/bin/bash
#SBATCH --partition=cca
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --constraint=genoa
#SBATCH --mem=0

#SBATCH --time=24:00:00
#SBATCH --job-name=ar_trace_templates
#SBATCH --output=slurm_logs/%x_%j.out
# ------------------------------------------------------------------------------
# Rebuild the per-epoch trace-template calibration from the 2D quartz flats of a reduction
# (README.md in this folder). Runs the same way locally and under Slurm:
#
#   local : bash  scripts/cal/trace_templates/run_trace_templates.sh <root> <workdir> <version> <epochs.csv> [parallel] [nmax]
#   slurm : sbatch scripts/cal/trace_templates/run_trace_templates.sh <root> <workdir> <version> <epochs.csv> [parallel] [nmax]
#
#   <root>       reduction output directory (apred/<mjd>/ar2Dcal_*_quartzflat.h5); read-only
#   <workdir>    where the unit builds (units/), the packed file <version>.h5 and the QA (qa/) go
#   <epochs.csv> CONFIRMED epoch table (tele,epoch,mjd_start,mjd_end); propose_trace_epochs.jl
#                writes a draft from the data, a person confirms it (README, "Epochs")
#   [parallel]   units built at once, one single-threaded process each (default: all units,
#                i.e. 42 for 14 epochs x 3 chips; set it to the cores you may use locally)
#   [nmax]       nights per unit (default 0 = every night of the epoch)
#
# Units already built (units/template_<tele>_<chip>_<epoch>.h5 present) are skipped, so the
# script is resumable. Cost of the first full build (pass-1 DR21, 42 units, NMAX = 0): see README.
set -e
set -o pipefail

if [ -n "${SLURM_JOB_ID:-}" ] ; then
    script_path=$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}' | awk '{print $1}')
else
    script_path=$(realpath "$0")
fi
base_dir="$(dirname "$(dirname "$(dirname "$(dirname "$script_path")")")")"
echo "base_dir: $base_dir  host: $(hostname)"

julia_version="1.11.6"
juliaup add $julia_version

root=$1
workdir=$2
version=$3
epochs=$4
parallel=${5:-0}
nmax=${6:-0}
[ -n "$root" ] && [ -n "$workdir" ] && [ -n "$version" ] && [ -f "$epochs" ] || {
    echo "usage: run_trace_templates.sh <root> <workdir> <version> <epochs.csv> [parallel] [nmax]"; exit 1; }

export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 JULIA_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p "$workdir/units" "$workdir/logs" "$workdir/qa"
cp "$epochs" "$workdir/trace_epochs_used.csv"
jl="julia +$julia_version --project=$base_dir"
here="$base_dir/scripts/cal/trace_templates"
$jl -e 'using ApogeeReduction'      # precompile once before the parallel builds

# 1. one unit per (epoch, chip), in parallel
jobs=()
while IFS=, read -r tele epoch m0 m1; do
    for chip in R G B; do
        [ -f "$workdir/units/template_${tele}_${chip}_${epoch}.h5" ] && continue
        jobs+=("$tele $chip $epoch $m0 $m1")
    done
done < <(tail -n +2 "$epochs" | tr -d '\r')
echo "units to build: ${#jobs[@]}"
[ "$parallel" -gt 0 ] || parallel=${#jobs[@]}
[ "$parallel" -gt 0 ] || parallel=1
failed=0
running=0
for j in "${jobs[@]}"; do
    set -- $j
    ( $jl "$here/build_trace_template.jl" --root "$root" --tele "$1" --chip "$2" --epoch "$3" \
        --mjd_start "$4" --mjd_end "$5" --nmax "$nmax" --outdir "$workdir/units" \
        > "$workdir/logs/build_${1}_${2}_${3}.log" 2>&1 || echo "FAILED $1 $2 $3" >> "$workdir/logs/build_${1}_${2}_${3}.log" ) &
    running=$((running + 1))
    if [ "$running" -ge "$parallel" ]; then wait -n || true; running=$((running - 1)); fi
done
wait
grep -l "^FAILED\|ERROR" "$workdir"/logs/build_*.log 2>/dev/null && failed=1 || true

# 2. pack (refuses an incomplete set unless the units are all there) and 3. check
$jl "$here/pack_trace_templates.jl" --epochs "$epochs" --units_dir "$workdir/units" \
    --out "$workdir/$version.h5" --version "$version" 2>&1 | tee "$workdir/logs/pack.log"
$jl "$here/check_trace_templates.jl" --cal "$workdir/$version.h5" --units_dir "$workdir/units" \
    --outdir "$workdir/qa" --require_complete true 2>&1 | tee "$workdir/logs/check.log"
echo "packed: $workdir/$version.h5 (install by copying into cal_ref/trace_templates/ and pointing --trace_template_cal at it)"
exit $failed
