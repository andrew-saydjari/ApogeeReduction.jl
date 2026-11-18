#!/bin/bash
#SBATCH --partition=cca
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem=128G
#SBATCH --constraint="[genoa|icelake|rome]"

#SBATCH --time=4:00:00
#SBATCH --job-name=ar_chip_transform
#SBATCH --output=slurm_logs/%x_%j.out
# ------------------------------------------------------------------------------
# Rebuild the comb-measured chip-transform calibration from a whole reduction
# (see README.md in this folder). Runs the same way locally and under Slurm:
#
#   local : bash scripts/cal/chip_transform/run_chip_transform.sh <root> <workdir> <version> [workers] [extra args]
#   slurm : sbatch scripts/cal/chip_transform/run_chip_transform.sh <root> <workdir> <version> [ignored] [extra args]
#
# Under Slurm every task becomes a single-threaded Distributed worker (SlurmManager);
# locally [workers] single-threaded workers are started (default 6).
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
workers=${4:-6}
shift 4 || shift $#

export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 JULIA_NUM_THREADS=1
# --cover_from: the first block is extended back to the start of the release so the plate era
# (no comb) uses the earliest FPS block -- validated by the stage-4 plate-era check (README).
julia +$julia_version --project=$base_dir $base_dir/scripts/cal/chip_transform/build_chip_transform.jl \
    --stage all --root "$root" --workdir "$workdir" --version "$version" --workers "$workers" \
    --tele apo,lco "$@"
