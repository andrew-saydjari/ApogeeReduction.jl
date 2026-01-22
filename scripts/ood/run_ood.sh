#!/bin/bash
#SBATCH --partition=cca,gen,genx
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --mem=0 #requesting all of the memory on the node

#SBATCH --time=8:00:00
#SBATCH --job-name=ar_ood
#SBATCH --output=slurm_logs/%x_%j.out
# ------------------------------------------------------------------------------
set -e # exit immediately if any of the steps returns a non-zero exit code
set -o pipefail

if [ -n "$SLURM_JOB_NODELIST" ]; then
    hostname
    echo $SLURM_JOB_NODELIST
else
    hostname
fi
echo "running from $(pwd)"

if [ -n "${SLURM_JOB_ID:-}" ] ; then
    script_path=$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')
else
    script_path=$(realpath "$0")
fi
base_dir="$(dirname "$(dirname "$(dirname "$script_path")")")"
echo "base_dir: $base_dir"

julia_version="1.11.6"
juliaup add $julia_version

# hardcode the mjd and expid for now
outdir=$1
mjd_start=$2
mjd_end=$3

ood_dir="ood/"
runname="allobs_${mjd_start}_${mjd_end}"
almanac_file=${outdir}almanac/${runname}.h5
runlist=${outdir}almanac/runlist_${runname}.h5
tele_list=("apo" "lco")

# nice function to print the elapsed time at each step
print_elapsed_time() {
    formatted_time=$(printf '%dd %dh:%dm:%ds\n' $(($SECONDS/86400)) $(($SECONDS%86400/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60)))
    echo
    echo "Elapsed time: $formatted_time"
    echo
    echo "--------- $1 ---------"
    echo
}

# make the ood exposure lists
print_elapsed_time "Making OOD Exposure Lists"
for tele in ${tele_list[@]}
do
    julia +$julia_version --project=$base_dir $base_dir/scripts/ood/gen_train_cov.jl --tele $tele --ood_dir $ood_dir --outdir $outdir --runlist $runlist --runname $runname
done


# Clean up logs and Report Timing
print_elapsed_time "Job completed"
