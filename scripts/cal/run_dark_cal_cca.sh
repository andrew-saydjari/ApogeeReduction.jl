#!/bin/bash
#SBATCH --partition=cca,gen,genx
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --mem=0 #requesting all of the memory on the node

#SBATCH --time=8:00:00
#SBATCH --job-name=ar_dark_cal
#SBATCH --output=../slurm_logs/%x_%j.out
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

julia_version="1.11.6"
juliaup add $julia_version

# hardcode the mjd and expid for now
tele=$1
mjd_start=$2
mjd_end=$3
dropfirstn=${4:-0}

runname="dark_cal_${tele}_${mjd_start}_${mjd_end}"
outdir="../outdir/"
doutdir=$outdir
almanac_file=${outdir}almanac/${runname}.h5
runlist=${outdir}almanac/runlist_${runname}.h5
gain_read_cal_dir="/mnt/ceph/users/asaydjari/working/2025_07_31/pass_clean/"

# get the runlist file (julia projects seem to refer to where your cmd prompt is when you call the shell. Here I imaging sitting at ApogeeReduction.jl level)
julia +$julia_version --project="./" scripts/cal/make_runlist_darks.jl --tele $tele --almanac_file $almanac_file --output $runlist

# run the reduction pipeline (all cals like dark sub/flats that would be a problem, should be post 3D->2D extraction)
## 100 read darks require fewer workers per node (a problem for the ultra darks)
julia +$julia_version --project="./" pipeline.jl --workers_per_node 28 --tele $tele --runlist $runlist --outdir $doutdir --runname $runname --chips "RGB" --doCal2d false --cluster cca --gain_read_cal_dir $gain_read_cal_dir

# make the stacked darks
mkdir -p ${outdir}darks
julia +$julia_version --project="./" scripts/cal/make_stack_darks.jl --mjd-start $mjd_start --mjd-end $mjd_end --tele $tele --chip "R" --dropfirstn $dropfirstn --dark_dir ${doutdir} --runlist $runlist
julia +$julia_version --project="./" scripts/cal/make_stack_darks.jl --mjd-start $mjd_start --mjd-end $mjd_end --tele $tele --chip "G" --dropfirstn $dropfirstn --dark_dir ${doutdir} --runlist $runlist
julia +$julia_version --project="./" scripts/cal/make_stack_darks.jl --mjd-start $mjd_start --mjd-end $mjd_end --tele $tele --chip "B" --dropfirstn $dropfirstn --dark_dir ${doutdir} --runlist $runlist

# Clean up logs and Report Timing
formatted_time=$(printf '%dd %dh:%dm:%ds\n' $(($SECONDS/86400)) $(($SECONDS%86400/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60)))
echo "Job completed in $formatted_time"
