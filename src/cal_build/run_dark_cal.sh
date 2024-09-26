#!/bin/bash
## This is a local test script, but would be replaced by a SLURM script for a cluster

# Exit the script immedatiely if any command returns a non-zero status.
#set -Eeuo pipefail
set -eo pipefail

tele=$1
mjd_start=$2
mjd_end=$3

echo "Running dark calibration for $tele from $mjd_start to $mjd_end"

# load all of the modules to talk to the database (need to be on Utah)
# should turn this off as an option for users once the MJD summaries are generated
module load sdssdb/main almanac sdsstools postgresql ffmpeg

export SLURM_TASKS_PER_NODE="64(x2)"
export SLURM_CLUSTERS="notchpeak"
export SLURM_CLUSTER="notchpeak.peaks"
export SLURM_NNODES="2"
export SLURM_NPROCS="128"

# Not in original
export SLURM_JOB_NUM_NODES="2"
export SLURM_JOB_QOS="sdss-np"
export SLURM_JOB_ACCOUNT="sdss-np"
export SLURM_STEP_NUM_TASKS=128
export SLURM_STEP_NUM_NODES=2
export SLURM_NTASKS=32
export SLURM_JOB_NODELIST=$SLURM_NODELIST

runname="dark_cal_${mjd_start}_${mjd_end}"
outdir="../outdir/"
doutdir=$outdir
almanac_file=${outdir}almanac/${runname}.h5
runlist=${outdir}almanac/runlist_${runname}.jld2

# set up the output directory (if does not exist)
mkdir -p ${outdir}almanac

# get the data summary file for the MJD
almanac -v -p 12 --mjd-start $mjd_start --mjd-end $mjd_end --${tele} --output $almanac_file

# get the runlist file (julia projects seem to refer to where your cmd prompt is when you call the shell. Here I imaging sitting at ApogeeReduction.jl level)
julia +1.10.0 --project="./" src/cal_build/make_runlist_darks.jl --tele $tele --almanac_file $almanac_file --output $runlist

# run the reduction pipeline (all cals like dark sub/flats that would be a problem, should be post 3D->2D extraction)
echo " before A prep $tele with $runlist and $doutdir and $runname"
julia +1.10.0 --project="./" pipeline.jl --tele $tele --runlist $runlist --outdir $doutdir --runname $runname --chip "a"
echo " before B prep $tele with $runlist and $doutdir and $runname"
julia +1.10.0 --project="./" pipeline.jl --tele $tele --runlist $runlist --outdir $doutdir --runname $runname --chip "b"
echo " before C prep $tele with $runlist and $doutdir and $runname"
julia +1.10.0 --project="./" pipeline.jl --tele $tele --runlist $runlist --outdir $doutdir --runname $runname --chip "c"
echo " done prep $tele with $runlist and $doutdir and $runname"

# make the stacked darks
mkdir -p ${outdir}darks
echo " before stack A prep $tele with $mjd_start and $mjd_end and $doutdir"
julia +1.10.0 --project="./" src/cal_build/make_stack_darks.jl --mjd-start $mjd_start --mjd-end $mjd_end --tele $tele --chip "a" --dark_dir $doutdir
echo " before stack B prep $tele with $mjd_start and $mjd_end and $doutdir"
julia +1.10.0 --project="./" src/cal_build/make_stack_darks.jl --mjd-start $mjd_start --mjd-end $mjd_end --tele $tele --chip "b" --dark_dir $doutdir
echo " before stack C prep $tele with $mjd_start and $mjd_end and $doutdir"
julia +1.10.0 --project="./" src/cal_build/make_stack_darks.jl --mjd-start $mjd_start --mjd-end $mjd_end --tele $tele --chip "c" --dark_dir $doutdir
echo " done stack things prep $tele with $mjd_start and $mjd_end and $doutdir"

# Clean up logs and Report Timing
formatted_time=$(printf '%dd %dh:%dm:%ds\n' $(($SECONDS/86400)) $(($SECONDS%86400/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60)))
echo "Job completed in $formatted_time"