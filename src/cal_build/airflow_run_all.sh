#!/bin/bash
# Exit the script immedatiely if any command returns a non-zero status.
set -eo pipefail

tele=$1
mjd_start=$2
mjd_end=$3

echo `printenv`

echo "AFTER:"

# TODO: Move this to airflow's `srun`
export SLURM_TASKS_PER_NODE="64(x2)"
export SLURM_CLUSTERS="notchpeak"
export SLURM_CLUSTER="notchpeak.peaks"
export SLURM_NNODES="2"
export SLURM_NPROCS="128"
export SLURM_JOB_NUM_NODES="2"
export SLURM_JOB_QOS="sdss-np"
export SLURM_JOB_ACCOUNT="sdss-np"
export SLURM_STEP_NUM_TASKS=128
export SLURM_STEP_NUM_NODES=2
export SLURM_NTASKS=32
export SLURM_JOB_NODELIST=$SLURM_NODELIST

echo `printenv`

runname="all_${mjd_start}_${mjd_end}"
outdir="../outdir/${mjd_start}_${mjd_end}/"
doutdir=$outdir
almanac_file=${outdir}almanac/${runname}.h5
runlist=${outdir}almanac/runlist_${runname}.jld2

module load sdssdb/main sdsstools postgresql ffmpeg

echo `module list`
julia +1.10.0 --project="./" src/cal_build/make_runlist.jl --tele $tele --almanac_file $almanac_file --output $runlist

# run the reduction pipeline (all cals like dark sub/flats that would be a problem, should be post 3D->2D extraction)
julia +1.10.0 --project="./" pipeline.jl --tele $tele --runlist $runlist --outdir $doutdir --runname $runname --chip "a"
julia +1.10.0 --project="./" pipeline.jl --tele $tele --runlist $runlist --outdir $doutdir --runname $runname --chip "b"
julia +1.10.0 --project="./" pipeline.jl --tele $tele --runlist $runlist --outdir $doutdir --runname $runname --chip "c"

# make the stacked darks
julia +1.10.0 --project="./" src/cal_build/make_stack_darks.jl --mjd-start $mjd_start --mjd-end $mjd_end --tele $tele --chip "a" --dark_dir $doutdir
julia +1.10.0 --project="./" src/cal_build/make_stack_darks.jl --mjd-start $mjd_start --mjd-end $mjd_end --tele $tele --chip "b" --dark_dir $doutdir
julia +1.10.0 --project="./" src/cal_build/make_stack_darks.jl --mjd-start $mjd_start --mjd-end $mjd_end --tele $tele --chip "c" --dark_dir $doutdir

# Clean up logs and Report Timing
formatted_time=$(printf '%dd %dh:%dm:%ds\n' $(($SECONDS/86400)) $(($SECONDS%86400/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60)))
echo "Job completed in $formatted_time"