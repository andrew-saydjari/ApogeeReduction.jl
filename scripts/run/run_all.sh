#!/bin/bash
# Run all the data for a given night and telescope.

# Arguments documented below, but for example (from the repo root dir):
# sbatch ./scripts/run/run_all.sh apo 60855 --mail-type=NONE

# ------------------------------------------------------------------------------
#SBATCH --account=sdss-np
#SBATCH --partition=sdss-shared-np
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64

#SBATCH --mem=0 #requesting all of the memory on the node

#SBATCH --time=8:00:00
#SBATCH --job-name=ar_all
#SBATCH --output=slurm_logs/%x_%j.out

#SBATCH --mail-type=ALL
#SBATCH --mail-user=7155301634@vtext.com
# ------------------------------------------------------------------------------

# exit the script immediately if a command fails
set -e
set -o pipefail

# load all of the modules to talk to the database (need to be on Utah)
# should turn this off as an option for users once the MJD summaries are generated
# TODO switch to almanac/default once that's working.
module load sdssdb/main almanac/main sdsstools postgresql ffmpeg
if [ -n "$SLURM_JOB_NODELIST" ]; then
    echo $SLURM_JOB_NODELIST
else
    hostname
fi
echo "running from $(pwd)"

juliaup add 1.11.0

# ARGUMENTS
# hardcode the mjd and expid for now
tele=$1
mjd=$2
run_2d_only=${3:-false}  # Third argument, defaults to false if not provided
outdir=${4:-"../outdir/"}  # Fourth argument, defaults to "../../outdir/" if not provided
caldir_darks=${5:-"/uufs/chpc.utah.edu/common/home/sdss42/sdsswork/users/u6039752-1/working/2025_06_16/outdir/"}
caldir_flats=${6:-"/uufs/chpc.utah.edu/common/home/sdss42/sdsswork/users/u6039752-1/working/2025_06_16/outdir/"}
path2arMADGICS=${7:-"../arMADGICS.jl/"}

runname="objects_${mjd}"
almanac_file=${outdir}/almanac/${runname}.h5
runlist=${outdir}/almanac/runlist_${runname}.h5

# set up the output directory (if does not exist)
mkdir -p ${outdir}/almanac

# nice function to print the elapsed time at each step
print_elapsed_time() {
    formatted_time=$(printf '%dd %dh:%dm:%ds\n' $(($SECONDS/86400)) $(($SECONDS%86400/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60)))
    echo
    echo "Elapsed time: $formatted_time"
    echo
    echo "--------- $1 ---------"
    echo
}

# get the data summary file for the MJD
print_elapsed_time "Running Almanac"
# switch to almanac -vvv for true verbosity (but only after upgrading to almanac 0.1.5)
almanac -v -p 12 --mjd-start $mjd --mjd-end $mjd --${tele} --output $almanac_file --fibers

# get the runlist file (julia projects seem to refer to where your cmd prompt is when you call the shell. Here I imagine sitting at ApogeeReduction.jl level)
print_elapsed_time "Building Runlist"
julia +1.11.0 --project="./" scripts/run/make_runlist_all.jl --tele $tele --almanac_file $almanac_file --output $runlist

print_elapsed_time "Running 3D->2D/2Dcal Pipeline"
# --workers_per_node 28 ## sometimes have to adjust this, could programmatically set based on the average or max read number in the exposures for that night
julia +1.11.0 --project="./" pipeline.jl --tele $tele --runlist $runlist --outdir $outdir --runname $runname --chips "RGB" --caldir_darks $caldir_darks --caldir_flats $caldir_flats --workers_per_node 50

# Only continue if run_2d_only is false
if [ "$run_2d_only" != "true" ]; then
    print_elapsed_time "Extracting Traces from Dome and Quartz Flats"
    ./scripts/cal/run_trace_cal.sh $tele $mjd $mjd $caldir_darks $caldir_flats

    print_elapsed_time "Running 2D->1D Pipeline"
    julia +1.11.0 --project="./" pipeline_2d_1d.jl --tele $tele --runlist $runlist --outdir $outdir --runname $runname --workers_per_node 32

    print_elapsed_time "Making Plots"
    julia +1.11.0 --project="./" scripts/run/plot_all.jl --tele $tele --runlist $runlist --outdir $outdir --runname $runname --chips "RGB"

    print_elapsed_time "Generating plot page for web viewing"
    julia +1.11.0 --project="./" scripts/run/generate_dashboard.jl --mjd $mjd --outdir $outdir
fi

## arMADGICS
if [ -d ${path2arMADGICS} ]; then
    print_elapsed_time "Running arMADGICS"
    julia +1.11.0 --project=${path2arMADGICS} ${path2arMADGICS}pipeline.jl --redux_base $outdir --almanac_file $almanac_file

    print_elapsed_time "Running arMADGICS Workup"
    julia +1.11.0 --project=${path2arMADGICS} ${path2arMADGICS}workup.jl --outdir ${outdir}arMADGICS/raw/
fi

print_elapsed_time "Job Completed"