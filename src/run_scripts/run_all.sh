#!/bin/bash

# run all the data for a given night. Arguments documented below

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

# load all of the modules to talk to the database (need to be on Utah)
# should turn this off as an option for users once the MJD summaries are generated
# TODO switch to almanac/default once that's working.
module load sdssdb/main almanac/main sdsstools postgresql ffmpeg
if [ -n "$SLURM_JOB_NODELIST" ]; then
    echo $SLURM_JOB_NODELIST
else
    hostname
fi
juliaup add 1.11.0

# ARGUMENTS
# hardcode the mjd and expid for now
tele=$1
mjd=$2
run_2d_only=${3:-false}  # Third argument, defaults to false if not provided
outdir=${4:-"../outdir/"}  # Fourth argument, defaults to "../../outdir/" if not provided
caldir_darks=${5:-"/uufs/chpc.utah.edu/common/home/sdss42/sdsswork/users/u6039752-1/working/2025_02_24/outdir/"}
caldir_flats=${6:-"/uufs/chpc.utah.edu/common/home/sdss42/sdsswork/users/u6039752-1/working/2025_02_24/outdir/"}

runname="objects_${mjd}"
almanac_file=${outdir}/almanac/${runname}.h5
runlist=${outdir}/almanac/runlist_${runname}.jld2

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
julia +1.11.0 --project="./" src/run_scripts/make_runlist_all.jl --tele $tele --almanac_file $almanac_file --output $runlist

# Run the reduction pipeline to 2D/2Dcal and stop
print_elapsed_time "Running 3D->2D/2Dcal Pipeline"
# --workers_per_node 28 ## sometimes have to adjust this, could programmatically set based on the average or max read number in the exposures for that night
julia +1.11.0 --project="./" pipeline.jl --tele $tele --runlist $runlist --outdir $outdir --runname $runname --chips "abc" --caldir_darks $caldir_darks --caldir_flats $caldir_flats --workers_per_node 28

# Only continue if run_2d_only is false
if [ "$run_2d_only" != "true" ]; then
    # Extract traces from dome flats
    print_elapsed_time "Extracting Traces from Dome and Quartz Flats"
    ./src/cal_build/run_trace_cal.sh $tele $mjd $mjd $caldir_darks $caldir_flats

    # Run pipeline 1D only
    print_elapsed_time "Running 2D->1D Pipeline"
    julia +1.11.0 --project="./" pipeline_2d_1d.jl --tele $tele --runlist $runlist --outdir $outdir --runname $runname --workers_per_node 32

    # End of night plotting script
    print_elapsed_time "Making Plots"
    julia +1.11.0 --project="./" src/run_scripts/plot_all.jl --tele $tele --runlist $runlist --outdir $outdir --runname $runname --chips "abc"
fi

# Clean up logs and Report Timing
print_elapsed_time "Job Completed"
