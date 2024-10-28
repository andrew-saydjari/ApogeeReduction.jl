#!/bin/bash
## This is a local test script, but would be replaced by a SLURM script for a cluster

# load all of the modules to talk to the database (need to be on Utah)
# should turn this off as an option for users once the MJD summaries are generated
module load sdssdb/main almanac sdsstools postgresql ffmpeg
juliaup add 1.11.0

# hardcode the mjd and expid for now
tele=$1
mjd=$2

runname="objects_${mjd}"
outdir="../outdir/"
doutdir=$outdir
almanac_file=${outdir}almanac/${runname}.h5
runlist=${outdir}almanac/runlist_${runname}.jld2

# set up the output directory (if does not exist)
mkdir -p ${outdir}almanac

# get the data summary file for the MJD
almanac -v -p 12 --mjd-start $mjd --mjd-end $mjd --${tele} --output $almanac_file

# get the runlist file (julia projects seem to refer to where your cmd prompt is when you call the shell. Here I imaging sitting at ApogeeReduction.jl level)
julia +1.11.0 --project="./" src/cal_build/make_runlist_all.jl --tele $tele --almanac_file $almanac_file --output $runlist

# Run the reduction pipeline to 2D/2Dcal and stop
julia +1.11.0 --project="./" pipeline.jl --tele $tele --runlist $runlist --outdir $doutdir --runname $runname --chips "abc" --caldir_darks "/uufs/chpc.utah.edu/common/home/sdss42/sdsswork/users/u6039752-1/working/2024_09_21/outdir/" --caldir_flats "/uufs/chpc.utah.edu/common/home/sdss42/sdsswork/users/u6039752-1/working/2024_10_03/outdir/" 

# Extract traces from dome flats
./src/cal_build/run_trace_cal.sh $tele $mjd $mjd

# Run pipeline 1D only
julia +1.11.0 --project="./" pipeline_2d_1d.jl --tele $tele --runlist $runlist --outdir $doutdir --runname $runname --chips "abc" --caldir_darks "/uufs/chpc.utah.edu/common/home/sdss42/sdsswork/users/u6039752-1/working/2024_09_21/outdir/" --caldir_flats "/uufs/chpc.utah.edu/common/home/sdss42/sdsswork/users/u6039752-1/working/2024_10_03/outdir/"

# End of night plotting script
julia +1.11.0 --project="./" src/run_scripts/plot_all.jl --tele $tele --runlist $runlist --outdir $doutdir --runname $runname --chip "a" 

# Clean up logs and Report Timing
formatted_time=$(printf '%dd %dh:%dm:%ds\n' $(($SECONDS/86400)) $(($SECONDS%86400/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60)))
echo "Job completed in $formatted_time"