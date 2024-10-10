#!/bin/bash
## This is a local test script, but would be replaced by a SLURM script for a cluster

# load all of the modules to talk to the database (need to be on Utah)
# should turn this off as an option for users once the MJD summaries are generated
module load sdssdb/main almanac sdsstools postgresql ffmpeg
juliaup add 1.11.0

# hardcode the mjd and expid for now
tele=$1
mjd_start=$2
mjd_end=$3

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
julia +1.11.0 --project="./" src/cal_build/make_runlist_darks.jl --tele $tele --almanac_file $almanac_file --output $runlist

# run the reduction pipeline (all cals like dark sub/flats that would be a problem, should be post 3D->2D extraction)
julia +1.11.0 --project="./" pipeline.jl --tele $tele --runlist $runlist --outdir $doutdir --runname $runname --chips "abc"

# make the stacked darks
mkdir -p ${outdir}darks
julia +1.11.0 --project="./" src/cal_build/make_stack_darks.jl --mjd-start $mjd_start --mjd-end $mjd_end --tele $tele --chip "a" --dark_dir ${doutdir} --runlist $runlist
julia +1.11.0 --project="./" src/cal_build/make_stack_darks.jl --mjd-start $mjd_start --mjd-end $mjd_end --tele $tele --chip "b" --dark_dir ${doutdir} --runlist $runlist
julia +1.11.0 --project="./" src/cal_build/make_stack_darks.jl --mjd-start $mjd_start --mjd-end $mjd_end --tele $tele --chip "c" --dark_dir ${doutdir} --runlist $runlist

# Clean up logs and Report Timing
formatted_time=$(printf '%dd %dh:%dm:%ds\n' $(($SECONDS/86400)) $(($SECONDS%86400/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60)))
echo "Job completed in $formatted_time"