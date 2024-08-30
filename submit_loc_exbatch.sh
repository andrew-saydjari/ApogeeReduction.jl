#!/bin/bash
## This is a local test script, but would be replaced by a SLURM script for a cluster

# load all of the modules to talk to the database (need to be on Utah)
# should turn this off as an option for users once the MJD summaries are generated
module load sdssdb/main almanac sdsstools postgresql

# hardcode the mjd and expid for now
runname="testrun"
runlist="/uufs/chpc.utah.edu/common/home/u6039752/scratch1/working/2024_08_14/testrun.jld2"
outdir="../outdir/"
tele="apo"

# set up the output directory (if does not exist)
mkdir -p ${outdir}almanac

# Extract MJD min and max values using Julia
mjd_values=$(julia --project="./" -e "
    using JLD2
    file = jldopen(\"$runlist\", \"r\")
    mjd_data = read(file, \"mjd\")
    close(file)
    println(minimum(mjd_data))
    println(maximum(mjd_data))
")

# Extract the start and end MJD from the Julia output
mjd_start=$(echo "$mjd_values" | sed -n '1p')
mjd_end=$(echo "$mjd_values" | sed -n '2p')

# get the data summary file for the MJD (for a specific telescope... do we need to handle options for both? but cals would be hard)
## should I choose to invoke almanac with parallism? This was buggy before.
almanac --mjd-start $mjd_start --mjd-end $mjd_end --${tele} --output ${outdir}almanac/${runname}.h5

# run the reduction pipeline
julia +1.10.0 --project="./" pipeline.jl --tele $tele --runlist $runlist --outdir $outdir --runname $runname

# Clean up logs and Report Timing
formatted_time=$(printf '%dd %dh:%dm:%ds\n' $(($SECONDS/86400)) $(($SECONDS%86400/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60)))
echo "Job completed in $formatted_time"