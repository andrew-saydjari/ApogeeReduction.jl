#!/bin/bash
## This is a local test script, but would be replaced by a SLURM script for a cluster

# load all of the modules to talk to the database (need to be on Utah)
# should turn this off as an option for users once the MJD summaries are generated
module load sdssdb/main almanac sdsstools postgresql

# hardcode the mjd and expid for now
mjd=60356
expid=8
outdir="../outdir/"

# set up the output directory (if does not exist)
mkdir -p ${outdir}almanac

# get the data summary file for the MJD (APO for now hardcoded)
almanac --mjd $mjd --apo --fibers --output ${outdir}almanac/${mjd}.h5

# run the reduction pipeline
julia +1.10.0 pipeline.jl --mjd $mjd --expid $expid --outdir $outdir

# Clean up logs and Report Timing
formatted_time=$(printf '%dd %dh:%dm:%ds\n' $(($SECONDS/86400)) $(($SECONDS%86400/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60)))
echo "Job completed in $formatted_time"