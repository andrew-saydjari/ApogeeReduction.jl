#!/bin/bash
## This is a local test script, but would be replaced by a SLURM script for a cluster

# load all of the modules to talk to the database (need to be on Utah)
# should turn this off as an option for users once the MJD summaries are generated
module load sdssdb/main almanac sdsstools postgresql ffmpeg

# hardcode the mjd and expid for now
# runname="arctest"
# mjd=60356
# expid=8
# outdir="../outdir/"

# runname="darktest"
# mjd=60399
# expid=79
# outdir="../outdir/"

# runname="darktest2"
# tele="apo"
# mjd=60350
# expid=16
# outdir="../outdir/"

# runname="lcodometest"
# tele="lco"
# mjd=60555
# expid=10
# outdir="../outdir/"

runname="apodometest"
tele="apo"
mjd=60546
expid=13
outdir="../outdir/"

# set up the output directory (if does not exist)
mkdir -p ${outdir}almanac

# get the data summary file for the MJD (APO for now hardcoded)
almanac --mjd $mjd --${tele} --output ${outdir}almanac/${runname}.h5

# run the reduction pipeline
julia +1.10.0 --project="./" pipeline.jl --tele $tele --mjd $mjd --expid $expid --outdir $outdir --runname $runname --chip "a"
julia +1.10.0 --project="./" pipeline.jl --tele $tele --mjd $mjd --expid $expid --outdir $outdir --runname $runname --chip "b"
julia +1.10.0 --project="./" pipeline.jl --tele $tele --mjd $mjd --expid $expid --outdir $outdir --runname $runname --chip "c"

# Clean up logs and Report Timing
formatted_time=$(printf '%dd %dh:%dm:%ds\n' $(($SECONDS/86400)) $(($SECONDS%86400/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60)))
echo "Job completed in $formatted_time"