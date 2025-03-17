#!/bin/bash
## This is a local test script, but would be replaced by a SLURM script for a cluster

# load all of the modules to talk to the database (need to be on Utah)
# should turn this off as an option for users once the MJD summaries are generated
module load sdssdb/main almanac/0.1.4 sdsstools postgresql ffmpeg
juliaup add 1.11.0

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

# runname="apodometest"
# tele="apo"
# mjd=60546
# expid=13
# outdir="../outdir/"

# runname="apointernalflattest"
# tele="apo"
# mjd=60575
# expid=96
# outdir="../outdir/"

# runname="apointernalflattest"
# tele="apo"
# mjd=60574
# expid=49
# outdir="../outdir/"

# runname="lcointernalflattest"
# tele="lco"
# mjd=59236
# expid=14
# outdir="../outdir/"

# runname="singleCalTest"
# tele="apo"
# mjd=60606
# expid=19
# outdir="../outdir/"

# runname="singleCalTest"
# tele="lco"
# mjd=60606
# expid=45
# outdir="../outdir/"

runname="ap1DTest"
tele="apo"
mjd=60606
expid=19
outdir="../../2024_12_02/outdir/"

# set up the output directory (if does not exist)
mkdir -p ${outdir}almanac

# get the data summary file for the MJD (APO for now hardcoded)
almanac --mjd $mjd --${tele} --output ${outdir}almanac/${runname}.h5

# run the reduction pipeline
julia +1.11.0 --project="./" pipeline.jl --tele $tele --mjd $mjd --expid $expid --outdir $outdir --runname $runname --chip "abc" --caldir_darks "/uufs/chpc.utah.edu/common/home/sdss42/sdsswork/users/u6039752-1/working/2024_09_21/outdir/" --caldir_flats "/uufs/chpc.utah.edu/common/home/sdss42/sdsswork/users/u6039752-1/working/2024_10_03/outdir/"

# this part only works if you have the traces extracted, which takes like 10 min per chip per domeflat exposure... (removed b and c for this reason in my tests)
julia +1.11.0 --project="./" pipeline_2d_1d.jl --tele $tele --mjd $mjd --expid $expid --outdir $outdir --runname $runname --chip "a"

julia +1.11.0 --project="./" src/run_scripts/plot_all.jl --chip "a" --tele $tele --mjd $mjd --expid $expid --outdir $outdir --runname $runname

# Clean up logs and Report Timing
formatted_time=$(printf '%dd %dh:%dm:%ds\n' $(($SECONDS/86400)) $(($SECONDS%86400/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60)))
echo "Job completed in $formatted_time"
