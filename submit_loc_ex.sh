#!/bin/bash
## This is a local test script, but would be replaced by a SLURM script for a cluster

# load all of the modules to talk to the database (need to be on Utah)
# should turn this off as an option for users once the MJD summaries are generated
module load sdssdb/main almanac sdsstools postgresql ffmpeg
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

runname="apointernalflattest"
tele="apo"
mjd=59550
expid=22
outdir="../outdir/"


# set up the output directory (if does not exist)
mkdir -p ${outdir}almanac

# get the data summary file for the MJD (APO for now hardcoded)
almanac --mjd $mjd --${tele} --output ${outdir}almanac/${runname}.h5

# run the reduction pipeline
julia +1.11.0 --project="./" pipeline.jl --tele $tele --mjd $mjd --expid $expid --outdir $outdir --runname $runname --chip "a" --caldir_darks "/uufs/chpc.utah.edu/common/home/sdss42/sdsswork/users/u6039752-1/working/2024_09_21/outdir/" --caldir_flats "/uufs/chpc.utah.edu/common/home/sdss42/sdsswork/users/u6039752-1/working/2024_10_03/outdir/"
julia +1.11.0 --project="./" pipeline.jl --tele $tele --mjd $mjd --expid $expid --outdir $outdir --runname $runname --chip "b" --caldir_darks "/uufs/chpc.utah.edu/common/home/sdss42/sdsswork/users/u6039752-1/working/2024_09_21/outdir/" --caldir_flats "/uufs/chpc.utah.edu/common/home/sdss42/sdsswork/users/u6039752-1/working/2024_10_03/outdir/"
julia +1.11.0 --project="./" pipeline.jl --tele $tele --mjd $mjd --expid $expid --outdir $outdir --runname $runname --chip "c" --caldir_darks "/uufs/chpc.utah.edu/common/home/sdss42/sdsswork/users/u6039752-1/working/2024_09_21/outdir/" --caldir_flats "/uufs/chpc.utah.edu/common/home/sdss42/sdsswork/users/u6039752-1/working/2024_10_03/outdir/"

# Clean up logs and Report Timing
formatted_time=$(printf '%dd %dh:%dm:%ds\n' $(($SECONDS/86400)) $(($SECONDS%86400/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60)))
echo "Job completed in $formatted_time"