#!/bin/bash
#SBATCH --account=sdss-np
#SBATCH --partition=sdss-shared-np
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64

#SBATCH --mem=0 #requesting all of the memory on the node

#SBATCH --time=8:00:00
#SBATCH --job-name=ApogeeReduction_trace_cal
#SBATCH --output=slurm_logs/%x_%j.out

#SBATCH --mail-type=ALL
#SBATCH --mail-user=7155301634@vtext.com
# ------------------------------------------------------------------------------

# load all of the modules to talk to the database (need to be on Utah)
# should turn this off as an option for users once the MJD summaries are generated
module load sdssdb/main almanac/main sdsstools postgresql ffmpeg
juliaup add 1.11.0

# hardcode the mjd and expid for now
tele=$1
mjd_start=$2
mjd_end=$3

runname="trace_cal_${tele}_${mjd_start}_${mjd_end}"
outdir="../outdir/"
doutdir=$outdir
almanac_file=${outdir}almanac/${runname}.h5
runlist=${outdir}almanac/runlist_${runname}.jld2
caldir_darks=${4:-"/uufs/chpc.utah.edu/common/home/sdss42/sdsswork/users/u6039752-1/working/2025_02_24/outdir/"}
caldir_flats=${5:-"/uufs/chpc.utah.edu/common/home/sdss42/sdsswork/users/u6039752-1/working/2025_02_24/outdir/"}

# set up the output directory (if does not exist)
mkdir -p ${outdir}almanac

# get the data summary file for the MJD
almanac -v -p 12 --mjd-start $mjd_start --mjd-end $mjd_end --${tele} --output $almanac_file

# get the runlist file (julia projects seem to refer to where your cmd prompt is when you call the shell. Here I imaging sitting at ApogeeReduction.jl level)
echo "Getting domeflat runlist..."
julia +1.11.0 --project="./" src/cal_build/make_runlist_dome_flats.jl --tele $tele --almanac_file $almanac_file --output $runlist

# run the reduction pipeline (all cals like dark sub/flats that would be a problem, should be post 3D->2D extraction)
echo "Running 3D->2D..."
julia +1.11.0 --project="./" pipeline.jl --tele $tele --runlist $runlist --outdir $outdir --runname $runname --chips "abc" --caldir_darks $caldir_darks --caldir_flats $caldir_flats

mkdir -p ${outdir}dome_flats
julia +1.11.0 --project="./" src/cal_build/make_traces_domeflats.jl --mjd-start $mjd_start --mjd-end $mjd_end --tele $tele --trace_dir ${doutdir} --runlist $runlist

echo "Running 2D->1D Pipeline for DomeFlats"
julia +1.11.0 --project="./" pipeline_2d_1d.jl --tele $tele --runlist $runlist --outdir $outdir --runname $runname

runlist_quartz=${outdir}almanac/runlist_${runname}.jld2
# get the runlist file (julia projects seem to refer to where your cmd prompt is when you call the shell. Here I imaging sitting at ApogeeReduction.jl level)
echo "Getting quartzflat runlist..."
julia +1.11.0 --project="./" src/cal_build/make_runlist_quartz_flats.jl --tele $tele --almanac_file $almanac_file --output $runlist_quartz

# run the reduction pipeline (all cals like dark sub/flats that would be a problem, should be post 3D->2D extraction)
echo "Running 3D->2D..."
julia +1.11.0 --project="./" pipeline.jl --tele $tele --runlist $runlist_quartz --outdir $outdir --runname $runname --chips "abc" --caldir_darks $caldir_darks --caldir_flats $caldir_flats

mkdir -p ${outdir}quartz_flats
julia +1.11.0 --project="./" src/cal_build/make_traces_quartzflats.jl --mjd-start $mjd_start --mjd-end $mjd_end --tele $tele --trace_dir ${doutdir} --runlist $runlist_quartz

echo "Running 2D->1D Pipeline for QuartzFlats"
julia +1.11.0 --project="./" pipeline_2d_1d.jl --tele $tele --runlist $runlist_quartz --outdir $outdir --runname $runname

# Clean up logs and Report Timing
formatted_time=$(printf '%dd %dh:%dm:%ds\n' $(($SECONDS/86400)) $(($SECONDS%86400/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60)))
echo "Job completed in $formatted_time"
