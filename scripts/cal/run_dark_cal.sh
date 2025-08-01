#!/bin/bash
#SBATCH --account=sdss-np
#SBATCH --partition=sdss-shared-np
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64

#SBATCH --mem=0 #requesting all of the memory on the node

#SBATCH --time=8:00:00
#SBATCH --job-name=ar_dark_cal
#SBATCH --output=slurm_logs/%x_%j.out

#SBATCH --mail-type=ALL
#SBATCH --mail-user=7155301634@vtext.com
# ------------------------------------------------------------------------------

# exit immediately if any of the steps returns a non-zero exit code
set -e
set -o pipefail

# load all of the modules to talk to the database (need to be on Utah)
# should turn this off as an option for users once the MJD summaries are generated
module load sdssdb/main almanac/main sdsstools postgresql ffmpeg
juliaup add 1.11.0

# hardcode the mjd and expid for now
tele=$1
mjd_start=$2
mjd_end=$3
dropfirstn=${4:-0}

runname="dark_cal_${tele}_${mjd_start}_${mjd_end}"
outdir="../outdir/"
doutdir=$outdir
almanac_file=${outdir}almanac/${runname}.h5
runlist=${outdir}almanac/runlist_${runname}.h5

# set up the output directory (if does not exist)
mkdir -p ${outdir}almanac

# get the data summary file for the MJD
almanac -v -p 12 --mjd-start $mjd_start --mjd-end $mjd_end --${tele} --output $almanac_file

# get the runlist file (julia projects seem to refer to where your cmd prompt is when you call the shell. Here I imaging sitting at ApogeeReduction.jl level)
julia +1.11.0 --project="./" scripts/cal/make_runlist_darks.jl --tele $tele --almanac_file $almanac_file --output $runlist

# run the reduction pipeline (all cals like dark sub/flats that would be a problem, should be post 3D->2D extraction)
## 100 read darks require fewer workers per node (a problem for the ultra darks)
julia +1.11.0 --project="./" pipeline.jl --workers_per_node 28 --tele $tele --runlist $runlist --outdir $doutdir --runname $runname --chips "RGB" --doCal2d false

# make the stacked darks
mkdir -p ${outdir}darks
julia +1.11.0 --project="./" scripts/cal/make_stack_darks.jl --mjd-start $mjd_start --mjd-end $mjd_end --tele $tele --chip "R" --dropfirstn $dropfirstn --dark_dir ${doutdir} --runlist $runlist
julia +1.11.0 --project="./" scripts/cal/make_stack_darks.jl --mjd-start $mjd_start --mjd-end $mjd_end --tele $tele --chip "G" --dropfirstn $dropfirstn --dark_dir ${doutdir} --runlist $runlist
julia +1.11.0 --project="./" scripts/cal/make_stack_darks.jl --mjd-start $mjd_start --mjd-end $mjd_end --tele $tele --chip "B" --dropfirstn $dropfirstn --dark_dir ${doutdir} --runlist $runlist

# Clean up logs and Report Timing
formatted_time=$(printf '%dd %dh:%dm:%ds\n' $(($SECONDS/86400)) $(($SECONDS%86400/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60)))
echo "Job completed in $formatted_time"
