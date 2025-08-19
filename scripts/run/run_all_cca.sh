#!/bin/bash
# Run all the data for a given night and telescope.

# Ideally you run this from a parent directory, which contains
# both the output directory and the slurm logs directory.
# (Both will be created if they don't exist.)
#
# Arguments documented below, but for example:
# sbatch path/to/ApogeeReduction/scripts/run/run_all_cca.sh apo 60855

# ------------------------------------------------------------------------------
#SBATCH --partition=cca,gen,genx
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --mem=0 #requesting all of the memory on the node

#SBATCH --time=8:00:00
#SBATCH --job-name=ar_all
#SBATCH --output=slurm_logs/%x_%j.out
# ------------------------------------------------------------------------------
set -e # exit immediately if any of the steps returns a non-zero exit code
set -o pipefail

if [ -n "$SLURM_JOB_NODELIST" ]; then
    hostname
    echo $SLURM_JOB_NODELIST
else
    hostname
fi
echo "running from $(pwd)"

if [ -n "${SLURM_JOB_ID:-}" ] ; then
    script_path=$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')
else
    script_path=$(realpath "$0")
fi
base_dir="$(dirname "$(dirname "$(dirname "$script_path")")")"
echo "base_dir: $base_dir"

git -C $base_dir status

julia_version="1.11.0" # 1.11.6
juliaup add $julia_version

# ARGUMENTS
tele=$1
mjd=$2
run_2d_only=${3:-false}  # Third argument, defaults to false if not provided
outdir=${4:-"outdir/"}  # Fourth argument, defaults to "outdir/" if not provided
outdir_madgics=${outdir}arMADGICS/raw/ # this is where the arMADGICS output will be saved
caldir_darks=${5:-"/mnt/ceph/users/asaydjari/working/2025_07_31/outdir_ref/"}
caldir_flats=${6:-"/mnt/ceph/users/asaydjari/working/2025_07_31/outdir_ref/"}
gain_read_cal_dir=${7:-"/mnt/ceph/users/asaydjari/working/2025_07_31/pass_clean/"}
path2arMADGICS=${8:-"$(dirname "$base_dir")/arMADGICS.jl/"}

runname="objects_${mjd}"
#almanac_file=${outdir}/almanac/${runname}.h5
almanac_file=${outdir}/almanac/${runname}.h5
runlist=${outdir}/almanac/runlist_${runname}.h5

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
# print_elapsed_time "Running Almanac"
# # switch to almanac -vvv for true verbosity (but only after upgrading to almanac 0.1.5)
# almanac -v -p 12 --mjd-start $mjd --mjd-end $mjd --${tele} --output $almanac_file --fibers

#print_elapsed_time "Building Runlist"
#julia +$julia_version --project=$base_dir $base_dir/scripts/run/make_runlist_all.jl --tele $tele --almanac_file $almanac_file --output $runlist
#
#print_elapsed_time "Running 3D->2D/2Dcal Pipeline"
## --workers_per_node 28 ## sometimes have to adjust this, could programmatically set based on the average or max read number in the exposures for that night
#julia +$julia_version --project=$base_dir $base_dir/pipeline.jl --tele $tele --runlist $runlist --outdir $outdir --runname $runname --chips "RGB" --caldir_darks $caldir_darks --caldir_flats $caldir_flats --workers_per_node 50 --cluster cca --gain_read_cal_dir $gain_read_cal_dir

# Only continue if run_2d_only is false
if [ "$run_2d_only" != "true" ]; then
    print_elapsed_time "Extracting Traces from Dome and Quartz Flats"
    $base_dir/scripts/cal/run_trace_cal_cca.sh $tele $mjd $mjd $outdir $caldir_darks $caldir_flats $gain_read_cal_dir

    print_elapsed_time "Running 2D->1D Pipeline"
    julia +$julia_version --project=$base_dir $base_dir/pipeline_2d_1d.jl --tele $tele --runlist $runlist --outdir $outdir --runname $runname --workers_per_node 32

    print_elapsed_time "Making Plots"
    julia +$julia_version --project=$base_dir $base_dir/scripts/run/plot_all.jl --tele $tele --runlist $runlist --outdir $outdir --runname $runname --chips "RGB"

    print_elapsed_time "Generating plot page for web viewing"
    julia +$julia_version --project=$base_dir $base_dir/scripts/run/generate_dashboard.jl --mjd $mjd --outdir $outdir
fi

## arMADGICS
if [ -d ${path2arMADGICS} ]; then
    print_elapsed_time "Running arMADGICS"
    julia +$julia_version --project=${path2arMADGICS} ${path2arMADGICS}pipeline_cca.jl --redux_base $outdir --almanac_file $almanac_file --outdir $outdir_madgics

    print_elapsed_time "Running arMADGICS Workup"
    julia +$julia_version --project=${path2arMADGICS} ${path2arMADGICS}workup.jl --outdir $outdir_madgics
fi

print_elapsed_time "Job Completed"
