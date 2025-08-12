#!/bin/bash
# Run all the data for a given night and telescope.

# Arguments documented below, but for example (from the repo root dir):
# sbatch ./scripts/bulk/run_bulk.sh 60584 60591
#60584 60614, 60796 60826

# ------------------------------------------------------------------------------
#SBATCH --partition=cca,gen
#SBATCH --constraint="[genoa|icelake|rome]"
#SBATCH --nodes=1 #3

#SBATCH --time=8:00:00
#SBATCH --job-name=ar_bulk
#SBATCH --output=slurm_logs/%x_%j.out
# ------------------------------------------------------------------------------
SLURM_NTASKS=$(($SLURM_CPUS_ON_NODE * $SLURM_NNODES))
export SLURM_NTASKS
# Print all SLURM environment variables
env | grep SLURM | while read -r line; do
    echo "$line"
done

# exit the script immediately if a command fails
set -e
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

julia_version="1.11.0" # 1.11.6
juliaup add $julia_version

# ARGUMENTS
# hardcode the mjd and expid for now
mjd_start=$1
mjd_end=$2
run_2d_only=${3:-false}  # Third argument, defaults to false if not provided
outdir=${4:-"outdir/"}  # Fourth argument, defaults to "outdir/" if not provided
caldir_darks=${5:-"/mnt/ceph/users/asaydjari/working/2025_07_31/outdir_ref/"}
caldir_flats=${6:-"/mnt/ceph/users/asaydjari/working/2025_07_31/outdir_ref/"}
gain_read_cal_dir=${7:-"/mnt/ceph/users/asaydjari/working/2025_07_31/pass_clean/"}
path2arMADGICS=${8:-"$(dirname "$base_dir")/arMADGICS.jl/"}

runname="allobs_${mjd_start}_${mjd_end}"
almanac_file=${outdir}almanac/${runname}.h5
runlist=${outdir}almanac/runlist_${runname}.h5
tele_list=("apo" "lco")
flat_types=("quartz" "dome")
# set up the output directory (if does not exist)
mkdir -p ${outdir}/almanac

# nice function to print the elapsed time at each step
print_elapsed_time() {
    current_seconds=$SECONDS
    if [ -n "$LAST_TIME" ]; then
        diff_seconds=$((current_seconds - LAST_TIME))
        diff_time=$(printf '%dd %dh:%dm:%ds\n' $(($diff_seconds/86400)) $(($diff_seconds%86400/3600)) $(($diff_seconds%3600/60)) $(($diff_seconds%60)))
        echo
        echo "Time since last step: $diff_time"
    fi
    formatted_time=$(printf '%dd %dh:%dm:%ds\n' $(($current_seconds/86400)) $(($current_seconds%86400/3600)) $(($current_seconds%3600/60)) $(($current_seconds%60)))
    echo "Elapsed time: $formatted_time"
    echo
    echo "--------- $1 ---------"
    echo
    LAST_TIME=$current_seconds
}

# # get the data summary file for the MJD
# print_elapsed_time "Running Almanac"
# # would like to switch back to -p 12 once the multinode context is fixed
# almanac -v --mjd-start $mjd_start --mjd-end $mjd_end --output $almanac_file --fibers

print_elapsed_time "Building Runlist"
julia +$julia_version --project=$base_dir $base_dir/scripts/run/make_runlist_all.jl --almanac_file $almanac_file --output $runlist

# for tele in ${tele_list[@]}
# do
#     print_elapsed_time "Running 3D->2D/2Dcal Pipeline for $tele"
#     ## sometimes have to adjust workers_per_node based on nreads, could programmatically set based on the average or max read number in the exposures for that night
#     julia +$julia_version --project=$base_dir $base_dir/pipeline.jl --tele $tele --runlist $runlist --outdir $outdir --runname $runname --chips "RGB" --caldir_darks $caldir_darks --caldir_flats $caldir_flats --cluster cca --gain_read_cal_dir $gain_read_cal_dir
# done

# Only continue if run_2d_only is false
if [ "$run_2d_only" != "true" ]; then
    # ### Traces and Refluxing
    # # would really like to combine the different telescopes here
    # for flat_type in ${flat_types[@]}
    # do
    #     flatrunlist=${outdir}almanac/runlist_${flat_type}_${runname}.h5
    #     print_elapsed_time "Making runlist for $flat_type Flats"
    #     julia +$julia_version --project=$base_dir $base_dir/scripts/cal/make_runlist_fiber_flats.jl --almanac_file $almanac_file --output $flatrunlist --flat_type $flat_type

    #     for tele in ${tele_list[@]}
    #     do
    #         print_elapsed_time "Extracting Traces from $flat_type Flats for $tele"
    #         mkdir -p ${outdir}${flat_type}_flats
    #         julia +$julia_version --project=$base_dir $base_dir/scripts/cal/make_traces_from_flats.jl --tele $tele --trace_dir ${outdir} --runlist $flatrunlist --flat_type $flat_type --slack_quiet true

    #         print_elapsed_time "Running 2D->1D Pipeline without relFlux for $flat_type Flats for $tele"
    #         mkdir -p ${outdir}/apredrelflux
    #         julia +$julia_version --project=$base_dir $base_dir/pipeline_2d_1d.jl --tele $tele --runlist $flatrunlist --outdir $outdir --runname $runname --relFlux false --waveSoln false
    #     done

    #     print_elapsed_time "Making relFlux for $flat_type Flats"
    #     julia +$julia_version --project=$base_dir $base_dir/scripts/cal/make_relFlux.jl --trace_dir ${outdir} --runlist $flatrunlist --runname $runname   
    # done

    ### 2D->1D
    for tele in ${tele_list[@]}
    do
        print_elapsed_time "Running 2D->1D Pipeline for $tele"
        julia +$julia_version --project=$base_dir $base_dir/pipeline_2d_1d.jl --tele $tele --runlist $runlist --outdir $outdir --runname $runname
    done
    ### Plots
    #     print_elapsed_time "Making Plots"
    #     julia +$julia_version --project=$base_dir $base_dir/scripts/run/plot_all.jl --tele $tele --runlist $runlist --outdir $outdir --runname $runname --chips "RGB"

    # ## Dashboard
    # print_elapsed_time "Generating plot page for web viewing"
    # julia +$julia_version --project=$base_dir $base_dir/scripts/run/generate_dashboard.jl --mjd $mjd --outdir $outdir

    ## arMADGICS
    if [ -d ${path2arMADGICS} ]; then
        print_elapsed_time "Running arMADGICS"
        julia +$julia_version --project=${path2arMADGICS} ${path2arMADGICS}pipeline.jl --redux_base $outdir --almanac_file $almanac_file

        print_elapsed_time "Running arMADGICS Workup"
        julia +$julia_version --project=${path2arMADGICS} ${path2arMADGICS}workup.jl --outdir ${outdir}arMADGICS/raw/
    fi
fi

print_elapsed_time "Job Completed"