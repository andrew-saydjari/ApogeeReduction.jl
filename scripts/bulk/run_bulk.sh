#!/bin/bash
# Run all the data for a given night and telescope.

# Arguments documented below, but for example (from the repo root dir):
# sbatch ./scripts/bulk/run_bulk.sh 60584 60591
#60584 60614, 60796 60826

#SBATCH --account=sdss-np
#SBATCH --partition=sdss-shared-np
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=64

#SBATCH --mem=0 #requesting all of the memory on the node

#SBATCH --time=8:00:00
#SBATCH --job-name=ar_bulk
#SBATCH --output=slurm_logs/%x_%j.out

# ------------------------------------------------------------------------------

# exit the script immediately if a command fails
set -e
set -o pipefail

# load all of the modules to talk to the database (need to be on Utah)
# should turn this off as an option for users once the MJD summaries are generated
# TODO switch to almanac/default once that's working.
module load sdssdb/main almanac/main sdsstools postgresql ffmpeg
if [ -n "$SLURM_JOB_NODELIST" ]; then
    hostname
    echo $SLURM_JOB_NODELIST
else
    hostname
fi
echo "running from $(pwd)"

juliaup add 1.11.0

# ARGUMENTS
# hardcode the mjd and expid for now
mjd_start=$1
mjd_end=$2
run_2d_only=${3:-false}  # Third argument, defaults to false if not provided
outdir=${4:-"../outdir/"}  # Fourth argument, defaults to "../outdir/" if not provided
caldir_darks=${5:-"/uufs/chpc.utah.edu/common/home/sdss42/sdsswork/users/u6039752-1/working/2025_06_16/outdir/"}
caldir_flats=${6:-"/uufs/chpc.utah.edu/common/home/sdss42/sdsswork/users/u6039752-1/working/2025_06_16/outdir/"}
path2arMADGICS=${7:-"../arMADGICS.jl/"}

runname="allobs_${mjd_start}_${mjd_end}"
almanac_file=${outdir}almanac/${runname}.h5
runlist=${outdir}almanac/runlist_${runname}.h5
tele_list=("apo" "lco")
flat_types=("quartz" "dome")
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

### ANDREW YOU REALLY NEED TO THINK ABOUT WORKERS PER NODE BEFORE LAUNCHING
### REMOVING SOME OF THE UNCAL SAVES WOULD ALSO BE A GOOD IDEA

# # get the data summary file for the MJD
# print_elapsed_time "Running Almanac"
# # would like to switch back to -p 12 once the multinode context is fixed
# almanac -v --mjd-start $mjd_start --mjd-end $mjd_end --output $almanac_file --fibers

# print_elapsed_time "Building Runlist"
# julia +1.11.0 --project="./" scripts/run/make_runlist_all.jl --almanac_file $almanac_file --output $runlist

# for tele in ${tele_list[@]}
# do
#     print_elapsed_time "Running 3D->2D/2Dcal Pipeline for $tele"
#     # --workers_per_node 28 ## sometimes have to adjust this, could programmatically set based on the average or max read number in the exposures for that night
#     julia +1.11.0 --project="./" pipeline.jl --tele $tele --runlist $runlist --outdir $outdir --runname $runname --chips "RGB" --caldir_darks $caldir_darks --caldir_flats $caldir_flats --workers_per_node 50
# done

# Only continue if run_2d_only is false
# if [ "$run_2d_only" != "true" ]; then
#     ### Traces and Refluxing
#     # would really like to combine the different telescopes here
#     for flat_type in ${flat_types[@]}
#     do
        # flatrunlist=${outdir}almanac/runlist_${flat_type}_${runname}.h5
        # print_elapsed_time "Making runlist for $flat_type Flats"
        # julia +1.11.0 --project="./" scripts/cal/make_runlist_fiber_flats.jl --almanac_file $almanac_file --output $flatrunlist --flat_type $flat_type

        # for tele in ${tele_list[@]}
        # do
        #     print_elapsed_time "Extracting Traces from $flat_type Flats for $tele"
        #     mkdir -p ${outdir}${flat_type}_flats
        #     julia +1.11.0 --project="./" scripts/cal/make_traces_from_flats.jl --tele $tele --trace_dir ${outdir} --runlist $flatrunlist --flat_type $flat_type

        #     print_elapsed_time "Running 2D->1D Pipeline without relFlux for $flat_type Flats for $tele"
        #     mkdir -p ${outdir}/apredrelflux
        #     julia +1.11.0 --project="./" pipeline_2d_1d.jl --tele $tele --runlist $flatrunlist --outdir $outdir --runname $runname --relFlux false --waveSoln false
        # done

        # print_elapsed_time "Making relFlux for $flat_type Flats"
        # julia +1.11.0 --project="./" scripts/cal/make_relFlux.jl --trace_dir ${outdir} --runlist $flatrunlist --runname $runname   
    # done

    ### 2D->1D
    for tele in ${tele_list[@]}
    do
        print_elapsed_time "Running 2D->1D Pipeline"
        julia +1.11.0 --project="./" pipeline_2d_1d.jl --tele $tele --runlist $runlist --outdir $outdir --runname $runname --workers_per_node 64
    done
    ### Plots
    #     print_elapsed_time "Making Plots"
    #     julia +1.11.0 --project="./" scripts/run/plot_all.jl --tele $tele --runlist $runlist --outdir $outdir --runname $runname --chips "RGB"

    # ## Dashboard
    # print_elapsed_time "Generating plot page for web viewing"
    # julia +1.11.0 --project="./" scripts/run/generate_dashboard.jl --mjd $mjd --outdir $outdir

    ### arMADGICS
    # if [ -d ${path2arMADGICS} ]; then
    #     print_elapsed_time "Running arMADGICS"
    #     julia +1.11.0 --project=${path2arMADGICS} ${path2arMADGICS}pipeline.jl --redux_base $outdir --almanac_file $almanac_file

    #     print_elapsed_time "Running arMADGICS Workup"
    #     julia +1.11.0 --project=${path2arMADGICS} ${path2arMADGICS}workup.jl --outdir ${outdir}arMADGICS/raw/
    # fi
fi

print_elapsed_time "Job Completed"