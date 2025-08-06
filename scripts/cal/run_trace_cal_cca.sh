#!/bin/bash
#SBATCH --partition=cca,gen,genx
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --mem=0 #requesting all of the memory on the node

#SBATCH --time=8:00:00
#SBATCH --job-name=ar_dark_cal
#SBATCH --output=../slurm_logs/%x_%j.out
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

julia_version="1.11.6"
juliaup add $julia_version

tele=$1
mjd_start=$2
mjd_end=$3

outdir="/mnt/ceph/users/asaydjari/working/2025_07_31/outdir/"
doutdir=$outdir
runname="objects_${mjd_start}"
almanac_file=${outdir}almanac/${runname}.h5
domerunlist=${outdir}almanac/runlist_dome_${runname}.h5
quartzrunlist=${outdir}almanac/runlist_quartz_${runname}.h5
caldir_darks=${4:-"/mnt/ceph/users/asaydjari/working/2025_07_31/outdir_ref/"}
caldir_flats=${5:-"/mnt/ceph/users/asaydjari/working/2025_07_31/outdir_ref/"}
gain_read_cal_dir=${6:-"/mnt/ceph/users/asaydjari/working/2025_07_31/pass_clean/"}

# nice function to print the elapsed time at each step
print_elapsed_time() {
    formatted_time=$(printf '%dd %dh:%dm:%ds\n' $(($SECONDS/86400)) $(($SECONDS%86400/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60)))
    echo
    echo "Elapsed time: $formatted_time"
    echo
    echo "--------- $1 ---------"
    echo
}

# get the runlist file (julia projects seem to refer to where your cmd prompt is when you call the shell. Here I imaging sitting at ApogeeReduction.jl level)
print_elapsed_time "Getting domeflat runlist..."
julia +$julia_version --project=$base_dir $base_dir/scripts/cal/make_runlist_dome_flats.jl --tele $tele --almanac_file $almanac_file --output $domerunlist

# run the reduction pipeline (all cals like dark sub/flats that would be a problem, should be post 3D->2D extraction)
print_elapsed_time "Running 3D->2D..."
julia +$julia_version --project=$base_dir $base_dir/pipeline.jl --tele $tele --runlist $domerunlist --outdir $outdir --runname $runname --chips "RGB" --caldir_darks $caldir_darks --caldir_flats $caldir_flats --cluster cca --gain_read_cal_dir $gain_read_cal_dir

mkdir -p ${outdir}dome_flats
julia +$julia_version --project=$base_dir $base_dir/scripts/cal/make_traces_domeflats.jl --mjd-start $mjd_start --mjd-end $mjd_end --tele $tele --trace_dir ${doutdir} --runlist $domerunlist

#### Parallel reduction for quartz flats ####
print_elapsed_time "Getting quartzflat runlist..."
julia +$julia_version --project=$base_dir $base_dir/scripts/cal/make_runlist_quartz_flats.jl --tele $tele --almanac_file $almanac_file --output $quartzrunlist

# run the reduction pipeline (all cals like dark sub/flats that would be a problem, should be post 3D->2D extraction)
print_elapsed_time "Running 3D->2D..."
julia +$julia_version --project=$base_dir $base_dir/pipeline.jl --tele $tele --runlist $quartzrunlist --outdir $outdir --runname $runname --chips "RGB" --caldir_darks $caldir_darks --caldir_flats $caldir_flats --cluster cca --gain_read_cal_dir $gain_read_cal_dir

mkdir -p ${outdir}quartz_flats
julia +$julia_version --project=$base_dir $base_dir/scripts/cal/make_traces_quartzflats.jl --mjd-start $mjd_start --mjd-end $mjd_end --tele $tele --trace_dir ${doutdir} --runlist $quartzrunlist

### 2D->1D and relFlux ####
# DomeFlats
print_elapsed_time "Running 2D->1D Pipeline for DomeFlats"
julia +$julia_version --project=$base_dir $base_dir/pipeline_2d_1d.jl --tele $tele --runlist $domerunlist --outdir $outdir --runname $runname --relFlux false

julia +$julia_version --project=$base_dir $base_dir/scripts/cal/make_relFlux.jl --mjd-start $mjd_start --mjd-end $mjd_end --tele $tele --trace_dir ${doutdir} --runlist $domerunlist --runname $runname

# QuartzFlats
print_elapsed_time "Running 2D->1D Pipeline for QuartzFlats"
julia +$julia_version --project=$base_dir $base_dir/pipeline_2d_1d.jl --tele $tele --runlist $quartzrunlist --outdir $outdir --runname $runname --relFlux false

julia +$julia_version --project=$base_dir $base_dir/scripts/cal/make_relFlux.jl --mjd-start $mjd_start --mjd-end $mjd_end --tele $tele --trace_dir ${doutdir} --runlist $quartzrunlist --runname $runname

# Clean up logs and Report Timing
print_elapsed_time "Job completed"
