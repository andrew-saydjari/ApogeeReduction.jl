#!/bin/bash
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

# Define sbatchCustom function
function sbatchCustom() {
    export SLURM_ID=$(sbatch --parsable "$@")
    mkdir -p slurm_logs
    scontrol write batch_script $SLURM_ID slurm_logs/slurm-$SLURM_ID.sh
}

data_dir=$(realpath "outdir/")

# set up the output directory (if does not exist)
mkdir -p ${data_dir}/plots

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

echo "Processing back-to-back flat observations from special_cal_obs.txt..."

# Arrays to store job information
declare -a job_ids=()
declare -a job_info=()

# First pass: Submit all sbatch jobs
while read -r line; do
    # Skip comments and empty lines
    [[ "$line" =~ ^#.*$ || -z "$line" ]] && continue

    # Read fields from each line
    read -r tele mjd expid_range cal_type description <<< "$line"

    # Only process b2b calibration types
    if [ "$cal_type" = "b2b" ]; then
        # Extract start and end expids from range and remove first 4 digits
        expid_start=$(echo $expid_range | cut -d'-' -f1 | sed 's/^[0-9]\{4\}//')
        expid_end=$(echo $expid_range | cut -d'-' -f2 | sed 's/^[0-9]\{4\}//')

        print_elapsed_time "Submitting Back2Back Flats job for $tele $mjd"
        sbatchCustom --job-name=b2b_${tele}_${mjd} $base_dir/scripts/daily/run_all.sh ${tele} ${mjd} true ${data_dir}

        # Store job ID and associated information
        job_ids+=($SLURM_ID)
        job_info+=("$tele|$mjd|$expid_start|$expid_end")
    fi
done < "metadata/special_cal_obs.txt"

#Wait for all jobs to complete
print_elapsed_time "Waiting for all reduction jobs to complete"
for job_id in "${job_ids[@]}"; do
    # Check if job exists and starts with b2b
    squeue -h -j $job_id -n "b2b*" > /dev/null 2>&1
    while [ $? -eq 0 ]; do
        sleep 1
        squeue -h -j $job_id -n "b2b*" > /dev/null 2>&1
    done
done

# Run all back2back_flats.jl scripts
print_elapsed_time "Running all Back2Back plots"
for i in "${!job_info[@]}"; do
    # Extract stored information
    IFS='|' read -r tele mjd expid_start expid_end <<< "${job_info[$i]}"

    print_elapsed_time "Making Back2Back Plots $tele $mjd"
    julia +$julia_version --project=$base_dir $base_dir/scripts/verify/back2back_flats.jl \
        --tele $tele \
        --mjd $mjd \
        --expid-start $expid_start \
        --expid-end $expid_end \
        --data_dir $data_dir
done

# Clean up logs and Report Timing
print_elapsed_time "Job Completed"
