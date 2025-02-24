#!/bin/bash
# ------------------------------------------------------------------------------

# load all of the modules to talk to the database (need to be on Utah)
# should turn this off as an option for users once the MJD summaries are generated
module load sdssdb/main almanac sdsstools postgresql ffmpeg
if [ -n "$SLURM_JOB_NODELIST" ]; then
    echo $SLURM_JOB_NODELIST
else
    hostname
fi
juliaup add 1.11.0

# Define sbatchCustom function
function sbatchCustom() {
    export SLURM_ID=$(sbatch --parsable "$@")
    mkdir -p slurm_logs
    scontrol write batch_script $SLURM_ID slurm_logs/slurm-$SLURM_ID.sh
}

data_dir=$(realpath "../outdir/")

# set up the output directory (if does not exist)
mkdir -p ${data_dir}/plots

# nice function to print the elapsed time at each step
print_elapsed_time() {
    formatted_time=$(printf '%dd %dh:%dm:%ds\n' $(($SECONDS/86400)) $(($SECONDS%86400/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60)))
    echo
    echo "Elapsed time: $formatted_time"
    echo
    echo "--------- $1 ---------"
    echo
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
     
        # print_elapsed_time "Submitting Back2Back Flats job for $tele $mjd"
        sbatchCustom --job-name=b2b_${tele}_${mjd} ./src/run_scripts/run_all.sh ${tele} ${mjd} true ${data_dir}
        
        # Store job ID and associated information
        job_ids+=($SLURM_ID)
        job_info+=("$tele|$mjd|$expid_start|$expid_end")
    fi
done < "metadata/special_cal_obs.txt"

# Wait for all jobs to complete
print_elapsed_time "Waiting for all reduction jobs to complete"
for job_id in "${job_ids[@]}"; do
    squeue -h -j $job_id > /dev/null 2>&1
    while [ $? -eq 0 ]; do
        sleep 1
        squeue -h -j $job_id > /dev/null 2>&1
    done
done

# Run all back2back_flats.jl scripts
print_elapsed_time "Running all Back2Back plots"
for i in "${!job_info[@]}"; do
    # Extract stored information
    IFS='|' read -r tele mjd expid_start expid_end <<< "${job_info[$i]}"
    
    print_elapsed_time "Making Back2Back Plots $tele $mjd"
    julia +1.11.0 --project="./" src/verify_scripts/back2back_flats.jl \
        --tele $tele \
        --mjd $mjd \
        --expid-start $expid_start \
        --expid-end $expid_end \
        --data_dir $data_dir
done

# Clean up logs and Report Timing
print_elapsed_time "Job Completed"