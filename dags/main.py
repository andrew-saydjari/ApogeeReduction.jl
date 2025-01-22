import os
import re
import subprocess
import time
from datetime import datetime, timedelta
from astropy.time import Time
from airflow import DAG
from airflow.sensors.filesystem import FileSensor
from airflow.operators.bash import BashOperator
from airflow.operators.python import PythonOperator
from airflow.utils.task_group import TaskGroup
from airflow.exceptions import AirflowSkipException
from airflow.providers.slack.notifications.slack import send_slack_notification

REPO_DIR = "/uufs/chpc.utah.edu/common/home/sdss51/sdsswork/mwm/sandbox/airflow-ApogeeReduction.jl/ApogeeReduction.jl"
REPO_BRANCH = "airflow"

def send_slack_notification_partial(text):
    return send_slack_notification(text=text, channel="#apogee-reduction-jl")

# Add this function to check SLURM job status
def wait_for_slurm(job_id):
    while True:
        result = subprocess.run(['squeue', '-j', str(job_id)], capture_output=True, text=True)
        if "Invalid job id specified" in result.stderr or result.stdout.count('\n') <= 1:
            return  # Job is done
        time.sleep(5)  # Check every minute

# Modify your BashOperator to capture and wait for the job ID
def submit_and_wait(bash_command, **context):
    # Now bash_command comes directly from the arguments
    print(f"Submitting command: {bash_command}")
    result = subprocess.run(bash_command.split(), capture_output=True, text=True)
    
    if result.returncode != 0:
        raise Exception(f"SLURM submission failed: {result.stderr}")
        
    job_id = result.stdout.strip().split()[-1]  # Get job ID from "Submitted batch job XXXXX"
    print(f"Submitted SLURM job {job_id}")
    
    wait_for_slurm(job_id)
    return job_id

observatories = ("apo", "lco")

sbatch_prefix = re.sub(r"\s+", " ", f"""
    sbatch 
    -vvv
    -D {REPO_DIR}
""") 
# text=f"ApogeeReduction-main DAG failed on {{{{ ds }}}}: {DAG_URL}",
with DAG(
    "ApogeeReduction-main", 
    start_date=datetime(2024, 10, 10), # datetime(2014, 7, 18), 
    schedule="0 5 * * *", # 8 am ET
    max_active_runs=2,
    default_args=dict(retries=0),
    catchup=False,
    on_failure_callback=[
        send_slack_notification_partial(
            text=f"ApogeeReduction-main DAG failed on {{{{ ds }}}}",
        )
    ]    
) as dag:

    with TaskGroup(group_id="update") as group_git:
        (
            BashOperator(
                task_id="repo",
                bash_command=(
                    # assumes you have run 'git checkout airflow' to be on the airflow branch
                    f"cd {REPO_DIR}; "
                    "git add -A; "  # Stage all changes, including deletions
                    "git commit -m 'Auto-commit local changes'; "  # Commit changes with a message
                    f"git push origin {REPO_BRANCH}; "  # Push local changes
                    # create a PR against main with these local changes ('|| true' prevents failure if PR already exists)
                    f"gh pr create --title 'Automated updates from airflow pipeline' --body 'This PR was automatically created by the airflow pipeline.' --base main --head {REPO_BRANCH} || true; "
                    # auto-merge the PR
                    "gh pr merge --auto --merge --delete-branch=false; "
                    # get main and use it to merge into current branch
                    "git fetch origin main; "  # Get latest main
                    "git merge origin/main --no-edit; "  # Merge main into current branch
                    f"git pull origin {REPO_BRANCH}"  # Pull latest changes
                ),
                # Add error handling for git commands
                bash_command_opts="-e",  # Exit on error
            )
        ) >> (
            BashOperator(
                task_id="dag",
                bash_command="airflow dags report",
            )
        )

    with TaskGroup(group_id="setup") as group_setup:

        BashOperator(
            task_id="julia",
            bash_command=(
                f'cd {REPO_DIR}; '
                'juliaup add +1.11.0; '
                'julia +1.11.0 --project="./" -e \''
                    'using Pkg; '
                    'Pkg.add(url = "https://github.com/andrew-saydjari/SlackThreads.jl.git"); '
                    'Pkg.add(url = "https://github.com/nasa/SIRS.git"); '
                    'Pkg.resolve(); '
                    'Pkg.instantiate(); '
                '\''
            )
        )
        PythonOperator(
            task_id="mjd",
            python_callable=lambda data_interval_start, **_: int(Time(data_interval_start).mjd) + 1 # +1 offset to get the most recent day
        )

    observatory_groups = []
    for observatory in observatories:
        with TaskGroup(group_id=observatory) as group:
            transfer = FileSensor(
                task_id="transfer",
                filepath=f"/uufs/chpc.utah.edu/common/home/sdss50/sdsswork/data/staging/{observatory}/log/mos/{{{{ ti.xcom_pull(task_ids='setup.mjd') }}}}/transfer-{{{{ ti.xcom_pull(task_ids='setup.mjd') }}}}.done",
                on_success_callback=[
                    send_slack_notification_partial(
                        text=f"{observatory.upper()} data transfer complete for SJD {{{{ ti.xcom_pull(task_ids='setup.mjd') }}}} (night of {{{{ macros.ds_add(ds, -1) }}}})."
                             f"Transfer file last modified at {{{{ macros.datetime.fromtimestamp(os.path.getmtime(filepath)).astimezone(macros.pendulum.timezone('America/New_York')).strftime('%Y-%m-%d %I:%M %p %Z') }}}}."
                             f"Starting reduction pipeline.",
                    )
                ],
                on_failure_callback=[
                    send_slack_notification_partial(
                        text=f"{observatory.upper()} data transfer on SJD {{{{ ti.xcom_pull(task_ids='setup.mjd') }}}} (night of {{{{ macros.ds_add(ds, -1) }}}}) is incomplete. Please check https://data.sdss5.org/sas/sdsswork/data/staging/{observatory}/log/mos/ and investigate.",
                    )
                ],
                timeout=60*60*12, # 12 hours: ~8pm ET
                poke_interval=600, # 10 minutes
                mode="poke",
            )
            
            # Could set this up to take in directories for cleaner airflow running
            darks = PythonOperator(
                task_id="darks",
                python_callable=submit_and_wait,
                op_kwargs={'bash_command': f"{sbatch_prefix} src/cal_build/run_dark_cal.sh {observatory} {{{{ ti.xcom_pull(task_ids='setup.mjd') }}}} {{{{ ti.xcom_pull(task_ids='setup.mjd') }}}}"},
            )

            flats = PythonOperator(
                task_id="flats",
                python_callable=submit_and_wait,
                op_kwargs={'bash_command': f"{sbatch_prefix} src/cal_build/run_flat_cal.sh {observatory} {{{{ ti.xcom_pull(task_ids='setup.mjd') }}}} {{{{ ti.xcom_pull(task_ids='setup.mjd') }}}}"},
            )
            
            science = PythonOperator(
                task_id="science",
                python_callable=submit_and_wait,
                op_kwargs={
                    'bash_command': f"{sbatch_prefix} src/run_scripts/run_all.sh {observatory} {{{{ ti.xcom_pull(task_ids='setup.mjd') }}}}"
                },
                on_success_callback=[
                    send_slack_notification_partial(
                        text=f"{observatory.upper()} science frames reduced for SJD {{{{ ti.xcom_pull(task_ids='setup.mjd') }}}} (night of {{{{ macros.ds_add(ds, -1) }}}}).",
                    )
                ],        
                on_failure_callback=[
                    send_slack_notification_partial(
                        text=f"{observatory.upper()} science frame reduction failed for SJD {{{{ ti.xcom_pull(task_ids='setup.mjd') }}}} (night of {{{{ macros.ds_add(ds, -1) }}}}). :picard_facepalm:",
                    )
                ]               
            )   
            transfer >> darks >> flats >> science
        
        observatory_groups.append(group)

        # Add final notification task
    
    final_notification = PythonOperator(
        task_id="completion_notification",
        python_callable=lambda **context: send_slack_notification_partial(
            text=f"ApogeeReduction pipeline completed successfully for SJD {{{{ ti.xcom_pull(task_ids='setup.mjd') }}}} (night of {{{{ macros.ds_add(ds, -1) }}}}). Both observatories processed."
        )(**context),
        dag=dag
    )

    group_git >> group_setup >> observatory_groups >> final_notification
    


