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

REPO_DIR = "/mnt/home/sdssv/gitcode/ApogeeReduction.jl"
REPO_BRANCH = "airflow"
OUT_DIR = "/mnt/ceph/users/sdssv/work/daily/"

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
    # Set environment variable for the subprocess
    env = os.environ.copy()
    env["SLACK_CHANNEL"] = "C08B7FKMP16" # apogee-reduction-jl
    
    # Now bash_command comes directly from the arguments
    print(f"Submitting command: {bash_command}")
    result = subprocess.run(bash_command.split(), capture_output=True, text=True, env=env)
    
    if result.returncode != 0:
        raise Exception(f"SLURM submission failed: {result.stderr}")
        
    job_id = result.stdout.strip().split()[-1]  # Get job ID from "Submitted batch job XXXXX"
    print(f"Submitted SLURM job {job_id}")
    
    wait_for_slurm(job_id)
    return job_id

class TransferFileSensor(FileSensor):
    def poke(self, context) -> bool:
        if Time(context["data_interval_start"]).mjd < 59148:
            raise AirflowSkipException("Data interval range is ancient, assuming transfer is complete.")
        return super().poke(context)
        

observatories = ("apo", "lco")

# Precompute yesterday's date for use in notifications
yesterday_date = "{{ (ds | string | to_datetime - timedelta(days=1)).strftime('%Y-%m-%d') }}"

# mails to Andrew Saydjari (can turn off by commenting out the mail user)
sbatch_prefix = re.sub(r"\s+", " ", f"""
    sbatch 
    -vvv
    -D {OUT_DIR}
    --mail-type=ALL
    --mail-user=slurm_notify-aaaaq7zc7ou7enlnbsexygd3he@sdss5.slack.com
""") 

with DAG(
    "ApogeeReduction-airflow", 
    start_date=datetime(2024, 10, 10), # could be start of survey, but not needed bc of bulk
    schedule="0 7 * * *", # 7 am ET
    max_active_runs=2,
    default_args=dict(retries=0),
    catchup=False,
            on_failure_callback=[
            send_slack_notification_partial(
                text=f"ApogeeReduction-airflow DAG failed on {{{{ yesterday_date }}}}",
            )
        ]    
) as dag:

    with TaskGroup(group_id="update") as group_git:
        BashOperator(
            task_id="repo",
            bash_command=(
                f"cd {REPO_DIR}\n"
                "echo '=== Git Status ==='\n"
                "git status\n"
                "echo '=== Git Log (last 3 commits) ==='\n"
                "git log --oneline -3\n"
                "echo '=== Git Remote Status ==='\n"
                "git remote -v\n"
            ),
        )

    with TaskGroup(group_id="setup") as group_setup:
        sdsscore = BashOperator(
            task_id="sdsscore",
            bash_command=(
                "ORIG_PWD=$(pwd)\n"
                "cd /mnt/ceph/users/sdssv/raw/APOGEE/sdsscore/\n"
                "./update.sh\n" 
                "cd $ORIG_PWD"
            )
        )

        mjd = PythonOperator(
            task_id="mjd",
            python_callable=lambda data_interval_start, **_: int(Time(data_interval_start).mjd) -1 # +1 offset to get the most recent day
        )

        # Set order: sdsscore must complete before mjd
        sdsscore >> mjd

    observatory_groups = []
    for observatory in ["lco", "apo"]:  # Changed order to LCO first to allow serial case
        with TaskGroup(group_id=observatory) as group:
            initial_notification = PythonOperator(
                task_id="initial_notification",
                python_callable=lambda **_: None,  # Simple no-op function
                on_success_callback=[
                    send_slack_notification_partial(
                        text=f"Starting reduction for {observatory.upper()} for SJD {{{{ task_instance.xcom_pull(task_ids='setup.mjd') }}}} (night of {{{{ yesterday_date }}}})."
                    )
                ]
            )

            # could add in Globus CLI transfer triggered on public https://data.sdss5.org/sas/sdsswork/data/staging/apo/log/mos/60918/transfer-60918.done existence
            
            science = PythonOperator(
                task_id="science",
                python_callable=submit_and_wait,
                op_kwargs={
                    'bash_command': f"{sbatch_prefix} --job-name=ar_all_{observatory}_{{{{ ti.xcom_pull(task_ids='setup.mjd') }}}} /mnt/home/sdssv/gitcode/ApogeeReduction.jl/scripts/daily/run_all.sh {observatory} {{{{ ti.xcom_pull(task_ids='setup.mjd') }}}} {OUT_DIR} /mnt/home/sdssv/gitcode/arMADGICS.jl false"
                },
                on_success_callback=[
                    send_slack_notification_partial(
                        text=f"{observatory.upper()} science frames reduced for SJD {{{{ ti.xcom_pull(task_ids='setup.mjd') }}}} (night of {{{{ yesterday_date }}}}).",
                    )
                ],        
                on_failure_callback=[
                    send_slack_notification_partial(
                        text=f"{observatory.upper()} science frame reduction failed for SJD {{{{ ti.xcom_pull(task_ids='setup.mjd') }}}} (night of {{{{ yesterday_date }}}}). :picard_facepalm:",
                    )
                ]               
            )   
            initial_notification >> science
            
        observatory_groups.append(group)

    # Add final notification task
    final_notification = PythonOperator(
        task_id="completion_notification",
        python_callable=lambda **_: None,  # dummy function that does nothing
        on_success_callback=[
            send_slack_notification_partial(
                text=f"ApogeeReduction pipeline completed successfully for SJD {{{{ ti.xcom_pull(task_ids='setup.mjd') }}}} (night of {{{{ yesterday_date }}}}). Both observatories processed."
            )
        ],
        dag=dag
    )
    
    # Modify the final dependencies to chain the observatories sequentially
    group_git >> group_setup >> observatory_groups[0] >> observatory_groups[1] >> final_notification
    # group_git >> group_setup >> observatory_groups >> final_notification ## this is the original we want to use during parallel reductions
