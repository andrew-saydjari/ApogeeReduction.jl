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

# I struggled to move this to sandbox, because the fileloc would not update.
REPO_DIR = f"{os.path.expandvars('$MWM_SANDBOX')}/airflow/ApogeeReduction.jl/"
REPO_BRANCH = "airflow"
DAG_NAME = "ApogeeReduction-prod"

def send_slack_notification_partial(text):
    return send_slack_notification(text=f"[prod] {text}", channel="#apogee-reduction-jl")

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

sbatch_prefix = re.sub(r"\s+", " ", f"""
    sbatch 
    -vvv
    -D {REPO_DIR}
""") 


with DAG(
    DAG_NAME,
    start_date=datetime(2014, 7, 18), 
    schedule_interval=timedelta(days=8),
    max_active_runs=1,
    default_args=dict(retries=0),
    catchup=True,
    on_failure_callback=[
        send_slack_notification_partial(
            text=f"{DAG_NAME} DAG failed on {{{{ ds }}}}",
        )
    ]    
) as dag:

    with TaskGroup(group_id="update") as group_git:
        (
            BashOperator(
                task_id="repo",
                bash_command=(
                    f"cd {REPO_DIR}\n"
                    "set -e\n"
                    "git add -A\n"
                    # Check if there are changes to commit, if there are PR and merge
                    'if [[ -n "$(git status --porcelain)" ]]; then\n'
                    "    git commit -m 'Auto-commit local changes'\n"
                    f"    git push origin {REPO_BRANCH}\n"
                    f'    gh pr create --title "Automated updates" --body "Auto-created by airflow" --base main --head {REPO_BRANCH}\n'
                    # '    PR_NUM=$(echo "$PR_OUT" | grep -o \'/[0-9]*$\' | tr -d \'/\')\n'
                    '    gh pr merge "$PR_NUM" --admin --merge --delete-branch=false\n'
                    '    echo "Merged PR #$PR_NUM"\n'
                    "    git fetch origin main\n"
                    "    git merge origin/main --no-edit\n"
                    f"    git pull origin {REPO_BRANCH}\n"
                    "else\n"
                    '    echo "No changes to commit"\n'
                    '    git fetch origin main\n'
                    '    git merge origin/main --no-edit\n'
                    "fi\n"
                    # --- Separate check for unpushed commits (ahead of origin) ---
                    'read behind ahead < <(git rev-list --left-right --count origin/' + f'{REPO_BRANCH}' + '...HEAD)\n'
                    'if [[ "$ahead" -gt 0 ]]; then\n'
                    f'    echo "Branch is ahead of origin/{REPO_BRANCH} by $ahead commits. Pushing..."\n'
                    f'    git push origin HEAD:{REPO_BRANCH}\n'
                    "else\n"
                    '    echo "Branch is up to date with origin/{REPO_BRANCH}"\n'
                    "fi\n"
                ),
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
                'juliaup add 1.11.0; '
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

    # Replace the observatory_groups section with this new implementation
    # Process LCO first, then APO (temporaily, while we are still doing serial reductinos)
    observatory_groups = []
    for observatory in ["lco", "apo"]:  # Changed order to LCO first
        with TaskGroup(group_id=observatory) as group:
            initial_notification = PythonOperator(
                task_id="initial_notification",
                python_callable=lambda **_: None,  # Simple no-op function
                on_success_callback=[
                    send_slack_notification_partial(
                        text=f"Waiting for {observatory.upper()} data transfer for SJD {{{{ task_instance.xcom_pull(task_ids='setup.mjd') }}}} "
                            f"(night of {{{{ ds }}}}). "
                            f"Check here for transfer status: https://data.sdss5.org/sas/sdsswork/data/staging/{observatory}/log/mos/"
                    )
                ]
            )

            filepath = f"/uufs/chpc.utah.edu/common/home/sdss50/sdsswork/data/staging/{observatory}/log/mos/{{{{ task_instance.xcom_pull(task_ids='setup.mjd') }}}}/transfer-{{{{ task_instance.xcom_pull(task_ids='setup.mjd') }}}}.done"
            transfer = TransferFileSensor(
                task_id="transfer",
                filepath=filepath,
                on_success_callback=[
                    send_slack_notification_partial(
                        text=f"{observatory.upper()} data transfer complete for SJD {{{{ task_instance.xcom_pull(task_ids='setup.mjd') }}}} "
                             f"(night of {{{{ ds }}}}). Starting reduction pipeline. Exposure list available at https://data.sdss5.org/sas/sdsswork/data/apogee/{observatory}/{{{{ task_instance.xcom_pull(task_ids='setup.mjd') }}}}/{{{{ task_instance.xcom_pull(task_ids='setup.mjd') }}}}.log.html"
                    )
                ],
                on_failure_callback=[
                    send_slack_notification_partial(
                        text=f"{observatory.upper()} data transfer on SJD {{{{ task_instance.xcom_pull(task_ids='setup.mjd') }}}} "
                             f"(night of {{{{ ds }}}}) is incomplete. "
                             f"Please check https://data.sdss5.org/sas/sdsswork/data/staging/{observatory}/log/mos/ and investigate."
                    )
                ],
                on_skipped_callback=[
                    send_slack_notification_partial(
                        text=f"{observatory.upper()} data on SJD {{{{ task_instance.xcom_pull(task_ids='setup.mjd') }}}} "
                             f"(night of {{{{ ds }}}}) assumed to already exist: not awaiting `done` file. Starting reduction pipeline."
                    )
                ],
                timeout=60*60*18, # 18 hours: midnight ET
                poke_interval=600, # 10 minutes
                mode="poke",
            )
            
            # Could set this up to take in directories for cleaner airflow running
            darks = PythonOperator(
                task_id="darks",
                python_callable=submit_and_wait,
                op_kwargs={'bash_command': f"{sbatch_prefix} --job-name=ar_dark_cal_{observatory}_{{{{ ti.xcom_pull(task_ids='setup.mjd') }}}} src/cal_build/run_dark_cal.sh {observatory} {{{{ ti.xcom_pull(task_ids='setup.mjd') }}}} {{{{ ti.xcom_pull(task_ids='setup.mjd') }}}}"},
                trigger_rule="none_failed" # so this doesn't get skipped if the FileSensor is skipped
            )

            flats = PythonOperator(
                task_id="flats",
                python_callable=submit_and_wait,
                op_kwargs={'bash_command': f"{sbatch_prefix} --job-name=ar_flat_cal_{observatory}_{{{{ ti.xcom_pull(task_ids='setup.mjd') }}}} src/cal_build/run_flat_cal.sh {observatory} {{{{ ti.xcom_pull(task_ids='setup.mjd') }}}} {{{{ ti.xcom_pull(task_ids='setup.mjd') }}}}"},
            )
            
            science = PythonOperator(
                task_id="science",
                python_callable=submit_and_wait,
                op_kwargs={
                    'bash_command': f"{sbatch_prefix} --job-name=ar_all_{observatory}_{{{{ ti.xcom_pull(task_ids='setup.mjd') }}}} src/run_scripts/run_all.sh {observatory} {{{{ ti.xcom_pull(task_ids='setup.mjd') }}}}"
                },
                on_success_callback=[
                    send_slack_notification_partial(
                        text=f"{observatory.upper()} science frames reduced for SJD {{{{ ti.xcom_pull(task_ids='setup.mjd') }}}} (night of {{{{ ds }}}}).",
                    )
                ],        
                on_failure_callback=[
                    send_slack_notification_partial(
                        text=f"{observatory.upper()} science frame reduction failed for SJD {{{{ ti.xcom_pull(task_ids='setup.mjd') }}}} (night of {{{{ ds }}}}). :picard_facepalm:",
                    )
                ]               
            )   
            initial_notification >> transfer >> darks >> flats >> science
            
        observatory_groups.append(group)

    # Add final notification task
    final_notification = PythonOperator(
        task_id="completion_notification",
        python_callable=lambda **_: None,  # dummy function that does nothing
        on_success_callback=[
            send_slack_notification_partial(
                text=f"ApogeeReduction pipeline completed successfully for SJD {{{{ ti.xcom_pull(task_ids='setup.mjd') }}}} (night of {{{{ ds }}}}). Both observatories processed."
            )
        ],
        dag=dag
    )
    
    # Modify the final dependencies to chain the observatories sequentially
    group_git >> group_setup >> observatory_groups[0] >> observatory_groups[1] >> final_notification
    # group_git >> group_setup >> observatory_groups >> final_notification ## this is the original we want to use during parallel reductions
