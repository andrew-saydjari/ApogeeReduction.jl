import os
import re
from datetime import datetime, timedelta
from astropy.time import Time
from airflow import DAG
from airflow.sensors.filesystem import FileSensor
from airflow.operators.bash import BashOperator
from airflow.operators.python import PythonOperator
from airflow.utils.task_group import TaskGroup
from airflow.exceptions import AirflowSkipException
from airflow.providers.slack.notifications.slack import send_slack_notification

SLURM_CLUSTER = nk, *_ = "notchpeak"
REPO_DIR = "/uufs/chpc.utah.edu/common/home/u6039752/scratch1/working/2025_01_21/ApogeeReduction.jl"
REPO_BRANCH = "2025_01_21"

def send_slack_notification_partial(text):
    return send_slack_notification(text=text, channel="#apogee-reduction-jl")

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
    schedule="0 8 * * *", # 8 am ET
    max_active_runs=1,
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
                    f"cd {REPO_DIR}; "
                    f"git checkout {REPO_BRANCH}; "
                    "git add -A; "  # Stage all changes, including deletions
                    "git commit -m 'Auto-commit local changes'; "  # Commit changes with a message
                    "git push; "  # Push local changes
                    "git pull"  # Pull latest changes
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
                        text=f"{observatory.upper()} data transfer complete for SJD {{{{ ti.xcom_pull(task_ids='setup.mjd') }}}} ({{{{ ds }}}}). Starting reduction pipeline.",
                    )
                ],
                on_failure_callback=[
                    send_slack_notification_partial(
                        text=f"{observatory.upper()} data transfer on SJD {{{{ ti.xcom_pull(task_ids='setup.mjd') }}}} ({{{{ ds }}}}) is incomplete. Please investigate.",
                    )
                ],
                timeout=60*60*12, # 12 hours: ~8pm ET
                poke_interval=600, # 10 minutes
                mode="poke",
            )
            
            # Darks and Flats and Science are not taking in directories and almanac paths as arguments {OUTPUT_DIR_TEMPLATE} {ALMANAC_PATH}
            darks = BashOperator(
                task_id="darks",
                bash_command=f"{sbatch_prefix} src/cal_build/run_dark_cal.sh {observatory} {{{{ ti.xcom_pull(task_ids='setup.mjd') }}}} {{{{ ti.xcom_pull(task_ids='setup.mjd') }}}}",
            )

            flats = BashOperator(
                task_id="flats",
                bash_command=f"{sbatch_prefix} src/cal_build/run_flat_cal.sh {observatory} {{{{ ti.xcom_pull(task_ids='setup.mjd') }}}} {{{{ ti.xcom_pull(task_ids='setup.mjd') }}}}",
            )
            
            science = BashOperator(
                task_id="science",
                bash_command=f"{sbatch_prefix} src/run_scripts/run_all.sh {observatory} {{{{ ti.xcom_pull(task_ids='setup.mjd') }}}}",
                on_success_callback=[
                    send_slack_notification_partial(
                        text=f"{observatory.upper()} science frames reduced for SJD {{{{ ti.xcom_pull(task_ids='setup.mjd') }}}} ({{{{ ds }}}}).",
                    )
                ],        
                on_failure_callback=[
                    send_slack_notification_partial(
                        text=f"{observatory.upper()} science frame reduction failed for SJD {{{{ ti.xcom_pull(task_ids='setup.mjd') }}}} ({{{{ ds }}}}). :picard_facepalm:",
                    )
                ]               
            )        
            transfer >> darks >> flats >> science
        
        observatory_groups.append(group)
            
    group_git >> group_setup >> observatory_groups
    


