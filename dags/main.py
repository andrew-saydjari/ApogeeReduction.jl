import os
import re
from datetime import datetime, timedelta
from astropy.time import Time
from airflow import DAG
from airflow.operators.bash import BashOperator
from airflow.operators.python import PythonOperator
from airflow.utils.task_group import TaskGroup

from slurm_utils import SlurmJobStatusSensor

SLURM_CLUSTER = nk, *_ = "notchpeak"
REPO_DIR = "$MWM_SANDBOX/ApogeeReduction.jl"

sbatch_kwargs = re.sub(r"\s+", " ", 
    f"""
    -vvv
    --job-name=DreamTeam-{{{{ ds}}}}
    --account=sdss-{nk}p 
    --partition=sdss-{nk}p
    --clusters={SLURM_CLUSTER}
    --nodes=2
    --ntasks=128
    --exclusive
    --mem=0
    --time=01:00:00
    """.replace("\n", " ")
)
srun_prefix = re.sub(r"\s+", " ", f"""
    srun 
    -v
    --external-launcher
    --ntasks=1
    --ntasks-per-node=16
    --jobid={{{{ ti.xcom_pull(task_ids='slurm.submit') }}}} 
    -M {SLURM_CLUSTER}
    -D {REPO_DIR}
    --mem=400G
""") 
#         
#    --cpus-per-task=8
#    --ntasks-per-node=16




with DAG(
    "ApogeeReduction-main", 
    start_date=datetime(2021, 11, 30), 
    end_date=datetime(2021, 12, 1),
    schedule_interval="@daily",
    max_active_runs=1,
    catchup=True,
) as dag:
    with TaskGroup(group_id="update") as group_git:
        (
            BashOperator(
                task_id="repo",
                bash_command=f"cd {REPO_DIR}; git checkout main; git pull",
            )
        ) >> (
            BashOperator(
                task_id="dag",
                bash_command="airflow dags report",
            )
        )

    with TaskGroup(group_id="slurm") as group_slurm:
        (
            BashOperator(
                task_id="submit", 
                bash_command=f"sbatch {sbatch_kwargs} --wrap=\"printenv; sleep 365d & wait\"",
                output_processor=lambda o: re.match("Submitted batch job (\d+)", o).group(1)
            )
        ) >> (
            SlurmJobStatusSensor(
                task_id="wait", 
                job_id="{{ ti.xcom_pull(task_ids='slurm.submit') }}",
            )
        )

    with TaskGroup(group_id="mjd") as group_mjd:
        (
            PythonOperator(
                task_id="start",
                python_callable=lambda data_interval_start, **_: int(Time(data_interval_start).mjd)
            )
        ) >> (
            PythonOperator(
                task_id="end",
                python_callable=lambda data_interval_end, **_: int(Time(data_interval_end).mjd)
            )
        )


    group_observatories = []
    for observatory in ("apo", "lco"):
        with TaskGroup(group_id=observatory) as group:
            with TaskGroup(group_id="calibrations") as group_calibrations:
                darks = BashOperator(
                    task_id="darks",
                    bash_command=f"{srun_prefix} src/cal_build/run_dark_cal.sh {observatory} {{{{ ti.xcom_pull(task_ids='mjd.start') }}}} {{{{ ti.xcom_pull(task_ids='mjd.end') }}}}"
                )
        group_observatories.append(group)

    scancel = BashOperator(
        task_id="teardown",
        bash_command=f"scancel -M {SLURM_CLUSTER} {{{{ ti.xcom_pull(task_ids='slurm.submit') }}}}",
        trigger_rule="all_done",
    )

    group_mjd >> group_observatories
    group_git >> group_slurm >> group_observatories >> scancel
    
