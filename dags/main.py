import os
import re
from datetime import datetime, timedelta
from astropy.time import Time
from airflow import DAG
from airflow.operators.bash import BashOperator
from airflow.operators.python import PythonOperator, BranchPythonOperator
from airflow.utils.task_group import TaskGroup
from airflow.exceptions import AirflowSkipException

from slurm_utils import SlurmJobStatusSensor

SLURM_CLUSTER = nk, *_ = "notchpeak"
REPO_DIR = "$MWM_SANDBOX/ApogeeReduction.jl"
REPO_BRANCH = "main"

OUTPUT_DIR_TEMPLATE = (
    "$MWM_SANDBOX/"
    "outdir/"
    "{{ ti.xcom_pull(task_ids='setup.mjd_start') }}_{{ ti.xcom_pull(task_ids='setup.mjd_end') }}/"
)

ALMANAC_PATH = f"{OUTPUT_DIR_TEMPLATE}/almanac/all_{{{{ ti.xcom_pull(task_ids='setup.mjd_start') }}}}_{{{{ ti.xcom_pull(task_ids='setup.mjd_end') }}}}.h5"

sbatch_kwargs = re.sub(r"\s+", " ", 
    f"""
    -vvv
    --job-name=DreamTeam-{{{{ ds }}}}
    --account=sdss-{nk}p 
    --partition=sdss-{nk}p
    --clusters={SLURM_CLUSTER}
    --nodes=2
    -o {OUTPUT_DIR_TEMPLATE}/logs/vmstat_{{{{ ds }}}}_{{{{ ti.run_id }}}}.out
    -e {OUTPUT_DIR_TEMPLATE}/logs/vmstat_{{{{ ds }}}}-{{{{ ti.run_id }}}}.err
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
    --export=ALL,SLURM_TASKS_PER_NODE=64\(x2\),SLURM_JOB_QOS=sdss-{nk}p,SLURM_JOB_ACCOUNT=sdss-{nk}k,SLURM_NPROCS=128,SLURM_STEP_NUM_TASKS=128,SLURM_NTASKS=32,SLURM_NNODES=2,SLURM_JOB_NUM_NODES=2,SLURM_STEP_NUM_NODES=2,SLURM_CLUSTERS={SLURM_CLUSTER},SLURM_CLUSTER={SLURM_CLUSTER}.peaks
    --mem=250G
""") 


def observatories_with_data(almanac_path, **kwargs):
    import h5py as h5
    with h5.File(os.path.expandvars(almanac_path)) as fp:
        if len(fp.keys()) == 0:
            raise AirflowSkipException(f"No data in {almanac_path}")
        return list(fp.keys())

with DAG(
    "ApogeeReduction-main", 
    #start_date=datetime(2014, 7, 18), 
    start_date=datetime(2023, 9, 1),
    schedule_interval="@monthly",
    max_active_runs=1,
    default_args=dict(
        retries=3,
    ),
    catchup=True,
) as dag:
    with TaskGroup(group_id="update") as group_git:
        (
            BashOperator(
                task_id="repo",
                bash_command=(
                    f"cd {REPO_DIR}; "
                    f"git checkout {REPO_BRANCH}; "
                    "git pull"
                )
            )
        ) >> (
            BashOperator(
                task_id="dag",
                bash_command="airflow dags report",
            )
        )

    with TaskGroup(group_id="setup") as group_setup:
        [
            PythonOperator(
                task_id="mjd_start",
                python_callable=lambda data_interval_start, **_: int(Time(data_interval_start).mjd)
            ),
            #PythonOperator(
            #    task_id="mjd_end",
            #    python_callable=lambda data_interval_end, **_: int(Time(data_interval_end).mjd)
            #)
            PythonOperator(
                task_id="mjd_end",
                python_callable=lambda data_interval_start, **_: int(Time(data_interval_start).mjd) # one night at a time
            )
        ] >> (
            BashOperator(
                task_id="mkdirs",
                bash_command=(
                    f"mkdir -p {OUTPUT_DIR_TEMPLATE};"
                    f"mkdir -p {OUTPUT_DIR_TEMPLATE}/almanac;"
                    f"mkdir -p {OUTPUT_DIR_TEMPLATE}/logs;"
                    f"mkdir -p {OUTPUT_DIR_TEMPLATE}/darks;"
                    f"mkdir -p {OUTPUT_DIR_TEMPLATE}/flats;"
                )
            )        
        ) >> (
            BashOperator(
                task_id="almanac",
                bash_command=(
                    "module load almanac sdssdb/main sdsstools postgresql; "
                    "almanac -v -p 12 --apo --lco "
                    "--mjd-start {{ ti.xcom_pull(task_ids='setup.mjd_start') }} "
                    "--mjd-end {{ ti.xcom_pull(task_ids='setup.mjd_end') }} "
                    f"--output {ALMANAC_PATH}"
                )
            )
        ) >> (            
            PythonOperator(
                task_id="any_data",
                python_callable=observatories_with_data,
                op_kwargs=dict(almanac_path=f"{ALMANAC_PATH}")
            )
        )

    with TaskGroup(group_id="slurm") as group_slurm:
        (
            BashOperator(
                task_id="submit", 
                bash_command=f"sbatch {sbatch_kwargs} --wrap=\"vmstat -t 1\"",
                output_processor=lambda o: re.match("Submitted batch job (\d+)", o).group(1)
            )
        ) >> (
            SlurmJobStatusSensor(
                task_id="wait", 
                job_id="{{ ti.xcom_pull(task_ids='slurm.submit') }}",
            )
        )

    branch = BranchPythonOperator(
        task_id="branch",
        python_callable=observatories_with_data,
        op_kwargs=dict(almanac_path=f"{ALMANAC_PATH}")
    )

    observatories = []
    for observatory in ("apo", "lco"):
        with TaskGroup(group_id=observatory) as group:                
            with TaskGroup(group_id="calibrations") as g:
                BashOperator(
                    task_id="3d",
                    bash_command=f"{srun_prefix} src/cal_build/airflow_run_all.sh {observatory} {{{{ ti.xcom_pull(task_ids='setup.mjd_start') }}}} {{{{ ti.xcom_pull(task_ids='setup.mjd_end') }}}} {OUTPUT_DIR_TEMPLATE} {ALMANAC_PATH}",
                )
            
        observatories.append(g)

    scancel = BashOperator(
        task_id="teardown",
        bash_command=f"scancel -M {SLURM_CLUSTER} {{{{ ti.xcom_pull(task_ids='slurm.submit') }}}}",
        trigger_rule="all_done",
    )

    group_git >> group_setup >> group_slurm >> branch >> observatories >> scancel
    
