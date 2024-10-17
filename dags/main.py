import os
import re
from datetime import datetime, timedelta
from astropy.time import Time
from airflow import DAG
from airflow.sensors.filesystem import FileSensor
from airflow.operators.bash import BashOperator
from airflow.operators.python import PythonOperator, BranchPythonOperator
from airflow.utils.task_group import TaskGroup
from airflow.exceptions import AirflowSkipException
from airflow.providers.slack.notifications.slack import send_slack_notification

from slurm_utils import SlurmJobStatusSensor, plot_vmstats

SLURM_CLUSTER = nk, *_ = "notchpeak"
REPO_DIR = "$MWM_SANDBOX/ApogeeReduction.jl"
REPO_BRANCH = "main"
DAG_URL = "http://astra-vm.sdss.utah.edu/dags/ApogeeReduction-main/grid"


OUTPUT_DIR_TEMPLATE = (
    "$MWM_SANDBOX/"
    "outdir/"
    "{{ ti.xcom_pull(task_ids='setup.mjd') }}_{{ ti.xcom_pull(task_ids='setup.mjd') }}/"
)

ALMANAC_PATH = f"{OUTPUT_DIR_TEMPLATE}/almanac/all_{{{{ ti.xcom_pull(task_ids='setup.mjd') }}}}_{{{{ ti.xcom_pull(task_ids='setup.mjd') }}}}.h5"
def send_slack_notification_partial(text):
    return send_slack_notification(text=text, channel="#apogee-reduction-jl")

observatories = ("apo", "lco")

sbatch_kwargs = re.sub(r"\s+", " ", 
    f"""
    -vvv
    --account=sdss-{nk}p 
    --partition=sdss-{nk}p
    --clusters={SLURM_CLUSTER}
    --nodes=2
    --exclusive
    --mem=0
    --time=01:00:00
    -o {OUTPUT_DIR_TEMPLATE}/logs/vmstat_{{{{ ds }}}}_{{{{ ti.run_id }}}}.out
    -e {OUTPUT_DIR_TEMPLATE}/logs/vmstat_{{{{ ds }}}}_{{{{ ti.run_id }}}}.err
    """.replace("\n", " ")
)
#    --ntasks=128


srun_prefix = re.sub(r"\s+", " ", f"""
    srun 
    -v
    --external-launcher
    --ntasks=1
    -M {SLURM_CLUSTER}
    -D {REPO_DIR}
    --export=ALL,SLURM_TASKS_PER_NODE=64\(x2\),SLURM_JOB_QOS=sdss-{nk}p,SLURM_JOB_ACCOUNT=sdss-{nk}k,SLURM_NPROCS=128,SLURM_STEP_NUM_TASKS=128,SLURM_NTASKS=32,SLURM_NNODES=2,SLURM_JOB_NUM_NODES=2,SLURM_STEP_NUM_NODES=2,SLURM_CLUSTERS={SLURM_CLUSTER},SLURM_CLUSTER={SLURM_CLUSTER}.peaks
""") 
#    --ntasks-per-node=64


def observatories_with_data(almanac_path, **kwargs):
    import h5py as h5
    with h5.File(os.path.expandvars(almanac_path)) as fp:
        if len(fp.keys()) == 0:
            raise AirflowSkipException(f"No data in {almanac_path}")
        return list(fp.keys())



with DAG(
    "ApogeeReduction-main", 
    #start_date=datetime(2014, 7, 18), 
    start_date=datetime(2024, 10, 10),
    schedule="0 10 * * *", # 6 am Eastern
    max_active_runs=1,
    default_args=dict(retries=2),
    catchup=False,
    on_failure_callback=[
        send_slack_notification_partial(
            text=f"ApogeeReduction-main DAG failed on {{{{ ds }}}}: {DAG_URL}",
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
                    "git pull"
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
                'julia +1.11.0 --project="./" -e \'using Pkg; Pkg.add(url = "https://github.com/andrew-saydjari/SlackThreads.jl.git"); Pkg.add(url = "https://github.com/nasa/SIRS.git"); Pkg.resolve(); Pkg.instantiate()\''
            )
        )
        PythonOperator(
            task_id="mjd",
            python_callable=lambda data_interval_start, **_: int(Time(data_interval_start).mjd)
        ) >> (
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
                timeout=60*60*6, # 6 hours
            )
            almanac = (
                BashOperator(
                    task_id="almanac",
                    bash_command=(
                        "module load almanac sdssdb/main sdsstools postgresql; "
                        f"almanac -v -p 12 --{observatory} "
                        "--mjd-start {{ ti.xcom_pull(task_ids='setup.mjd') }} "
                        "--mjd-end {{ ti.xcom_pull(task_ids='setup.mjd') }} "
                        f"--output {ALMANAC_PATH}"
                    )
                )
            ) 
            
            with TaskGroup(group_id="slurm") as group_slurm:
                (
                    BashOperator(
                        task_id="submit", 
                        bash_command=f"sbatch --job-name=DreamTeam-{{{{ ds }}}}-{observatory} {sbatch_kwargs} --wrap=\"vmstat -t 1\"",
                        output_processor=lambda o: re.match("Submitted batch job (\d+)", o).group(1)
                    )
                ) >> (
                    SlurmJobStatusSensor(
                        task_id="wait", 
                        job_id=f"{{{{ ti.xcom_pull(task_ids='{observatory}.slurm.submit') }}}}",
                    )
                )

            pipeline = BashOperator(
                task_id="pipeline",
                bash_command=f"{srun_prefix} --jobid={{{{ ti.xcom_pull(task_ids='{observatory}.slurm.submit') }}}} src/cal_build/airflow_run_all.sh {observatory} {{{{ ti.xcom_pull(task_ids='setup.mjd') }}}} {{{{ ti.xcom_pull(task_ids='setup.mjd') }}}} {OUTPUT_DIR_TEMPLATE} {ALMANAC_PATH}",
                on_success_callback=[
                    send_slack_notification_partial(
                        text=f"{observatory.upper()} reductions complete for SJD {{{{ ti.xcom_pull(task_ids='setup.mjd') }}}} ({{{{ ds }}}}).",
                    )
                ],
                on_failure_callback=[
                    send_slack_notification_partial(
                        text=f"{observatory.upper()} reductions failed for SJD {{{{ ti.xcom_pull(task_ids='setup.mjd') }}}} ({{{{ ds }}}}). Please investigate: {DAG_URL}/",
                    )
                ]                    
            )

            scancel = BashOperator(
                task_id="teardown",
                bash_command=f"scancel -M {SLURM_CLUSTER} {{{{ ti.xcom_pull(task_ids='{observatory}.slurm.submit') }}}}",
                trigger_rule="all_done",
            )            
            transfer >> almanac >> group_slurm >> pipeline >> scancel
        
        observatory_groups.append(group)
            
    '''
    op_plot_vmstats = PythonOperator(
        task_id="plot_vmstats",
        python_callable=plot_vmstats,
        op_kwargs=dict(
            vmstat_path=f"{OUTPUT_DIR_TEMPLATE}/logs/vmstat_{{{{ ds }}}}_{{{{ ti.run_id }}}}.out",
            figure_path=f"{OUTPUT_DIR_TEMPLATE}/plots/vmstat_{{{{ ds }}}}_{{{{ ti.run_id }}}}.png",
        ),
    )
    '''
    group_git >> group_setup >> observatory_groups
    

    