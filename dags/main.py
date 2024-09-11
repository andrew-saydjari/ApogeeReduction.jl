import os
import re
from airflow import DAG
from airflow.operators.bash import BashOperator
from airflow.operators.python import PythonOperator
from airflow.exceptions import AirflowSkipException
from datetime import datetime, timedelta
from airflow.utils.task_group import TaskGroup

from airflow.sensors.base_sensor_operator import BaseSensorOperator
from subprocess import check_output, Popen, PIPE

slurm_partition = "sdss-np"
slurm_job_name_template = "DreamTeam-{{ ds }}"

SLURM_ALLOC_KWARGS = re.sub(r"\s+", " ", 
    f"""
    -vvv
    --job-name={slurm_job_name_template}
    --account={slurm_partition} 
    --partition={slurm_partition}
    --nodes=1
    --ntasks=1
    --mem-per-cpu=100M
    --cpus-per-task=1
    --time=01:00:00
    """.replace("\n", " ")
)

# TODO: Put this slurm stuff elsewhere, and make it not require --account or a full path
def get_slurm_queue():
    pattern = (
        "(?P<id>\d+)+\s+(?P<name>[-\w\d_\.]+)\s+(?P<user>[\w\d]+)\s+(?P<group>\w+)"
        "\s+(?P<account>[-\w]+)\s+(?P<partition>[-\w]+)\s+(?P<time_limit>[-\d\:]+)\s+"
        "(?P<time_left>[-\d\:]+)\s+(?P<status>\w*)\s+(?P<nodelist>[\w\d\(\)]+)"
    )
    process = Popen(
        [
            "/uufs/notchpeak.peaks/sys/installdir/slurm/std/bin/squeue",
            "--account=sdss-np,sdss-kp,notchpeak-gpu,sdss-np-fast",
            '--format="%14i %j %10u %10g %13a %13P %11l %11L %2t %R"',
        ],
        stdin=PIPE,
        stdout=PIPE,
        stderr=PIPE,
        universal_newlines=True,
    )
    output, error = process.communicate()
    yield from (match.groupdict() for match in re.finditer(pattern, output))
    

class SenseSlurmJobReadiness(BaseSensorOperator):

    template_fields = ("job_id", "job_name")

    def __init__(
        self,
        job_id: str | int = None,
        job_name: str = None,
        poke_interval: int = 60,
        exponential_backoff: bool = True,
        timeout: int = 86_400,
        soft_fail: bool = True,
        **kwargs
    ) -> None:
        super().__init__(
            mode="poke", 
            poke_interval=poke_interval,
            exponential_backoff=exponential_backoff,
            timeout=timeout,
            soft_fail=soft_fail,
            **kwargs
        )
        if job_id is None and job_name is None:
            raise ValueError("Either job_id or job_name must be provided")
        self.job_name = job_name
        self.job_id = job_id
        return None


    def poke(self, context) -> bool:
        for job in get_slurm_queue():
            if (
                (self.job_id is not None and f"{self.job_id}" == job["id"])
            or  (self.job_name is not None and self.job_name == job["name"])
            ):
                return (job["status"] == "R")
        # If we get to this point, the job is not even in the list, so we should fail or skip
        raise AirflowSkipException(f"Job not found in queue (id={self.job_id}, name={self.job_name})")


with DAG(
    "ApogeeReduction-main", 
    start_date=datetime(2023, 11, 13), 
    schedule_interval="@daily",
    max_active_runs=1,
    catchup=False,
) as dag:

    with TaskGroup(group_id="update") as group_git:
        p = BashOperator(
            task_id="git_pull",
            bash_command=f"cd $MWM_SANDBOX/ApogeeReduction.jl; git checkout main; git pull",
        )
        r = BashOperator(
            task_id="refresh_dag",
            bash_command="airflow dags report"
        )
        p >> r

    with TaskGroup(group_id="slurm") as group_slurm:

        salloc = BashOperator(task_id="allocate", bash_command=f"module load slurm/notchpeak; salloc {SLURM_ALLOC_KWARGS} /usr/bin/sleep 60 </dev/null &")
        wait = SenseSlurmJobReadiness(task_id="wait", job_name=slurm_job_name_template)
        salloc >> wait


    # Need to wait for the job to be allocated
    test_srun = BashOperator(
        task_id="test_srun",
        bash_command=f"srun --job-name={slurm_job_name_template} --mem=10mb bash -c \"echo 'Hello world ' `hostname`\"",
    )

    test_failure = BashOperator(
        task_id="test_failure",
        bash_command=f"srun --job-name={slurm_job_name_template} --mem=10mb python -c \"1/0\""
    )


    group_git >> group_slurm >> (test_srun, test_failure)
