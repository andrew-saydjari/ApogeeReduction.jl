import os
import re
import subprocess
import time
from datetime import datetime, timedelta
from astropy.time import Time
from airflow import DAG
from airflow.sensors.filesystem import FileSensor
from airflow.operators.bash import BashOperator
from airflow.operators.python import PythonOperator, BranchPythonOperator
from airflow.utils.task_group import TaskGroup
from airflow.exceptions import AirflowSkipException
from airflow.providers.slack.notifications.slack import send_slack_notification

# I struggled to move this to sandbox, because the fileloc would not update.
REPO_DIR = f"{os.path.expandvars('$MWM_SANDBOX')}/airflow/ApogeeReduction.jl/"
REPO_BRANCH = "airflow-prod"
DAG_NAME = "ApogeeReduction-prod"

SLACK_CHANNEL = "#apogee-reduction-jl"
SLACK_CHANNEL_KEYS = {
    "#apogee-reduction-jl": "C08B7FKMP16",
    "#apogee-reduction-jl-dev": "C07KQ7BJY5P"
}

n_days_recency = 3
nth_day_verbose = 10

# ~3736 days since 2025-03-01
# ~1.5 node-hr per calendar day


# Want to complete within 10 days, using 2 active runs (so 4 nodes per day).
# -> 2 nodes * 24 hrs * 10 days = 480 node-hrs
# -> 480 node-hrs / 1.5 node-hr per calendar day = 320 calendar days
# -> 3736 days / 320 days = ~12 days

# now let's scale it to actually do 3 active runs and we will use what's available
schedule_interval = timedelta(days=12)
max_active_runs = 3 # reduce to 3 when LCO is up and running

observatories = ("apo", "lco")
sbatch_prefix = re.sub(r"\s+", " ", f"""
    sbatch 
    -vvv
    -D {REPO_DIR}
""") 


class TransferFileSensor(FileSensor):
    def poke(self, context) -> bool:
        if (datetime.now(context["data_interval_start"].tzinfo) - context["data_interval_start"]).days > n_days_recency:
            raise AirflowSkipException("Data interval range is ancient, assuming transfer is complete.")
        return super().poke(context)

def to_sloan_modified_date(data_interval_start):
    return int(Time(data_interval_start).mjd) + 1 # +1 offset to get the most recent day

def send_slack_notification_partial(text):
    def f(context, **k):
        silenced = bool(context["ti"].xcom_pull(task_ids="silenced"))
        print(f"Text (silenced={silenced}): {text}")
        if not silenced:
            return send_slack_notification(text=text, channel=SLACK_CHANNEL)(context)
    return f
    
# Add this function to check SLURM job status
def wait_for_slurm(job_id, min_rows=0):
    t_init = time.time()
    state = None
    while True:
        result = subprocess.run(
            ['squeue', '-j', str(job_id), '--noheader'], 
            capture_output=True, 
            text=True
        )
        if "Invalid job id specified" in result.stderr or result.stdout.count('\n') <= min_rows:
            print(f"Job {job_id} not in queue after {time.time() - t_init:.0f} seconds")
            return  # Job is done
        try:
            this_state = result.strip().split()[4]
        except:
            None
        else:
            if this_state != state:
                if state is not None:
                    print(
                        f"Job {job_id} state changed from {state} to {this_state} "
                        f"after {time.time() - t_init:.0f} seconds"
                    )
                state = this_state
        time.sleep(5)

def submit_and_wait(bash_command, **context):
    silenced = bool(context["ti"].xcom_pull(task_ids="silenced"))

    # Set environment variable for the subprocess
    env = os.environ.copy()
    if silenced:
        env.pop("SLACK_CHANNEL", None)
    else:
        env["SLACK_CHANNEL"] = SLACK_CHANNEL_KEYS[SLACK_CHANNEL]
    
    # Now bash_command comes directly from the arguments
    print(f"Submitting command (silenced={silenced}): {bash_command}")
    result = subprocess.run(bash_command.split(), capture_output=True, text=True, env=env)
    
    if result.returncode != 0:
        raise Exception(f"SLURM submission failed: {result.stderr}")
    
    # Depending on your slurm settings, sbatch might send back something like
    # - "Submitted batch job XXXXX"
    # - "Submitted batch job XXXXX on cluster YYYY"
    job_id = result.stdout.strip().split("job")[1].split()[0]
    print(f"Submitted SLURM job {job_id}")
    wait_for_slurm(job_id)
    return job_id

def is_recent(data_interval_start, **kwargs):
    return (datetime.now(data_interval_start.tzinfo) - data_interval_start).days <= n_days_recency

def data_dir_exists(observatory, data_interval_start, **kwargs):
    return os.path.exists(f"/uufs/chpc.utah.edu/common/home/sdss50/sdsswork/data/apogee/{observatory}/{to_sloan_modified_date(data_interval_start)}/")

def skip_if_not_true(criteria):
    if not criteria:
        raise AirflowSkipException("Condition not met for skipping.")
    return criteria

# Notable dates: https://sdss-wiki.atlassian.net/wiki/spaces/MWM/pages/14659365/Notable+dates
# - 2024-07-18: New APO/APOGEE blue chip.


with DAG(
    DAG_NAME,
    start_date=datetime(2014, 7, 18), 
    schedule_interval=schedule_interval,
    max_active_runs=max_active_runs,
    default_args=dict(retries=1, retry_delay=timedelta(minutes=5)),
    catchup=True,
    on_failure_callback=[
        send_slack_notification_partial(
            text=f"{DAG_NAME} DAG failed on {{{{ ds }}}}",
        )
    ]    
) as dag:

    with TaskGroup(group_id="check") as group_check:
        task_is_recent = PythonOperator(
            task_id="is_recent",
            python_callable=lambda **k: skip_if_not_true(is_recent(**k)),
        )
        task_apo_data_dir_exists = PythonOperator(
            task_id="apo_data_dir_exists",
            python_callable=lambda **k: skip_if_not_true(data_dir_exists("apo", **k)),
        )
        task_lco_data_dir_exists = PythonOperator(
            task_id="lco_data_dir_exists",
            python_callable=lambda **k: skip_if_not_true(data_dir_exists("lco", **k)),
        )

    sjd = PythonOperator(
        task_id="sjd",
        python_callable=lambda data_interval_start, **_: to_sloan_modified_date(data_interval_start),
        trigger_rule="one_success"
    )
    silenced = PythonOperator(
        task_id="silenced",
        # Every 10 days, notify. Otherwise silenced.
        python_callable=lambda ti, **k: skip_if_not_true((ti.xcom_pull(task_ids='sjd') % nth_day_verbose) > 0),
    )

    with TaskGroup(group_id="update") as group_update:
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
        ) >> (
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
        )

    branch = BranchPythonOperator(
        task_id="branch",
        python_callable=lambda **k: [o for o in observatories if data_dir_exists(o, **k)],
        trigger_rule="none_failed_min_one_success"
    )

    group_observatories = []
    for observatory in observatories:  # Changed order to LCO first
        with TaskGroup(group_id=observatory) as group:

            initial_notification = PythonOperator(
                task_id="initial_notification",
                python_callable=lambda **_: None,  # Simple no-op function
                on_success_callback=[
                    send_slack_notification_partial(
                        text=f"Waiting for {observatory.upper()} data transfer for SJD {{{{ task_instance.xcom_pull(task_ids='sjd') }}}} "
                            f"(night of {{{{ ds }}}}). "
                            f"Check here for transfer status: https://data.sdss5.org/sas/sdsswork/data/staging/{observatory}/log/mos/",
                    )
                ],
            )

            filepath = f"/uufs/chpc.utah.edu/common/home/sdss50/sdsswork/data/staging/{observatory}/log/mos/{{{{ task_instance.xcom_pull(task_ids='sjd') }}}}/transfer-{{{{ task_instance.xcom_pull(task_ids='sjd') }}}}.done"
            transfer = TransferFileSensor(
                task_id="transfer",
                filepath=filepath,
                on_success_callback=[
                    send_slack_notification_partial(
                        text=f"{observatory.upper()} data transfer complete for SJD {{{{ task_instance.xcom_pull(task_ids='sjd') }}}} "
                             f"(night of {{{{ ds }}}}). Starting reduction pipeline. Exposure list available at https://data.sdss5.org/sas/sdsswork/data/apogee/{observatory}/{{{{ task_instance.xcom_pull(task_ids='sjd') }}}}/{{{{ task_instance.xcom_pull(task_ids='sjd') }}}}.log.html",
                    )
                ],
                on_failure_callback=[
                    send_slack_notification_partial(
                        text=f"{observatory.upper()} data transfer on SJD {{{{ task_instance.xcom_pull(task_ids='sjd') }}}} "
                             f"(night of {{{{ ds }}}}) is incomplete. "
                             f"Please check https://data.sdss5.org/sas/sdsswork/data/staging/{observatory}/log/mos/ and investigate.",
                    )
                ],
                on_skipped_callback=[
                    send_slack_notification_partial(
                        text=f"{observatory.upper()} data on SJD {{{{ task_instance.xcom_pull(task_ids='sjd') }}}} "
                             f"(night of {{{{ ds }}}}) assumed to already exist: not awaiting `done` file. Starting reduction pipeline.",
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
                op_kwargs=dict(
                    bash_command=f"{sbatch_prefix} --job-name=ar_dark_cal_{observatory}_{{{{ ti.xcom_pull(task_ids='sjd') }}}} src/cal_build/run_dark_cal.sh {observatory} {{{{ ti.xcom_pull(task_ids='sjd') }}}} {{{{ ti.xcom_pull(task_ids='sjd') }}}}",
                ),
                trigger_rule="none_failed_min_one_success" # requires one success from initial notification or file sensor
            )

            flats = PythonOperator(
                task_id="flats",
                python_callable=submit_and_wait,
                op_kwargs=dict(
                    bash_command=f"{sbatch_prefix} --job-name=ar_flat_cal_{observatory}_{{{{ ti.xcom_pull(task_ids='sjd') }}}} src/cal_build/run_flat_cal.sh {observatory} {{{{ ti.xcom_pull(task_ids='sjd') }}}} {{{{ ti.xcom_pull(task_ids='sjd') }}}}",
                ),
            )
            
            science = PythonOperator(
                task_id="science",
                python_callable=submit_and_wait,
                op_kwargs=dict(
                    bash_command=f"{sbatch_prefix} --job-name=ar_all_{observatory}_{{{{ ti.xcom_pull(task_ids='sjd') }}}} src/run_scripts/run_all.sh {observatory} {{{{ ti.xcom_pull(task_ids='sjd') }}}}",
                ),
                on_success_callback=[
                    send_slack_notification_partial(
                        text=f"{observatory.upper()} science frames reduced for SJD {{{{ ti.xcom_pull(task_ids='sjd') }}}} (night of {{{{ ds }}}}).",
                    )
                ],        
                on_failure_callback=[
                    send_slack_notification_partial(
                        text=f"{observatory.upper()} science frame reduction failed for SJD {{{{ ti.xcom_pull(task_ids='sjd') }}}} (night of {{{{ ds }}}}). :picard_facepalm:",
                    )
                ]               
            )   
            initial_notification >> darks # if initial notification is skipped, darks will skip too
            initial_notification >> transfer >> darks >> flats >> science
            #elif observatory == "lco":
            #    task_lco_data_dir_exists >> initial_notification >> transfer >> darks >> flats >> science
            #else:
            #    raise ValueError(f"Unknown observatory: {observatory}")
            
        group_observatories.append(group)

    # Add final notification task
    final_notification = PythonOperator(
        task_id="completion_notification",
        python_callable=lambda **_: None,  # dummy function that does nothing
        on_success_callback=[
            send_slack_notification_partial(
                text=f"ApogeeReduction pipeline completed successfully for SJD {{{{ ti.xcom_pull(task_ids='sjd') }}}} (night of {{{{ ds }}}}). Both observatories processed.",
            )
        ],
        trigger_rule="none_failed_min_one_success"
    )
    
    group_check >> sjd >> (silenced, group_update) >> branch >> group_observatories >> final_notification
