"""apogee_heartbeat — liveness beacon for the Airflow instance (card O3/N3).

Every 10 minutes, touch a timestamp file on ceph. Any other machine can then
tell whether the Airflow scheduler on ccalin051 is alive by checking the file
age — that check is airflow/scripts/check_airflow_heartbeat.sh, intended to
run from scrontab on rusty (cross-machine watchdog; see airflow/INSTALL.md).

The same script also implements the SLA-miss substitute (sla_miss_callback
was removed in Airflow 3.0): it alerts when the daily metrics table has not
gained a fresh row recently.
"""

from __future__ import annotations

import json
import os
import socket
import sys
from datetime import datetime, timezone

# make sibling ar_common importable regardless of how the bundle is parsed
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from airflow import DAG
from airflow.providers.standard.operators.python import PythonOperator

import ar_common as C


def beat(**context):
    path = C.AR_HEARTBEAT_FILE
    os.makedirs(os.path.dirname(path), exist_ok=True)
    payload = {
        "utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "host": socket.gethostname(),
        "dag_run": getattr(context.get("dag_run"), "run_id", "?"),
    }
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        json.dump(payload, f)
        f.write("\n")
    os.replace(tmp, path)
    print(f"heartbeat -> {path}: {payload}")


with DAG(
    dag_id="apogee_heartbeat",
    description="10-min liveness beacon consumed by the cross-machine "
                "watchdog script",
    start_date=datetime(2026, 9, 1),
    schedule="*/10 * * * *",
    catchup=False,
    max_active_runs=1,
    tags=["apogee", "watchdog"],
) as dag:
    PythonOperator(task_id="beat", python_callable=beat)
