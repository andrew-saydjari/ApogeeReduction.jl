# Installing the APOGEE Airflow layer (card O3) — steps requiring AKS action

Everything in `airflow/` is developed and tested from the repo; nothing here
self-installs. The steps below are the deploy checklist. Items marked **AKS**
need the sdssv service account, secrets, or a policy decision.

## 0. State found 2026-09-02 (why some of this is needed)

- `AIRFLOW_HOME=/mnt/home/sdssv/airflow` exists (owned `asaydjari:sdssv`,
  group-writable) with an Airflow **3.0.6** sqlite db and the previous
  PRODUCTION DAG, `dags/ar_main.py` — CCA-native (sbatch of run_all.sh from
  `/mnt/home/sdssv/gitcode/ApogeeReduction.jl` into
  `/mnt/ceph/users/sdssv/work/daily/`; zero Utah references), the DAG that
  ran the dailies until 2026-05. Its imports are **Airflow-3-compatible**
  (verified by importing it under 3.0.6: `airflow.operators.bash`,
  `airflow.sensors.filesystem`, `airflow.operators.python` are provider
  shims, and astropy/pytz/`apache-airflow-providers-slack` are all present
  in the shared env). It is superseded feature-for-feature by the new DAGs
  (mapping table in README.md) and should be retired to `dags/attic/` —
  moved aside, not deleted — so the old and new daily DAGs never schedule
  side by side.
- The uv env `/mnt/home/sdssv/uv_env/airflow_env` is **broken**: it was
  built against `/usr/bin/python3.11` (`pyvenv.cfg: home = /usr/bin`),
  which no longer exists after the OS upgrade → every entry point fails
  with `bad interpreter` (`airflow version` itself cannot run; the 3.0.6
  version is read from the env's `apache_airflow-3.0.6.dist-info`). It must
  be rebuilt (step 2).
- `sla_miss_callback` was **removed in Airflow 3.0**; the SLA layer here is
  `dagrun_timeout` (16 h) + the cross-machine watchdog (step 6).

## 1. Deploy the DAGs (symlink from the repo checkout)

As on Utah, DAGs live in the repo and are symlinked into the sandbox:

```bash
cd /mnt/home/sdssv/airflow/dags
# retire the previous production DAG (superseded feature-for-feature by the
# new ones — see the mapping table in README.md; keep it for reference)
mkdir -p attic && mv ar_main.py attic/
REPO=/mnt/home/sdssv/gitcode/ApogeeReduction.jl  # once o3/airflow-dags merges
ln -s $REPO/airflow/dags/ar_common.py .
ln -s $REPO/airflow/dags/apogee_daily.py .
ln -s $REPO/airflow/dags/apogee_heartbeat.py .
```

## 2. Rebuild the uv env (**AKS**, or anyone in group sdssv)

```bash
cd /mnt/home/sdssv/uv_env
mv airflow_env airflow_env.broken_2026_09       # keep for forensics
uv python install 3.11
uv venv --python 3.11 airflow_env
uv pip install -p airflow_env/bin/python \
    "apache-airflow==3.0.6" apache-airflow-providers-standard \
    --constraint https://raw.githubusercontent.com/apache/airflow/constraints-3.0.6/constraints-3.11.txt
```

Notes: pin to a **uv-managed** python (not `/usr/bin/...`) so the next OS
upgrade doesn't break it again. The DAGs use only stdlib + the standard
provider (Slack goes through `urllib`, MJD math avoids astropy), so no other
packages are required.

## 3. Secrets / environment (**AKS**)

```bash
cp $REPO/airflow/scripts/airflow_env.sh.example /mnt/home/sdssv/airflow/airflow_env.sh
chmod 600 /mnt/home/sdssv/airflow/airflow_env.sh
# then fill in: SLACK_TOKEN (bot apogeereductionjl), PUBLIC_URL_SLUG,
# AR_SLACK_CHANNEL (dev C07KQ7BJY5P until trusted, then prod C08B7FKMP16)
```

The almanac DB identity is NOT in this file: it comes from
`~/.almanac/config.yaml` + `~/.pgpass` of the account running Airflow, and
the DB is reached through localhost:63333 forwarded over the persistent
`mwm` ssh ControlMaster. The `tunnel_check` task verifies the master
(`ssh -O check mwm`) and (re)adds the port-forward, but **only AKS can
reopen the master itself** (`sshpass -f ~/pass_file ssh -fN mwm`) — the task
fails with exactly that message when the master is down. The account running
Airflow therefore needs the `mwm` ssh alias + ControlMaster config, and
`~/.almanac/config.yaml` + `~/.pgpass` (currently on asaydjari, not sdssv —
decide which account hosts the service, step 7).

## 4. systemd user service (**AKS**; replaces `startAirflow()`/tmux)

On ccalin051, as the account chosen in step 7:

```bash
loginctl enable-linger $USER        # keeps user services alive w/o login
mkdir -p ~/.config/systemd/user
cp $REPO/airflow/systemd/airflow-standalone.service \
   $REPO/airflow/systemd/airflow-failure-notify.service \
   ~/.config/systemd/user/
systemctl --user daemon-reload
systemctl --user enable --now airflow-standalone.service
systemctl --user status airflow-standalone.service
```

`Restart=on-failure` (60 s backoff, 5 tries/30 min) restarts crashes; when
the restart budget is exhausted the `OnFailure=` unit curls a Slack alert.
The interim `startAirflow()` tmux helper keeps working for manual runs — do
not run both at once (sqlite + two schedulers).

If `loginctl enable-linger` is refused on the node, ask SCC to enable
lingering for the service account (this is the one hard admin dependency).

## 5. Unpause + first supervised run

```bash
airflow dags list
airflow dags unpause apogee_heartbeat
airflow dags unpause apogee_daily_apo
airflow dags unpause apogee_daily_lco
# supervised first night (dev Slack channel is the default):
airflow dags trigger apogee_daily_apo --conf '{"mjd": <recent>, "run_kind": "test"}'
```

## 6. Cross-machine watchdog (**AKS**: scrontab on rusty)

The heartbeat DAG touches
`/mnt/ceph/users/sdssv/work/daily/metrics/airflow_heartbeat.txt` every
10 min. From a rusty login node (different failure domain than ccalin051):

```
scrontab -e   # add:
*/30 * * * * SLACK_TOKEN=... /mnt/home/sdssv/gitcode/ApogeeReduction.jl/airflow/scripts/check_airflow_heartbeat.sh >> /mnt/home/sdssv/airflow/logs/watchdog.log 2>&1
```

(or source the token from a chmod-600 file inside a tiny wrapper rather than
inlining it in scrontab). It alerts on: stale heartbeat (>45 min) and no new
daily-metrics row in >30 h (the SLA-miss substitute), rate-limited to one
alert/6 h per condition.

## 7. Open decisions (**AKS**)

1. **Which account hosts the service?** AIRFLOW_HOME is sdssv's, but the
   mwm ControlMaster, `~/.almanac` identity, juliaup, and public_www slug
   currently live on asaydjari. Either move those onto sdssv (cleaner) or
   run the systemd unit as asaydjari with AIRFLOW_HOME=/mnt/home/sdssv/airflow.
2. **local vs slurm mode** — `local` (default) runs the chain under
   `nice -n 10` on ccalin051; recommended: the node is dedicated, the chain
   fits (~1h50m at 8 workers, less at 16), and no scheduler round-trips.
   `slurm` mode (one sbatch of run_all.sh per night, 5-min squeue polls) is
   fully written but OFF: automating sbatch is your call per cluster
   etiquette. Flip with `AR_AIRFLOW_MODE=slurm` when/if desired.
3. **Slack channel promotion** dev → prod (`AR_SLACK_CHANNEL=C08B7FKMP16`)
   once a few supervised nights look right.
4. **sqlite + LocalExecutor**: fine at this scale (2 DAG runs/day + a
   heartbeat); if concurrency is ever raised, move to postgres.
