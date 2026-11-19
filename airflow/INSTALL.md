# Installing the APOGEE Airflow layer (card O3) — steps requiring AKS action

Everything in `airflow/` is developed and tested from the repo; nothing here
self-installs. The steps below are the deploy checklist. Items marked **AKS**
need the sdssv service account, secrets, or a policy decision.

## 0. State found 2026-09-02 (why some of this is needed)

- `AIRFLOW_HOME=/mnt/home/sdssv/airflow` exists (owned `asaydjari:sdssv`,
  group-writable) with a sqlite db from the broken env's Airflow **3.0.6**
  (rebuild target below is **3.3.1**, AKS 2026-09-03) and the previous
  PRODUCTION DAG, `dags/ar_main.py` — CCA-native (sbatch of run_all.sh from
  `/mnt/home/sdssv/gitcode/ApogeeReduction.jl` into
  `/mnt/ceph/users/sdssv/work/daily/`; zero Utah references), the DAG that
  ran the dailies until 2026-05. Its imports are **Airflow-3-compatible**
  (verified by importing it under both 3.0.6 and 3.3.1:
  `airflow.operators.bash`,
  `airflow.sensors.filesystem`, `airflow.operators.python` are provider
  shims, and astropy/pytz/`apache-airflow-providers-slack` are all present
  in the shared env). It is superseded feature-for-feature by the new
  `apogee_daily` DAG (mapping table in README.md) and should be retired to
  `dags/attic/` — moved aside, not deleted — so the old and new daily DAGs
  never schedule side by side.
- The uv env `/mnt/home/sdssv/uv_env/airflow_env` is **broken**: it was
  built against `/usr/bin/python3.11` (`pyvenv.cfg: home = /usr/bin`),
  which no longer exists after the OS upgrade → every entry point fails
  with `bad interpreter` (`airflow version` itself cannot run; the 3.0.6
  version is read from the env's `apache_airflow-3.0.6.dist-info`). It must
  be rebuilt (step 2).
- `sla_miss_callback` was **removed in Airflow 3.0**; the SLA layer here is
  `dagrun_timeout` (20 h) + the cross-machine watchdog (step 6).

## 1a. ONE-TIME: migrate the production checkout to main (**AKS — deploy
## decision, review before running**)

The production clone `/mnt/home/sdssv/gitcode/ApogeeReduction.jl` (what
`AR_REPO` defaults to and what `update.repo` pulls) is currently on the old
`airflow` branch, whose `run_all.sh` takes the OLD 5-arg signature. The DAG
drives the NEW 12-arg run_all.sh (R1 layout), so the clone must move to the
consolidated `main` of andrew-saydjari/ApogeeReduction.jl (already its
`origin`):

```bash
cd /mnt/home/sdssv/gitcode/ApogeeReduction.jl
git fetch origin
git checkout main
git pull --ff-only
```

After this, the DAG's `update.repo` task keeps it current with
`git pull --ff-only` every morning (deliberate NEW behavior — the old
`update.repo` was status-only, so pulling was effectively off).

Context (AKS 2026-09-03, confirmed): the old `airflow` branch existed to
insulate production from comment-toggle edits during bulk runs — a churn
class eliminated by the R1 `run_madgics`/env flags, so dailies now run off
consolidated `main`. Fallback if main churn ever causes issues: cut a
`production-airflow` branch and update it deliberately.

## 1b. Deploy the DAGs (symlink from the repo checkout)

As on Utah, DAGs live in the repo and are symlinked into the sandbox:

```bash
cd /mnt/home/sdssv/airflow/dags
# retire the previous production DAG (superseded feature-for-feature by the
# new one — see the mapping table in README.md; keep it for reference)
mkdir -p attic && mv ar_main.py attic/
REPO=/mnt/home/sdssv/gitcode/ApogeeReduction.jl  # once o3/airflow-dags merges
ln -s $REPO/airflow/dags/ar_common.py .
ln -s $REPO/airflow/dags/apogee_daily.py .
ln -s $REPO/airflow/dags/apogee_heartbeat.py .
```

## 2. Rebuild the uv env (**AKS**, or anyone in group sdssv)

Target (AKS 2026-09-03): **apache-airflow==3.3.1 on Python 3.12** — the
version pair the DAGs were validated against in staging (parse + full
`dags test`).

```bash
cd /mnt/home/sdssv/uv_env
mv airflow_env airflow_env.broken_2026_09       # keep for forensics
uv python install 3.12
uv venv --python 3.12 airflow_env
uv pip install -p airflow_env/bin/python \
    "apache-airflow==3.3.1" apache-airflow-providers-standard \
    --constraint https://raw.githubusercontent.com/apache/airflow/constraints-3.3.1/constraints-3.12.txt
```

Migrating the existing sqlite db (created by 3.0.6) is handled by
`airflow db migrate` on first start; the on-disk DAG-run history survives.

Notes: pin to a **uv-managed** python (not `/usr/bin/...`) so the next OS
upgrade doesn't break it again. The DAGs use only stdlib + the standard
provider (Slack goes through `urllib`, MJD math avoids astropy), so no other
packages are required.

## 2b. ONE-TIME: almanac local clone + env (done 2026-09-04)

almanac runs from a local git clone with an editable install, like the
AR/arM pipelines — `update.repo` pulls it daily, every invocation logs
`git rev-parse HEAD`, and the uvx wheel-cache trap cannot occur:

```bash
cd /mnt/home/sdssv/gitcode
git clone https://github.com/andrew-saydjari/almanac.git
cd /mnt/home/sdssv/uv_env
uv venv --python 3.12 almanac_env
uv pip install -p almanac_env/bin/python -e /mnt/home/sdssv/gitcode/almanac
```

After each daily pull, `update.repo` re-runs that `uv pip install -e` through
`airflow/scripts/sync_almanac_env.sh`, so a dependency added to almanac's
`pyproject.toml` (e.g. healpy) is installed the next morning; if the install
or `import almanac.cli` fails, the task fails (and retries, then notifies).
It does not read almanac's `uv.lock`: installed versions that still satisfy
`pyproject.toml` are kept. Manual runs of `scripts/bulk/run_bulk.sh` use the
same env and get no re-sync of their own.

Override paths via `ALMANAC_DIR` / `ALMANAC_BIN` (read by ar_common.py,
run_all.sh, run_bulk.sh). Smoke-tested 2026-09-04 with a live DB query
through the tunnel (raw/ layout + missing_exposures verified).

## 3. Secrets / environment (**AKS**)

```bash
cp $REPO/airflow/scripts/airflow_env.sh.example /mnt/home/sdssv/airflow/airflow_env.sh
chmod 600 /mnt/home/sdssv/airflow/airflow_env.sh
# then fill in: SLACK_TOKEN (bot apogeereductionjl), PUBLIC_URL_SLUG.
# AR_SLACK_CHANNEL ships as prod C08B7FKMP16 (promoted 2026-09-06); leave it
# as-is. Demote to dev C07KQ7BJY5P only for testing (see step 7.3).
```

The almanac DB identity is NOT in this file: it comes from
`~/.almanac/config.yaml` + `~/.pgpass` of the account running Airflow, and
the DB is reached through localhost:63333, forwarded by the
`mwm-tunnel.service` systemd user unit (step 4), which owns the `mwm` ssh
ControlMaster in the foreground and establishes the forward itself. The
`tunnel_check` task checks the master (`ssh -O check mwm`) and runs a real
`select 1` through the forward; when either fails, it names the unit to
check (`systemctl --user status mwm-tunnel.service`, `journalctl --user -u
mwm-tunnel.service`). The account running Airflow therefore needs the `mwm`
ssh alias (with a `ControlPath`), `~/.almanac/config.yaml` + `~/.pgpass`,
and `psql` on the PATH of airflow_env.sh (currently on asaydjari, not sdssv —
step 7). `tunnel_check`'s DB user/name default to `sdss`/`sdss5db` (the
current `~/.almanac/config.yaml`); override with `AR_DB_USER`/`AR_DB_NAME`
if that file changes.

## 4. systemd user units (**AKS**; replaces `startAirflow()`/tmux)

Three user units, kept in the repo as templates (`airflow/systemd/*.service.in`)
whose paths are parameters, not `%h` (the account running the units is not
the one whose home holds AIRFLOW_HOME — `%h` made the unit loop 7,760 times
on a missing EnvironmentFile on 2026-11-17):

| unit | role |
|---|---|
| `mwm-tunnel.service` | the `mwm` ssh ControlMaster **and** the DB forward `localhost:63333 → operations.sdss.org:5432`, one foreground ssh (`ControlMaster=yes`, `ControlPersist=no`, keep-alives 30 s × 3, `ExitOnForwardFailure=yes`). Pre-start asks any existing master to exit and removes a stale socket. Restart cap 3 starts/h, so a changed password cannot lock the Utah account; after that it stays failed until reset by hand |
| `airflow-standalone.service` | `airflow standalone`; ordered `After=`/`Wants=` the tunnel; restart cap 5 starts/30 min, then `OnFailure=` → Slack |
| `airflow-failure-notify.service` | the Slack outage alert |

On ccalin051, as the account chosen in step 7 (today: asaydjari). Render
into a scratch directory, compare with what is installed, then copy:

```bash
loginctl enable-linger $USER        # keeps user services alive w/o login
REPO=/mnt/home/sdssv/gitcode/ApogeeReduction.jl
OUT=$(mktemp -d)
$REPO/airflow/scripts/render_systemd_units.sh \
    --sshpass /mnt/home/asaydjari/local/bin/sshpass \
    --pass-file /mnt/home/asaydjari/pass_file \
    "$OUT"
# defaults: --airflow-home /mnt/home/sdssv/airflow,
#   --env-file <airflow-home>/airflow_env.sh,
#   --airflow-bin /mnt/home/sdssv/uv_env/airflow_env/bin/airflow,
#   --tunnel-host mwm, --db-forward 63333:operations.sdss.org:5432
for u in mwm-tunnel airflow-standalone airflow-failure-notify; do
  diff -u ~/.config/systemd/user/$u.service "$OUT/$u.service"
done
systemd-analyze --user verify "$OUT"/*.service
mkdir -p ~/.config/systemd/user
for u in mwm-tunnel airflow-standalone airflow-failure-notify; do
  cp ~/.config/systemd/user/$u.service ~/.config/systemd/user/$u.service.bak-$(date +%Y%m%d) 2>/dev/null
  cp "$OUT/$u.service" ~/.config/systemd/user/
done
systemctl --user daemon-reload
systemctl --user enable mwm-tunnel.service airflow-standalone.service
```

Against the units installed on 2026-11-17 the rendered `airflow-*` units
differ only in comments; `mwm-tunnel.service` gains the DB forward
(`-L 63333:operations.sdss.org:5432`, `ExitOnForwardFailure=yes`), keep-alive
30 s instead of 60 s, and an explicit password-file path. **The new tunnel
unit takes effect only when it restarts** — pick a moment when no daily is
in the almanac/tunnel stage:

```bash
systemctl --user restart mwm-tunnel.service   # drops + rebuilds the master
systemctl --user status mwm-tunnel.service
ssh -O check mwm
PGCONNECT_TIMEOUT=10 psql -X -w -h localhost -p 63333 -U sdss -d sdss5db -Atc 'select 1'
```

Restarting `airflow-standalone.service` is not needed for these changes
(its unit file changes are comments only; DAG changes are picked up by the
DAG processor after promotion).

If the tunnel unit hit its start limit (journal shows "Start request
repeated too quickly"), fix the cause (e.g. the password file after a
password change), then
`systemctl --user reset-failed mwm-tunnel.service && systemctl --user start mwm-tunnel.service`.

If `loginctl enable-linger` is refused on the node, ask SCC to enable
lingering for the service account (this is the one hard admin dependency).
Do not run the interim `startAirflow()` tmux helper alongside the unit
(sqlite + two schedulers).

## 5. Unpause + first supervised run

There is ONE daily DAG (`apogee_daily`, LCO then APO serially) plus the
heartbeat:

```bash
airflow dags list                      # expect: apogee_daily, apogee_heartbeat
airflow dags unpause apogee_heartbeat
airflow dags unpause apogee_daily
# Routine posts go to the PROD channel (promoted 2026-09-06; airflow_env.sh
# sets AR_SLACK_CHANNEL=C08B7FKMP16 explicitly). A one-off test run can still
# be routed to dev with conf, without touching the env file. Runs in slurm
# mode (the default) unless overridden:
airflow dags trigger apogee_daily --conf \
  '{"mjd": <recent>, "run_kind": "test", "slack_channel": "C07KQ7BJY5P"}'
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
alert/6 h per condition. It posts to prod `C08B7FKMP16` unless
`AR_SLACK_ALERT_CHANNEL` says otherwise — this scrontab environment does not
see airflow_env.sh, so no channel needs to be passed here.

## 7. Open decisions (**AKS**)

1. **Which account hosts the service?** AIRFLOW_HOME is sdssv's, but the
   mwm ControlMaster, `~/.almanac` identity, juliaup, and public_www slug
   currently live on asaydjari. Either move those onto sdssv (cleaner) or
   run the systemd unit as asaydjari with AIRFLOW_HOME=/mnt/home/sdssv/airflow.
2. **local vs slurm mode** — `slurm` is the DEFAULT (production dailies as
   run historically: one sbatch of run_all.sh per night, 5-s squeue polls;
   AKS launches the Airflow instance himself, so submission is
   human-initiated per cluster etiquette). `local` runs every step under
   `nice -n 10` on ccalin051 (~1h50m at 8 workers) — the testing mode, or a
   deliberate alternative on the dedicated node. Flip with
   `AR_AIRFLOW_MODE=local` or conf `{"mode": "local"}`.
3. **Slack channel** — SETTLED 2026-09-06: promoted to prod. The supervised
   nights on dev are over, and `airflow_env.sh` (step 3) now sets
   `AR_SLACK_CHANNEL=C08B7FKMP16` **explicitly**, which is also the code
   default (`ar_common.py` `SLACK_CHANNEL_PROD`, ar_main.py's behavior).
   Explicit rather than simply deleting the override, because the
   deployed channel should be greppable in one place, and because
   "unset" used to be the state that misrouted the outage alert (below).
   To demote for testing: swap to `C07KQ7BJY5P` in airflow_env.sh (the
   commented line is right there), or, better for a one-off, pass conf
   `{"slack_channel": "C07KQ7BJY5P"}` and leave the deployment alone.

   **Outage alerts are a separate knob.** `airflow-failure-notify.service`
   and `check_airflow_heartbeat.sh` read `AR_SLACK_ALERT_CHANNEL` and
   default to prod `C08B7FKMP16`. They used to default to *dev*, so an
   unset (or dev-demoted) `AR_SLACK_CHANNEL` sent "the APOGEE daily
   orchestrator is down" to a channel nobody watches. Demoting the DAG's
   routine posts no longer moves these alerts; set
   `AR_SLACK_ALERT_CHANNEL` explicitly only to test the alerting path.
4. **sqlite + LocalExecutor**: fine at this scale (2 DAG runs/day + a
   heartbeat); if concurrency is ever raised, move to postgres.

## 8. Re-running a failed night

Use a **fresh trigger** of the night with `almanac_clobber`:

```bash
airflow dags list-runs apogee_daily --state failed     # find the night
MJD=61360
NOW=$(date -u +%Y-%m-%dT%H:%M:%S+00:00)
airflow dags trigger apogee_daily \
    --logical-date "$NOW" \
    --run-id "manual__rerun_${MJD}_$(date -u +%Y%m%dT%H%M%SZ)" \
    --conf "{\"mjd\": $MJD, \"almanac_clobber\": true}"
```

(Checked against `--help` of the deployed Airflow 3.3.1: `dags trigger`
takes `-c/--conf`, `-r/--run-id`, `-l/--logical-date`; `dags list-runs`
takes `--state`.) Notes:

- **Pass `--logical-date`.** In Airflow 3 a trigger without it creates a
  run with no data interval; before this change `setup.date_mjd` then
  failed with `KeyError: 'data_interval_start'` (seen 2026-11-17 on the
  61360 re-run). `setup.mjd`/`setup.date_mjd` now fall back to
  `logical_date`, then `run_after`, then now, and with `mjd` in conf the
  date label is derived from the MJD itself (it used to come from the
  trigger time, i.e. name the wrong night for a past-day re-run), so the
  command also works without `--logical-date` once this is promoted. Keep
  passing it anyway: production on the pre-promotion DAG needs it. The
  logical date only has to be unique; the night comes from `mjd`.

- A manual run is never "stale"; it waits while another run is active
  (`max_active_runs=1`), and stays queued while the DAG is paused.
- It repeats the cheap update/setup/resolve/wait/tunnel tasks and runs
  both observatories; for the one that already succeeded, almanac is
  rebuilt (clobber) and run_all.sh's `checkpoint_mode=commit_exists` skips
  existing products.
- A killed almanac used to leave a truncated 96-byte HDF5 stub, and
  **`almanac_clobber` did not recover from it** on the pre-change DAG:
  almanac opens its output in append mode and failed on the truncated
  file (`OSError: Unable to synchronously open file (truncated file:
  eof = 96 ...)`; the earlier note that append mode would write into the
  stub was wrong). Before re-running a night on the **pre-promotion** DAG,
  move such a stub aside by hand (e.g. `mv allobs_lco_61360.h5
  allobs_lco_61360.h5.stub-YYYYMMDD` in `/mnt/ceph/users/sdssv/work/daily/almanac/`).
  With this change almanac never touches an existing output: it writes a
  fresh `<file>.partial` that replaces the file only on success, an
  existing file is reused only if it opens and holds the night, and no
  stub is left by a killed run in the first place. `almanac_clobber` then
  just forces a fresh almanac, which is what a re-run after a DB/tunnel
  problem should have.

**Why not clear the failed run's tasks.** On 2026-11-17 clearing tasks of a
failed run did not restart it (the exact cause could not be established
afterwards), and even when clearing does re-queue a run it is the wrong
tool here:

- A cleared run keeps its **original conf**: you cannot add
  `almanac_clobber` (or change `mjd`, `mode`, ...) to it. Before this
  change that meant almanac skipped on the 96-byte stub and runlist read an
  empty file.
- `airflow tasks clear` selects by **DAG, task regex and date window, not
  by run**: without `-s/-e` it clears the matching tasks in *every* run of
  the DAG, including old failed ones (which would then be failed by
  `stale_run_guard` at best).
- Clearing only the failed task (`-f`, or the UI default) without
  downstream (`-d`) leaves the rest of the chain `upstream_failed`; and the
  `all_done` tail tasks (metrics, summary) that already ran are not redone,
  so the metrics/Slack record of the night is incomplete.
- A cleared run competes for the single active-run slot like any other.
