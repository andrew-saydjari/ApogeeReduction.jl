# APOGEE daily pipeline on Airflow (cards O3 / N1 / N3)

Airflow layer for the daily reduction: ONE serial DAG (`apogee_daily`, LCO
then APO) whose per-observatory task groups split at
`scripts/daily/run_all.sh`'s natural step boundaries, Slack notifications, a
heartbeat/watchdog layer, and the N1 per-night metrics seed.

```
airflow/
├── dags/
│   ├── ar_common.py          shared config/helpers (Slack, paths, metrics)
│   ├── apogee_daily.py       apogee_daily (ONE DAG, lco→apo serial)
│   └── apogee_heartbeat.py   10-min liveness beacon
├── systemd/                  user unit TEMPLATES (*.service.in)
├── scripts/
│   ├── check_airflow_heartbeat.sh   cross-machine watchdog (scrontab)
│   ├── check_almanac_output.py      does an almanac file hold a night?
│   ├── render_systemd_units.sh      fill the unit templates (never installs)
│   └── airflow_env.sh.example       env template (secrets NOT committed)
├── tests/                    unittest suite (see "Tests" below)
├── INSTALL.md                deploy checklist (AKS actions)
└── README.md                 this file
```

## DAG structure (`apogee_daily` — ONE DAG, LCO then APO serially)

One DAG, observatories chained serially (LCO first), restoring the
production ar_main.py template per AKS review 2026-09-03: the Airflow
metadata DB is SQLite (single-writer), so parallel per-observatory runs
contend. A skipped observatory (no data) does not block the other
(`none_failed` on each group's root); a FAILED LCO still blocks APO
(ar_main.py's serial semantics).

```
# Regenerated from the DAG's actual edge wiring (65 tasks, 75 edges;
# edge-set sha256[:16] b6d9abefd75ee370 after the 2026-11-17 stale-run
# guard; sha = sha256 of sorted "up>down" lines):
stale_run_guard >> update.repo
update.repo >> update.sdsscore >> update.sync_logs >> setup.{mjd, date_mjd}
setup.{mjd, date_mjd} >> lco.resolve_night

# per observatory group (lco shown; apo identical):
lco.resolve_night >> lco.wait_for_raw >> lco.tunnel_check >> lco.almanac
  >> lco.runlist >> lco.start_notification >> lco.select_mode
lco.select_mode >> {lco.local.p3d2d | lco.slurm_submit}          # branch
lco.local.p3d2d >> lco.local.quartz_flats.(runlist >> traces >> extract
  >> relflux) >> lco.local.dome_flats.(runlist >> traces >> extract
  >> relflux) >> lco.local.p2d1d_full >> lco.local.plots
  >> lco.local.dashboard
lco.local.dashboard >> {lco.local.madgics_gate, lco.join, lco.chain_status}
lco.local.madgics_gate >> lco.local.madgics_pipeline
  >> lco.local.spectra_workup >> lco.local.workup_cleanup
  # single workup (legacy workup.jl retired 2026-09-04); workup_cleanup =
  # W3 validator gate, then rm of the raw batch fiber dirs
lco.local.workup_cleanup >> {lco.join, lco.chain_status}
lco.slurm_submit >> {lco.join, lco.chain_status}
lco.join >> lco.metrics_append >> lco.daily_summary

# serial edge + tail:
{lco.chain_status, lco.daily_summary} >> apo.resolve_night
{apo.chain_status, apo.daily_summary} >> completion_notification  # only leaf

# non-default trigger rules:
#   *.resolve_night, completion_notification         none_failed
#   *.join, *.chain_status              none_failed_min_one_success
#   *.metrics_append, *.daily_summary                     all_done
```

- **Night selection** (`setup.mjd` / `setup.date_mjd`): ar_main.py's exact
  arithmetic — `int(Time(data_interval_start).mjd) - 1` and the ET date
  − 2 days. The Globus Utah→CCA sync is a fixed 4:00–4:30 am daily task,
  so the 7 am ET run processes the data that landed that morning. On this
  Airflow 3 instance a cron schedule is a CronTriggerTimetable
  (`data_interval_start` == trigger time), and the 300 historical daily
  almanac files confirm the expression reproduces exactly what production
  processed (see the compute_mjd docstring). A manual trigger with conf
  `{"mjd": N}` overrides cleanly through the whole chain (past-day reruns).
- **wait_for_raw** gates reduction on raw-transfer COMPLETION (the
  `<mjd>.md5sum` marker) — ported from ar_main.py's `TransferFileSensor`
  semantics, including "ancient night (mjd < 59148) → assume transferred".
  Cheap insurance only: normally an instant pass (first poke true);
  reschedule pokes every 30 min, 12 h timeout skips the observatory
  cleanly (`soft_fail`); conf `raw_nowait` skips immediately (smoke tests).
  No-data fast path: a frameless night skips on the first poke once the
  sync round that would have carried it has demonstrably run (the other
  observatory's same-night marker, or a later night's marker for this one)
  — so an observatory that took no APOGEE data (weather/engineering, e.g.
  lco 61284/61285) does not hold the serial chain for the 12 h timeout;
  "incomplete" (frames present, md5sum pending) still waits.
- **update.repo**: reports the production checkout state AND
  `git pull --ff-only` — for the AR clone, the arMADGICS clone, AND the
  almanac clone (DELIBERATE new behavior — ar_main.py's repo task was
  status-only; the arM pull was added after the 2026-09-04 first night ran
  a months-stale arM main that predated the raw/-layout almanac handling;
  almanac runs from a local clone + editable env instead of a uvx commit
  pin, per AKS 2026-09-04 — hash logged at every invocation, and the
  update task re-syncs the env for dependency changes);
  **update.sync_logs** rsyncs both observatories'
  `.log.html` pages into the public `apogee_logs` mirror.
- **update.sdsscore is default ON** (O1 punch list #2 — a fresh night after
  any gap dies on missing confSummary otherwise), with new-submodule
  bootstrap (--init + SSH→HTTPS rewrite).
- **tunnel_check** verifies the `mwm` ControlMaster (`ssh -O check`) and
  runs a real `select 1` through the localhost:63333 forward (psql, 10 s
  connect timeout, 15 s statement timeout, 30 s wall clock). The master and
  the forward belong to the `mwm-tunnel.service` systemd unit; on failure
  the task (and its Slack hint) names that unit and the `systemctl`/
  `journalctl` commands to check it. Stopgap: if the master is up but the
  port is not forwarded (a tunnel unit installed before it carried `-L`),
  it adds the forward over the master once and logs a warning to reinstall
  the unit.
- **almanac** (see "Robustness" below): reuses an existing output only if
  it really holds the night, writes to a `.partial` file moved into place
  on success, and runs under a statement timeout, an execution timeout and
  one retry.
- **runlist** maps `make_runlist_all.jl` exit 16 ("no exposures") to a clean
  Airflow skip.
- **Execution modes**: `slurm` (DEFAULT — production dailies as run
  historically; AKS launches the Airflow instance himself, so submission is
  human-initiated) submits the whole run_all.sh chain as one sbatch job with
  ar_main.py's 5-second `squeue -j` polls + a final `sacct` verdict;
  `local` runs every step directly on the
  host (ccalin051) under `nice -n 10` with
  `AR_LOCAL_WORKERS`/`--workers_per_node` = `workers` param and all
  `SLURM_*` env scrubbed — used for testing (`AR_AIRFLOW_MODE=local` or
  conf `{"mode": "local"}`).
- **arMADGICS** is gated by the `run_madgics` param only (no env fallback),
  default TRUE — production ran arMADGICS + its workup unconditionally
  every night (AKS 2026-09-03).
- **spectra_workup** — the SINGLE workup (AKS 2026-09-04: the legacy
  `workup.jl` is retired from the chain — same aggregation job but no row
  contract, no integrity checks, and it deleted the raw batches out from
  under the new workup, the apo 61286 first-night failure): calls the
  arMADGICS-side entrypoint `workup/run_workup.sh <rawdir> <redux> <outdir>`
  (PR #26, MPI tier default) under the same madgics gate. Local mode calls
  it explicitly; slurm mode inherits the identical step from run_all.sh.
  Existence-guarded (skips with a warning until PR #26 is in the checkout).
- **workup_cleanup** — validated deletion, replacing workup.jl's rm role:
  the W3 validator (`workup/validate_workup.jl --K 500`) must PASS, then
  the raw batch fiber dirs are removed (batch_info.txt / full_list_info.h5
  kept as provenance). `AR_KEEP_MADGICS_BATCHES=true` skips deletion
  (dailies). Bulk keeps batches by DEFAULT per the W5 decision
  (`AR_CLEAN_MADGICS_BATCHES=true` opts in to validated cleanup there).
  Cal-only nights (0 batches) skip both validator and rm.
- **dagrun_timeout=20 h** (two serial observatories) stands in for SLAs
  (removed in Airflow 3.0), together with the metrics-freshness check in
  the watchdog script. Run-state correctness: metrics/summary run
  `all_done`, so each group's `chain_status` gate (which also carries the
  serial edge to the next observatory) propagates a chain failure into
  `completion_notification` — the DAG's only leaf, `none_failed` — and the
  DagRun is marked FAILED (verified live: an `all_done` leaf here masked a
  failed chain in the structural test).
- **Skipped-observatory hygiene**: an observatory skipped before any chain
  step (no data / no exposures) writes NO metrics row, and
  `both_observatories_done` requires `n_steps > 0`, so
  `completion_notification` cannot post a false "both observatories
  processed" (also a structural-test catch).

## Relationship to the previous production DAG (`dags/ar_main.py`)

`ar_main.py` in the Airflow sandbox (`/mnt/home/sdssv/airflow/dags/`) is the
CCA-NATIVE production DAG (`ApogeeReduction-airflow`) that ran the dailies
until 2026-05 — sbatch of run_all.sh from the sdssv checkout into
`/mnt/ceph/users/sdssv/work/daily/`. Its imports are Airflow-3-compatible
(verified under both 3.0.6 and the 3.3.1 deploy target — the
`airflow.operators.*`/`airflow.sensors.*` paths are provider shims, and
astropy/pytz/the slack provider are in the shared env). It is retired to
`dags/attic/` ONLY because this DAG supersedes it feature-for-feature:

| ar_main.py | apogee_daily |
|---|---|
| one DAG, lco→apo serial, `0 7 * * *` ET, max_active_runs=2 | one DAG, lco→apo serial, same schedule, **max_active_runs=1** since 2026-11-17 (see "Robustness") |
| `update.repo` (git status/log/remote report; status-only) | `update.repo` (same report + `git pull --ff-only` — deliberate NEW behavior, AKS 2026-09-03) |
| `update.sdsscore` (`update.sh`) | `update.sdsscore` (param-gated, default ON; + new-submodule bootstrap) |
| `update.sync_logs` (rsync both observatories' `.log.html`) | `update.sync_logs` (same, both observatories) |
| `setup.mjd` (`int(Time(data_interval_start).mjd) - 1`) | `setup.mjd` — VERBATIM arithmetic (stdlib mjd; conf `{"mjd": N}` override for past-day reruns) |
| `setup.date_mjd` (ET date − 2 d, pytz) | `setup.date_mjd` — VERBATIM arithmetic (zoneinfo) |
| `TransferFileSensor` (defined but never wired; ancient-date skip + transfer-done poke) | `<obs>.wait_for_raw` (wired as cheap insurance; md5sum marker; mjd < 59148 assumed complete; timeout → clean skip) |
| `<obs>.initial_notification` (public link interpolated **SLACK_TOKEN**) | `<obs>.start_notification` (link from **PUBLIC_URL_SLUG**; token/slug decoupled) |
| `<obs>.science` (`submit_and_wait`: sbatch `-vvv -D outdir --mail-*`, 5 s `squeue -j` polls, SLACK_CHANNEL exported) | `<obs>.slurm_submit` — same mechanics restored (sbatch prefix, 5 s polls, SLACK_CHANNEL env) + a final `sacct` verdict (new: a FAILED job fails the task); local mode replaces it with the split step chain |
| arMADGICS + workup ran unconditionally every night | `run_madgics` param default TRUE (only knob; no env fallback); + `spectra_workup` |
| per-obs success/failure Slack (`:picard_facepalm:`) | `daily_summary` census per obs + per-task `on_failure_callback` (`:picard_facepalm:` + per-task hints) |
| `completion_notification` ("Both observatories processed") | `completion_notification` (single, end of DAG; implemented via the N1 metrics table) |
| always posted to prod `#apogee-reduction-jl` | SAME default (prod channel); dev override via `AR_SLACK_CHANNEL` env / `slack_channel` conf |

New relative to ar_main.py: tunnel_check, split almanac/runlist with
exit-16 skip, the whole local mode, spectra_workup, N1 metrics,
heartbeat/watchdog/systemd, checkpoint-aware reruns, skipped-observatory
pass-through (no-data LCO no longer blocks APO).

## Robustness (2026-11-17 restart)

The dailies were down 2026-09-26 → 2026-11-17. The restart surfaced the
failure modes below; each now has a guard. All knobs are env variables
(set them in `airflow_env.sh`); a value that does not parse falls back to
the default with a warning in the DAG-processor log rather than breaking
the DAG import.

| guard | default | env knob | why this value |
|---|---|---|---|
| almanac `statement_timeout` (PGOPTIONS) | 5 h | `AR_ALMANAC_STATEMENT_TIMEOUT` | the slow step is one statement, the APOGEE-ID lookup: ~40 s warm, 11–40 min on September nights, ~3.5–4 h fully cold (after the 2026-11-16 DB restart). 5 h = worst cold case + 25 %, so a legitimate cold run finishes and only a truly stuck statement errors (with a PostgreSQL message in `almanac.log`) |
| almanac `execution_timeout` | 6 h | `AR_ALMANAC_EXECUTION_TIMEOUT` | statement budget + 1 h for the rest of almanac. Must exceed the statement timeout (enforced). Lower both (e.g. 45 min / 1 h) once the DB-side covering index or precomputed lookup lands |
| almanac retries | 1, after 10 min | `AR_ALMANAC_RETRIES`, `AR_ALMANAC_RETRY_DELAY` | covers a transient tunnel drop; the first attempt has partly warmed the DB cache |
| startup retries (`stale_run_guard`, `update.*`, `setup.*`) | 3, 60 s apart | `AR_STARTUP_RETRIES`, `AR_STARTUP_RETRY_DELAY` | runs launched ~20 s after service start failed after "Pre Execute", before any command ran (most likely before the API server was up). These tasks are idempotent |
| stale scheduled runs | 36 h | `AR_STALE_RUN_MAX_AGE_H` (conf `stale_run_max_age_h`; 0 disables) | the 2026-09-27 scheduled run resumed 51 days later next to the day's run. 36 h still lets a run delayed by up to a day go |
| `max_active_runs` | 1 (was 2) | — | serial chain + SQLite; two active runs doubled the slow DB lookup load. A run still going at the next 7 am makes the new run queue (and the stale guard fails it if it waits > 36 h) |

Durations take a unit (`ms`, `s`, `min`, `h`, `d`); bare numbers are
rejected because PostgreSQL reads them as ms and Airflow as s.

- **almanac output handling.** Skip only if the existing file holds
  `raw/<obs>/<mjd>/exposures` and `.../fibers`
  (`scripts/check_almanac_output.py`, run with the almanac env's python,
  which has h5py) — not merely because the file exists: the killed
  2026-11-17 run left a 96-byte stub that made every retry skip. almanac
  writes to `<file>.partial`, moved into place only on exit 0, so a killed
  or failed run never leaves a stub at the real path, and a failed
  `almanac_clobber` re-run keeps the previous good file. almanac is never
  pointed at an existing output: it opens its output in append mode and
  fails on a truncated stub, so on the old DAG even `almanac_clobber` could
  not recover a killed night. Old-layout files
  (pre-`raw/`, before 2026-09) count as absent and are rebuilt, which is
  what the current pipeline needs anyway.
- **"0 % CPU is not a hang".** Before almanac starts, the task log says
  that the APOGEE-ID lookup can take hours at 0 % CPU on a cold DB cache,
  with the active limits and the path of the full step log (step output
  only reaches the Airflow log as a tail at the end). Check
  `pg_stat_activity` for `DataFileRead` before killing anything. A line
  from *inside* almanac right before the lookup (worker side,
  `apogee.py`, before `create_apogee_id_lookup()`) belongs to almanac
  itself and is left to the almanac PR.
- **Manual runs without a data interval.** A trigger without
  `--logical-date` has no `data_interval_start` in Airflow 3;
  `setup.mjd`/`setup.date_mjd` fall back to `logical_date`, then
  `run_after`, then now, and with conf `mjd` the date label is derived
  from the MJD (`ar_common.date_label_from_mjd`, the same convention as the
  scheduled arithmetic, checked for every day of a year in the tests).
- **stale_run_guard** is the DAG's first task. Limitation: a run that had
  already passed it before an outage resumes where it stopped.
- **Re-running a failed night**: a fresh trigger, not clearing tasks —
  INSTALL.md, "Re-running a failed night".

## Tests

`airflow/tests/` (stdlib `unittest`; no pytest needed):

- `test_robustness_helpers.py` — duration parsing, timeout consistency,
  stale-run check, step announce, failure hints (any python3).
- `test_render_units.py` — unit templates + `render_systemd_units.sh`
  (any python3; uses `systemd-analyze verify` if present).
- `test_check_almanac_output.py` — the almanac output check, including
  the 96-byte stub (needs h5py: run with the almanac env's python).
- `test_dag.py` — DagBag import with the deployed Airflow, structure
  (guard is the only root, retries, timeouts, max_active_runs), every
  BashOperator template rendered and `bash -n`-checked, and the almanac
  step run end to end against a stand-in almanac (stub → rebuild, failure
  → no stub, good file reused) — needs Airflow, and h5py for the e2e part.

```bash
cd <repo>
# DAG + helper tests with the deployed Airflow, on a throw-away AIRFLOW_HOME
# (nothing touches the production metadata DB):
AIRFLOW_HOME=$(mktemp -d) AIRFLOW__CORE__LOAD_EXAMPLES=False \
PYTHONDONTWRITEBYTECODE=1 \
  /mnt/home/sdssv/uv_env/airflow_env/bin/python -m unittest discover -s airflow/tests -v
# almanac output check (h5py lives in the almanac env):
PYTHONDONTWRITEBYTECODE=1 /mnt/home/sdssv/uv_env/almanac_env/bin/python \
  -m unittest discover -s airflow/tests -p 'test_check_almanac_output.py' -v
```

CI runs the same suite (`.github/workflows/airflow-dags.yml`) with
apache-airflow 3.3.1 on Python 3.12 plus h5py.

## Notifications (N3)

- Every task failure fires `on_failure_callback` → Slack `chat.postMessage`
  via `SLACK_TOKEN` (stdlib urllib; no Airflow connection needed), with
  per-task hints (tunnel down, sdsscore stale, ...). DAG-level callback
  covers run failure/timeout.
- `daily_summary` posts a one-line census at DAG end: steps run/failed, wall
  time, warning counts by class, product counts, dashboard status.
- `slack_mode` param: `full` (production: step posts from the Julia code +
  alerts + summary), `summary_only` (exactly one message — used by smoke
  tests; SLACK_TOKEN is scrubbed from the Julia steps), `off`.
- Channel: production `C08B7FKMP16` by default (ar_main.py always posted to
  #apogee-reduction-jl); dev `C07KQ7BJY5P` via `AR_SLACK_CHANNEL` env or
  `slack_channel` conf when testing. The deployed airflow_env.sh has set prod
  explicitly since the 2026-09-06 promotion (INSTALL.md step 7.3).
- Outage alerts (`airflow-failure-notify.service`,
  `check_airflow_heartbeat.sh`) are on a **separate** knob,
  `AR_SLACK_ALERT_CHANNEL`, also defaulting to prod: demoting the routine
  posts to dev must not silently demote "the orchestrator is down".

## Metrics (N1 seed)

`metrics_append` (trigger_rule=all_done, so failures are recorded too)
appends to CSVs under `/mnt/ceph/users/sdssv/work/daily/metrics/`:

- `daily_metrics.csv` — one row per DAG run:
  `ts_utc, run_kind(daily|test|backfill), dag_id, run_id, tele, mjd, mode,
  workers, n_steps, n_steps_failed, total_wall_s, warn_no_wavecal,
  warn_no_relflux, warn_other, n_ar2D, n_ar2Dcal, n_ar1Dcal, n_ar1Dunical,
  n_plots, dashboard_ok`
- `step_metrics.csv` — one row per chain step:
  `ts_utc, run_id, tele, mjd, step, exit_code, wall_s`

Raw material: each local-mode step appends `step,rc,wall_s` to
`<outdir>/logs/<tele>_<mjd>/steps.csv` and keeps its full log next to it;
warning counts are parsed from those logs (classes match the O1/bulk census:
no-wavecal-fallback, no-relflux; everything else tallied as `warn_other`).
CSV over HDF5: append-safe, human-readable, trivial to ingest for the N2
drift dashboards; promote to HDF5/DB when N2 lands if needed.

## Heartbeat + watchdog

`apogee_heartbeat` touches
`/mnt/ceph/users/sdssv/work/daily/metrics/airflow_heartbeat.txt` every
10 min. `scripts/check_airflow_heartbeat.sh` (run from rusty scrontab —
different failure domain) alerts to Slack when the heartbeat or the daily
metrics table goes stale. systemd user units restart `airflow standalone` on
failure and alert when restarts are exhausted. Installation requires AKS
action — see INSTALL.md.

## Smoke-testing without spamming Slack

```bash
airflow dags trigger apogee_daily --conf '{
  "mjd": 61284, "mode": "local", "raw_nowait": true,
  "outdir": "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_02/o3_dag_test/",
  "slack_mode": "summary_only", "test_label": true,
  "slack_channel": "C07KQ7BJY5P",
  "run_kind": "test", "run_madgics": false, "workers": 12}'
```

produces at most ONE Slack message (the labeled summary) on the dev channel.
