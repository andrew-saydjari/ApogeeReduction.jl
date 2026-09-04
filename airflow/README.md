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
├── systemd/                  user units (prepared, NOT installed)
├── scripts/
│   ├── check_airflow_heartbeat.sh   cross-machine watchdog (scrontab)
│   └── airflow_env.sh.example       env template (secrets NOT committed)
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
# Regenerated from the DAG's actual edge wiring (62 tasks; edge-set
# sha256 275ec83d2857e1c7, identical under Airflow 3.0.6 and 3.3.1):
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
  >> lco.local.madgics_workup >> lco.local.spectra_workup
lco.local.spectra_workup >> {lco.join, lco.chain_status}
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
- **update.repo**: reports the production checkout state AND
  `git pull --ff-only` (DELIBERATE new behavior — ar_main.py's repo task
  was status-only); **update.sync_logs** rsyncs both observatories'
  `.log.html` pages into the public `apogee_logs` mirror.
- **update.sdsscore is default ON** (O1 punch list #2 — a fresh night after
  any gap dies on missing confSummary otherwise), with new-submodule
  bootstrap (--init + SSH→HTTPS rewrite).
- **tunnel_check** verifies the `mwm` ControlMaster (`ssh -O check`) and
  (re)adds the 63333→operations.sdss.org:5432 forward; if the master is down
  it fails with a notification saying only AKS can reopen it.
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
- **spectra_workup** (first-class chain step per AKS 2026-09-03; was a
  manual afterstep invoked by neither driver historically): calls the
  arMADGICS-side entrypoint `workup/run_workup.sh <rawdir> <redux> <outdir>`
  (PR #26, MPI tier default) under the same madgics gate. Local mode calls
  it explicitly; slurm mode inherits the identical step from run_all.sh.
  Existence-guarded (skips with a warning until PR #26 is in the checkout).
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
| one DAG, lco→apo serial, `0 7 * * *` ET, max_active_runs=2 | SAME: one DAG, lco→apo serial, same schedule, max_active_runs=2 (SQLite single-writer) |
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
  `slack_channel` conf when testing.

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
