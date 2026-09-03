# APOGEE daily pipeline on Airflow (cards O3 / N1 / N3)

Airflow layer for the daily reduction: per-observatory DAGs split at
`scripts/daily/run_all.sh`'s natural step boundaries, Slack notifications, a
heartbeat/watchdog layer, and the N1 per-night metrics seed.

```
airflow/
├── dags/
│   ├── ar_common.py          shared config/helpers (Slack, paths, metrics)
│   ├── apogee_daily.py       apogee_daily_apo + apogee_daily_lco
│   └── apogee_heartbeat.py   10-min liveness beacon
├── systemd/                  user units (prepared, NOT installed)
├── scripts/
│   ├── check_airflow_heartbeat.sh   cross-machine watchdog (scrontab)
│   └── airflow_env.sh.example       env template (secrets NOT committed)
├── INSTALL.md                deploy checklist (AKS actions)
└── README.md                 this file
```

## DAG structure (`apogee_daily_{apo,lco}`)

```
resolve_night → wait_for_raw → update_repo → update_sdsscore → sync_logs
   → tunnel_check → almanac → runlist → start_notification → select_mode
                                                                        │
              ┌──────────────── mode="local" (testing) ─────────────────┤
              │ p3d2d → quartz_flats(runlist→traces→extract→relflux)    │
              │       → dome_flats(...) → p2d1d_full → plots            │
              │       → dashboard → madgics_gate → madgics_pipeline     │
              │       → madgics_workup → spectra_workup                 │
              └─ mode="slurm" (DEFAULT): slurm_submit (sbatch run_all.sh) ┘
                                                                        │
                 join → metrics_append (N1) → daily_summary (N3)
                      → completion_rollup            + chain_status leaf
```

- **Night selection**: `mjd` param (`-1` = auto: `int(current UTC MJD) - 1`,
  i.e. the night that ended this morning).
- **wait_for_raw** gates reduction on raw-transfer COMPLETION (the
  `<mjd>.md5sum` marker), not just the date — ported from ar_main.py's
  `TransferFileSensor` semantics, including "ancient night (mjd < 59148) →
  assume transferred"; reschedule-mode pokes every 30 min, and a 12 h
  timeout skips the run cleanly (`soft_fail`).
- **update_repo / sync_logs**: ported from ar_main.py's `update` group —
  repo-state report (plus opt-in `pull_repo` fast-forward) and the rsync of
  the telescope's nightly `.log.html` pages into the public `apogee_logs`
  mirror.
- **update_sdsscore is default ON** (O1 punch list #2 — a fresh night after
  any gap dies on missing confSummary otherwise).
- **tunnel_check** verifies the `mwm` ControlMaster (`ssh -O check`) and
  (re)adds the 63333→operations.sdss.org:5432 forward; if the master is down
  it fails with a notification saying only AKS can reopen it.
- **runlist** maps `make_runlist_all.jl` exit 16 ("no exposures") to a clean
  Airflow skip.
- **Execution modes**: `slurm` (DEFAULT — production dailies as run
  historically; AKS launches the Airflow instance himself, so submission is
  human-initiated) submits the whole run_all.sh chain as one sbatch job with
  gentle 5-min single-job polling; `local` runs every step directly on the
  host (ccalin051) under `nice -n 10` with
  `AR_LOCAL_WORKERS`/`--workers_per_node` = `workers` param and all
  `SLURM_*` env scrubbed — used for testing (`AR_AIRFLOW_MODE=local` or
  conf `{"mode": "local"}`).
- **arMADGICS** is gated by the `run_madgics` param / `RUN_MADGICS` env
  (default true), mirroring run_all.sh arg 12.
- **spectra_workup** (first-class chain step per AKS 2026-09-03; was a
  manual afterstep invoked by neither driver historically): calls the
  arMADGICS-side entrypoint `workup/run_workup.sh <rawdir> <redux> <outdir>`
  (PR #26, MPI tier default) under the same madgics gate. Local mode calls
  it explicitly; slurm mode inherits the identical step from run_all.sh.
  Existence-guarded (skips with a warning until PR #26 is in the checkout).
- **dagrun_timeout=16 h** stands in for SLAs (removed in Airflow 3.0),
  together with the metrics-freshness check in the watchdog script; the
  `chain_status` leaf makes the DagRun state track the chain outcome even
  though metrics/summary always run.

## Relationship to the previous production DAG (`dags/ar_main.py`)

`ar_main.py` in the Airflow sandbox (`/mnt/home/sdssv/airflow/dags/`) is the
CCA-NATIVE production DAG (`ApogeeReduction-airflow`) that ran the dailies
until 2026-05 — sbatch of run_all.sh from the sdssv checkout into
`/mnt/ceph/users/sdssv/work/daily/`. Its imports are Airflow-3-compatible
(verified under 3.0.6 — the `airflow.operators.*`/`airflow.sensors.*` paths
are provider shims, and astropy/pytz/the slack provider are in the shared
env). It is retired to `dags/attic/` ONLY because these DAGs supersede it
feature-for-feature:

| ar_main.py | new DAGs |
|---|---|
| `update.repo` (git status/log/remote report) | `update_repo` (same report; `pull_repo=true` adds `git pull --ff-only`) |
| `update.sdsscore` (`update.sh`) | `update_sdsscore` (param-gated, default ON) |
| `update.sync_logs` (rsync both observatories' `.log.html`) | `sync_logs` (per-observatory half in each DAG) |
| `setup.mjd` (`int(Time(interval_start).mjd) - 1`) | `resolve_night` auto-mjd (`int(current MJD) - 1`, `mjd` param override) |
| `setup.date_mjd` (ET date − 2 d) | `resolve_night` `night_date` (MJD → civil date) |
| `TransferFileSensor` (defined but never wired; ancient-date skip + transfer-done poke) | `wait_for_raw` (actually wired; md5sum marker; mjd < 59148 assumed complete; timeout → clean skip) |
| `<obs>.initial_notification` (public link interpolated **SLACK_TOKEN**) | `start_notification` (link from **PUBLIC_URL_SLUG**; token/slug decoupled) |
| `<obs>.science` (`submit_and_wait`: sbatch + `squeue -j` poll every 5 s) | `slurm_submit` (sbatch `--parsable` + 300 s polls + `sacct` verdict); in local mode the split step chain replaces it |
| per-obs success/failure Slack (`:picard_facepalm:`) | `daily_summary` census per obs + per-task `on_failure_callback` (with `:picard_facepalm:` + per-task hints) |
| `completion_notification` ("Both observatories processed") | `completion_rollup` (posted by whichever obs DAG finishes second, judged from the N1 metrics table) |
| one DAG, lco→apo sequential, 7 am ET | two per-obs DAGs (9/10 am ET), independently parallel |

New relative to ar_main.py: tunnel_check, split almanac/runlist with
exit-16 skip, the whole local mode, N1 metrics, heartbeat/watchdog/systemd,
checkpoint-aware reruns.

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
- Channel: dev `C07KQ7BJY5P` by default; production `C08B7FKMP16` via
  `AR_SLACK_CHANNEL` (see INSTALL.md).

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
airflow dags trigger apogee_daily_apo --conf '{
  "mjd": 61284,
  "outdir": "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_02/o3_dag_test/",
  "slack_mode": "summary_only", "test_label": true,
  "run_kind": "test", "run_madgics": false, "workers": 12}'
```

produces at most ONE Slack message (the labeled summary) on the dev channel.
