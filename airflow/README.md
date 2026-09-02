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
resolve_night → update_sdsscore → tunnel_check → almanac → runlist → select_mode
                                                                        │
              ┌─────────────────── mode="local" (default) ──────────────┤
              │ p3d2d → quartz_flats(runlist→traces→extract→relflux)    │
              │       → dome_flats(...) → p2d1d_full → plots            │
              │       → dashboard → madgics_gate → madgics_pipeline     │
              │       → madgics_workup                                  │
              └── mode="slurm": slurm_submit (one sbatch of run_all.sh) ┘
                                                                        │
                                join → metrics_append (N1) → daily_summary (N3)
```

- **Night selection**: `mjd` param (`-1` = auto: `int(current UTC MJD) - 1`,
  i.e. the night that ended this morning). `resolve_night` skips the run
  cleanly when the raw mirror has no frames yet or the transfer lacks its
  `.md5sum` completion marker.
- **update_sdsscore is default ON** (O1 punch list #2 — a fresh night after
  any gap dies on missing confSummary otherwise).
- **tunnel_check** verifies the `mwm` ControlMaster (`ssh -O check`) and
  (re)adds the 63333→operations.sdss.org:5432 forward; if the master is down
  it fails with a notification saying only AKS can reopen it.
- **runlist** maps `make_runlist_all.jl` exit 16 ("no exposures") to a clean
  Airflow skip.
- **Execution modes**: `local` runs every step directly on the host
  (ccalin051) under `nice -n 10` with `AR_LOCAL_WORKERS`/`--workers_per_node`
  = `workers` param and all `SLURM_*` env scrubbed; `slurm` submits the whole
  run_all.sh chain as one job (written, disabled by default — AKS's call).
- **arMADGICS** is gated by the `run_madgics` param / `RUN_MADGICS` env
  (default true), mirroring run_all.sh arg 12.
- **dagrun_timeout=16 h** stands in for SLAs (removed in Airflow 3.0),
  together with the metrics-freshness check in the watchdog script.

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
