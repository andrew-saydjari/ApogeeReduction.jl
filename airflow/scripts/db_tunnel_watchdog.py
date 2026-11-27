#!/usr/bin/env python3
"""DB tunnel watchdog: one tick, run every 10 min by ar-db-tunnel-watchdog.timer
(airflow/systemd/ar-db-tunnel-watchdog.*.in; install: airflow/INSTALL.md
step 4).

Each tick runs the same `select 1` through localhost:<AR_DB_FORWARD_PORT> as
the DAG's tunnel_check (ar_common.check_tunnel), restarts the tunnel unit
under ar_common.watchdog_tick's bounds (at most one restart per failure
episode and one per hour; never when the unit is failed/start-limited; never
reset-failed) and posts to Slack only on state changes, through
ar_common.slack_post (SLACK_TOKEN from the environment; channel
AR_SLACK_ALERT_CHANNEL, default prod). The environment is the airflow
EnvironmentFile, so psql sees the same PATH, PGPASSFILE and DB identity as
the DAG. Stdlib only.

State: a small JSON file, by default $XDG_RUNTIME_DIR/ar-db-tunnel-watchdog.json
(host-local; override with AR_TUNNEL_WATCHDOG_STATE or --state).

    db_tunnel_watchdog.py               # one tick (what the timer runs)
    db_tunnel_watchdog.py --probe-only  # health check: probe + saved state,
                                        # no restart, no Slack; exit 0 iff OK
"""

import argparse
import json
import os
import socket
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "dags"))

import ar_common as C  # noqa: E402


def default_state_path() -> str:
    explicit = os.environ.get("AR_TUNNEL_WATCHDOG_STATE")
    if explicit:
        return explicit
    base = (os.environ.get("XDG_RUNTIME_DIR") or
            os.environ.get("AIRFLOW_HOME") or os.path.expanduser("~"))
    return os.path.join(base, "ar-db-tunnel-watchdog.json")


def load_state(path: str) -> dict:
    try:
        with open(path) as f:
            st = json.load(f)
        return st if isinstance(st, dict) else {}
    except (OSError, ValueError):
        return {}


def save_state(path: str, state: dict) -> None:
    tmp = f"{path}.tmp{os.getpid()}"
    with open(tmp, "w") as f:
        json.dump(state, f, indent=1, sort_keys=True)
        f.write("\n")
    os.replace(tmp, path)


def main(argv=None, env=None, post=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--state", default=None,
                    help="state file [$XDG_RUNTIME_DIR/ar-db-tunnel-watchdog.json]")
    ap.add_argument("--probe-only", action="store_true",
                    help="probe and print the saved state; no restart, no Slack")
    args = ap.parse_args(argv)
    path = args.state or default_state_path()
    env = env or C.TunnelEnv()
    post = post or (lambda text: C.slack_post(
        text, channel=C.AR_SLACK_ALERT_CHANNEL))

    if args.probe_only:
        verdict, detail = C.probe_db(env)
        print(f"probe: {verdict} ({detail}); unit {C.TUNNEL_UNIT}: "
              f"{C.unit_state(env)}")
        print(f"watchdog state ({path}): {load_state(path) or 'none yet'}")
        return 0 if verdict == C.TV_OK else 1

    old = load_state(path)
    new, alerts = C.watchdog_tick(old, env)
    try:
        save_state(path, new)
    except OSError as e:
        print(f"WARNING: cannot save state to {path}: {e}")
    host = socket.gethostname().split(".")[0]
    for text in alerts:
        post(f"[{host}] APOGEE DB tunnel watchdog: {text}")
    print(f"{old.get('status', 'healthy')} -> {new['status']} "
          f"[{new['verdict']}], restarts this episode "
          f"{new['episode_restarts']}, alerts {len(alerts)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
