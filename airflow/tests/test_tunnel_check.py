"""Tests for the DB tunnel check, its bounded self-heal, and the watchdog
state machine (ar_common.check_tunnel / watchdog_tick,
scripts/db_tunnel_watchdog.py). Stdlib only, fakes only: no psql, ssh,
systemctl or Slack is ever run."""

import json
import os
import shutil
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "dags"))
sys.path.insert(0, os.path.join(HERE, "..", "scripts"))

import ar_common as C  # noqa: E402
import db_tunnel_watchdog as W  # noqa: E402

ACTIVE = {"LoadState": "loaded", "ActiveState": "active",
          "SubState": "running", "Result": "success"}
FAILED_LIMIT = {"LoadState": "loaded", "ActiveState": "failed",
                "SubState": "failed", "Result": "start-limit-hit"}
CYCLING = {"LoadState": "loaded", "ActiveState": "activating",
           "SubState": "auto-restart", "Result": "exit-code"}
INACTIVE = {"LoadState": "loaded", "ActiveState": "inactive",
            "SubState": "dead", "Result": "success"}

OK = (0, "1")
REFUSED = (2, "psql: error: connection to server at \"localhost\" (::1), "
              "port 63333 failed: Connection refused")
CLOSED_BY_FAR_END = (2, "psql: error: server closed the connection "
                        "unexpectedly")
AUTH = (2, "psql: error: FATAL:  password authentication failed for user "
           "\"sdss\"")


class FakeEnv:
    """Scripted world. `db` is the psql answer, `port` whether something
    listens, `unit` the systemctl show answer; `after_restart` (if given)
    is applied to those when `systemctl --user restart` is called."""

    def __init__(self, db=OK, port=True, unit=ACTIVE, psql=True,
                 after_restart=None, t0=1_000_000.0, last_restart=None):
        self.db, self.port, self.unit, self.psql = db, port, unit, psql
        self.after_restart = after_restart
        self.t = t0
        self._last_restart = last_restart
        self.calls = []

    # TunnelEnv interface
    def which(self, name):
        return "/usr/bin/psql" if self.psql else None

    def run(self, argv, timeout, extra_env=None):
        self.calls.append(list(argv))
        if argv[0] == C.PSQL_BIN:
            assert extra_env["PGCONNECT_TIMEOUT"] == "10"
            assert "statement_timeout" in extra_env["PGOPTIONS"]
            return self.db
        if argv[:3] == ["systemctl", "--user", "show"]:
            if isinstance(self.unit, tuple):
                return self.unit
            return 0, "\n".join(f"{k}={v}" for k, v in self.unit.items())
        if argv[:3] == ["systemctl", "--user", "restart"]:
            if self.after_restart:
                for k, v in self.after_restart.items():
                    setattr(self, k, v)
            return 0, ""
        raise AssertionError(f"unexpected command {argv}")

    def port_open(self, port):
        return self.port

    def sleep(self, s):
        self.t += s

    def now(self):
        return self.t

    def last_restart(self):
        return self._last_restart

    def mark_restart(self):
        self._last_restart = self.t

    # helpers
    def restarts(self):
        return [c for c in self.calls if c[:3] == ["systemctl", "--user",
                                                    "restart"]]

    def ssh_or_sshpass(self):
        return [c for c in self.calls
                if any("ssh" in a.split("/")[-1] for a in c)]


def quiet(*_a, **_k):
    pass


class CheckTunnel(unittest.TestCase):
    def check(self, env, **kw):
        r = C.check_tunnel(env, log=quiet, **kw)
        self.assertEqual(env.ssh_or_sshpass(), [])  # never ssh / sshpass
        self.assertNotIn("reset-failed", json.dumps(env.calls))
        return r

    def test_query_ok_passes_without_touching_ssh_or_systemd(self):
        env = FakeEnv()
        r = self.check(env)
        self.assertTrue(r.ok)
        self.assertEqual(r.verdict, C.TV_OK)
        self.assertEqual([c[0] for c in env.calls], [C.PSQL_BIN])

    def test_query_fails_one_restart_then_pass(self):
        env = FakeEnv(db=REFUSED, port=False,
                      after_restart={"db": OK, "port": True})
        r = self.check(env)
        self.assertTrue(r.ok, r.message)
        self.assertTrue(r.restarted)
        self.assertEqual(env.restarts(),
                         [["systemctl", "--user", "restart",
                           C.TUNNEL_UNIT]])
        self.assertIn("after restarting", r.message)

    def test_never_more_than_one_restart(self):
        env = FakeEnv(db=REFUSED, port=False)  # restart does not help
        r = self.check(env)
        self.assertFalse(r.ok)
        self.assertEqual(len(env.restarts()), 1)
        self.assertEqual(r.verdict, C.TV_PORT_CLOSED)
        self.assertTrue(r.retryable)
        self.assertIn("yes, without success", r.message)
        self.assertIn("ss -ltnp", r.message)

    def test_start_limited_unit_is_not_restarted(self):
        env = FakeEnv(db=REFUSED, port=False, unit=FAILED_LIMIT)
        r = self.check(env)
        self.assertFalse(r.ok)
        self.assertEqual(env.restarts(), [])
        self.assertEqual(r.verdict, C.TV_UNIT_FAILED)
        self.assertFalse(r.retryable)
        self.assertIn("start-limit-hit", r.message)
        self.assertIn("password", r.message)
        self.assertIn(f"systemctl --user reset-failed {C.TUNNEL_UNIT}",
                      r.message)  # told to the human, never run

    def test_cycling_unit_is_left_to_systemd(self):
        env = FakeEnv(db=REFUSED, port=False, unit=CYCLING)
        r = self.check(env)
        self.assertEqual(env.restarts(), [])
        self.assertEqual(r.verdict, C.TV_UNIT_CYCLING)
        self.assertTrue(r.retryable)

    def test_port_open_but_db_rejects_login(self):
        env = FakeEnv(db=AUTH, port=True)
        r = self.check(env)
        self.assertFalse(r.ok)
        self.assertEqual(r.verdict, C.TV_DB_REJECTED)
        self.assertEqual(env.restarts(), [])
        self.assertFalse(r.retryable)
        self.assertIn("The tunnel is fine", r.message)
        self.assertIn("PGPASSFILE", r.message)
        self.assertIn("password authentication failed", r.message)

    def test_port_open_but_query_failing(self):
        env = FakeEnv(db=CLOSED_BY_FAR_END, port=True)
        r = self.check(env)
        self.assertFalse(r.ok)
        self.assertEqual(r.verdict, C.TV_QUERY_FAILED)
        self.assertIn("accepts connections but select 1 fails", r.message)
        self.assertEqual(len(env.restarts()), 1)  # the channel may be dead
        self.assertTrue(r.retryable)

    def test_psql_missing(self):
        env = FakeEnv(psql=False)
        r = self.check(env)
        self.assertEqual(r.verdict, C.TV_PSQL_MISSING)
        self.assertEqual(env.calls, [])
        self.assertFalse(r.retryable)
        self.assertIn("AR_PSQL", r.message)

    def test_systemctl_unusable(self):
        env = FakeEnv(db=REFUSED, port=False,
                      unit=(1, "Failed to connect to bus"))
        r = self.check(env)
        self.assertEqual(r.verdict, C.TV_UNIT_UNKNOWN)
        self.assertEqual(env.restarts(), [])
        self.assertIn("Failed to connect to bus", r.message)

    def test_unit_not_installed(self):
        env = FakeEnv(db=REFUSED, port=False,
                      unit={"LoadState": "not-found",
                            "ActiveState": "inactive", "SubState": "dead",
                            "Result": "success"})
        r = self.check(env)
        self.assertEqual(r.verdict, C.TV_UNIT_UNKNOWN)
        self.assertEqual(env.restarts(), [])

    def test_inactive_unit_is_started_once(self):
        env = FakeEnv(db=REFUSED, port=False, unit=INACTIVE)
        r = self.check(env)
        self.assertEqual(len(env.restarts()), 1)
        self.assertEqual(r.verdict, C.TV_UNIT_INACTIVE)

    def test_restart_disallowed(self):
        env = FakeEnv(db=REFUSED, port=False)
        r = self.check(env, allow_restart=False)
        self.assertEqual(env.restarts(), [])
        self.assertFalse(r.ok)
        self.assertIn("Restarted", r.message)

    def test_shared_cooldown(self):
        # restarted 5 min ago by the task or the watchdog: wait, don't restart
        env = FakeEnv(db=REFUSED, port=False, t0=10_000.0,
                      last_restart=10_000.0 - 300)
        self.check(env)
        self.assertEqual(env.restarts(), [])
        env = FakeEnv(db=REFUSED, port=False, t0=10_000.0,
                      last_restart=10_000.0 - C.TUNNEL_RESTART_COOLDOWN_S)
        self.check(env)
        self.assertEqual(len(env.restarts()), 1)

    def test_unit_fails_during_wait(self):
        env = FakeEnv(db=REFUSED, port=False,
                      after_restart={"unit": FAILED_LIMIT})
        r = self.check(env)
        self.assertEqual(r.verdict, C.TV_UNIT_FAILED)
        self.assertEqual(len(env.restarts()), 1)

    def test_wait_is_bounded(self):
        env = FakeEnv(db=REFUSED, port=False)
        t0 = env.t
        self.check(env)
        self.assertLessEqual(env.t - t0, C.TUNNEL_RESTART_WAIT_S)


class Watchdog(unittest.TestCase):
    def tick(self, state, env):
        new, alerts = C.watchdog_tick(state, env, log=quiet)
        self.assertEqual(env.ssh_or_sshpass(), [])
        self.assertNotIn("reset-failed", json.dumps(env.calls))
        return new, alerts

    def test_healthy_stays_silent(self):
        st, alerts = self.tick({}, FakeEnv())
        self.assertEqual(st["status"], "healthy")
        self.assertEqual(alerts, [])
        st, alerts = self.tick(st, FakeEnv())
        self.assertEqual(alerts, [])

    def test_self_heal_in_one_tick_is_one_message(self):
        env = FakeEnv(db=REFUSED, port=False,
                      after_restart={"db": OK, "port": True})
        st, alerts = self.tick({"status": "healthy"}, env)
        self.assertEqual(st["status"], "healthy")
        self.assertEqual(len(alerts), 1)
        self.assertIn("restarted", alerts[0])
        self.assertEqual(st["last_watchdog_restart"], env.t - 5)

    def test_broken_episode_bounded_and_alerts_only_on_transitions(self):
        t = 1_000_000.0
        # tick 1: breaks; one restart that does not help -> broken, 1 alert
        env = FakeEnv(db=REFUSED, port=False, t0=t)
        st, alerts = self.tick({"status": "healthy"}, env)
        self.assertEqual(st["status"], "broken")
        self.assertEqual(len(env.restarts()), 1)
        self.assertEqual(len(alerts), 1)
        self.assertIn("broken", alerts[0])
        # tick 2 (10 min later): still broken -> no second restart, give up
        env = FakeEnv(db=REFUSED, port=False, t0=t + 600)
        st, alerts = self.tick(st, env)
        self.assertEqual(env.restarts(), [])
        self.assertEqual(st["status"], "gave_up")
        self.assertEqual(len(alerts), 1)
        self.assertIn("needs a human", alerts[0])
        # ticks 3..: silent, never restarts
        for k in range(3, 12):
            env = FakeEnv(db=REFUSED, port=False, t0=t + 600 * k)
            st, alerts = self.tick(st, env)
            self.assertEqual((env.restarts(), alerts), ([], []))
            self.assertEqual(st["status"], "gave_up")
        # human fixes it -> one recovery message
        env = FakeEnv(t0=t + 9000)
        st, alerts = self.tick(st, env)
        self.assertEqual(st["status"], "healthy")
        self.assertEqual(len(alerts), 1)
        self.assertIn("healthy again", alerts[0])

    def test_start_limited_gives_up_at_once_without_restart(self):
        env = FakeEnv(db=REFUSED, port=False, unit=FAILED_LIMIT)
        st, alerts = self.tick({"status": "healthy"}, env)
        self.assertEqual(st["status"], "gave_up")
        self.assertEqual(env.restarts(), [])
        self.assertEqual(len(alerts), 1)
        self.assertIn("needs a human", alerts[0])
        self.assertIn("start-limit-hit", alerts[0])

    def test_cycling_stays_broken_until_failed(self):
        env = FakeEnv(db=REFUSED, port=False, unit=CYCLING)
        st, alerts = self.tick({"status": "healthy"}, env)
        self.assertEqual((st["status"], len(alerts)), ("broken", 1))
        self.assertEqual(env.restarts(), [])
        env = FakeEnv(db=REFUSED, port=False, unit=CYCLING, t0=env.t + 600)
        st, alerts = self.tick(st, env)
        self.assertEqual((st["status"], alerts), ("broken", []))
        env = FakeEnv(db=REFUSED, port=False, unit=FAILED_LIMIT,
                      t0=env.t + 600)
        st, alerts = self.tick(st, env)
        self.assertEqual(st["status"], "gave_up")
        self.assertEqual(len(alerts), 1)
        self.assertEqual(env.restarts(), [])

    def test_hourly_cooldown_across_episodes(self):
        t = 1_000_000.0
        env = FakeEnv(db=REFUSED, port=False, t0=t,
                      after_restart={"db": OK, "port": True})
        st, _ = self.tick({}, env)                    # episode 1: healed
        self.assertEqual(len(env.restarts()), 1)
        # episode 2, 20 min later: within the hour -> no restart, broken
        env = FakeEnv(db=REFUSED, port=False, t0=t + 1200,
                      after_restart={"db": OK, "port": True})
        st, alerts = self.tick(st, env)
        self.assertEqual(env.restarts(), [])
        self.assertEqual(st["status"], "broken")
        self.assertEqual(len(alerts), 1)
        # an hour after the first restart: allowed again (still episode 2,
        # whose restart budget is unused)
        env = FakeEnv(db=REFUSED, port=False, t0=t + 3700,
                      after_restart={"db": OK, "port": True})
        st, alerts = self.tick(st, env)
        self.assertEqual(len(env.restarts()), 1)
        self.assertEqual(st["status"], "healthy")
        self.assertIn("healthy again", alerts[0])

    def test_db_rejection_needs_a_human(self):
        env = FakeEnv(db=AUTH, port=True)
        st, alerts = self.tick({}, env)
        self.assertEqual(st["status"], "gave_up")
        self.assertEqual(env.restarts(), [])
        self.assertIn("DB rejected", alerts[0])


class WatchdogScript(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.state = os.path.join(self.d, "state.json")
        self.posted = []

    def tearDown(self):
        shutil.rmtree(self.d, ignore_errors=True)

    def run_main(self, env, *args):
        return W.main(["--state", self.state, *args], env=env,
                      post=self.posted.append)

    def test_state_persisted_and_alerts_posted_once(self):
        self.assertEqual(self.run_main(FakeEnv(db=REFUSED, port=False)), 0)
        with open(self.state) as f:
            self.assertEqual(json.load(f)["status"], "broken")
        self.assertEqual(len(self.posted), 1)
        self.assertIn("APOGEE DB tunnel watchdog", self.posted[0])
        self.run_main(FakeEnv(db=REFUSED, port=False, unit=CYCLING))
        self.assertEqual(len(self.posted), 1)
        self.run_main(FakeEnv())
        self.assertEqual(len(self.posted), 2)
        with open(self.state) as f:
            self.assertEqual(json.load(f)["status"], "healthy")

    def test_probe_only_never_restarts_or_posts(self):
        env = FakeEnv(db=REFUSED, port=False)
        self.assertEqual(self.run_main(env, "--probe-only"), 1)
        self.assertEqual(env.restarts(), [])
        self.assertEqual(self.posted, [])
        self.assertFalse(os.path.exists(self.state))
        self.assertEqual(self.run_main(FakeEnv(), "--probe-only"), 0)

    def test_corrupt_state_is_treated_as_fresh(self):
        with open(self.state, "w") as f:
            f.write("{not json")
        self.assertEqual(self.run_main(FakeEnv()), 0)
        self.assertEqual(self.posted, [])

    def test_default_state_path_is_host_local(self):
        old = dict(os.environ)
        try:
            os.environ.pop("AR_TUNNEL_WATCHDOG_STATE", None)
            os.environ["XDG_RUNTIME_DIR"] = "/run/user/1234"
            self.assertEqual(W.default_state_path(),
                             "/run/user/1234/ar-db-tunnel-watchdog.json")
        finally:
            os.environ.clear()
            os.environ.update(old)


class RealEnvStamp(unittest.TestCase):
    def test_stamp_roundtrip(self):
        d = tempfile.mkdtemp()
        try:
            env = C.TunnelEnv(stamp_path=os.path.join(d, "stamp"))
            self.assertIsNone(env.last_restart())
            env.mark_restart()
            self.assertAlmostEqual(env.last_restart(), env.now(), delta=5)
        finally:
            shutil.rmtree(d, ignore_errors=True)


if __name__ == "__main__":
    unittest.main()
