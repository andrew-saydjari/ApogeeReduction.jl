"""Unit tests for the pure-Python robustness helpers in ar_common.py.

Stdlib only (ar_common imports nothing from Airflow), so this runs under any
python3:  python3 -m unittest discover -s airflow/tests -v
"""

import os
import sys
import unittest
from datetime import datetime, timedelta, timezone
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "dags"))

import ar_common as C  # noqa: E402


class ParseDuration(unittest.TestCase):
    def test_units(self):
        self.assertEqual(C.parse_duration_s("5h"), 18000)
        self.assertEqual(C.parse_duration_s("300min"), 18000)
        self.assertEqual(C.parse_duration_s("90s"), 90)
        self.assertEqual(C.parse_duration_s("1d"), 86400)
        self.assertEqual(C.parse_duration_s("1500ms"), 2)
        self.assertEqual(C.parse_duration_s(" 2.5 h "), 9000)

    def test_rejects_bare_numbers_zero_and_junk(self):
        # bare numbers are ambiguous (PostgreSQL: ms, Airflow: s)
        for bad in ("3600", "0h", "0", "", "5 hours", "-1h", "h", "1e3s",
                    "500ms", None):
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                C.parse_duration_s(bad)

    def test_pgoptions_fragment_is_milliseconds(self):
        self.assertEqual(C.pg_statement_timeout_option("5h"),
                         "-c statement_timeout=18000000")
        self.assertEqual(C.pg_statement_timeout_option("15s"),
                         "-c statement_timeout=15000")


class AlmanacTimeouts(unittest.TestCase):
    def test_defaults_are_consistent_and_match_doc(self):
        self.assertEqual(C.ALMANAC_STATEMENT_TIMEOUT_S, 5 * 3600)
        self.assertEqual(C.ALMANAC_EXECUTION_TIMEOUT_S, 6 * 3600)
        self.assertEqual(C.ALMANAC_PGOPTIONS,
                         "-c statement_timeout=18000000")
        # above the measured cold lookup (~4 h)
        self.assertGreater(C.ALMANAC_STATEMENT_TIMEOUT_S, 4 * 3600)
        C.check_almanac_timeouts(C.ALMANAC_STATEMENT_TIMEOUT_S,
                                 C.ALMANAC_EXECUTION_TIMEOUT_S)

    def test_statement_must_be_shorter_than_execution(self):
        with self.assertRaises(ValueError):
            C.check_almanac_timeouts(3600, 3600)
        with self.assertRaises(ValueError):
            C.check_almanac_timeouts(7200, 3600)
        C.check_almanac_timeouts(3599, 3600)

    def test_env_overrides_and_fallbacks(self):
        env = {"AR_ALMANAC_STATEMENT_TIMEOUT": "45min",
               "AR_ALMANAC_EXECUTION_TIMEOUT": "1h"}
        with mock.patch.dict(os.environ, env):
            self.assertEqual(C._almanac_timeouts(), (2700, 3600))
        # inconsistent pair -> both defaults (never a broken DAG import)
        env = {"AR_ALMANAC_STATEMENT_TIMEOUT": "2h",
               "AR_ALMANAC_EXECUTION_TIMEOUT": "1h"}
        with mock.patch.dict(os.environ, env):
            self.assertEqual(C._almanac_timeouts(), (18000, 21600))
        # unparsable value -> that value's default
        env = {"AR_ALMANAC_STATEMENT_TIMEOUT": "3600",
               "AR_ALMANAC_EXECUTION_TIMEOUT": "8h"}
        with mock.patch.dict(os.environ, env):
            self.assertEqual(C._almanac_timeouts(), (18000, 28800))

    def test_env_int(self):
        with mock.patch.dict(os.environ, {"X_TEST_INT": "4"}):
            self.assertEqual(C._env_int("X_TEST_INT", 1), 4)
        with mock.patch.dict(os.environ, {"X_TEST_INT": "-2"}):
            self.assertEqual(C._env_int("X_TEST_INT", 1), 1)
        with mock.patch.dict(os.environ, {"X_TEST_INT": "two"}):
            self.assertEqual(C._env_int("X_TEST_INT", 1), 1)


class StaleRun(unittest.TestCase):
    NOW = datetime(2026, 11, 17, 21, 32, tzinfo=timezone.utc)

    def test_the_2026_09_27_run_is_stale(self):
        logical = datetime(2026, 9, 27, 11, 0, tzinfo=timezone.utc)
        reason = C.stale_run_reason("scheduled", logical, self.NOW, 36)
        self.assertIsNotNone(reason)
        self.assertIn("stale scheduled run", reason)
        self.assertIn("almanac_clobber", reason)

    def test_todays_late_run_is_not_stale(self):
        logical = datetime(2026, 11, 17, 12, 0, tzinfo=timezone.utc)
        self.assertIsNone(
            C.stale_run_reason("scheduled", logical, self.NOW, 36))

    def test_boundary(self):
        self.assertIsNone(C.stale_run_reason(
            "scheduled", self.NOW - timedelta(hours=36), self.NOW, 36))
        self.assertIsNotNone(C.stale_run_reason(
            "scheduled", self.NOW - timedelta(hours=36, seconds=1),
            self.NOW, 36))

    def test_manual_and_backfill_runs_are_never_stale(self):
        old = datetime(2025, 1, 1, tzinfo=timezone.utc)
        for rt in ("manual", "backfill", "asset_triggered", None, ""):
            with self.subTest(rt=rt):
                self.assertIsNone(C.stale_run_reason(rt, old, self.NOW, 36))

    def test_enum_like_run_type(self):
        class RT:  # stands in for airflow's DagRunType(str, Enum)
            value = "scheduled"
        old = datetime(2025, 1, 1, tzinfo=timezone.utc)
        self.assertIsNotNone(C.stale_run_reason(RT(), old, self.NOW, 36))

    def test_disabled_or_no_date(self):
        old = datetime(2025, 1, 1, tzinfo=timezone.utc)
        self.assertIsNone(C.stale_run_reason("scheduled", old, self.NOW, 0))
        self.assertIsNone(C.stale_run_reason("scheduled", None, self.NOW, 36))


class ManualRunDates(unittest.TestCase):
    """setup.mjd / setup.date_mjd for runs without a data interval
    (2026-11-17: a trigger without --logical-date -> KeyError
    'data_interval_start' in setup.date_mjd)."""

    T = datetime(2026, 11, 17, 12, 0, tzinfo=timezone.utc)  # 7 am ET

    def test_reference_time_fallbacks(self):
        class Run:
            run_after = datetime(2026, 11, 18, 1, tzinfo=timezone.utc)
            logical_date = None
        self.assertEqual(C.run_reference_time(
            {"data_interval_start": self.T, "logical_date": None}),
            (self.T, "data_interval_start"))
        self.assertEqual(C.run_reference_time({"logical_date": self.T}),
                         (self.T, "logical_date"))
        self.assertEqual(C.run_reference_time(
            {"data_interval_start": None, "dag_run": Run()}),
            (Run.run_after, "dag_run.run_after"))
        now = datetime(2026, 11, 18, 2, tzinfo=timezone.utc)
        self.assertEqual(C.run_reference_time({}, now=now), (now, "now"))

    def test_label_from_mjd_matches_the_scheduled_arithmetic(self):
        # the run of 2026-11-17 processed mjd 61360, labelled 2026-11-15
        label, _ = C.compute_date_label({"data_interval_start": self.T,
                                         "params": {"mjd": -1}})
        self.assertEqual(label, "2026-11-15")
        self.assertEqual(C.date_label_from_mjd(61360), "2026-11-15")
        # and for every 7 am ET trigger over a year (DST both ways)
        from zoneinfo import ZoneInfo
        et = ZoneInfo("America/New_York")
        day = datetime(2026, 1, 1, 7, tzinfo=et)
        for i in range(366):
            t = (day + timedelta(days=i)).replace(hour=7)
            with self.subTest(t=t):
                sched, _ = C.compute_date_label({"data_interval_start": t,
                                                 "params": {}})
                mjd = int(C.mjd_of(t)) - 1
                self.assertEqual(sched, C.date_label_from_mjd(mjd))

    def test_manual_run_with_mjd_and_no_data_interval(self):
        ctx = {"params": {"mjd": 61360, "almanac_clobber": True}}
        label, src = C.compute_date_label(ctx)
        self.assertEqual(label, "2026-11-15")
        self.assertIn("param mjd=61360", src)

    def test_manual_run_without_mjd_or_data_interval(self):
        now = datetime(2026, 11, 18, 2, tzinfo=timezone.utc)
        label, src = C.compute_date_label({"params": {"mjd": -1}}, now=now)
        self.assertEqual(label, "2026-11-15")  # 2026-11-17 21:00 ET - 2 d
        self.assertTrue(src.startswith("now="))


class StepCmdAnnounce(unittest.TestCase):
    def test_announce_goes_before_the_redirected_command(self):
        cmd = C.step_cmd("almanac", "true", "lco", announce="hello")
        self.assertLess(cmd.index("almanac: hello"),
                        cmd.index('true > "$LOGDIR/almanac.log"'))
        self.assertNotIn("hello", C.step_cmd("almanac", "true", "lco"))

    def test_announce_rejects_double_quotes(self):
        with self.assertRaises(ValueError):
            C.step_cmd("almanac", "true", "lco", announce='say "hi"')


class Hints(unittest.TestCase):
    def test_tunnel_hint_names_the_unit_not_a_manual_command(self):
        hint = C.FAILURE_HINTS["tunnel_check"]
        self.assertIn(C.TUNNEL_UNIT, hint)
        self.assertNotIn("sshpass", hint)
        self.assertNotIn("pass_file", hint)


if __name__ == "__main__":
    unittest.main()
