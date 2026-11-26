"""Unit tests for the Slurm job verdict helpers in ar_common.py
(parse_sacct_state / wait_for_sacct_terminal / slurm_job_verdict).

Stdlib only; no Slurm command is ever run (the runner and sleep are fakes):
  python3 -m unittest discover -s airflow/tests -v
"""

import os
import sys
import unittest
from types import SimpleNamespace

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "dags"))

import ar_common as C  # noqa: E402


class FakeSacct:
    """Returns the scripted sacct outputs in order (the last one repeats)."""

    def __init__(self, outputs, returncodes=None):
        self.outputs = list(outputs)
        self.returncodes = list(returncodes or [0] * len(self.outputs))
        self.calls = []

    def __call__(self, argv):
        i = min(len(self.calls), len(self.outputs) - 1)
        self.calls.append(list(argv))
        return SimpleNamespace(returncode=self.returncodes[i],
                               stdout=self.outputs[i], stderr="boom")


class FakeSleep:
    def __init__(self):
        self.calls = []

    def __call__(self, s):
        self.calls.append(s)


def _quiet(*_a, **_k):
    pass


def verdict(outputs, **kw):
    run, sleep = FakeSacct(outputs, kw.pop("returncodes", None)), FakeSleep()
    try:
        return C.slurm_job_verdict("7548417", runner=run, sleep=sleep,
                                   log=_quiet, **kw), run, sleep, None
    except RuntimeError as e:
        return None, run, sleep, e


class ParseSacctState(unittest.TestCase):
    def test_forms(self):
        cases = {
            "COMPLETED|0:0\n": ("COMPLETED", "0:0"),
            "  COMPLETED      0:0 \n": ("COMPLETED", "0:0"),
            "CANCELLED by 1234|0:15": ("CANCELLED", "0:15"),
            "CANCELLED by 1234   0:15": ("CANCELLED", "0:15"),
            "CANCELLED+   0:15": ("CANCELLED", "0:15"),
            "FAILED|1:0": ("FAILED", "1:0"),
            "OUT_OF_MEMORY|0:125": ("OUT_OF_MEMORY", "0:125"),
            "RUNNING|0:0": ("RUNNING", "0:0"),
            "": (None, None),
            "\n  \n": (None, None),
            "REQUEUED|0:0\nCOMPLETED|0:0\n": ("COMPLETED", "0:0"),
        }
        for out, want in cases.items():
            with self.subTest(out=out):
                self.assertEqual(C.parse_sacct_state(out), want)


class SlurmJobVerdict(unittest.TestCase):
    def test_lagging_accounting_then_completed_succeeds(self):
        jid, run, sleep, err = verdict(
            ["RUNNING|0:0", "RUNNING|0:0", "COMPLETED|0:0"])
        self.assertIsNone(err)
        self.assertEqual(jid, "7548417")
        self.assertEqual(len(run.calls), 3)
        self.assertEqual(len(sleep.calls), 2)

    def test_immediately_completed_needs_no_sleep(self):
        jid, run, sleep, err = verdict(["COMPLETED|0:0"])
        self.assertEqual(jid, "7548417")
        self.assertEqual(sleep.calls, [])
        self.assertEqual(run.calls, [["sacct", "-j", "7548417", "-n", "-X",
                                      "-P", "-o", "State,ExitCode"]])

    def test_real_failure_after_lag_names_failed(self):
        _, run, _, err = verdict(["RUNNING|0:0", "FAILED|1:0"])
        self.assertIsInstance(err, RuntimeError)
        self.assertNotIsInstance(err, C.SlurmStateUnsettled)
        self.assertIn("FAILED", str(err))
        self.assertIn("1:0", str(err))
        self.assertNotIn("RUNNING", str(err))
        self.assertEqual(len(run.calls), 2)

    def test_never_settles_raises_unsettled(self):
        _, run, sleep, err = verdict(["RUNNING|0:0"], timeout_s=600,
                                     poll_s=30)
        self.assertIsInstance(err, C.SlurmStateUnsettled)
        self.assertIn("unsettled", str(err))
        self.assertIn("UNKNOWN", str(err))
        self.assertNotIn("ended in state RUNNING", str(err))
        # bounded: total sleep within the timeout, 1 + 600/30 queries
        self.assertLessEqual(sum(sleep.calls), 600)
        self.assertEqual(len(run.calls), 21)

    def test_empty_output_then_completed(self):
        jid, run, sleep, err = verdict(["", "\n", "COMPLETED|0:0"])
        self.assertIsNone(err)
        self.assertEqual(jid, "7548417")
        self.assertEqual(len(run.calls), 3)

    def test_failed_sacct_call_is_not_settled(self):
        jid, run, _, err = verdict(["", "COMPLETED|0:0"],
                                   returncodes=[1, 0])
        self.assertIsNone(err)
        self.assertEqual(len(run.calls), 2)

    def test_cancelled_by_user(self):
        for out in ("CANCELLED by 1234|0:15", "CANCELLED by 1234   0:15"):
            with self.subTest(out=out):
                _, run, sleep, err = verdict(["RUNNING|0:0", out])
                self.assertIsInstance(err, RuntimeError)
                self.assertNotIsInstance(err, C.SlurmStateUnsettled)
                self.assertIn("ended in state CANCELLED", str(err))

    def test_other_not_settled_states_keep_waiting(self):
        for st in ("PENDING", "COMPLETING", "REQUEUED", "RESIZING",
                   "SUSPENDED"):
            with self.subTest(state=st):
                jid, run, _, err = verdict([f"{st}|0:0", "COMPLETED|0:0"])
                self.assertIsNone(err)
                self.assertEqual(len(run.calls), 2)

    def test_terminal_non_completed_states_fail(self):
        for st in ("TIMEOUT", "OUT_OF_MEMORY", "NODE_FAIL", "PREEMPTED",
                   "BOOT_FAIL", "DEADLINE"):
            with self.subTest(state=st):
                _, run, sleep, err = verdict([f"{st}|0:1"])
                self.assertIn(f"ended in state {st}", str(err))
                self.assertEqual(sleep.calls, [])

    def test_completed_requires_exit_0_0(self):
        _, _, _, err = verdict(["COMPLETED|1:0"])
        self.assertIsInstance(err, RuntimeError)
        self.assertIn("COMPLETED", str(err))

    def test_sleeps_between_sacct_calls_at_least_30s(self):
        # even when a caller asks for a faster poll, the floor holds
        for poll in (None, 0, 5, 29.9, 30, 45):
            with self.subTest(poll=poll):
                _, _, sleep, err = verdict(["RUNNING|0:0"], timeout_s=300,
                                           poll_s=poll)
                self.assertIsInstance(err, C.SlurmStateUnsettled)
                self.assertTrue(sleep.calls)
                self.assertTrue(all(s >= 30 for s in sleep.calls),
                                sleep.calls)
        self.assertGreaterEqual(C.SACCT_SETTLE_POLL_S, 30)
        self.assertEqual(C.SACCT_MIN_POLL_S, 30)

    def test_default_timeout_is_about_ten_minutes(self):
        self.assertEqual(C.SACCT_SETTLE_TIMEOUT_S, 600)


class DagUsesSharedHelper(unittest.TestCase):
    """No single immediate sacct call may creep back into the DAGs."""

    def test_no_inline_sacct_in_dags(self):
        dags = os.path.join(HERE, "..", "dags")
        for name in sorted(os.listdir(dags)):
            if not name.endswith(".py") or name == "ar_common.py":
                continue
            with open(os.path.join(dags, name)) as f:
                src = f.read()
            with self.subTest(dag=name):
                self.assertNotIn('"sacct"', src)
        with open(os.path.join(dags, "apogee_daily.py")) as f:
            self.assertIn("C.slurm_job_verdict(job_id)", f.read())


if __name__ == "__main__":
    unittest.main()
