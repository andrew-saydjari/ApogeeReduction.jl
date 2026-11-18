"""DAG-level tests for apogee_daily (needs Airflow; skipped without it).

Run with the deployed Airflow env against a throw-away AIRFLOW_HOME (nothing
touches the production metadata DB):

    AIRFLOW_HOME=$(mktemp -d) AIRFLOW__CORE__LOAD_EXAMPLES=False \
    PYTHONDONTWRITEBYTECODE=1 \
    /mnt/home/sdssv/uv_env/airflow_env/bin/python -m unittest \
        discover -s airflow/tests -v

The almanac end-to-end test also needs a python with h5py for the stand-in
almanac (AR_TEST_H5_PYTHON, default: the production almanac env's python);
it is skipped when none is found.
"""

import os
import shutil
import subprocess
import sys
import tempfile
import textwrap
import unittest
import warnings
from datetime import datetime, timezone

HERE = os.path.dirname(os.path.abspath(__file__))
DAGS = os.path.abspath(os.path.join(HERE, "..", "dags"))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
RENDER = os.path.join(HERE, "_render_task.py")
sys.path.insert(0, DAGS)
warnings.simplefilter("ignore")

try:
    import airflow  # noqa: F401
    try:
        from airflow.sdk.exceptions import AirflowFailException
    except ImportError:  # pragma: no cover
        from airflow.exceptions import AirflowFailException
    HAVE_AIRFLOW = True
except ImportError:  # pragma: no cover
    HAVE_AIRFLOW = False

H5_PYTHON = os.environ.get(
    "AR_TEST_H5_PYTHON", "/mnt/home/sdssv/uv_env/almanac_env/bin/python")


def have_h5_python():
    if not os.path.exists(H5_PYTHON):
        return False
    r = subprocess.run([H5_PYTHON, "-c", "import h5py"],
                       capture_output=True)
    return r.returncode == 0


def render_task(task_id, outdir, mjd=61360, env=None, **params):
    args = [sys.executable, RENDER, task_id, outdir, str(mjd)]
    args += [f"{k}={v}" for k, v in params.items()]
    r = subprocess.run(args, capture_output=True, text=True,
                       env={**os.environ, **(env or {})})
    if r.returncode != 0:
        raise AssertionError(f"render {task_id} failed:\n{r.stderr[-3000:]}")
    return r.stdout


@unittest.skipUnless(HAVE_AIRFLOW, "airflow not importable")
class DagImport(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        try:
            from airflow.dag_processing.dagbag import DagBag
        except ImportError:  # older layout
            from airflow.models.dagbag import DagBag
        cls.bag = DagBag(dag_folder=DAGS)
        cls.dag = cls.bag.dags.get("apogee_daily")
        import ar_common
        cls.C = ar_common

    def test_no_import_errors(self):
        self.assertEqual(dict(self.bag.import_errors), {})
        self.assertIn("apogee_daily", self.bag.dags)
        self.assertIn("apogee_heartbeat", self.bag.dags)

    def test_stale_guard_is_the_only_root(self):
        self.assertEqual([t.task_id for t in self.dag.roots],
                         ["stale_run_guard"])
        self.assertIn("update.repo",
                      self.dag.get_task("stale_run_guard").downstream_task_ids)

    def test_max_active_runs(self):
        self.assertEqual(self.dag.max_active_runs, 1)

    def test_startup_tasks_retry(self):
        for tid in ("stale_run_guard", "update.repo", "update.sdsscore",
                    "update.sync_logs", "setup.mjd", "setup.date_mjd"):
            with self.subTest(tid=tid):
                t = self.dag.get_task(tid)
                self.assertEqual(t.retries, self.C.STARTUP_RETRIES)
                self.assertEqual(t.retry_delay.total_seconds(),
                                 self.C.STARTUP_RETRY_DELAY_S)
                self.assertGreaterEqual(t.retries, 1)

    def test_almanac_limits(self):
        for tele in ("lco", "apo"):
            t = self.dag.get_task(f"{tele}.almanac")
            self.assertEqual(t.execution_timeout.total_seconds(),
                             self.C.ALMANAC_EXECUTION_TIMEOUT_S)
            self.assertEqual(t.retries, self.C.ALMANAC_RETRIES)
            self.assertEqual(t.retry_delay.total_seconds(),
                             self.C.ALMANAC_RETRY_DELAY_S)
            self.assertIn(self.C.ALMANAC_PGOPTIONS, t.bash_command)
            self.assertIn('--output "$TMP"', t.bash_command)
            self.assertIn("check_almanac_output.py", t.bash_command)

    def test_tunnel_check_runs_a_query(self):
        cmd = self.dag.get_task("lco.tunnel_check").bash_command
        self.assertIn("select 1", cmd)
        self.assertIn(self.C.TUNNEL_UNIT, cmd)
        self.assertNotIn("sshpass", cmd)
        self.assertNotIn("pass_file", cmd)

    def test_stale_run_guard_callable(self):
        import apogee_daily as D

        class Run:
            def __init__(self, run_type):
                self.run_type = run_type

        old = datetime(2026, 9, 27, 11, tzinfo=timezone.utc)
        ctx = {"params": {"stale_run_max_age_h": 36},
               "dag_run": Run("scheduled"), "logical_date": old}
        with self.assertRaises(AirflowFailException):
            D.stale_run_guard(**ctx)
        ctx["dag_run"] = Run("manual")
        D.stale_run_guard(**ctx)  # manual re-runs of old nights pass
        ctx["dag_run"] = Run("scheduled")
        ctx["logical_date"] = datetime.now(timezone.utc)
        D.stale_run_guard(**ctx)

    def test_every_bash_command_renders_and_parses(self):
        """Render each BashOperator template with stand-in context and run
        `bash -n` on it (catches quoting mistakes in the nested bash -c)."""
        from airflow.providers.standard.operators.bash import BashOperator
        sys.path.insert(0, HERE)
        import _render_task as R
        bash_tasks = [t for t in self.dag.tasks if isinstance(t, BashOperator)]
        self.assertGreater(len(bash_tasks), 10)
        with tempfile.TemporaryDirectory() as d:
            for t in bash_tasks:
                with self.subTest(task=t.task_id):
                    text = R.render(t, d + "/", 61360)
                    r = subprocess.run(["bash", "-n"], input=text,
                                       capture_output=True, text=True)
                    self.assertEqual(r.returncode, 0, r.stderr)


FAKE_ALMANAC = textwrap.dedent("""\
    #!/bin/bash
    # stand-in for almanac: --output PATH, behaviour from FAKE_MODE
    out=""; tele=""; mjd=""
    while [ $# -gt 0 ]; do
      case "$1" in
        --output) out="$2"; shift;;
        --lco|--apo) tele="${1#--}";;
        --mjd-start) mjd="$2"; shift;;
      esac
      shift
    done
    echo "fake almanac PGOPTIONS=$PGOPTIONS" >> "$FAKE_LOG"
    echo "$out" >> "$FAKE_LOG"
    case "$FAKE_MODE" in
      ok) exec "$H5_PYTHON" -c "import h5py,sys; f=h5py.File(sys.argv[1],'a'); \
    f.require_group('raw/%s/%s/exposures'%(sys.argv[2],sys.argv[3])); \
    f.require_group('raw/%s/%s/fibers'%(sys.argv[2],sys.argv[3])); f.close()" \
            "$out" "$tele" "$mjd";;
      fail) printf 'partial' > "$out"; exit 3;;
    esac
    """)


@unittest.skipUnless(HAVE_AIRFLOW, "airflow not importable")
@unittest.skipUnless(have_h5_python(), "no python with h5py for the fake")
class AlmanacStepEndToEnd(unittest.TestCase):
    """Run the rendered almanac step with a stand-in almanac binary."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.out = self.tmp + "/out/"
        os.makedirs(self.out + "almanac")
        self.fake = os.path.join(self.tmp, "almanac")
        with open(self.fake, "w") as f:
            f.write(FAKE_ALMANAC)
        os.chmod(self.fake, 0o755)
        self.log = os.path.join(self.tmp, "fake.log")
        self.af = self.out + "almanac/allobs_lco_61360.h5"
        self.env = {"ALMANAC_BIN": self.fake, "ALMANAC_PYTHON": H5_PYTHON,
                    "ALMANAC_DIR": self.tmp, "AR_REPO": REPO}

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def run_step(self, mode, **params):
        cmd = render_task("lco.almanac", self.out, env=self.env,
                          slack_mode="off", **params)
        env = {**os.environ, "FAKE_MODE": mode, "FAKE_LOG": self.log,
               "H5_PYTHON": H5_PYTHON}
        return subprocess.run(["bash", "-c", cmd], capture_output=True,
                              text=True, env=env)

    def calls(self):
        if not os.path.exists(self.log):
            return 0
        with open(self.log) as f:
            return sum(1 for line in f if line.startswith("fake almanac"))

    def test_stub_is_rebuilt_then_reused(self):
        with open(self.af, "wb") as f:
            f.write(b"\x89HDF\r\n\x1a\n" + b"\0" * 88)  # 96-byte stub
        r = self.run_step("ok")
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertIn("does not hold this night", r.stdout)
        self.assertGreater(os.path.getsize(self.af), 96)
        self.assertFalse(os.path.exists(self.af + ".partial"))
        self.assertEqual(self.calls(), 1)
        with open(self.log) as f:
            log = f.read()
        self.assertIn("statement_timeout=", log)  # PGOPTIONS reached almanac
        self.assertIn(self.af + ".partial", log)  # wrote to the temp path
        # second run: valid output -> skip without calling almanac
        r = self.run_step("ok")
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertIn("already holds this night", r.stdout)
        self.assertEqual(self.calls(), 1)
        # clobber forces a rebuild even of a valid file
        r = self.run_step("ok", almanac_clobber="true")
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertEqual(self.calls(), 2)

    def test_failure_leaves_no_stub(self):
        r = self.run_step("fail")
        self.assertEqual(r.returncode, 3, r.stdout + r.stderr)
        self.assertFalse(os.path.exists(self.af))
        self.assertFalse(os.path.exists(self.af + ".partial"))

    def test_failure_keeps_previous_good_output(self):
        self.assertEqual(self.run_step("ok").returncode, 0)
        before = os.path.getsize(self.af)
        r = self.run_step("fail", almanac_clobber="true")
        self.assertEqual(r.returncode, 3)
        self.assertEqual(os.path.getsize(self.af), before)
        self.assertFalse(os.path.exists(self.af + ".partial"))

    def test_announce_reaches_the_task_log(self):
        r = self.run_step("ok")
        self.assertIn("NOT a hang", r.stdout)
        self.assertIn("full log", r.stdout)


if __name__ == "__main__":
    unittest.main()
