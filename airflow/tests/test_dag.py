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

    def test_setup_tasks_without_data_interval(self):
        """A manual trigger without --logical-date has no data interval
        (2026-11-17: KeyError 'data_interval_start' in setup.date_mjd)."""
        import apogee_daily as D

        class Run:
            run_type = "manual"
            run_after = datetime(2026, 11, 18, 2, tzinfo=timezone.utc)
            logical_date = None

        ctx = {"params": {"mjd": 61360, "almanac_clobber": True},
               "dag_run": Run(), "logical_date": None}
        self.assertEqual(D.compute_mjd(**ctx), 61360)
        self.assertEqual(D.compute_date_mjd(**ctx), "2026-11-15")
        ctx["params"] = {"mjd": -1}
        self.assertEqual(D.compute_mjd(**ctx), 61361)  # from run_after
        self.assertEqual(D.compute_date_mjd(**ctx), "2026-11-15")
        D.stale_run_guard(**{**ctx, "params": {"stale_run_max_age_h": 36}})

    def test_p2d1d_passes_wavecal_order_and_history_roots(self):
        """The inline p2d1d_full step passes --wavecal_order and
        --arcfirst_history_roots (default "": no history, the arc_first
        cross-night step reports "unchecked"), like run_all.sh's env route."""
        sys.path.insert(0, HERE)
        import _render_task as R
        self.assertEqual(self.C.ARCFIRST_HISTORY_ROOTS, "")
        tasks = [t for t in self.dag.tasks if t.task_id.endswith("p2d1d_full")]
        self.assertEqual(len(tasks), 2)
        with tempfile.TemporaryDirectory() as d:
            for t in tasks:
                with self.subTest(task=t.task_id):
                    text = R.render(t, d + "/", 61360)
                    self.assertIn(f"--wavecal_order {self.C.WAVECAL_ORDER} ", text)
                    self.assertIn("--arcfirst_history_roots '' ", text)

    def test_slurm_submit_exports_wavecal_env(self):
        """slurm_submit_and_wait exports AR_WAVECAL_ORDER and
        AR_ARCFIRST_HISTORY_ROOTS into the sbatch'd run_all.sh (subprocess
        is stubbed: no Slurm command runs)."""
        from unittest import mock
        import apogee_daily as D
        seen = {}

        class Res:
            def __init__(self, out):
                self.returncode, self.stdout, self.stderr = 0, out, ""

        def fake_run(cmd, **kw):
            if cmd[0] == "sbatch":
                seen["env"] = kw.get("env") or {}
                return Res("4242\n")
            if cmd[0] == "squeue":
                return Res("")
            return Res("COMPLETED 0:0\n")

        class TI:
            def __init__(self, outdir):
                self.outdir = outdir

            def xcom_pull(self, task_ids=None, **_):
                return {"outdir": self.outdir, "mjd": 61360}

        with tempfile.TemporaryDirectory() as d:
            for roots in ("", "/a/wavecal_arcfirst,/b/wavecal_arcfirst"):
                with self.subTest(roots=roots), \
                        mock.patch.object(D.C, "ARCFIRST_HISTORY_ROOTS", roots), \
                        mock.patch.object(D.subprocess, "run", fake_run):
                    seen.clear()
                    D.slurm_submit_and_wait("apo", params={}, ti=TI(d))
                    self.assertEqual(seen["env"]["AR_WAVECAL_ORDER"], self.C.WAVECAL_ORDER)
                    self.assertEqual(seen["env"]["AR_ARCFIRST_HISTORY_ROOTS"], roots)

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


class RunAllWavecalOptions(unittest.TestCase):
    """scripts/daily/run_all.sh turns AR_WAVECAL_ORDER / AR_ARCFIRST_HISTORY_ROOTS
    into pipeline_2d_1d.jl flags (unset or "" -> no flag). Runs only the
    option blocks of the script, extracted from it (no Airflow needed)."""
    SCRIPT = os.path.join(REPO, "scripts", "daily", "run_all.sh")

    def _block(self, var):
        lines = open(self.SCRIPT).read().splitlines()
        i = next(k for k, l in enumerate(lines)
                 if l.startswith("if ") and "${" + var + "+set}" in l)
        j = next(k for k in range(i, len(lines)) if lines[k] == "fi")
        return "\n".join(lines[i:j + 1])

    def _opts(self, var, arr, value):
        env = {k: v for k, v in os.environ.items() if k != var}
        if value is not None:
            env[var] = value
        script = self._block(var) + '\nprintf "%s\\n" "${' + arr + '[@]}"\n'
        r = subprocess.run(["bash", "-c", script], capture_output=True,
                           text=True, env=env)
        self.assertEqual(r.returncode, 0, r.stderr)
        return [l for l in r.stdout.splitlines() if l]

    def test_history_roots(self):
        v, a = "AR_ARCFIRST_HISTORY_ROOTS", "arcfirst_hist_opt"
        self.assertEqual(self._opts(v, a, None), [])
        self.assertEqual(self._opts(v, a, ""), [])
        self.assertEqual(self._opts(v, a, "/x/wavecal_arcfirst,/y/wavecal_arcfirst"),
                         ["--arcfirst_history_roots",
                          "/x/wavecal_arcfirst,/y/wavecal_arcfirst"])

    def test_wavecal_order(self):
        v, a = "AR_WAVECAL_ORDER", "wavecal_order_opt"
        self.assertEqual(self._opts(v, a, None), [])
        self.assertEqual(self._opts(v, a, "arc_first"),
                         ["--wavecal_order", "arc_first"])

    def test_flags_reach_pipeline_2d_1d(self):
        # the science 2D->1D call (the flats' relFlux pass runs --waveSoln false)
        call = [l for l in open(self.SCRIPT).read().splitlines()
                if "pipeline_2d_1d.jl" in l and l.lstrip().startswith("julia")
                and "--waveSoln false" not in l]
        self.assertEqual(len(call), 1)
        self.assertIn('"${wavecal_order_opt[@]}"', call[0])
        self.assertIn('"${arcfirst_hist_opt[@]}"', call[0])


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

    def test_truncated_file_with_clobber_is_replaced(self):
        """2026-11-17: with almanac_clobber, almanac itself failed opening
        the killed run's truncated 96-byte file to append. almanac must
        never be pointed at an existing output: it writes a fresh .partial
        that replaces the file on success."""
        with open(self.af, "wb") as f:
            f.write(b"\x89HDF\r\n\x1a\n" + b"\0" * 88)  # 96-byte stub
        r = self.run_step("ok", almanac_clobber="true")
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        with open(self.log) as f:
            outputs = [ln.strip() for ln in f if not ln.startswith("fake")]
        self.assertEqual(outputs, [self.af + ".partial"])
        chk = subprocess.run(
            [H5_PYTHON, os.path.join(REPO, "airflow", "scripts",
                                     "check_almanac_output.py"),
             self.af, "lco", "61360", "--fibers"],
            capture_output=True, text=True)
        self.assertEqual(chk.returncode, 0, chk.stdout)
        self.assertFalse(os.path.exists(self.af + ".partial"))

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
