"""airflow/scripts/sync_almanac_env.sh (update.repo's almanac env re-sync).

Stdlib only, with stand-in `uv` and python executables, so this runs under any
python3:  python3 -m unittest discover -s airflow/tests -v
Nothing touches a real environment.
"""

import os
import stat
import subprocess
import tempfile
import textwrap
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.abspath(os.path.join(HERE, "..", "scripts", "sync_almanac_env.sh"))
DAG = os.path.abspath(os.path.join(HERE, "..", "dags", "apogee_daily.py"))


def _exe(path, body):
    with open(path, "w") as fp:
        fp.write("#!/bin/bash\n" + textwrap.dedent(body))
    os.chmod(path, os.stat(path).st_mode | stat.S_IXUSR)
    return path


class SyncAlmanacEnv(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.d = self.tmp.name
        self.log = os.path.join(self.d, "uv_calls.txt")

    def tearDown(self):
        self.tmp.cleanup()

    def run_script(self, uv_body, py_body, extra_env=None):
        uv = _exe(os.path.join(self.d, "uv"), f'echo "$@" >> {self.log}\n' + uv_body)
        py = _exe(os.path.join(self.d, "python"), py_body)
        env = {**os.environ, "UV": uv, **(extra_env or {})}
        return subprocess.run(["bash", SCRIPT, py, "/clone/almanac"], capture_output=True,
                              text=True, env=env)

    def test_success(self):
        r = self.run_script("echo ' + healpy==1.20.1' >&2; exit 0\n",
                            "echo '0.4.3 /clone/almanac/src/almanac/__init__.py'\n")
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertIn("+ healpy==1.20.1", r.stdout)  # what changed is logged
        self.assertIn("almanac env OK", r.stdout)
        with open(self.log) as fp:
            call = fp.read().split()
        self.assertEqual(call[:2], ["pip", "install"])
        self.assertIn("-e", call)
        self.assertEqual(call[call.index("-e") + 1], "/clone/almanac")
        self.assertEqual(call[call.index("-p") + 1], os.path.join(self.d, "python"))

    def test_failed_install_fails(self):
        # The old inline command (`... | tail -1 || echo WARNING`) exited 0 here.
        r = self.run_script("echo 'error: No solution found' >&2; exit 2\n", "exit 0\n")
        self.assertEqual(r.returncode, 1)
        self.assertIn("No solution found", r.stdout)
        self.assertIn("FATAL: almanac env re-sync failed", r.stdout)

    def test_failed_import_fails(self):
        r = self.run_script("exit 0\n",
                            "echo \"ModuleNotFoundError: No module named 'healpy'\"; exit 1\n")
        self.assertEqual(r.returncode, 1)
        self.assertIn("healpy", r.stdout)
        self.assertIn("cannot import almanac.cli", r.stdout)

    def test_missing_uv_fails(self):
        py = _exe(os.path.join(self.d, "python"), "exit 0\n")
        r = subprocess.run(["bash", SCRIPT, py, "/clone/almanac"], capture_output=True, text=True,
                           env={**os.environ, "UV": os.path.join(self.d, "no-such-uv")})
        self.assertEqual(r.returncode, 1)
        self.assertIn("not found", r.stdout)

    def test_old_inline_command_masked_failures(self):
        """Documents the bug being fixed: a pipeline takes tail's exit status."""
        r = subprocess.run(["bash", "-c", "false 2>&1 | tail -1 || echo WARNING"],
                           capture_output=True, text=True)
        self.assertEqual(r.returncode, 0)
        self.assertNotIn("WARNING", r.stdout)

    def test_dag_update_repo_uses_the_script_and_fails_on_it(self):
        with open(DAG) as fp:
            src = fp.read()
        self.assertIn("/airflow/scripts/sync_almanac_env.sh ", src)
        self.assertIn('{C.ALMANAC_PYTHON} {C.ALMANAC_DIR} || exit 1', src)
        self.assertNotIn("--quiet 2>&1 | tail -1 ||", src)

    def test_script_is_executable(self):
        self.assertTrue(os.access(SCRIPT, os.X_OK))


if __name__ == "__main__":
    unittest.main()
