"""Tests for airflow/scripts/render_systemd_units.sh and the unit templates.
Stdlib only (needs bash + sed). Optional: `systemd-analyze --user verify`
on the rendered units when systemd-analyze is available."""

import os
import re
import shutil
import subprocess
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(HERE, "..", "scripts", "render_systemd_units.sh")
UNITS = ("airflow-standalone.service", "airflow-failure-notify.service",
         "mwm-tunnel.service")
REQ = ["--sshpass", "/opt/bin/sshpass", "--pass-file", "/home/u/pass_file"]


def render(outdir, *args):
    return subprocess.run(["bash", SCRIPT, *args, outdir],
                          capture_output=True, text=True)


def sections(path):
    """Unit file -> {section: {key: [values]}} (keys may repeat)."""
    out, cur = {}, None
    with open(path) as f:
        text = f.read().replace("\\\n", " ")
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        m = re.fullmatch(r"\[(\w+)\]", line)
        if m:
            cur = out.setdefault(m.group(1), {})
            continue
        k, _, v = line.partition("=")
        cur.setdefault(k, []).append(v)
    return out


class RenderUnits(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.d, ignore_errors=True)

    def test_render_defaults(self):
        r = render(self.d, *REQ)
        self.assertEqual(r.returncode, 0, r.stderr)
        for u in UNITS:
            with open(os.path.join(self.d, u)) as f:
                text = f.read()
            self.assertNotRegex(text, r"@[A-Z_]+@", u)
            # no %h in any directive (comments may mention it)
            directives = [ln for ln in text.splitlines()
                          if ln.strip() and not ln.lstrip().startswith("#")]
            self.assertFalse([ln for ln in directives if "%h" in ln], u)

        a = sections(os.path.join(self.d, "airflow-standalone.service"))
        self.assertEqual(a["Service"]["EnvironmentFile"],
                         ["/mnt/home/sdssv/airflow/airflow_env.sh"])
        self.assertEqual(a["Service"]["WorkingDirectory"],
                         ["/mnt/home/sdssv/airflow"])
        self.assertIn("mwm-tunnel.service", a["Unit"]["After"][0])
        self.assertEqual(a["Unit"]["Wants"], ["mwm-tunnel.service"])
        self.assertIn("StartLimitBurst", a["Unit"])
        self.assertNotIn("StartLimitBurst", a["Service"])
        self.assertNotIn("StartLimitIntervalSec", a["Service"])

        n = sections(os.path.join(self.d, "airflow-failure-notify.service"))
        self.assertEqual(n["Service"]["EnvironmentFile"],
                         ["/mnt/home/sdssv/airflow/airflow_env.sh"])

        t = sections(os.path.join(self.d, "mwm-tunnel.service"))
        self.assertEqual(t["Unit"]["StartLimitBurst"], ["3"])
        self.assertEqual(t["Unit"]["StartLimitIntervalSec"], ["3600"])
        self.assertNotIn("StartLimitBurst", t["Service"])
        start = t["Service"]["ExecStart"][0]
        for frag in ("/opt/bin/sshpass -f /home/u/pass_file",
                     "ControlMaster=yes", "ControlPersist=no",
                     "ServerAliveInterval=30", "ServerAliveCountMax=3",
                     "ExitOnForwardFailure=yes",
                     "-L 63333:operations.sdss.org:5432", " -N "):
            self.assertIn(frag, start)
        self.assertTrue(start.endswith(" mwm"))
        self.assertIn("rm -f", t["Service"]["ExecStartPre"][0])
        self.assertEqual(t["Service"]["Restart"], ["always"])

    def test_overrides(self):
        r = render(self.d, *REQ, "--airflow-home", "/srv/af",
                   "--tunnel-host", "utah", "--db-forward",
                   "5433:db.example.org:5432")
        self.assertEqual(r.returncode, 0, r.stderr)
        a = sections(os.path.join(self.d, "airflow-standalone.service"))
        self.assertEqual(a["Service"]["EnvironmentFile"],
                         ["/srv/af/airflow_env.sh"])
        t = sections(os.path.join(self.d, "mwm-tunnel.service"))
        self.assertIn("-L 5433:db.example.org:5432 utah",
                      t["Service"]["ExecStart"][0])
        self.assertIn("ssh -G utah", t["Service"]["ExecStartPre"][0])

    def test_refusals(self):
        cases = [
            [],                                        # no --sshpass etc.
            ["--sshpass", "/opt/sshpass"],             # no --pass-file
            REQ + ["--airflow-home", "relative/dir"],  # not absolute
            REQ + ["--env-file", "/a b/env"],          # whitespace
            REQ + ["--db-forward", "63333"],           # malformed -L
            REQ + ["--bogus", "x"],
        ]
        for args in cases:
            with self.subTest(args=args):
                self.assertNotEqual(render(self.d, *args).returncode, 0)

    @unittest.skipUnless(shutil.which("systemd-analyze"),
                         "systemd-analyze not available")
    def test_systemd_analyze_verify(self):
        self.assertEqual(render(self.d, *REQ).returncode, 0)
        r = subprocess.run(
            ["systemd-analyze", "--user", "verify",
             *[os.path.join(self.d, u) for u in UNITS]],
            capture_output=True, text=True, timeout=60)
        # executables named by the units need not exist on a CI runner;
        # only syntax/semantic complaints about the files themselves count
        problems = [ln for ln in r.stderr.splitlines()
                    if "is not executable" not in ln
                    and "No such file" not in ln
                    and "Failed to prepare filename" not in ln
                    and "Command " not in ln and ln.strip()]
        self.assertEqual(problems, [], r.stderr)


if __name__ == "__main__":
    unittest.main()
