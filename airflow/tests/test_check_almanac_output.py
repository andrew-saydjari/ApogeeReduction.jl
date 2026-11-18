"""Tests for airflow/scripts/check_almanac_output.py (needs h5py; skipped
without it). In production h5py lives in the almanac env:

    /mnt/home/sdssv/uv_env/almanac_env/bin/python -m unittest \
        discover -s airflow/tests -p 'test_check_almanac_output.py' -v
"""

import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(HERE, "..", "scripts", "check_almanac_output.py")
sys.path.insert(0, os.path.join(HERE, "..", "scripts"))

import check_almanac_output as K  # noqa: E402

try:
    import h5py
except ImportError:  # pragma: no cover
    h5py = None


def make_almanac(path, groups):
    with h5py.File(path, "w") as fp:
        for g in groups:
            fp.require_group(g)
        fp.create_dataset("missing_exposures/mjd", data=[0])


@unittest.skipIf(h5py is None, "h5py not importable in this interpreter")
class AlmanacNightStatus(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.d = self.tmp.name

    def tearDown(self):
        self.tmp.cleanup()

    def p(self, name):
        return os.path.join(self.d, name)

    def test_complete_night_is_present(self):
        f = self.p("allobs_lco_61360.h5")
        make_almanac(f, ["raw/lco/61360/exposures", "raw/lco/61360/fibers",
                         "raw/lco/61360/sequences"])
        ok, why = K.almanac_night_status(f, "lco", 61360, fibers=True)
        self.assertTrue(ok, why)

    def test_missing_file(self):
        ok, why = K.almanac_night_status(self.p("nope.h5"), "lco", 61360)
        self.assertFalse(ok)
        self.assertIn("does not exist", why)

    def test_96_byte_stub_from_a_killed_run(self):
        # what the killed 2026-11-17 run left: an HDF5 superblock only
        f = self.p("stub.h5")
        make_almanac(f, ["raw/lco/61360/exposures"])
        with open(f, "r+b") as fh:
            fh.truncate(96)
        ok, why = K.almanac_night_status(f, "lco", 61360, fibers=True)
        self.assertFalse(ok)
        self.assertIn("not a readable HDF5", why)

    def test_empty_and_garbage_files(self):
        for name, data in (("empty.h5", b""), ("junk.h5", b"not hdf5\n")):
            with self.subTest(name=name):
                f = self.p(name)
                with open(f, "wb") as fh:
                    fh.write(data)
                self.assertFalse(K.almanac_night_status(f, "lco", 61360)[0])

    def test_other_night_or_observatory_is_absent(self):
        f = self.p("other.h5")
        make_almanac(f, ["raw/apo/61360/exposures", "raw/apo/61360/fibers",
                         "raw/lco/61359/exposures", "raw/lco/61359/fibers"])
        self.assertFalse(K.almanac_night_status(f, "lco", 61360, True)[0])
        self.assertTrue(K.almanac_night_status(f, "apo", 61360, True)[0])

    def test_fibers_required_only_with_flag(self):
        f = self.p("nofibers.h5")
        make_almanac(f, ["raw/lco/61360/exposures"])
        self.assertTrue(K.almanac_night_status(f, "lco", 61360, False)[0])
        ok, why = K.almanac_night_status(f, "lco", 61360, True)
        self.assertFalse(ok)
        self.assertIn("fibers", why)

    def test_night_group_without_exposures(self):
        f = self.p("noexp.h5")
        make_almanac(f, ["raw/lco/61360/fibers"])
        self.assertFalse(K.almanac_night_status(f, "lco", 61360, True)[0])

    def test_cli_exit_codes(self):
        good = self.p("good.h5")
        make_almanac(good, ["raw/lco/61360/exposures", "raw/lco/61360/fibers"])
        run = lambda *a: subprocess.run(  # noqa: E731
            [sys.executable, SCRIPT, *a], capture_output=True, text=True)
        r = run(good, "lco", "61360", "--fibers")
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertTrue(r.stdout.startswith("present:"))
        r = run(good, "apo", "61360", "--fibers")
        self.assertEqual(r.returncode, 1)
        self.assertTrue(r.stdout.startswith("absent:"))
        r = run(good, "xyz", "61360")
        self.assertEqual(r.returncode, 2)  # usage error -> caller rebuilds


if __name__ == "__main__":
    unittest.main()
