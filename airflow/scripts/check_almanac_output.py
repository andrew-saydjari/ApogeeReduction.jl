#!/usr/bin/env python3
"""Does an almanac HDF5 file actually hold a given observatory/MJD night?

Used by the apogee_daily `almanac` task to decide whether an existing output
may be reused (skip) or must be rebuilt. Before 2026-11-17 the task skipped
whenever the file merely EXISTED, so a killed run's 96-byte stub (not even a
readable HDF5 file) made every retry skip almanac and hand runlist an empty
file.

A night is present when the file opens and contains the group
``raw/<observatory>/<mjd>/exposures`` (written by almanac.io.update), plus
``raw/<observatory>/<mjd>/fibers`` when --fibers is given (the daily runs
almanac with --fibers, and almanac creates that group for every night it
writes in fiber mode).

Run with the almanac env's python (it has h5py):

    check_almanac_output.py FILE OBSERVATORY MJD [--fibers]

Exit status: 0 = night present (reuse), 1 = absent or unreadable (rebuild),
2 = usage error. One line on stdout says why.
"""

from __future__ import annotations

import argparse
import os
import sys


def night_group(observatory: str, mjd: int) -> str:
    return f"raw/{observatory}/{int(mjd)}"


def almanac_night_status(path: str, observatory: str, mjd: int,
                         fibers: bool = False) -> tuple[bool, str]:
    """(present, reason) for one night in an almanac file. Never raises for
    a missing or corrupt file: those are simply 'not present'."""
    if not os.path.exists(path):
        return False, f"{path} does not exist"
    size = os.path.getsize(path)
    try:
        import h5py
    except ImportError as e:  # wrong interpreter: refuse to guess "present"
        return False, f"cannot check {path}: h5py not importable ({e})"
    group = night_group(observatory, mjd)
    try:
        with h5py.File(path, "r") as fp:
            needed = [f"{group}/exposures"]
            if fibers:
                needed.append(f"{group}/fibers")
            for name in needed:
                obj = fp.get(name)
                if obj is None or not isinstance(obj, h5py.Group):
                    return False, (f"{path} ({size} B) has no group "
                                   f"{name}")
    except OSError as e:
        return False, f"{path} ({size} B) is not a readable HDF5 file: {e}"
    return True, f"{path} ({size} B) contains {group}"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("path")
    ap.add_argument("observatory", choices=["apo", "lco"])
    ap.add_argument("mjd", type=int)
    ap.add_argument("--fibers", action="store_true",
                    help="also require the night's fibers group")
    args = ap.parse_args(argv)
    ok, reason = almanac_night_status(args.path, args.observatory, args.mjd,
                                      fibers=args.fibers)
    print(("present: " if ok else "absent: ") + reason)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
