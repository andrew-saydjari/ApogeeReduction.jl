# Submitting the golden-baseline generation job

One multinode Slurm job generates all T1 golden baselines (6 runs: apo+lco
60000, apo 58588, 59337, 58011, 59429), mirroring the run_bulk.sh
SlurmClusterManager pattern. Preparation, validation, and almanac staging are
already done — **only the submission itself is left, and it is manual by
design.**

## Submit

```bash
cd /mnt/home/asaydjari/gitcode/worktrees/AR-T1     # the T1 branch checkout
mkdir -p slurm_logs                                # sbatch will not create it
sbatchAKS test/regression/submit_goldens.sh
```

Then **paste the job ID back into the Claude session** so progress can be
followed from the filesystem (no Slurm polling needed).

Allocation requested by the script: `-p cca --qos=cca -N 2
-C "[genoa|icelake|rome]" --mem=900G -t 1-12:00 -J ar_goldens`.
Expected wall time: ~1–2 h per full day (90–101 exposures) + ~20 min for the
small days — **roughly 6–10 h total**; the 36 h limit is generous margin.

Before submitting, optional sanity checks (no Slurm needed):

```bash
bash -n test/regression/submit_goldens.sh
test/regression/submit_goldens.sh --dry-run   # prints the exact 6-run sequence
```

The job refuses to start from a dirty git tree (goldens must be attributable
to a pushed sha) and verifies all staged almanac inputs exist before running
anything.

## Inputs (pre-staged, zero DB/tunnel dependency)

- Per-day almanacs: `/mnt/ceph/users/sdssv/work/asaydjari/2026_08_31/golden/almanac_inputs/allobs_<tele>_<mjd>.h5`
  (extracted 2026-08-31 from the 2026_05_01 bulk almanac's `raw/` layout).
- Raw `.apz`: cca mirror `/mnt/ceph/users/sdssv/raw/APOGEE` (via `--cluster cca`).
- Cals: darks/flats `2025_07_31/outdir_ref/`, gain/read `2025_07_31/pass_clean/`.
- No exposure-classifier model (decision 2026-08-31: goldens match run_all.sh
  production behavior; the classifier enters later via R1/O5).

## Monitoring while it runs (filesystem only — no squeue/sacct polling)

Let `GROOT=/mnt/ceph/users/sdssv/work/asaydjari/2026_08_31/golden/ApogeeReduction.jl@<shortsha>`
(the shortsha of the submitted HEAD; the dry-run prints it).

- **sbatch stdout log** (predictable path): `slurm_logs/ar_goldens_<jobid>.out`
  under the submission cwd above. Per-day sentinels inside:
  `=== GOLDEN BEGIN <tele> <mjd> <timestamp> ===` /
  `=== GOLDEN END <tele> <mjd> exit=N elapsed=... ===`.
- **`$GROOT/RUNSTATE`** — rewritten at every day boundary: which day is
  running, N of 6 done.
- **`$GROOT/<tele>_<mjd>/STATUS`** — written when that day finishes: exit
  code, elapsed, start/end timestamps.
- **`$GROOT/<tele>_<mjd>/logs/run_testday_<tele>_<mjd>.log`** — the live
  per-day step log (step banners + warnings census at the end).
- **`$GROOT/ALL_DONE`** — exists only when the whole job finished; contains
  the overall exit and day count.

## After completion

1. Check `$GROOT/MANIFEST.md`: branch/sha/config header + per-day exit-status
   table — all exit codes must be 0.
2. Read each day's warnings census (end of the per-day log) against the
   2026_05_01 baseline counts (REFACTOR_PLAN v1 §0).
3. **Sanity cross-check on apo 59429**: the ccalin051 smoke run of the same
   pipeline code exists at
   `/mnt/ceph/users/sdssv/work/asaydjari/2026_08_31/golden_smoke/apo_59429/`.
   Run:

   ```bash
   cd test/regression
   julia +1.11.0 --project=. h5diff_tree.jl \
       $GROOT/apo_59429 \
       /mnt/ceph/users/sdssv/work/asaydjari/2026_08_31/golden_smoke/apo_59429 \
       --out $GROOT/apo_59429_vs_smoke.md
   ```

   Expected: **all files identical** (modulo the provenance ignore list) —
   only `test/regression/` files changed since the smoke's sha, so the
   pipeline code is the same; this also cross-validates workstation-vs-Slurm
   reproducibility. Any DIFFERS here means environment-dependent behavior and
   must be investigated before the goldens are trusted.
4. The goldens then serve as the diff base for every fix/cleanup branch (see
   README.md, "expected-diff-statement workflow").

## The A3 gapped-almanac day (deferred)

A full scan of the 2026_05_01 bulk almanac (both telescopes, all ~6,200
mjd groups) found **zero days where the `exposure` column is non-dense** —
consistent with almanac's `organize_exposures` hole-filling contract (v1 §0).
The A3 test day therefore has to be **synthesized**: copy a staged per-day
almanac, delete one mid-sequence exposure row, and run it as an extra day.
Do this when card A3 starts (it needs no golden baseline from real data);
it would become the 7th run in this framework.
