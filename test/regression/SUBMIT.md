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

## Incident log: job 6969533 (2026-08-31) — missing wavecal/ dir

Day 1 (apo 60000) of the first submission died at "Skyline medwave/skyline
dither" (`pipeline_2d_1d.jl:427`): `skyline_medwavecal_skyline_dither`
(`src/wavecal.jl:1708`) does `h5open(<outdir>/wavecal/wavecalNightAve_*.h5,
"w")` **without mkpath-ing `<outdir>/wavecal/`**. `run_bulk.sh:159` works
around it with `mkdir -p ${outdir}/wavecal`; `run_all.sh` omits it and
survives only because production outdirs already contain `wavecal/`; the
proper code fix is on `fix-broken-fiber-fluxing` (lands via R1). The apo
59429 smoke never hit it (no object exposures → early return before the
h5open). The confusing "HDF5.API.H5Error: Error getting stack length" is
secondary: HDF5.jl cannot render a worker-remoted H5Error on the main
process; the true error is the `h5f_create` failure in the "caused by" block.

Mitigations applied:
- `run_testday.sh` now does `mkdir -p ${outdir}wavecal` before the final
  2D→1D stage (mirroring run_bulk.sh).
- The `wavecal/` dirs of all 6 day outdirs under
  `golden/ApogeeReduction.jl@f76194a/` were pre-created at 20:25 while job
  6969533 was in day 2's 3D stage, so days 2–6 of that job are expected to
  complete; only **apo 60000 needs a redo**.

### Redoing apo 60000 into the existing golden tree

Pipeline code is identical between `f76194a` and the fix commit (only
`test/regression/` changed), so the redo may legitimately write into the
`@f76194a` root — the MANIFEST appends a "Redo/partial run" section recording
the new sha:

```bash
cd /mnt/home/asaydjari/gitcode/worktrees/AR-T1
AR_DAYS="apo 60000" \
AR_OUTROOT=/mnt/ceph/users/sdssv/work/asaydjari/2026_08_31/golden/ApogeeReduction.jl@f76194a \
sbatchAKS test/regression/submit_goldens.sh
```

(`sbatchAKS` must forward the environment, as plain `sbatch` does by
default; if unsure use `sbatch --export=ALL,AR_DAYS="apo 60000",AR_OUTROOT=...`.)
With `commit_exists` checkpointing the completed 3D→2D/traces/relFlux
products of the failed attempt are reused; expect ~15–25 min. Then run the
post-completion checks above, including the apo 59429 vs smoke cross-diff.

## The A3 gapped-almanac day (deferred)

A full scan of the 2026_05_01 bulk almanac (both telescopes, all ~6,200
mjd groups) found **zero days where the `exposure` column is non-dense** —
consistent with almanac's `organize_exposures` hole-filling contract (v1 §0).
The A3 test day therefore has to be **synthesized**: copy a staged per-day
almanac, delete one mid-sequence exposure row, and run it as an extra day.
Do this when card A3 starts (it needs no golden baseline from real data);
it would become the 7th run in this framework.
