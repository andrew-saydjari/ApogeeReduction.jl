# Submitting the golden-baseline generation job

One multinode Slurm job generates all T1 golden baselines (6 test days:
apo+lco 60000, apo 58588, 59337, 58011, 59429) as **ONE true bulk run**: a
single combined runlist over every test-day exposure, then the exact
run_bulk.sh chain (per-tele pipeline.jl → flat-type cal loops → per-tele
pipeline_2d_1d.jl) — not per-day invocations, and deliberately NOT a job
array. Products land in one outdir whose `apred/<mjd>/` subdirs are per-day
as usual; a post-run bookkeeping step writes per-day STATUS files, the
MANIFEST table, and the ALL_DONE sentinel. Submission is **manual by
design.**

## Submit

```bash
cd /mnt/home/asaydjari/gitcode/worktrees/AR-T1     # the T1 branch checkout
mkdir -p slurm_logs                                # sbatch will not create it
sbatchAKS test/regression/submit_goldens.sh
```

Then **paste the job ID back into the Claude session** so progress can be
followed from the filesystem (no Slurm polling needed).

Allocation requested by the script: `-p cca -N 2
-C "[genoa|icelake|rome]" --mem=900G -t 1-12:00 -J ar_goldens`.
Expected wall time comparable to the old sequential mode (**roughly 6–10 h
total** for the full day set); the 36 h limit is generous margin.

Before submitting, optional sanity checks (no Slurm needed):

```bash
bash -n test/regression/submit_goldens.sh
test/regression/submit_goldens.sh --dry-run   # prints the bulk stage sequence
```

The job refuses to start from a dirty git tree (goldens must be attributable
to a pushed sha) and verifies the almanac file exists before running
anything.

## Inputs (zero DB/tunnel dependency)

- Almanac: the bulk `raw/`-layout file `AR_ALMANAC_SRC`
  (default `/mnt/ceph/users/sdssv/work/asaydjari/2026_05_01/outdir/almanac/allobs_57600_61160.h5`),
  consumed directly — see README.md. (The `@f76194a` goldens instead used
  per-day extracts under `golden/almanac_inputs/`; their MANIFEST records
  this. They also predate the 2026-08-31 main merge: A2/#369 changes
  domeflat `ivarimage`, which feeds trace fitting and legitimately cascades
  into most 1D products on flat-bearing days — the 2026-09-02 bulk-mode
  validation (`2026_09_02/t1_bulkmode_val/EXPECTED_DIFF.md`) proved this
  drift is merge/machine, not bulk-mode: bulk vs single-day on the same
  machine was bit-identical for both apo 59429 and apo 58011. Fresh bulk
  goldens from this script supersede `@f76194a`.)
- Raw `.apz`: cca mirror `/mnt/ceph/users/sdssv/raw/APOGEE` (via `--cluster cca`).
- Cals: darks/flats `2025_07_31/outdir_ref/`, gain/read `2025_07_31/pass_clean/`.
- No exposure-classifier model (decision 2026-08-31: goldens match run_all.sh
  production behavior; the classifier enters later via R1/O5).

## Monitoring while it runs (filesystem only — no squeue/sacct polling)

Let `GROOT=/mnt/ceph/users/sdssv/work/asaydjari/2026_08_31/golden/ApogeeReduction.jl@<shortsha>`
(the shortsha of the submitted HEAD; the dry-run prints it).

- **sbatch stdout log** (predictable path): `slurm_logs/ar_goldens_<jobid>.out`
  under the submission cwd above (mirrored to
  `$GROOT/logs/submit_goldens_<runname>.log`). Per-STAGE sentinels inside:
  `=== STAGE BEGIN <name> <timestamp> ===` /
  `=== STAGE END <name> elapsed=... ===`, and `=== STAGE FAILED <name> ===`
  if the chain dies (bulk mode is one chain: a failed stage aborts the rest —
  the last BEGIN without END/FAILED marks the crash point).
- **`$GROOT/RUNSTATE`** — rewritten at every stage boundary: which stage is
  running (or FAILED in which stage).
- **`$GROOT/status/STATUS_<tele>_<mjd>`** — written by the post-run
  bookkeeping step: per-day product counts (ar2D/ar2Dcal/ar1Dcal/ar1Dunical),
  the expected-vs-found missing-file check, and the per-day diff base paths.
- **`$GROOT/ALL_DONE`** — exists only when the whole job finished; contains
  the overall bookkeeping exit.

## After completion

1. Check `$GROOT/MANIFEST.md`: branch/sha/config header + the per-day product
   table — every verdict must be `ok` (product counts match the runlist's
   per-day exposure counts). The same numbers live in
   `$GROOT/status/STATUS_<tele>_<mjd>`.
2. Read the warnings census (end of the bulk log) against the 2026_05_01
   baseline counts (REFACTOR_PLAN v1 §0).
3. **Per-day diffs** point at the `apred/<mjd>` subtrees (see README.md):

   ```bash
   cd test/regression
   julia +1.11.0 --project=. h5diff_tree.jl \
       $GROOT/apred/<mjd> <candidateOutdir>/apred/<mjd> --out report.md
   ```

   At mjd 60000 both telescopes share `apred/60000/` — against an apo-only
   (lco-only) candidate add `--skip-files "_lco_"` (`"_apo_"`).
4. **Sanity cross-check on apo 59429**: the ccalin051 smoke run of the same
   pipeline-code lineage exists at
   `/mnt/ceph/users/sdssv/work/asaydjari/2026_08_31/golden_smoke/apo_59429/`.
   Run:

   ```bash
   cd test/regression
   julia +1.11.0 --project=. h5diff_tree.jl \
       $GROOT/apred/59429 \
       /mnt/ceph/users/sdssv/work/asaydjari/2026_08_31/golden_smoke/apo_59429/apred/59429 \
       --out $GROOT/apo_59429_vs_smoke.md
   ```

   Expected: **all files identical** (modulo the provenance ignore list) —
   only `test/regression/` files changed since the smoke's sha, so the
   pipeline code is the same; this also cross-validates workstation-vs-Slurm
   reproducibility. Any DIFFERS here means environment-dependent behavior and
   must be investigated before the goldens are trusted.

   **Finding (2026-09-02, T1-rework validation)**: this cross-check DIFFERS
   on the three `fpiPeaks_*` files (`fpi_line_mat`, `fpi_line_cov_mat`,
   `fpi_line_trace_centers`) — the FPI peak fitting is environment-dependent
   (Slurm compute node vs ccalin051), though deterministic per machine: two
   ccalin051 runs (pre-merge smoke vs post-merge-validation) were
   bit-identical there. Consequence: fpiPeaks golden diffs are only
   meaningful when golden and candidate ran on the same machine class; every
   other apo 59429 product is machine-independent. Reports:
   `2026_09_02/t1_rework_val/diff_golden_vs_smoke_apo_59429.md`. Candidate
   for its own refactor-plan card.
5. The goldens then serve as the diff base for every fix/cleanup branch (see
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
- `run_testday.sh` gained `mkdir -p ${outdir}wavecal` before the final
  2D→1D stage (mirroring run_bulk.sh). (Removed again after R1/#365 merged
  the proper `mkpath` into `skyline_medwavecal_skyline_dither`.)
- The `wavecal/` dirs of all 6 day outdirs under
  `golden/ApogeeReduction.jl@f76194a/` were pre-created at 20:25 while job
  6969533 was in day 2's 3D stage, so days 2–6 of that job are expected to
  complete; only **apo 60000 needs a redo**.

### Partial/redo runs (bulk mode)

`AR_DAYS="apo 60000;apo 58011"` restricts the ONE bulk run to a day subset
(the combined runlist then covers only those days); `AR_OUTROOT` points it
at an existing outdir — with `commit_exists` checkpointing, products that
already exist are reused, so a redo reprocesses only what is missing. The
MANIFEST appends a "Redo/partial run" section recording the new sha.
(`sbatchAKS` must forward the environment, as plain `sbatch` does by
default; if unsure use `sbatch --export=ALL,AR_DAYS=...,AR_OUTROOT=...`.)
NOTE: the old per-day-tree goldens under `@f76194a/<tele>_<mjd>/` predate
bulk mode and cannot be extended by it — a bulk redo belongs in a bulk-mode
outroot.

## The A3 gapped-almanac day (deferred)

A full scan of the 2026_05_01 bulk almanac (both telescopes, all ~6,200
mjd groups) found **zero days where the `exposure` column is non-dense** —
consistent with almanac's `organize_exposures` hole-filling contract (v1 §0).
The A3 test day therefore has to be **synthesized**: copy one day's group out
of the bulk almanac into a small `raw/`-layout file, delete one mid-sequence
exposure row, and run it as an extra day (pointing `AR_ALMANAC_SRC` at it).
Do this when card A3 starts (it needs no golden baseline from real data);
it would become the 7th run in this framework.
