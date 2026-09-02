# Golden-diff regression harness (task T1)

Tools to (1) reduce a single test day with any ApogeeReduction.jl checkout on a
workstation and (2) compare two output trees dataset-by-dataset, so that every
bug fix ships with a written *expected-diff statement* and every pure cleanup
proves it changed nothing. See
`/mnt/ceph/users/sdssv/work/asaydjari/2026_07_19/REFACTOR_PLAN.md` §4 and
`/mnt/ceph/users/sdssv/work/asaydjari/2026_08_31/REFACTOR_PLAN_v2.md` §2.

## Contents

| file | purpose |
|---|---|
| `run_testday.sh` | run the AR pipeline (runlist → 3D→2D → traces/relFlux → 2D→1D) for one `(tele, mjd)` from existing raw data + almanac, no Slurm, no Utah tunnel — the single-day tool for fix-branch runs |
| `h5diff_tree.jl` | walk two output trees, compare every HDF5 dataset + attribute, emit a markdown diff report |
| `submit_goldens.sh` | sbatch script generating ALL golden baselines as ONE true bulk run over the combined test-day exposure set (the exact run_bulk.sh chain), with per-day bookkeeping afterwards; see `SUBMIT.md` |
| `concat_runlists.jl` | merge per-telescope runlists into the one combined runlist the bulk chain uses (the test days are tele-specific pairs; run with `--project=<AR checkout>`) |
| `SUBMIT.md` | golden-job submission instructions + filesystem-only monitoring + post-run checks |
| `Project.toml` / `Manifest.toml` | Julia env for `h5diff_tree.jl` (HDF5 + ArgParse only — independent of the AR package env; `concat_runlists.jl` runs with the AR project instead) |

One-time setup: `julia +1.11.0 --project=. -e 'using Pkg; Pkg.instantiate()'`.

## Running a test day

```bash
cd test/regression
./run_testday.sh <AR_checkout_dir> <tele> <mjd> <outdir>
# e.g.
./run_testday.sh ~/gitcode/worktrees/AR-T1 apo 59429 \
    /mnt/ceph/users/sdssv/work/asaydjari/<date>/golden/<repo>@<sha>/apo_59429/
```

All paths are configurable by env var (see the header of `run_testday.sh`):
`AR_ALMANAC_SRC` (bulk almanac; default the 2026_05_01 run's
`allobs_57600_61160.h5`), `AR_RAW_CLUSTER` (`cca` → raw `.apz` mirror
`/mnt/ceph/users/sdssv/raw/APOGEE`, or an explicit base path),
`AR_CALDIR_DARKS` / `AR_CALDIR_FLATS` (default `2025_07_31/outdir_ref/`),
`AR_GAIN_READ_CAL_DIR` (default `2025_07_31/pass_clean/`), `AR_WORKERS`
(default 24), `AR_JULIA_VERSION` (default 1.11.0, matching run_all.sh),
`AR_CHECKPOINT_MODE`, `AR_CHIPS`, `AR_EXP_CLASS_MODEL`. Set
`AR_TESTDAY_CONFIG=<file>` to source a config file with those assignments.
`AR_SLURM` (auto/true/false) selects local Distributed workers vs the
run_bulk.sh SlurmClusterManager pattern; inside an sbatch allocation it
defaults to Slurm mode (`submit_goldens.sh` relies on this).

The full step log lands in `<outdir>/logs/run_testday_<tele>_<mjd>.log` and
ends with a **warnings census** (counts of the known warning classes from the
2026_05_01 bulk log — regression metrics per v1 §0).

Differences from `scripts/daily/run_all.sh` (deliberate):

- **No almanac invocation.** The bulk `raw/`-layout almanac file
  (`AR_ALMANAC_SRC`, one file covering all days) is consumed DIRECTLY:
  since R1 (#365) `read_almanac_exp_df` reads the `raw/` layout natively, so
  the harness just symlinks the bulk file to
  `<outdir>almanac/allobs_<tele>_<mjd>.h5` (the path every downstream stage
  derives from outdir+runname; all stages open it read-only) and selects the
  day with the runlist makers' `--mjds` flag (a comma-separated MJD list;
  `submit_goldens.sh` passes each telescope's full test-day list).
  *Historical note*: the `@f76194a` golden baselines predate this — they were
  generated via per-day extracts made by a since-deleted
  `extract_almanac_day.jl` (their MANIFEST records the staged
  `almanac_inputs/` directory). The extracted per-day groups were byte-copies
  of the bulk file's groups, so runlists and products are unaffected by the
  switch.
- **No sdsscore update, no plots/dashboard/arMADGICS** — the golden diff
  compares science products only. (arM regression is run separately.)
- **Slurm env is scrubbed** so `pipeline.jl` / `pipeline_2d_1d.jl` use local
  `addprocs` instead of SlurmClusterManager, and `--workers_per_node` is passed
  explicitly (their default `-1` means "all cores" only under Slurm; outside
  Slurm it would call `addprocs(-1)` and crash).

Known workstation caveats (pipeline code, not harness — candidates for later
cards): `make_traces_from_flats.jl` hardcodes `addprocs(16)` and
`make_relFlux.jl` hardcodes `addprocs(64)` when not under Slurm
(oversubscribes a 32-core node; harmless for light nights);
`make_relFlux.jl` unconditionally constructs a `SlackThread()` (no-op warning
without Slack credentials).

## Comparing two output trees

```bash
julia +1.11.0 --project=. h5diff_tree.jl <goldenDir> <newDir> --out report.md
```

Files are paired by relative path (symlinks followed — outdirs symlink
darkRate/flatFraction cals into `apred/<mjd>/`). The harness-staged almanac
INPUT (`almanac/allobs_<tele>_<mjd>.h5`) is excluded from the pairing by
default — since the bulk-almanac switch it is a symlink to the multi-GB
multi-day bulk file and an input, not a product (`--skip-files` overrides;
the runlist files, which ARE products of the runlist makers, are still
compared). Per-dataset verdicts:

- `identical` — bit-for-bit (elementwise `isequal`: NaN==NaN counts as equal;
  NaN-pattern *changes* are counted and always fail tolerance);
- `within-rtol` — every differing element satisfies
  `|a-b| <= atol + rtol*max(|a|,|b|)` (`--rtol/--atol`, default 0 = exact);
  max abs/rel diffs reported;
- `DIFFERS` — with localization: element count, max-diff location/values, and
  per-axis distinct-index summaries (axes of length 300 labeled `fiber`, 8700
  `unipix`, 2048 `xpix` per the product schemas in v1 §0);
- `missing-in-golden` / `missing-in-new` — for datasets, attributes, or whole
  files.

Provenance metadata that legitimately differs between runs is **ignored by
default** (still compared, reported as `IGNORED`, never fails the diff):
`metadata/{git_commit, git_branch, git_clean}`,
`metadata/trace_orig_param_fname` (embeds the checkout path), and the
top-level `trace_used_param_fname` (embeds the run's own outdir path;
present in ar1Dcal/ar2Dresidualscal). Extend with
`--ignore pat1,pat2` (path-suffix match, attrs as `path@attrname`); disable the
defaults with `--no-default-ignores`. NOTE: `metadata/mjd_mid_exposure*`,
`ndiff_used`, `nread_total` are exposure properties, NOT provenance — if they
differ between two runs of the same exposure something is genuinely wrong.

Exit code: 0 when everything is identical/within-rtol (modulo ignores), 1 when
anything DIFFERS or is missing — usable directly in CI/scripts.

**Per-day diffs against a bulk-mode tree**: the bulk goldens live in ONE
outdir whose `apred/<mjd>/` subdirs are per-day, while a single-day
run_testday tree nests the same layout under its own outdir — so point the
walker at the `apred/<mjd>` SUBTREES:

```bash
julia +1.11.0 --project=. h5diff_tree.jl \
    <bulkGoldenRoot>/apred/59429 <newSingleDayOutdir>/apred/59429 --out report.md
```

At mjd 60000 both telescopes share `apred/60000/`; when the other side is
apo-only (or lco-only), exclude the other telescope's files with
`--skip-files "_lco_"` (resp. `"_apo_"`) — all product filenames are
tele-tagged. Nightly wavecal solutions (`wavecal/wavecalNightAve_*`) and the
per-mjd `dome_flats/<mjd>`/`quartz_flats/<mjd>` trace dirs can be compared
the same way (subtree or `--skip-files` on the other days' `_<mjd>` tags).

## The expected-diff-statement workflow (v1 §1.3–1.4, §4.2)

1. **Baselines**: generate goldens ONCE from the pre-fix code
   (`exposure-type-classifier` HEAD for AR) on the test days, stored under
   `/mnt/ceph/users/sdssv/work/asaydjari/<date>/golden/<repo>@<sha>/<tele>_<mjd>/`.
2. Every **behavioral fix** ships with a written expected-diff statement of the
   form *"only datasets X for fibers satisfying Y change, nothing else"* —
   e.g. A1: "on apo 58588 only fiber 211's flux_1d/ivar_1d rows change, apo
   60000 is bit-identical." Run the fix branch on the test days, diff against
   the goldens, and check the report matches the statement (the per-axis
   `fiber` localization lines are what you check against Y).
3. Every **pure cleanup** requires an all-identical report on all test days;
   float-reassociation waivers must be explicit per dataset (use `--rtol` with
   a written justification).
4. When a fix legitimately changes outputs everywhere, regenerate the goldens
   (new `<repo>@<sha>` directory) and record the transition in the PR.

### Test days (v1 §4.2)

| tele/mjd | why |
|---|---|
| apo 60000 + lco 60000 | healthy full product set; the "nothing should change" day |
| apo 58588 | dead fiber 211/adjfib 90 (relthrpt 1.4e-4) → A1/M3 |
| apo 59337 | aborted exposure 55 → short-exposure/DCS paths, M1 snr |
| apo 58011 | negative relthrpt fiber → A1 sign handling |
| apo 59429 | single arclamp exposure only → missing-flats/wavecal fallbacks, A5 |
| gapped-almanac night | for A3 (synthesize by deleting a row in a copied almanac if none exists) |
