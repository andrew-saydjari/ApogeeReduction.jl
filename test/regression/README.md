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
`AR_GAIN_READ_CAL_DIR` (default `2026_09_06/pass_clean/`), `AR_WORKERS`
(default 24), `AR_JULIA_VERSION` (default 1.11.0, matching run_all.sh),
`AR_CHECKPOINT_MODE`, `AR_CHIPS`, `AR_EXP_CLASS_MODEL`. Set
`AR_TESTDAY_CONFIG=<file>` to source a config file with those assignments.
`AR_SLURM` (auto/true/false) selects local Distributed workers vs the
run_bulk.sh SlurmClusterManager pattern; inside an sbatch allocation it
defaults to Slurm mode (`submit_goldens.sh` relies on this).

The full step log lands in `<outdir>/logs/run_testday_<tele>_<mjd>.log` and
ends with a **warnings census** (counts of the known warning classes from the
2026_05_01 bulk log — regression metrics per v1 §0) followed by a
**warnings triage** stage.

### Warnings census vs. warnings triage

The census greps a fixed list of hand-named substrings. It is a convenience
view and it has two blind spots by construction: it cannot count a category
nobody named, and if a pattern drifts from the message text it reports `0`
rather than an error. Both bit us in job 6980442 — 88 "lamp turned off"
warnings were entirely uncounted, and `"no useful relfluxing files"` reported
`0x` against a real count of 3 because the message says *any* useful
relfluxing files.

`warnings_triage.sh` is the authoritative inventory. It enumerates every Julia
`@warn`/`@error` record in a log and groups them by **emit site**
(`Module src/file.jl:LINE`), which is exact, survives message-text edits, and
survives ProgressMeter re-renders being glued onto the message text.

```bash
# inventory for a run, with the per-(tele, mjd) distribution per site
test/regression/warnings_triage.sh --units <outdir>/logs/

# inventory + diff against the reference set (exit 1 iff a NEW site appeared)
test/regression/warnings_triage.sh -r test/regression/warnings_reference_6980442.tsv <log>

# start a fresh reference set (verdict/expect/note are filled in by hand)
test/regression/warnings_triage.sh --emit-reference <log> > new_reference.tsv
```

`warnings_reference_6980442.tsv` is the adjudicated reference set from the
200-MJD DR21 testbed (slurm 6980442, AR @ 824fd978): 454 records across 7 emit
sites, each with a verdict (`EXPECTED` / `ACTIONABLE` / `UNKNOWN`) and an
`expect` field saying how the count should move in the next run. Override the
file the harness diffs against with `AR_WARN_REFERENCE=<file>`.

**It is a reference set, not a baseline of acceptability.** Nobody had
adjudicated the previous run's warnings before 2026-09-07; the counts were
simply what happened. A category earns `EXPECTED` from a demonstrated
mechanism — the code path, the data condition that triggers it, and evidence
that the condition is legitimately present in this data — never from having
been there before. `UNKNOWN` is the correct verdict for anything whose
mechanism has not been established, and the diff re-prints those every run so
they do not quietly become furniture. Three of the seven sites in the current
reference set are `ACTIONABLE`; the adjudication and its evidence are at
`/mnt/ceph/users/sdssv/work/asaydjari/2026_09_07/warnings_baseline/WARNINGS.md`.

The triage stage never fails a run — a new emit site is for a human to read,
not a reason to discard ten hours of reduction.

### Cross-checking warnings against the exposure-type classifier

A warning tells you the pipeline objected. It does not tell you whether the
pipeline was *right*. The exposure-type classifier is the one other
per-exposure opinion available, and the goldens/testbed runs deliberately run
with `AR_EXP_CLASS_MODEL=""` (matching production `run_all.sh`), so its zero
warning count in a run log is structurally uninformative — it could not have
fired. That is exactly the kind of zero WARNINGS.md §7 warns against reading as
health.

Two scripts recover the signal *after* the run, without touching it:

```bash
# 1. classify every delivered exposure (read-only over <outdir>/apred/)
julia --project=. test/regression/classify_run.jl \
    --outdir  <outdir> \
    --model   <exposure_classifier_rf_v6.jld2> \
    --output  predictions.tsv \
    --nworkers 10

# 2. resolve the run's warnings to the exposure level
test/regression/warnings_triage.sh --exposures <outdir>/logs/ > warnings_by_exposure.tsv

# 3. join them
julia test/regression/classifier_crosscheck.jl \
    --predictions predictions.tsv --warnings warnings_by_exposure.tsv \
    --focus apo:60255:7
```

`classify_run.jl` calls the same `ApogeeReduction` functions the in-pipeline
check calls (`exposure_class_features`, `classify_exposure_type`,
`exposure_check_category`, `exposure_predicted_bad`) and reproduces
`pipeline.jl`'s `persistence_prior` post-step, so the post-hoc pass and the
in-pipeline check cannot silently diverge. Running it after the fact is not a
compromise: it keeps the classifier an *independent* second opinion on the
warnings rather than a co-author of the same log, and it leaves the run itself
comparable to previous runs.

Cost, MEASURED on ccalin051 over the 2026_09_03 testbed: ~2 s per exposure per
core, I/O-bound on three ar2D reads. 17 911 exposures took ~1 h at
`--nworkers 10`. No Slurm allocation is needed at this scale; `--append`
resumes a killed sweep.

`classifier_crosscheck.jl` prints a contingency table — CONFIRMED (warning and
classifier agree), WARN-ONLY (candidate false-positive warning), CLF-ONLY
(candidate missed detection), quiet — plus the per-emit-site breakdown, the
`predicted_bad` vs almanac `flagged_bad` divergence, and declared-vs-classified
type disagreements. It deliberately computes no accuracy or F-score: there is
no ground truth here, both signals are estimates, and scoring one against the
other would assume the classifier is right.

Read the quadrants with the classifier's blind spots in mind. It is a random
forest over 67 whole-chip summary statistics, so it is blind to per-fiber and
per-chip effects; a warning about one fiber or one chip can be entirely real
while the exposure classifies fine. WARN-ONLY is a *candidate* false positive,
never a refutation.

The model artifact is not in this repo (it is ~175 MB). The current one is
`exposure_classifier_rf_v6.jld2` under the 2026_07_14 scratch dir, with its
training pipeline and provenance in that directory's `README.md`.

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
