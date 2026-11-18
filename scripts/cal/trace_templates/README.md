# Per-epoch trace templates: rebuilding the calibration before a release

## What it calibrates

The robust-trace design (study and PI decisions: `2026_11_17/trace_robust/TRACE_ROBUST_V2.md`,
Sec. 3, 5 and 7) places every fiber of a night from a per-epoch **template** of the slit: for each
telescope, chip and instrument epoch, the fiber centers and Gauss-Hermite widths along x and the
relative fiber throughput. A night's trace is the template registered to the night by a few rigid
numbers (shift quadratic in x, scale linear in x), refined per fiber where there is light. The
templates are used for the template-seeded fits, the identification check, the unlit fibers and
the fallback levels (T1-T3); that code is in `src/` (slit model, PR B), not here.

All templates and the epoch table live in **one** versioned HDF5 file for both telescopes and all
three chips. Format, validation, lookup and identity: `src/trace_template_cal.jl` (format version 1).

| item | |
|---|---|
| epoch table | per telescope, contiguous inclusive MJD ranges; the last epoch may be open-ended (`--open_last_epoch`, then `mjd_end = 99999`) |
| template unit | per (telescope, chip, epoch): `center`, `width` [32 bands x 300 fibers] (px, detector rows on the scale of `trace_params`), `center_poly` / `width_poly` (per fiber, in u = (x - 1024)/1024), `height_rel` [300], `lit` [300], the input flats and their registration (d, s at x = 1024) |
| identity | per unit `trace_template_id = <tele>_<chip>_<epoch>:<fingerprint>` (FNV-1a 64 over the unit's arrays and epoch row, not the file bytes); per file `<calibration_version>:<content fingerprint>`; both checked on load |
| provenance | builder attributes per unit (inputs, all-lit counts, modal dark fibers, bootstrap flat, iterations, drift, code SHA-256), packing provenance (epoch table + SHA-256, unit directory, commit, choices, missing units) |

## Pipeline switch

`--trace_template_cal <path>|none` on **both** `scripts/cal/make_traces_from_flats.jl` and
`pipeline_2d_1d.jl` (use the same file for both).

- **OFF by default (`none`).** Traces are made exactly as in PR A (template-free). Whether it
  becomes ON by default is a PI decision; the switch is ready.
- Shell drivers: `AR_TRACE_TEMPLATE_CAL` in `scripts/bulk/run_bulk.sh` and
  `scripts/daily/run_all.sh` (unset = no flag = off; set = passed to the trace step and to every
  `pipeline_2d_1d.jl` call). Airflow: `TRACE_TEMPLATE_CAL` in `airflow/dags/ar_common.py`
  (`"none"`), passed by `apogee_daily.py` to every trace-making step and as
  `AR_TRACE_TEMPLATE_CAL` to `run_all.sh`.
- **Hard errors, never a silent fallback to the template-free traces:** a configured path that
  does not exist (at startup), a file that cannot be read, has an unknown format version, fails
  validation (epoch coverage, shapes, finite values, fibers ordered along the slit in every band)
  or whose content does not match its stored fingerprint.
- A file may be **incomplete** (some units missing, listed in its provenance and refused by
  `load_trace_template_cal(...; require_complete = true)`); lookup then returns `nothing` for
  those (telescope, chip, epoch), i.e. template-free behaviour there.
- **Checkpointing.** Every traceMain and every per-flat trace QA record carries
  `trace_template_id` (`"none"` when off), plus `trace_template_cal_id`, `trace_template_cal`,
  `trace_template_epoch` and `trace_template_provisional`. A checkpointed traceMain or per-flat
  trace whose recorded id differs from the configured one is recomputed, under `commit_same` and
  `commit_exists` alike. Changing ONE unit only recomputes the nights that unit serves.
  - Products downstream of the traceMain (ar1D, wavelength solutions) are checkpointed on the
    code commit only. After changing the template file at the same code commit, re-run
    `pipeline_2d_1d.jl` with `--checkpoint_mode clobber` for the affected nights.

## Bootstrap: rebuilding from nothing

The file is rebuilt from the data, never from a previous file (AKS: every learned calibration
quantity ships with scripts that rebuild it from the data before each release, and the pipeline
must run without the file).

1. **Reduce with the switch off**: `--trace_template_cal none` (`AR_TRACE_TEMPLATE_CAL=none`).
   This gives the 2D quartz flats (`apred/<mjd>/ar2Dcal_*_quartzflat.h5`) and template-free
   traceMains for every night. Nothing here needs a template.
2. **Epochs** (`propose_trace_epochs.jl --root <reduction> --outdir <dir>`): change points in the
   nightly registration series of the traceMains (shift, rotation, bow, scale, lit-fiber pattern;
   defaults k = 10 nights, thr = 8, floor = 0.15 px at the chip edge). Writes
   `epoch_candidates.csv` (evidence per boundary) and a **draft** `trace_epochs_draft.csv`.
   A person confirms or edits it: in the first table the two "pattern" boundaries caused by the
   2017 dark bundle were replaced by the shutdown boundary 57969, and five small no-gap steps
   (APO 59534, 59816, 60121, 60237; LCO 60249) are real (> 20 sigma) but their cause is unknown.
3. **Build and pack** (`run_trace_templates.sh <root> <workdir> <version> <epochs.csv>`, Slurm or
   local): one `build_trace_template.jl` per (epoch, chip), then `pack_trace_templates.jl`, then
   `check_trace_templates.jl`. Resumable (finished units are skipped).
4. **Install** by copying `<workdir>/<version>.h5` to `cal_ref/trace_templates/` and pointing
   `--trace_template_cal` / `TRACE_TEMPLATE_CAL` at it.

`test/trace_templates_bootstrap.jl` runs steps 3's build and pack on synthetic flats with no
template file anywhere and checks the recovered slit (< 0.05 px max, < 0.01 px rms after the gauge).

## The builder (`build_trace_template.jl`, one unit)

Algorithm of the study's `build_trace_templates.jl` (TRACE_ROBUST_V2 Sec. 7.2), with the fitter
vendored as `TemplateFit.jl`:

1. Candidates: the quartz flats of every night in the epoch (directory listing, or `--flat_list`),
   one per night (the last exposure), `--nmax` nights spread over the epoch (0 = all).
2. **All-lit flats only**, decided template-free: central-band peak census, gaps turned into
   missing slots; detected + missing = 300 and the missing set equal to the epoch's modal set
   (fibers dark on every flat are physically absent; any additional dark fiber refuses the flat).
3. Bootstrap from the highest-S/N all-lit flat (slots numbered 1..300 from the bottom of the chip).
4. **Drift within an epoch** (design Sec. 7.3): each flat is registered (shift quadratic, scale
   linear in x) and its free per-band centers are mapped back to the common frame before the
   per-fiber median, so rigid drift never enters the template; iterate to a p99 change < 0.003 px.
   Pattern drift (last third minus first third of the epoch, px rms) is measured and flagged above
   `--drift_warn` (0.01 px); the remedy (split the epoch) is a human decision. The per-flat
   registration is kept (`input_d_center`, `input_s_center`) for a later time-dependent T3 level.
5. **In-focus flats only** (`--fit_rms_max`, default 0.02 px = the trace gate's threshold): flats
   whose template fit is worse are dropped and the template rebuilt. In the first full build this
   would have removed 6 of 13,685 inputs, all on known problem nights (LCO 57803, 57813, 60702).
6. Validation: leave-nights-out and no-template mode (`--validate`).

Packing (`pack_trace_templates.jl`): the fibers dark in an epoch are never fitted by the builder;
their positions are re-placed between their lit neighbours when the builder never fitted them
(`height_rel` < 0.05) or they left their slot (> 0.5 spacing) (`--unlit_positions interpolate`,
`interpolate_unlit_template`). Measured on the first build: the APO dead fibers (24; 136, 276,
280) were frozen at their central-band position, 1-2.4 px off at the chip edges; on APO G apoE6
fibers 276/280 sat ~20 px (3 spacings) low, on lit neighbours. LCO fiber 270 (faint, lit,
`height_rel` 0.08-0.10) was fitted and is kept.

## QA (`check_trace_templates.jl`)

Per unit: coverage, rigid drift (NMAD of the per-flat shift), pattern drift, residual outliers,
leave-nights-out and no-template statistics; held-out flats above `--max_dc_flag` (0.3 px) are
listed; `--reference` compares with a previous file. The largest single-flat differences of the
first build (up to 2.5 px on APO R, 0.8 px on APO G, 1.6 px on LCO R) were diagnosed with
per-cell refits: they sit on ONE fiber at the chip edge where that fiber has little or no light
(APO R: the lowest lit fiber at x < 300; APO G: fibers 252/256/261 at x > 1890; LCO R: fiber 300
at x > 1800), recur on every flat of every epoch, and are where the template, the no-template fit
and production all disagree; not bad flats and not misplaced epoch boundaries. See the PR
description for the numbers.

## Cost

First full build (pass-1 DR21, 14 epochs x 3 chips = 42 units, every night, one exclusive genoa
node, one unit per core, with validation): 3.95 h wall time = ~4 node-hours, set by the largest
unit (APO 57969-59442, 1,146 candidate / 1,052 all-lit flats). Units of 150-600 flats take 15-130
min (~11 s per input flat including validation). Packing and QA take a few minutes each.

## Open items

- Switch the builder from the vendored `TemplateFit.jl` to the pipeline's slit model once that
  is validated to reproduce the templates (compare with `check_trace_templates.jl --reference`).
- Time-dependent epoch-median registration for T3 in drifting epochs (design Sec. 7.3 iii).
- The per-fiber residual fit at the chip edges where a fiber has no light (QA above).
