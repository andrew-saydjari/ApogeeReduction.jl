# ApogeeReduction  <img src="docs/src/assets/logo.png" alt="AR Logo" width="100" align="right"/>

[![Build Status](https://github.com/andrew-saydjari/ApogeeReduction.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/andrew-saydjari/ApogeeReduction.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![codecov](https://codecov.io/gh/andrew-saydjari/ApogeeReduction.jl/branch/main/graph/badge.svg?branch=main)](https://codecov.io/gh/andrew-saydjari/ApogeeReduction.jl)

## Files
The pipeline produces files at many stages of reduction.
- `ar3Dcal`: Raw 3D datacubes. These are zero point adjusted 3D inputs into the photoelectron rate extractions and are experimental (not always created).
- `ar2D`: 2D images after 3D→2D extraction, before calibration
- `ar2Dcal`: 2D calibrated images after dark subtraction and flat fielding
- `ar2Dresiduals`: Residuals from 2D extraction process
- `ar1D`: Extracted 1D spectra for each fiber, in detector pixel units (before wavelength calibration and resampling)
- `ar1Duni`: 1D spectra resampled onto a uniform wavelength grid
- `ar1Dunical`: Flux (relative) calibrated 1D spectra on the uniform wavelength grid

## Structure

There are four main types of files in this repository:
- scripts/\*/run_*.sh : wrapper scripts to run the pipeline (submit job, determine resources, etc.)
- scripts/\*/make_runlist_*.sh: interface to get the data to run the pipeline on
- pipeline_*.sh: how functions combine to process the data
- src/*.jl: core functions of the repository

### File Structure
```
├── src/ : core functions of the repository
├── scripts/ : scripts for running the pipeline
│   └── run/ : scripts general users will interact with to run the pipeline
│   └── cal/ : scripts to build the calibrations files
├── test/ : test files for the repository (name matched to the src/ files they test)
├── metadata/ : metadata files for the repository (mostly dates for instrument changes/special calibrations runs)
├── data/ : input data (e.g. sky line lists from HITRAN)
├── dags/ : dags to run the pipeline through Airflow automations
├── pipeline.jl : main pipeline function (3D → 2D)
└── pipeline_2d_1d.jl : 2D pipeline function (2D → 1D)
```

### Call Structure

Nightly Runs:
```
└── run_all.sh : run all the data for a given night
    ├── almanac: queries database containing targeting information and data transfer status
    ├── make_runlist_all.sh: convert almanac output into a runlist interpreted by the pipeline
    ├── pipeline.sh: reduces data from raw type (3D compressed) to 2D calibrated data
    ├── run_trace_cal.sh: extracts the traces from domeflats to define 1D extraction profiles
    │   ├── almanac
    │   ├── make_runlist_dome_flats.sh: scrape almanac outputs for dome flats
    │   ├── pipeline.sh
    │   └── make_traces_domeflats.jl: extracts/saves traces from dome flats via gaussian fits to the "y" direction
    ├── pipeline_2d_1d.sh: extracts and calibrates 1D spectra from 2D calibrated data
    └── plot_all.sh: makes end of night plots for validation/QA and posts them to Slack
```

Bulk reprocessing workflow is still TBD, but the massive parallelization we have designed even for nightly runs means it should be similar, with possible interruptions to build higher signal to noise calibrations based on combining many calibration exposures.

## Current Flag Bits

Certain pixels are entirely masked or have data of questionable quality. This pipeline bit gives insight into the root cause of why this (tiny fraction of the) data is unable to be processed.

| Bit   | Value     | Meaning     |
| ----- | --------- | ----------- |
| -     | 0         | No problems       |
| 0     | 1         | reference array pixels |
| 1     | 2         | reference pixels |
| 2     | 4         | bad reference pixels |
| 3     | 8         | pixels not dark corrected |
| 4     | 16        | pixels with negative dark current |
| 5     | 32        | pixels with large dark current |
| 6     | 64        | flat response too low |
| 7     | 128       | one diff was dropped because it is a likely cosmic ray |
| 8     | 256       | more than one diff was dropped because they were likely cosmic rays (sus) |
| 9     | 512       | bad linear SUTR chi2 |
| 10    | 1024      | failed 1D extraction |
| 11    | 2048      | no nearby good pixels in 1D extraction |
| 12    | 4096      | neff>10 in 1D extraction |
| 13    | 8192      | pixel partially saturated |
| 14    | 16384     | pixel fully saturated |

### Per-fiber relative-throughput bits (`bitmsk_relthrpt`)

A **different axis** from the per-pixel mask above and from the exposure-level
`exp_class_*` metadata: this one is per FIBER and per EXPOSURE. It is produced by
`get_relFlux` (`src/ar1D.jl`) from a dome flat, stored in the `relFlux_*` cal
products, copied into each chip's `ar1D*` product by `process_1D`, and stacked
into the resampled `ar1Duni*` products as an `(N_CHIPS, N_FIBERS)` array by
`reinterp_spectra`, alongside the `relthrpt` values themselves.

| Bit | Value | Constant | Meaning |
| --- | --- | --- | --- |
| -   | 0   | | Fiber throughput is normal |
| 0   | 1   | `RELTHRPT_WARN_BIT` | `relthrpt < 1 - sig_cut * IQR`: low, but still usable and still flux-scaled |
| 1   | 2   | `RELTHRPT_BROKEN_BIT` | `relthrpt < rel_val_cut` (0.07): fiber is dead/near-dead |
| 2   | 4   | `RELTHRPT_NOFILE_BIT` | No fluxing file was available; `relthrpt` forced to exactly 1 |
| 3   | 8   | `RELTHRPT_NOTFINITE_BIT` | `relthrpt` is NaN/Inf (fiber all-NaN or all-zero in the flat). Always accompanied by bit 1 |
| 4   | 16  | `RELTHRPT_LOWGOODPIX_BIT` | Fewer than `RELTHRPT_MIN_GOODPIX` pixels survived masking, so the median is noise, not throughput. Sets bits 3, 1 and 0 too |

`RELTHRPT_UNUSABLE_BITS = 2 | 8 | 16`. **A fiber carrying either of those bits is
deliberately NOT flux-scaled by `process_1D`.** Nothing is dropped from the
reduction: the spectrum is written normally, but its flux is left on an arbitrary
scale, so any chi2 computed against it downstream is meaningless. Consumers
should mask on `relthrpt_fiber_unusable(bitmsk_relthrpt)`;
`relthrpt_fiber_fluxable` gives the exact set of fibers that were scaled.

Bit 2 is deliberately excluded from `RELTHRPT_UNUSABLE_BITS`: in that case
`relthrpt` is forced to exactly 1, so the spectrum is simply unfluxed but
unscaled, a known and benign state rather than a broken fiber.

Fiber quality is per fiber-EPOCH, not per fiber: fibers break and fibers get
repaired. A static per-fiber blacklist is the wrong tool; use this per-exposure
flag.

### The throughput median honours the pixel mask

`get_relFlux` takes each fiber's median over pixels carrying no `bad_pix_bits`
(`use_pix_mask`, default **on**). Previously `mask_1d_good` was computed and then
ignored, so a bad-yet-finite-nonzero pixel — a cosmic ray, a saturated pixel —
contributed to the number the entire flux scale is built on. `nanzeromedian` only
ever dropped NaN and exact zeros.

MEASURED blast radius of turning it on, over 2,976 testbed domeflat products
(892,800 fiber measurements, all chips):

| | APO | LCO |
|---|---:|---:|
| median fractional change in `relthrpt` | 5.5e-4 | 1.0e-4 |
| p95 / p99 | 1.3e-2 / 5.7e-2 | 5.8e-4 / 1.9e-3 |
| fibers changing by >1% | 6.2% | 0.28% |
| fiber measurements changing FLAG state | 148 (2.4e-4) | 13 (4.5e-5) |
| fibers newly flagged broken | 12 | 2 |
| worst single flat | 20 fibers change flag, 4 newly broken | 1 / 1 |

So: typically sub-0.1%, with a small tail; the flag state essentially never moves.

`RELTHRPT_MIN_GOODPIX = 256` is the floor below which a fiber is flagged rather
than assigned a meaningless median. It is derived from the data, not picked
round: sub-sampling real domeflat fibers, 256 is the smallest `n` whose
95th-percentile median error (0.064) falls below `rel_val_cut` (0.07), the
sharpest cut this function makes. It cannot fire on healthy data — the minimum
observed good-pixel count across the testbed is 1805 of 2048, seven times the
floor, and zero fibers fell below it.

> **Interaction on the record.** `bad_pix_bits` (24566) does NOT include
> `pix_not_dark_corr_bits` (2^3 = 8), which is what the coherent APO chip-G
> defect block (columns ~512–531, rows ~1362–1377) reads. **Those pixels still
> contribute to throughput after this change.** That is expected: the
> defect-region bit is a separate change and lands outside `bad_pix_bits` first.

**Known limitation.** Only the last chip's (`B`) throughput solution is computed
and it is applied to R, G and B alike, even though `make_relFlux.jl` writes a
per-chip solution and the per-chip symlink is named as though it were the chip's
own. Chromatic throughput differences are therefore unmodelled by construction,
and a fiber dead on one chip only is not flagged unless it is also dead on B. The
chip actually used is now recorded in the product metadata as `relflux_chip`;
`process_1D(...; per_chip_relflux = true)` switches to genuine per-chip fluxing.
That is a survey-wide science change (it moves the flux scale of every R and G
spectrum) and is off by default.

### Per-fiber ingest bits (`bitmsk_ingest`, in `ar1Duni*`)

A second, independent per-fiber flag, computed by `reinterp_spectra` and — until
now — discarded at the end of the function instead of being written out.

| Bit | Value | Meaning |
| --- | --- | --- |
| 1   | 2   | Every pixel of the fiber's 1D flux is NaN or zero |
| 2   | 4   | Every pixel of the fiber is bad by `bad_pix_bits` or missing-chip |
| 3   | 8   | Every pixel of the fiber's ivar is NaN or zero |

Different axis from `bitmsk_relthrpt`: that one says the fiber's dome-flat
throughput is dead, this one says the extracted data itself is unusable. Read
both. Named `bitmsk_ingest`, not `ingestBit`, because arMADGICS has its own
per-spectrum `ingestBit` column with an entirely different bit table.

## Exposure-Level Flags (1D metadata)

The bit tables above are *per pixel* and *per fiber*. Separately, each exposure
is judged as a whole
between the 2D and 1D stages, and the verdict is carried into the `metadata`
group of the 1D data products (`ar1D*`, `ar1Dcal*`, and the reinterpolated
`ar1Duni*` / `ar1Dunical*`, which inherit it from the first chip's 1D file). A
consumer can read a 1D file and see whether the frame was judged bad, and why,
without re-deriving anything from the 2D products or the almanac.

| Field | Type | Meaning |
| ----- | ---- | ------- |
| `exposure_flags` | UInt8 | Bitmask; see the bit table and the reading rule below |
| `exp_class_status` | String | `ok`, `lamp_off_candidate`, `mislabel_candidate`, `faint_twilight`, `persistence_prior`, `unknown`, or `notrun` |
| `exp_class_pred` | String | Predicted content class, e.g. `quartzflat_q1t0u0` |
| `exp_class_labeled` | String | The commanded label it was compared against |
| `exp_class_prob` | Float64 | Max forest probability, `NaN` when there is no verdict |

`exposure_flags` bits (shared namespace — do not renumber):

| Bit | Value | Meaning |
| --- | ----- | ------- |
| 0   | 1     | `predicted_bad` — exposure-type classifier verdict |
| 1   | 2     | `engineering` — configuration's science fibers are an engineering carton |
| 2   | 4     | `notrun` — **no verdict was formed** for this exposure |

### Reading it correctly

`exposure_flags == 0` means **judged, and nothing wrong**. "Never judged" is a
*different value*: bit 2 set. The two are distinguishable from the byte alone —
that is what bit 2 is for.

- `flags == 0` → judged, fine
- `flags & 4` → no verdict; bit 0 is guaranteed clear (the two are mutually
  exclusive, and the writer throws rather than emit a byte asserting both)
- `flags & 1` → judged, and adverse

Read it with `ApogeeReduction.exposure_class_verdict(metadata)`, which returns
`:bad`, `:fine`, or `:unknown`. For "may I use this for science?", use
`exposure_ok_for_science(flags)` — note that bit 2 is deliberately **not** a
no-science bit: an unjudged exposure is not a known-bad one, and silently
dropping everything we failed to look at would turn a monitoring gap into
invisible data loss.

An exposure has no verdict when the check was deliberately disabled, its per-MJD
table is missing, or it errored on that exposure — a crashed check is a failure
to form an opinion, never an adverse one. `exp_class_status` distinguishes those
(`"notrun"` vs `"checkfail"`); they share bit 2 because no consumer would act on
them differently.

**Backward compatibility:** a product written before these fields exist carries
no `exposure_flags` at all. `exposure_class_verdict` reports that as `:unknown`
too, so both routes agree. That is a compatibility shim, not the design.

The almanac carries the same byte at
`exposure_class/<tele>/<mjd>/exposure_flags`, with every row that had no verdict
set to bit 2.

### It is advisory

No exposure is dropped from the reduction because of these flags: engineering
and known-bad frames are still reduced. The one place they exclude anything is
`make_runlist_fiber_flats.jl`, which drops flats with bit 0 set from the
trace/fluxing runlists and logs every exclusion.

### Configuration

The check runs by default. The classifier artifact is a **calibration input**,
configured by path exactly like `caldir_darks` / `caldir_flats` /
`gain_read_cal_dir`: `pipeline.jl --exp_class_model` holds the default and
`airflow/dags/ar_common.py` (`EXP_CLASS_MODEL`) sets it for production, so
swapping models is a config change beside the other calibration inputs. The
version is pinned (v6) — never a glob, never newest-wins.

To disable deliberately: `pipeline.jl --exp_class_model ""`, or
`AR_EXP_CLASS_MODEL="" ./run_all.sh ...`. Leaving `AR_EXP_CLASS_MODEL` unset
means on; setting it to the empty string means off. **A model path that does not
exist is a hard error at startup**, never a silent skip, so a moved or
cleaned-up artifact cannot quietly downgrade a run to "no classification".

Measured cost: **~1.75 s per exposure of worker time** (3 chips), of which
~1.34 s is re-reading the `ar2D` images and only ~0.3 ms is the forest itself,
plus a one-off ~3.6 s model load and ~183 MiB resident per worker process —
about **0.2%** of the reduction's total CPU.

### The engineering bit (a targeting check, not an image check)

Bit 1 is set by `exposure_engineering_from_almanac`, which reads the almanac's
own fiber table — no confSummary dependency. It is **ungated by
`--exp_class_model`**: the check needs only the almanac, so it runs even when the
image classifier is disabled. That independence is why bit 1 can co-occur with
bit 2 (`notrun`): an exposure the classifier never judged can still be known to
be an engineering frame.

It is set when any clause holds, against `ENGINEERING_CARTON_PREFIXES`
(currently only `manual_fps_position_stars`, which by prefix covers `_10`,
`_apogee_10`, and `_lco_apogee_10`):

1. the configuration HAS `category == "science"` fibers and **all** of them
   carry an engineering carton — purity, not majority
   (`ENGINEERING_CARTON_PURITY = 1.0`). Every configuration those cartons appear
   in is 100% that carton, and purity is what keeps the other 27 `manual_*`
   cartons out. A configuration that is *mostly but not purely* an engineering
   carton is NOT flagged and raises a loud warning — that has never happened in
   DR21, so it would be a real signal.
2. the configuration has **zero** science fibers and **any** fiber carries an
   engineering carton (`ENGINEERING_FLAG_SCIENCELESS_POSITION_STARS`). This
   catches the earliest APO FPS positioning configurations, which predate the
   science-category convention. It is keyed on the engineering carton itself,
   not a general "fall back to all fibers", so a science-less configuration
   carrying some other carton is untouched.
3. the configuration was built **without a robostrategy design**
   (`design_id == -999`, `ENGINEERING_FLAG_DESIGNLESS`). `-999` is genuine
   observatory output meaning "no design was generated", not an almanac
   sentinel (almanac's own missing value is `-1`). 17,311 rows carry it but
   17,296 are calibration frames, so this clause is applied **only** to
   `image_type == "object"` — see below.

Clauses 1 and 2 are mutually exclusive by construction (one requires science
fibers, the other requires none), so clause 2 cannot perturb clause 1. Clause 3
is applied per exposure row, outside the per-configuration cache, because
`design_id` lives on the exposure rather than the configuration.

Plate-era exposures have no configuration and no carton, so they are never
flagged engineering, and calibration frames are never flagged (only
`image_type == "object"` is checked — a dark taken while an engineering
configuration was loaded is still a good dark). That image-type restriction is
load-bearing rather than cosmetic: without it, clause 3 alone would discard
17,296 good FPS-era calibration frames.

Alongside the bit, `engineering_frac`, `engineering_carton` and
`engineering_basis` record the evidence for the verdict. `engineering_basis`
records which clause fired: `"science"` (clause 1),
`"scienceless_position_stars"` (clause 2), `"designless"` (clause 3), or
`"none"`.

The bit is computed in `pipeline.jl` right after the 2D stage and before the 1D
stage — the same point as the exposure-type classifier. It writes
`apred/<mjd>/exposureEngineering_<tele>_<mjd>.h5`.
`scripts/cal/decorate_almanac_exptype.jl` computes the same thing for a whole
almanac after the fact, writing `engineering`, `engineering_frac`,
`engineering_carton` and `engineering_basis` alongside `exposure_flags`.

## Testing

To test the pipeline, run the `run_all.sh` script with the desired tele and SJD. For example:

```bash
./src/run_scripts/run_all.sh apo 60639
```

This is good practice before asserting a PR with substantial changes is ready for a merge (in the absence of a CI pipeline, which is still in progress).


## Nomenclature
### SJD
SJD is an "SDSS Julian day," which is adjusted to roll-over earlier than the usual MJD (modified Julian day) so that the day roll-over does not collide with evening calibrations and preparations (defined in https://ui.adsabs.harvard.edu/abs/2015PASP..127..397W/abstract, updated for LCO see for example https://github.com/sdss/sdsstools/blob/main/src/sdsstools/time.py#L21).

The two APOGEE instruments are at two different observatories: APO (north) and LCO (south)

```
MJD = JD - 2400000.5
SJD = MJD + 0.3 # at APO
SJD = MJD + 0.4 # at LCO
```
- APO is MST/MDT. This means that a new SJD occurs at 10:48 AM MST (UTC-7), instead of 5:00 PM MST (UTC-7).
- LCO is CLT/CLST. This means that a new SJD occurs at 12:48 PM CLT (UTC-4), instead of 7:00 PM CLT (UTC-4).

SJD is only ever used for rough definitions of a "day" (taking only the integer part), used mostly for foldering and grouping nightly calibrations with observations. However, long daytime calibration runs can sometimes be broken up by the SJD switch. In call cases, when precise timing is necessary, we convert from TAI to JD, storing at Float64 precision.

## Contributing

All contributions are welcome! Please feel free to open a PR with any changes you would like to see. We will help you troubleshoot any test failures, so please feel free to open a PR with in progress code, or even code in another coding language that has the functionality you would like to see. If you don't have any code related to your idea, please feel free to open an issue and we will help you get started.

[Quick tips on Julia for Python programmers](https://docs.julialang.org/en/v1/manual/noteworthy-differences/#Noteworthy-differences-from-Python)

To enable the Slack Messaging functionality, you need the OAuth token for the bot to be in your bashrc. Please contact the current repo owner for that token. Please also change the the channel key `ENV["SLACK_CHANNEL"]` in the `src/utils.jl` file to the "dev" version during development to reduce noise on the daily processing channel.
