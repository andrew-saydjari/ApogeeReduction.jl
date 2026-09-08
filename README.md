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

## Exposure-Level Flags (1D metadata)

The bits above are *per pixel*. Separately, each exposure is judged as a whole
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
