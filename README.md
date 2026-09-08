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

## Exposure-Level Flag Bits

Separately from the per-pixel bits above, whole *exposures* can be flagged.
These bits are written by `scripts/cal/decorate_almanac_exptype.jl` into the
almanac's `exposure_class/<tele>/<mjd>/exposure_flags` dataset (`UInt8`, aligned
row-for-row with `raw/<tele>/<mjd>/exposures`).

| Bit   | Value     | Meaning     |
| ----- | --------- | ----------- |
| -     | 0         | No problems |
| 0     | 1         | `predicted_bad` — the image-content classifier says this exposure should not be used (policy: `exposure_predicted_bad`) |
| 1     | 2         | `engineering` — the configuration's science fibers carry an engineering carton, so the exposure exists to exercise the hardware, not to do science (policy: `exposure_is_engineering`) |

**These bits are advisory metadata, not a reduction filter.** Every exposure,
engineering ones included, is still reduced all the way to 1D; the bits exist so
that consumers who assemble *science* samples (prior builds, catalog
construction) can exclude them. `ApogeeReduction.exposure_ok_for_science(flags)`
is the single predicate for "safe to do science with", and
`ApogeeReduction.EXPFLAG_NO_SCIENCE` is the mask it applies.

The engineering bit is a targeting check, not an image check: it is set when
**all** of the configuration's `category == "science"` fibers have a
`firstcarton` matching one of `ENGINEERING_CARTON_PREFIXES` (currently only
`manual_fps_position_stars`, which by prefix covers `_10`, `_apogee_10`, and
`_lco_apogee_10`). The rule is purity, not majority
(`ENGINEERING_CARTON_PURITY = 1.0`): every configuration those cartons appear in
is 100% that carton, and purity is what keeps the other 27 `manual_*` cartons
out. A configuration that is *mostly but not purely* an engineering carton is
NOT flagged and raises a loud warning — that has never happened in DR21, so it
would be a real signal. Plate-era
exposures have no configuration and no carton, so they are never flagged
engineering, and calibration frames are never flagged (only `image_type ==
"object"` is checked — a dark taken while an engineering configuration was
loaded is still a good dark). Alongside the bit, `engineering_frac`,
`engineering_carton` and `engineering_basis` record the evidence for the verdict.

`ENGINEERING_FALLBACK_ALL_FIBERS` (default **off**) would extend the check to
every carton-bearing fiber for configurations with no `category == "science"`
fibers at all; five early-FPS configurations are engineering by carton but label
none of their fibers `science`. It is off because the rule as specified is
science fibers.

The bit is computed in `pipeline.jl` right after the 2D stage and before the 1D
stage — the same point as the exposure-type classifier, but ungated by
`--exp_class_model`, since the check needs only the almanac. It writes
`apred/<mjd>/exposureEngineering_<tele>_<mjd>.h5`.
`scripts/cal/decorate_almanac_exptype.jl` computes the same thing for a whole
almanac after the fact.


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
