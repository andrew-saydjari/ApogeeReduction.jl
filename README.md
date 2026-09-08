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

## Exposure-Level Classifier Fields (1D metadata)

The bits above are *per pixel*. Separately, the exposure-type classifier judges
each exposure as a whole from its 2D images, between the 2D and 1D stages, and
its verdict is carried into the `metadata` group of the 1D data products
(`ar1D*`, `ar1Dcal*`, and the reinterpolated `ar1Duni*` / `ar1Dunical*`, which
inherit it from the first chip's 1D file). A consumer can therefore read a 1D
file and see whether the frame was judged bad, and why, without re-deriving
anything from the 2D products or the almanac.

| Field | Type | Meaning |
| ----- | ---- | ------- |
| `exp_class_predicted_bad` | Int8 | **Tri-state.** `1` = classifier ran and judged this exposure bad; `0` = classifier ran and judged it fine; `-1` = **UNKNOWN**, no verdict exists |
| `exp_class_status` | String | Why: `ok`, `lamp_off_candidate`, `mislabel_candidate`, `faint_twilight`, `persistence_prior`, `unknown`, or `notrun` |
| `exp_class_pred` | String | Predicted content class, e.g. `quartzflat_q1t0u0`, `dark_q0t0u0` |
| `exp_class_labeled` | String | The commanded label it was compared against |
| `exp_class_prob` | Float64 | Max forest probability, `NaN` when there is no verdict |

**`-1` is a real value, not a filler.** The exposure-type check is **on by
default**: `pipeline.jl --exp_class_model` defaults to the pinned v6 artifact
`ApogeeReduction.DEFAULT_EXP_CLASS_MODEL`, and both DAGs inherit that. A product
can still legitimately carry no verdict — the check was deliberately disabled,
its per-MJD table is missing, it errored on that exposure, or the file predates
these fields — and in every such case the fields read `predicted_bad = -1`,
`status = "notrun"`, `pred = "unknown"`, `prob = NaN`. Never treat a missing or
`-1` value as a clean bill of health; only `== 1` means bad.

To turn the check off deliberately, pass an empty model path:
`pipeline.jl --exp_class_model ""`, or `AR_EXP_CLASS_MODEL="" ./run_all.sh ...`.
Leaving `AR_EXP_CLASS_MODEL` **unset** means on; setting it to the empty string
means off. A model path that does not exist is a hard error at startup rather
than a silent skip, so a moved or cleaned-up artifact can never quietly
downgrade a run to "no classification".

The artifact version is pinned deliberately in `src/exposureClassifier.jl`;
retraining means editing that constant, not dropping a newer file beside the old
one. Measured cost of the check on the testbed corpus: **~1.75 s per exposure of
worker time** (3 chips), of which ~1.34 s is re-reading the `ar2D` images and
only ~0.3 ms is the forest itself, plus a one-off ~3.6 s model load and ~183 MiB
resident per worker process — about **0.2%** of the reduction's total CPU.

This verdict is **advisory**. No exposure is dropped from the reduction because
of it: engineering and known-bad frames are still reduced. The one place it
excludes anything is `make_runlist_fiber_flats.jl`, which drops
`predicted_bad == 1` flats from the trace/fluxing runlists, and logs every
exclusion. Other per-exposure advisory flags get their own `exp_*` scalars
rather than bits inside `exp_class_predicted_bad`, so separate producers never
contend for one integer.

The same tri-state lands in the almanac as
`exposure_class/<tele>/<mjd>/predicted_bad`, written by
`scripts/cal/decorate_almanac_exptype.jl`, which is what the runlist builder
reads.

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
