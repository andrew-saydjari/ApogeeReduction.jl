# ApogeeReduction

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
- run_*.sh : wrapper scripts to run the pipeline (submit job, determine resources, etc.)
- make_runlist_*.sh: interface to get the data to run the pipeline on
- pipeline_*.sh: how functions combine to process the data
- src/*.jl: core functions of the repository

### File Structure
```
├── src/ : core functions of the repository
│   └── run_scripts/ : scripts general users will interact with to run the pipeline
│   └── cal_build/ : scripts to build the calibrations files
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

| Value         | Bit         | Meaning     |
| ----------- | ----------- | ----------- |
| 0     | -     | No problems       |
| 1     | 0     | reference array pixels |
| 2     | 1     | reference pixels |
| 4     | 2     | bad reference pixels |
| 8     | 3     | pixels not dark corrected |
| 16    | 4     | pixels with negative dark current |
| 32    | 5     | pixels with large dark current |
| 64    | 6     | flat response too low |
| 128   | 7     | reads dropped for CR rejection = 1 |
| 256   | 8     | reads dropped for CR rejection > 1 |
| 512   | 9     | bad linear SUTR chi2 |
| 1024  | 10    | failed 1D extraction |
| 2048  | 11    | no nearby good pixels in 1D extraction |
| 4096  | 12    | neff>10 in 1D extraction |
| 8192  | 13    | pixel likely saturated |


## Testing

To test the pipeline, run the `run_all.sh` script with the desired tele and MJD. For example:

```bash
./src/run_scripts/run_all.sh apo 60639
```

This is good practice before asserting a PR with substantial changes is ready for a merge (in the absence of a CI pipeline, which is still in progress).

## Contributing

All contributions are welcome! Please feel free to open a PR with any changes you would like to see. We will help you troubleshoot any test failures, so please feel free to open a PR with in progress code, or even code in another coding language that has the functionality you would like to see. If you don't have any code related to your idea, please feel free to open an issue and we will help you get started.

[Quick tips on Julia for Python programmers](https://docs.julialang.org/en/v1/manual/noteworthy-differences/#Noteworthy-differences-from-Python)

To enable the Slack Messaging functionality, you need the OAuth token for the bot to be in your bashrc. Please contact the current repo owner for that token. Please also change the the channel key `ENV["SLACK_CHANNEL"]` in the `src/utils.jl` file to the "dev" version during development to reduce noise on the daily processing channel.
