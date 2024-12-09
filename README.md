# ApogeeReduction

[![Build Status](https://github.com/andrew-saydjari/ApogeeReduction.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/andrew-saydjari/ApogeeReduction.jl/actions/workflows/CI.yml?query=branch%3Amain)

To enable the Slack Messaging functionality, you need the OAuth token for the bot to be in your bashrc. Please contact the current repo owner for that token.

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

## Testing

To test the pipeline, run the `run_all.sh` script with the desired tele and MJD. For example:

```bash
./src/run_scripts/run_all.sh apo 60639
```

This is good practice before asserting a PR with substantial changes is ready for a merge (in the absence of a CI pipeline, which is still in progress).

# Contributing

All contributions are welcome! Please feel free to open a PR with any changes you would like to see. We will help you troubleshoot any test failures, so please feel free to open a PR with in progress code, or even code in another coding language that has the functionality you would like to see. If you don't have any code related to your idea, please feel free to open an issue and we will help you get started.

[Quick tips on Julia for Python programmers](https://docs.julialang.org/en/v1/manual/noteworthy-differences/#Noteworthy-differences-from-Python)
