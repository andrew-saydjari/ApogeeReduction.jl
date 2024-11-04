# ApogeeReduction

[![Build Status](https://github.com/andrew-saydjari/ApogeeReduction.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/andrew-saydjari/ApogeeReduction.jl/actions/workflows/CI.yml?query=branch%3Amain)

To enable the Slack Messaging functionality, you need the OAuth token for the bot to be in your bashrc. Please contact the current repo owner for that token.

## Current Flag Bits

During ingestion, some of the exposure files may have issues that cause the spectrum to come through apMADGICS.jl as a vector of only NaNs. This pipeline bit gives insight into the root cause of why this (tiny fraction of the) data is unable to be processed.

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