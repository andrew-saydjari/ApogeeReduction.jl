#!/bin/bash
# Census of arMADGICS runtime diagnostics.
#
# WHY THIS EXISTS SEPARATELY: arM reports per-spectrum problems with bare
# `println`, not `@warn`. The warnings census and warnings_triage.sh in
# submit_goldens.sh both key on Julia `@warn`/`@error` records, so they see NONE
# of it. In job 7001233 that meant 1,708 silently dropped spectra were never
# counted anywhere -- 1,390 of them ingestBit=64, which the arM README then
# described as harmless.
#
# WHY IT IS NOT A STAGE INSIDE submit_goldens.sh: that script's census runs at
# the end of the REDUCTION stage, before arM has run at all. An arM census there
# could only ever report zero. Call this AFTER the arM stage.
#
# WHY IT NO LONGER GREPS THE LOG: `ingestBit` and `skyBit` are already
# per-spectrum COLUMNS in every arM batch product, and reading them there is
# strictly better than counting println lines:
#   - the sky-guard verdict is EXPOSURE-level but was printed once per TARGET
#     FIBER, so line counts overstated it ~130x (job 7001233: 479,570 lines,
#     3,676 unique exposures). arMADGICS now suppresses that print entirely, so
#     grepping for it would report ZERO from here on.
#   - `ingestBit` was visible only via "Skipping spectrum", i.e. only for FATAL
#     bits; the informational bits were invisible. The products carry all of it.
#   - a rotated, truncated or redirected log undercounts and looks exactly like a
#     clean run; an append-mode resume double-counts (that specific bug is fixed
#     in submit_goldens.sh, but the class of bug is inherent to logs).
#   - any reword of a println silently zeroes a grep-based census.
# The products are the record. The heavy lifting is in arm_census.jl, because
# 16k per-file HDF5 reads are not something bash can do.
#
# Usage:  arm_census.sh <arm-raw-dir> [arm-log ...]
#   <arm-raw-dir>  arM output directory containing batch_info.txt and the
#                  per-fiber NNN/ subdirectories of batch .h5 products.
#   [arm-log ...]  optional; only for diagnostics that have no product column
#                  yet (currently the prior-support guard).
# Env:    JULIA (default "julia"), ARM_CENSUS_THREADS (default 8)
# Exit:   0 always -- this is a reporting tool, not a gate.
set -uo pipefail

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
base_dir=$(cd "${here}/../.." && pwd)
JULIA=${JULIA:-julia}

if [ $# -lt 1 ]; then
    echo "arM census: no arM output directory given - arM diagnostics NOT counted (not zero, UNCOUNTED)"
    echo "usage: $(basename "$0") <arm-raw-dir> [arm-log ...]" >&2
    exit 0
fi

raw=$1
shift

# ---- flag columns, straight out of the batch products -----------------------
if [ ! -d "$raw" ]; then
    echo "arM census: $raw is not a directory - ingestBit/skyBit NOT counted (not zero, UNCOUNTED)"
elif ! command -v "$JULIA" >/dev/null 2>&1; then
    echo "arM census: no julia on PATH (set \$JULIA) - ingestBit/skyBit NOT counted (not zero, UNCOUNTED)"
else
    echo "arM census: reading flag columns from $raw"
    # arm_census.jl prints its own "arM census: ..." lines and exits non-zero only
    # when it could not count at all -- in which case it has already said NOT
    # COUNTED, so all that is left here is to not claim success.
    if ! "$JULIA" --project="$base_dir" -t "${ARM_CENSUS_THREADS:-8}" \
         "${here}/arm_census.jl" "$raw"; then
        echo "arM census: arm_census.jl did not complete - ingestBit/skyBit NOT counted (not zero, UNCOUNTED)"
    fi
fi

# ---- log-only diagnostics ---------------------------------------------------
# Everything below has no per-spectrum product column yet. Move each one to the
# products as its column appears; do not add new grep-based counters.
if [ $# -eq 0 ]; then
    echo "arM census: no log given - log-only diagnostics NOT counted (not zero, UNCOUNTED)"
fi
for log in "$@"; do
    if [ ! -r "$log" ]; then
        echo "arM census: $log not readable - log-only diagnostics NOT counted"
        continue
    fi
    echo "arM census: reading $log for log-only diagnostics"

    # prior-support guard (LCO-only and always 67 px in job 7001233).
    p=$(grep -c 'load_fiber_priors: prior-support guard' "$log" || true)
    echo "arM census: ${p}x prior-support guard"

    # The sky guard is deliberately NOT counted from the log any more: arMADGICS
    # no longer prints it, so a count here would be a silent zero rather than a
    # measurement. It comes from skyBit above. If the legacy lines DO appear, the
    # log predates that change -- say so rather than quietly mixing the two.
    stale=$(grep -c 'getSky4visit: sky guard flagged' "$log" || true)
    if [ "$stale" -gt 0 ]; then
        echo "arM census: NOTE ${stale}x legacy 'sky guard flagged' lines in this log - it"
        echo "arM census:   predates the arMADGICS change that suppressed them. The skyBit"
        echo "arM census:   tally above is the authoritative count; this line count is"
        echo "arM census:   inflated ~130x by per-target-fiber repetition."
    fi
done
exit 0
