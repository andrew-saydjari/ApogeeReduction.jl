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
# could only ever report zero. Call this AFTER the arM stage, on the log arM
# actually wrote to (usually the job's stdout).
#
# Usage:  arm_census.sh <arm-log> [more-logs...]
# Exit:   0 always -- this is a reporting tool, not a gate.
set -uo pipefail

if [ $# -lt 1 ]; then
    echo "arM census: no log given - arM diagnostics NOT counted (not zero, UNCOUNTED)"
    echo "usage: $(basename "$0") <arm-log> [more-logs...]" >&2
    exit 0
fi

for log in "$@"; do
    if [ ! -r "$log" ]; then
        echo "arM census: $log not readable - arM diagnostics NOT counted"
        continue
    fi
    echo "arM census: reading $log"

    # ---- ingestBit: a fatal bit means the spectrum was never fitted.
    # Bit meanings are in the arMADGICS README, "Ingest Module Flag Bits".
    # ingestBit != 0 is equivalent to RV_flag == 64 (verified on 1.6M spectra).
    n=$(grep -c 'Skipping spectrum (ingestBit=' "$log" || true)
    echo "arM census: ${n}x spectra skipped (ingestBit != 0)"
    grep -oE 'Skipping spectrum \(ingestBit=[0-9]+\)' "$log" 2>/dev/null \
        | grep -oE '=[0-9]+' | tr -d '=' | sort -n | uniq -c \
        | while read -r c bit; do echo "arM census:   ${c}x ingestBit=${bit}"; done

    # ---- sky guard. The verdict is EXPOSURE-level but printed once per target
    # fiber, so raw line counts overstate it by ~130x (job 7001233: 479,570
    # lines, 3,676 unique exposures). Report both; the unique count is the one
    # with physical meaning.
    sl=$(grep -c 'getSky4visit: sky guard flagged' "$log" || true)
    su=$(grep -oE 'getSky4visit: sky guard flagged tele=[a-z]+, mjd=[0-9]+, expnum=[0-9]+' \
         "$log" 2>/dev/null | sort -u | wc -l)
    echo "arM census: ${sl}x sky-guard lines (${su} unique exposures)"
    grep -oE 'skyBit=[0-9]+' "$log" 2>/dev/null | grep -oE '[0-9]+$' \
        | sort -n | uniq -c \
        | while read -r c bit; do echo "arM census:   ${c}x skyBit=${bit} (lines, not unique)"; done

    # ---- prior-support guard (LCO-only and always 67 px in job 7001233).
    p=$(grep -c 'load_fiber_priors: prior-support guard' "$log" || true)
    echo "arM census: ${p}x prior-support guard"
done
exit 0
