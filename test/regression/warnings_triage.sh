#!/bin/bash
# warnings_triage.sh — normalised warning/error inventory for a pipeline run,
# plus a diff against an adjudicated reference set.
#
# "Reference set", deliberately, not "baseline": the counts from a previous run
# are what happened, not what is acceptable. A category earns a benign verdict
# from a demonstrated mechanism — the code path, the data condition, and why
# that condition is legitimately present — never from having been there before.
#
# Motivation
# ----------
# submit_goldens.sh has a "warnings census" stage that greps a fixed list of
# hand-named substrings.  That only counts categories somebody thought to name,
# and it silently misses (a) categories nobody named and (b) categories whose
# grep pattern has drifted from the actual message text.  Both happened in job
# 6980442: 88 "lamp turned off" warnings were never censused, and the
# "no useful relfluxing files" pattern reported 0x against a real count of 3
# because the message reads "any useful relfluxing files".
#
# This script instead enumerates *every* Julia @warn/@error record in the log
# and groups them by EMIT SITE (module + file:line), which is exact, is immune
# to message-text edits, and is immune to ProgressMeter output being glued onto
# the message by carriage-return re-renders.
#
# Usage
# -----
#   warnings_triage.sh PATH [PATH...]
#       Print the inventory (TSV) for the given log files and/or directories.
#       Directories are scanned non-recursively for *.log and *.out.
#
#   warnings_triage.sh -r REFERENCE.tsv PATH [PATH...]
#       Print the inventory, then a diff against REFERENCE.tsv.
#       Exit status: 0 = no new sites, 1 = new emit site(s) appeared.
#       (-b/--baseline are accepted as aliases for -r.)
#
#   warnings_triage.sh --emit-reference PATH [PATH...] > reference.tsv
#       Print a reference skeleton; verdict/expect/note are left for a human,
#       and the honest starting verdict for an un-investigated site is UNKNOWN.
#
#   warnings_triage.sh --units PATH [PATH...]
#       Also print, per emit site, the distribution over (tele, mjd) units
#       parsed out of the message.  Concentration vs. spread is usually the
#       whole question: 260 warnings from one exposure is a data defect,
#       260 warnings spread over 200 nights is a threshold.
#
# Reference TSV columns (tab separated, '#' lines are comments):
#   site <TAB> count <TAB> verdict <TAB> expect <TAB> note
# where `verdict` is EXPECTED | ACTIONABLE | UNKNOWN and `expect` is one of
# stable | up | down | either, describing how the count should move in the next
# run given known upstream changes. UNKNOWN is the correct verdict for anything
# whose mechanism has not actually been established.
#
# Design note: this is bash+awk rather than Julia on purpose — it runs as the
# last stage inside the sbatch harness next to the existing census, must work
# with no Julia depot warm-up, and processes a multi-GB text stream.

set -u -o pipefail

emit_baseline=false
show_units=false
baseline=""
paths=()

while [ $# -gt 0 ]; do
    case "$1" in
        -r|--reference|-b|--baseline) baseline=$2; shift 2 ;;
        --emit-reference|--emit-baseline) emit_baseline=true; shift ;;
        --units)         show_units=true; shift ;;
        -h|--help)       sed -n '2,50p' "$0"; exit 0 ;;
        --)              shift; break ;;
        -*)              echo "unknown option: $1" >&2; exit 2 ;;
        *)               paths+=("$1"); shift ;;
    esac
done
paths+=("$@")

if [ ${#paths[@]} -eq 0 ]; then
    echo "usage: $(basename "$0") [-r REFERENCE.tsv] [--units] [--emit-reference] PATH [PATH...]" >&2
    exit 2
fi

# ---- collect input files -----------------------------------------------------
files=()
for p in "${paths[@]}"; do
    if [ -d "$p" ]; then
        for f in "$p"/*.log "$p"/*.out; do
            [ -f "$f" ] && files+=("$f")
        done
    elif [ -f "$p" ]; then
        files+=("$p")
    else
        echo "warning: no such file or directory: $p" >&2
    fi
done
if [ ${#files[@]} -eq 0 ]; then
    echo "error: no log files found" >&2
    exit 2
fi

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

# ---- flatten: split CR re-renders onto their own lines, strip ANSI ------------
# ProgressMeter writes with \r and ESC[K; without this the @warn records are
# buried inside multi-megabyte single "lines".
for f in "${files[@]}"; do
    tr '\r' '\n' < "$f"
done | sed -e 's/\x1b\[[0-9;?]*[a-zA-Z]//g' > "$tmp/flat.txt"

# ---- extract records ---------------------------------------------------------
# A Julia log record is:
#     ┌ Warning: <message>
#     │   <key> = <value>          (optional, repeated)
#     └ @ <Module> <path>:<line>
# We pair each ┌ with the next └ @.  Both can have progress-bar text glued to
# them, so we truncate defensively.
awk '
function clean(s) {
    sub(/\[K.*$/, "", s)                 # bare ESC[K residue and everything after
    sub(/ *[0-9]+%\|.*$/, "", s)         # a ProgressMeter re-render glued on
    sub(/[ \t]+$/, "", s)
    return s
}
function norm(s,   t) {
    t = s
    # order matters: most specific first
    gsub(/ar[0-9][A-Za-z]*_(apo|lco)_[0-9][0-9][0-9][0-9][0-9]_[0-9]+_[RGB]_[a-z]+\.h5/, "<PRODFILE>", t)
    gsub(/[^ ]*\.h5/,                       "<HDFFILE>", t)
    gsub(/\/[^ ]+/,                         "<PATH>", t)
    gsub(/(apo|lco) [0-9][0-9][0-9][0-9][0-9]/, "<TELE> <MJD>", t)
    gsub(/(apo|lco)/,                       "<TELE>", t)
    gsub(/ar[0-9][A-Za-z]*/,                "<FNAMETYPE>", t)
    gsub(/-?[0-9]+\.[0-9]+/,                "<F>", t)
    gsub(/[0-9]+/,                          "<N>", t)
    gsub(/  +/,                             " ", t)
    return t
}
# distinct (tele, mjd) units mentioned in the message
function units(s,   out, u, t) {
    out = ""
    t = s
    while (match(t, /(apo|lco)_[0-9][0-9][0-9][0-9][0-9]/)) {
        out = out " " substr(t, RSTART, RLENGTH)
        t = substr(t, RSTART + RLENGTH)
    }
    if (out != "") return out
    t = s
    while (match(t, /(apo|lco) (mjd )?[0-9][0-9][0-9][0-9][0-9]/)) {
        u = substr(t, RSTART, RLENGTH); sub(/ /, "_", u); sub(/mjd_/, "", u)
        out = out " " u
        t = substr(t, RSTART + RLENGTH)
    }
    return out
}
# "<Module> <path>:<line>" -> "<Module> src/foo.jl:123" (or "[depot] ...")
function shorten(site,   mod, path, i, tag) {
    i = index(site, " ")
    if (i == 0) return site
    mod  = substr(site, 1, i - 1)
    path = substr(site, i + 1)
    tag = (path ~ /\.julia\/packages\//) ? "[depot] " : ""
    if (match(path, /\/(src|scripts|test|ext)\//)) path = substr(path, RSTART + 1)
    return mod " " tag path
}
/┌ (Warning|Error):/ {
    match($0, /┌ (Warning|Error):/)
    rec = substr($0, RSTART)
    sev = (rec ~ /^┌ Error:/) ? "ERROR" : "WARN"
    sub(/^┌ (Warning|Error): */, "", rec)
    msg = clean(rec)
    pending = 1
    next
}
pending && /└ @ / {
    match($0, /└ @ /)
    site = substr($0, RSTART + RLENGTH)
    # "<Module> <path>:<line>" — cut anything glued after the line number
    if (match(site, /^[^ ]+ [^ ]*:[0-9]+/)) site = substr(site, RSTART, RLENGTH)
    site = shorten(site)

    n = norm(msg)
    key = sev "\t" site
    cnt[key]++
    # representative message: keep the SHORTEST normalisation seen. Records with
    # progress-bar text glued on are strictly longer, so this picks the clean one.
    if (!(key in rep) || length(n) < length(rep[key])) rep[key] = n
    u = units(msg)
    if (u != "") {
        split(u, arr, " ")
        for (j in arr) if (arr[j] != "") ukey[key SUBSEP arr[j]]++
    }
    pending = 0
    next
}
END {
    for (k in cnt) {
        nu = 0
        for (x in ukey) { split(x, p, SUBSEP); if (p[1] == k) nu++ }
        printf "%s\t%d\t%d\t%s\n", k, cnt[k], nu, rep[k]
    }
}
' "$tmp/flat.txt" | sort -t$'\t' -k3,3nr -k2,2 > "$tmp/inv.tsv"

# ---- per-site unit distribution ---------------------------------------------
if $show_units; then
    awk '
    function clean(s) { sub(/\[K.*$/, "", s); sub(/ *[0-9]+%\|.*$/, "", s); return s }
    function shorten(site,   mod, path, i, tag) {
        i = index(site, " "); if (i == 0) return site
        mod = substr(site, 1, i-1); path = substr(site, i+1)
        tag = (path ~ /\.julia\/packages\//) ? "[depot] " : ""
        if (match(path, /\/(src|scripts|test|ext)\//)) path = substr(path, RSTART+1)
        return mod " " tag path
    }
    /┌ (Warning|Error):/ { match($0, /┌ (Warning|Error):/); msg=clean(substr($0,RSTART)); pending=1; next }
    pending && /└ @ / {
        match($0, /└ @ /); site = substr($0, RSTART+RLENGTH)
        if (match(site, /^[^ ]+ [^ ]*:[0-9]+/)) site = substr(site, RSTART, RLENGTH)
        site = shorten(site)
        s = msg
        found = 0
        while (match(s, /(apo|lco)_[0-9][0-9][0-9][0-9][0-9]/)) {
            u = substr(s, RSTART, RLENGTH); gsub(/_/, " ", u)
            c[site SUBSEP u]++; found = 1
            s = substr(s, RSTART + RLENGTH)
        }
        if (!found) {
            s = msg
            if (match(s, /(tele|for) (apo|lco) (mjd )?[0-9][0-9][0-9][0-9][0-9]/)) {
                u = substr(s, RSTART, RLENGTH); sub(/^(tele|for) /, "", u); sub(/mjd /, "", u)
                c[site SUBSEP u]++
            } else {
                c[site SUBSEP "(no tele/mjd in message)"]++
            }
        }
        pending = 0; next
    }
    END { for (x in c) { split(x, p, SUBSEP); printf "%s\t%s\t%d\n", p[1], p[2], c[x] } }
    ' "$tmp/flat.txt" | sort -t$'\t' -k1,1 -k3,3nr > "$tmp/units.tsv"
fi

# ---- output ------------------------------------------------------------------
total=$(awk -F'\t' '{s+=$3} END{print s+0}' "$tmp/inv.tsv")
nsite=$(wc -l < "$tmp/inv.tsv")

if $emit_baseline; then
    echo "# warnings reference set — generated by test/regression/warnings_triage.sh"
    echo "# logs: ${files[*]}"
    echo "# generated: $(date -Is)"
    echo "# total records: ${total} across ${nsite} emit sites"
    echo "#"
    echo "# columns: site<TAB>count<TAB>verdict<TAB>expect<TAB>note"
    echo "#   verdict: EXPECTED | ACTIONABLE | UNKNOWN   (fill in by hand;"
    echo "#            UNKNOWN until a mechanism is actually demonstrated)"
    echo "#   expect:  stable | up | down | either       (how the count should move next run)"
    printf '#site\tcount\tverdict\texpect\tnote\n'
    awk -F'\t' '{printf "%s\t%d\tUNKNOWN\teither\t%s\n", $2, $3, $5}' "$tmp/inv.tsv"
    exit 0
fi

echo "=== warnings inventory ==="
echo "logs        : ${#files[@]} file(s)"
echo "records     : ${total}"
echo "emit sites  : ${nsite}"
echo
printf '%-8s %-58s %7s %6s  %s\n' SEVERITY SITE COUNT UNITS MESSAGE
awk -F'\t' '{printf "%-8s %-58s %7d %6d  %s\n", $1, $2, $3, $4, substr($5,1,110)}' "$tmp/inv.tsv"

if $show_units; then
    echo
    echo "=== per-site (tele, mjd) distribution ==="
    awk -F'\t' '
      $1 != last { print ""; print $1; last = $1 }
      { printf "    %-28s %6d\n", $2, $3 }
    ' "$tmp/units.tsv"
fi

rc=0
if [ -n "$baseline" ]; then
    if [ ! -f "$baseline" ]; then
        echo "error: reference set not found: $baseline" >&2
        exit 2
    fi
    echo
    echo "=== diff vs reference set: $baseline ==="
    awk -F'\t' -v OFS='\t' '
    FILENAME == ARGV[1] {
        if ($0 ~ /^#/ || NF < 2) next
        b[$1] = $2; verdict[$1] = $3; expect[$1] = $4; note[$1] = $5
        next
    }
    {
        cur[$2] = $3; msg[$2] = $5
    }
    END {
        newn = 0; gonen = 0; chg = 0; same = 0
        for (s in cur) {
            if (!(s in b)) {
                printf "NEW         %-58s %6d          %s\n", s, cur[s], substr(msg[s],1,90)
                newn++
            }
        }
        for (s in b) {
            if (!(s in cur)) {
                printf "RESOLVED    %-58s %6d -> 0    [%s] %s\n", s, b[s], verdict[s], note[s]
                gonen++
            } else if (cur[s] != b[s]) {
                d = cur[s] - b[s]
                dir = (d > 0) ? "up" : "down"
                ok = (expect[s] == "either" || expect[s] == dir) ? "as-expected" : \
                     ((expect[s] == "stable") ? "UNEXPECTED (reference says stable)" : \
                      "UNEXPECTED (reference says " expect[s] ")")
                printf "CHANGED     %-58s %6d -> %-6d (%+d, %s) [%s] %s\n", s, b[s], cur[s], d, ok, verdict[s], note[s]
                chg++
            } else {
                same++
            }
        }
        printf "\nsummary: %d new, %d resolved, %d changed, %d unchanged\n", newn, gonen, chg, same
        if (newn > 0) {
            print "NEW emit sites have never been adjudicated — triage them before accepting the run."
        }
        for (s in cur) if (s in b && verdict[s] == "UNKNOWN")
            printf "OPEN        %-58s %6d          verdict is still UNKNOWN — mechanism not established\n", s, cur[s]
        exit (newn > 0) ? 1 : 0
    }
    ' "$baseline" "$tmp/inv.tsv"
    rc=$?
fi

exit $rc
