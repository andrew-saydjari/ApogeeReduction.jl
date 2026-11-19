# classifier_crosscheck.jl — join the post-hoc exposure-type classifier sweep
# against the run's warnings and against the almanac's own bookkeeping.
#
# The question this exists to answer
# ---------------------------------
# The warnings adjudication (WARNINGS.md) established what the pipeline
# complained about.  It could not establish whether those complaints were RIGHT,
# because the only other per-exposure opinion available — the exposure-type
# classifier — was switched off for the run.  This script supplies that second
# opinion after the fact and lays the two side by side.
#
# Four quadrants, and what each one means:
#
#   CONFIRMED   warning + classifier flags it
#               Two independent signals agree.  The exposure is bad.
#
#   WARN-ONLY   warning + classifier says "ok"
#               Candidate FALSE-POSITIVE WARNING.  Not proof: the classifier is
#               a 67-feature random forest over whole-chip summary statistics
#               and is blind to per-fiber and per-chip effects, so a warning
#               about one chip or one fiber can be entirely real while the
#               exposure classifies fine.  Read the site before believing it.
#
#   CLF-ONLY    no warning + classifier flags it
#               Candidate MISSED DETECTION.  The most interesting quadrant:
#               nothing in the run's own logs objected to these exposures.
#
#   quiet       no warning + classifier ok.
#
# Deliberately NOT computed: an accuracy, a precision, or an F-score.  There is
# no ground truth here — both signals are estimates, and scoring one against the
# other would smuggle in the assumption that the classifier is right.  The
# output is a contingency table plus the exposure lists needed to go look.
#
# Usage
# -----
#   julia test/regression/classifier_crosscheck.jl \
#       --predictions exposure_class_predictions.tsv   (from classify_run.jl)
#       --warnings    warnings_by_exposure.tsv         (warnings_triage.sh --exposures)
#       [--focus apo:60255:7,lco:60078:44]             (print full detail rows)
#       [--list-limit 40]
#
# Reads two TSVs, writes a report to stdout.  Touches nothing else.

using Printf

function getopt(args, name, default)
    i = findfirst(==(name), args)
    i === nothing ? default : args[i + 1]
end
args = ARGS
predpath = getopt(args, "--predictions", "")
warnpath = getopt(args, "--warnings", "")
focusspec = getopt(args, "--focus", "")
listlimit = parse(Int, getopt(args, "--list-limit", "40"))
if predpath == "" || warnpath == ""
    println("usage: classifier_crosscheck.jl --predictions P.tsv --warnings W.tsv " *
            "[--focus tele:mjd:exp,...] [--list-limit N]")
    exit(2)
end

function readtsv(path)
    lines = readlines(path)
    isempty(lines) && error("empty file: $path")
    hdr = split(lines[1], '\t')
    rows = Vector{Dict{String, String}}()
    for ln in lines[2:end]
        isempty(strip(ln)) && continue
        p = split(ln, '\t')
        length(p) == length(hdr) || continue
        push!(rows, Dict(zip(hdr, p)))
    end
    rows
end

pred = readtsv(predpath)
warn = readtsv(warnpath)
key(r) = (r["tele"], parse(Int, r["mjd"]), parse(Int, r["exposure"]))

# ---- warnings, resolved and unresolved --------------------------------------
wsites = Dict{Tuple{String, Int, Int}, Dict{String, Int}}()
unresolved = Dict{String, Int}()   # site => record count with no exposure
nightonly = Dict{String, Int}()
for r in warn
    e = parse(Int, r["exposure"])
    n = parse(Int, r["count"])
    if e > 0
        get!(wsites, key(r), Dict{String, Int}())[r["site"]] = n
    elseif e == 0
        nightonly[r["site"]] = get(nightonly, r["site"], 0) + n
    else
        unresolved[r["site"]] = get(unresolved, r["site"], 0) + n
    end
end

# `status` values that mean "the classifier objects to this exposure".
# persistence_prior is informational in pipeline.jl (recorded, never warned), so
# it is NOT treated as an objection here either — keeping the two consistent.
const FLAGGED = Set(["unknown", "mislabel_candidate", "lamp_off_candidate",
    "faint_twilight"])
isflagged(r) = r["status"] in FLAGGED

println("=" ^ 100)
println("exposure-type classifier <-> warnings cross-check")
println("=" ^ 100)
println("predictions : ", predpath, "  (", length(pred), " exposures)")
println("warnings    : ", warnpath)
println()

# ---- classifier verdict distribution ----------------------------------------
statcount = Dict{String, Int}()
for r in pred
    statcount[r["status"]] = get(statcount, r["status"], 0) + 1
end
println("--- classifier verdicts over all classified exposures ---")
for (k, v) in sort(collect(statcount), by = x -> -x[2])
    @printf("  %-22s %7d  (%5.2f%%)\n", k, v, 100v / length(pred))
end
nflag = count(isflagged, pred)
@printf("  %-22s %7d  (%5.2f%%)\n", "TOTAL FLAGGED", nflag, 100nflag / length(pred))
println()

# ---- the contingency table ---------------------------------------------------
predkeys = Set(key(r) for r in pred)
confirmed = [r for r in pred if isflagged(r) && haskey(wsites, key(r))]
warnonly = [r for r in pred if !isflagged(r) && haskey(wsites, key(r))]
clfonly = [r for r in pred if isflagged(r) && !haskey(wsites, key(r))]
quiet = length(pred) - length(confirmed) - length(warnonly) - length(clfonly)
# warnings whose exposure never got classified (no ar2D triple, or no almanac row)
orphan = [k for k in keys(wsites) if !(k in predkeys)]

println("--- contingency table (exposures) ---")
@printf("%-28s %12s %12s\n", "", "clf FLAGS", "clf ok")
@printf("%-28s %12d %12d\n", "has warning(s)", length(confirmed), length(warnonly))
@printf("%-28s %12d %12d\n", "no warning", length(clfonly), quiet)
println()
println("  CONFIRMED : ", length(confirmed), " — warning and classifier agree")
println("  WARN-ONLY : ", length(warnonly),
    " — candidate false-positive warnings (see caveat in the header)")
println("  CLF-ONLY  : ", length(clfonly), " — candidate MISSED detections")
if !isempty(orphan)
    println("  orphan    : ", length(orphan),
        " warned exposure(s) absent from the prediction table (no ar2D triple / no almanac row)")
    for k in sort(orphan)[1:min(10, length(orphan))]
        println("              ", k)
    end
end
println()

if !isempty(unresolved) || !isempty(nightonly)
    println("--- warnings that CANNOT be cross-checked per exposure ---")
    println("(these are not failures of the join; the messages do not name an exposure)")
    for (s, n) in sort(collect(nightonly), by = x -> -x[2])
        @printf("  %-56s %6d records  (night only)\n", s, n)
    end
    for (s, n) in sort(collect(unresolved), by = x -> -x[2])
        @printf("  %-56s %6d records  (no tele/mjd/exposure in message)\n", s, n)
    end
    println()
end

# ---- per-emit-site agreement -------------------------------------------------
println("--- per emit site: does the classifier back the warning? ---")
# [n_exposures, n_flagged, n_ok, n_unclassified] — an exposure the sweep never
# saw is NOT evidence that the classifier is content with it, so it gets its own
# column rather than being folded into "clf ok".
sitestats = Dict{String, Vector{Int}}()
predindex = Dict(key(r) => r for r in pred)
for (k, sites) in wsites
    r = get(predindex, k, nothing)
    for s in keys(sites)
        v = get!(sitestats, s, [0, 0, 0, 0])
        v[1] += 1
        if r === nothing
            v[4] += 1
        elseif isflagged(r)
            v[2] += 1
        else
            v[3] += 1
        end
    end
end
@printf("%-56s %8s %8s %8s %10s\n", "SITE", "exps", "clfflag", "clf ok", "unclassif")
for (s, v) in sort(collect(sitestats), by = x -> -x[2][1])
    @printf("%-56s %8d %8d %8d %10d\n", s, v[1], v[2], v[3], v[4])
end
println()

# ---- detail listings ---------------------------------------------------------
sitesof(k) = join(sort(collect(keys(get(wsites, k, Dict{String, Int}())))), ", ")
function showrows(title, rows; withsites = true, limit = listlimit)
    println("--- ", title, " (", length(rows), ") ---")
    isempty(rows) && (println("  (none)"); println(); return)
    @printf("%-5s %6s %5s %-22s %-22s %6s %-20s %4s %4s\n",
        "tele", "mjd", "exp", "labeled", "pred", "p", "status", "fbad", "pbad")
    for r in rows[1:min(limit, length(rows))]
        @printf("%-5s %6s %5s %-22s %-22s %6s %-20s %4s %4s",
            r["tele"], r["mjd"], r["exposure"], r["labeled"], r["pred"],
            r["prob"], r["status"], r["flagged_bad"], r["predicted_bad"])
        withsites && print("  ", sitesof(key(r)))
        println()
    end
    length(rows) > limit && println("  ... ", length(rows) - limit, " more")
    println()
end

srt(rows) = sort(rows, by = r -> (r["tele"], parse(Int, r["mjd"]),
    parse(Int, r["exposure"])))
showrows("CONFIRMED — warning + classifier flags", srt(confirmed))
showrows("WARN-ONLY — warning, classifier says ok (candidate false-positive warnings)",
    srt(warnonly))

println("--- CLF-ONLY — classifier flags, NO warning (candidate missed detections) ---")
if isempty(clfonly)
    println("  (none)")
else
    grp = Dict{Tuple{String, String, String}, Int}()
    for r in clfonly
        k = (r["status"], r["labeled"], r["pred"])
        grp[k] = get(grp, k, 0) + 1
    end
    @printf("%-22s %-22s %-22s %8s\n", "status", "labeled", "pred", "n")
    for (k, v) in sort(collect(grp), by = x -> -x[2])
        @printf("%-22s %-22s %-22s %8d\n", k[1], k[2], k[3], v)
    end
    println()
    showrows("CLF-ONLY, highest-confidence first",
        sort(clfonly,
            by = r -> -(tryparse(Float64, r["prob"]) === nothing ? 0.0 :
                        parse(Float64, r["prob"]))); withsites = false)
end

# ---- almanac disagreements ---------------------------------------------------
println("--- disagreement with the almanac's own bookkeeping ---")
pb_fb = Dict{Tuple{String, String}, Int}()
for r in pred
    pb_fb[(r["predicted_bad"], r["flagged_bad"])] = get(pb_fb,
        (r["predicted_bad"], r["flagged_bad"]), 0) + 1
end
@printf("%-16s %-16s %8s\n", "predicted_bad", "flagged_bad", "n")
for (k, v) in sort(collect(pb_fb), by = x -> -x[2])
    @printf("%-16s %-16s %8d\n", k[1], k[2], v)
end
nmask = count(r -> r["predicted_bad"] == "1" && r["flagged_bad"] == "0", pred)
println()
println("  ", nmask,
    " exposure(s) the classifier would mask that the almanac still calls good ",
    "(flagged_bad = 0).")
println("  That divergence is the actionable output: either the classifier is ",
    "wrong or the almanac is stale.")
println()

println("--- declared type vs classified content (mislabel_candidate only) ---")
mis = [r for r in pred if r["status"] == "mislabel_candidate"]
if isempty(mis)
    println("  (none)")
else
    g = Dict{Tuple{String, String, String}, Int}()
    for r in mis
        g[(r["tele"], r["labeled"], r["pred"])] = get(g,
            (r["tele"], r["labeled"], r["pred"]), 0) + 1
    end
    @printf("%-5s %-24s %-24s %8s\n", "tele", "declared (almanac)", "content (classifier)",
        "n")
    for (k, v) in sort(collect(g), by = x -> -x[2])
        @printf("%-5s %-24s %-24s %8d\n", k[1], k[2], k[3], v)
    end
end
println()

# ---- focus rows --------------------------------------------------------------
if focusspec != ""
    println("--- focus exposures ---")
    for spec in split(focusspec, ",")
        p = split(strip(spec), ":")
        length(p) == 3 || (println("  bad --focus item: ", spec); continue)
        k = (String(p[1]), parse(Int, p[2]), parse(Int, p[3]))
        i = findfirst(r -> key(r) == k, pred)
        println()
        println("  ", k[1], " ", k[2], " exp ", k[3])
        if i === nothing
            println("    NOT in the prediction table")
        else
            r = pred[i]
            for c in ["image_type", "labeled", "pred", "prob", "status",
                "predicted_bad", "flagged_bad", "lamp_quartz", "lamp_thar",
                "lamp_une", "n_read", "q90max"]
                haskey(r, c) && @printf("    %-14s %s\n", c, r[c])
            end
        end
        ws = get(wsites, k, Dict{String, Int}())
        if isempty(ws)
            println("    warnings       (none)")
        else
            for (s, n) in sort(collect(ws), by = x -> -x[2])
                @printf("    warning        %-56s %6d records\n", s, n)
            end
        end
    end
    println()
end
