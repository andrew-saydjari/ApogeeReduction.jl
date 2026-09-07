# classify_run.jl — post-hoc exposure-type classifier sweep over a COMPLETED run.
#
# Why post-hoc
# ------------
# The in-pipeline check (pipeline.jl --exp_class_model) runs the classifier
# inside the reduction.  The goldens/testbed runs deliberately do NOT enable it
# (AR_EXP_CLASS_MODEL="" — matching production run_all.sh), so that the run stays
# comparable to previous runs.  That leaves the classifier's zero warning count
# uninformative: it could not have fired (see WARNINGS.md §7).
#
# This script recovers the signal without touching the run: it walks the
# delivered ar2D products, runs the SAME functions the pipeline would have run
# (ApogeeReduction.exposure_class_features / classify_exposure_type /
# exposure_check_category / exposure_predicted_bad), and writes its own table.
# It is strictly read-only over the run's products.
#
# Being post-hoc is not a compromise here — it is the stronger position for the
# question actually being asked.  The warnings already exist; running the
# classifier afterwards makes it an INDEPENDENT second opinion on them rather
# than a co-author of the same log.
#
# Usage
# -----
#   julia --project=<ApogeeReduction.jl> test/regression/classify_run.jl \
#       --outdir  /path/to/testbed_run \
#       --model   /path/to/exposure_classifier_rf_v6.jld2 \
#       --output  /path/to/exposure_class_predictions.tsv \
#       [--almanac PATH] [--nworkers 10] [--tele apo|lco] [--mjd 60255]
#       [--limit N] [--append]
#
# --almanac defaults to the single *.h5 in <outdir>/almanac/ that is not a
# runlist_/valid_ file (the bulk almanac the run consumed).
#
# Output: TSV, one row per (tele, mjd, exposure), columns
#   tele mjd exposure image_type labeled pred prob status predicted_bad
#   flagged_bad lamp_quartz lamp_thar lamp_une n_read q90max nchips
# `status` is the pipeline's flag vocabulary: ok | unknown | mislabel_candidate |
# lamp_off_candidate | faint_twilight | persistence_prior | checkfail.
#
# Cost (MEASURED on ccalin051, 2026_09_03 testbed): ~1.7 s/exposure/core, so
# 17_911 exposures is ~8.5 core-hours => ~55 min at --nworkers 10.  It is
# I/O-bound on three ~2048x2560 HDF5 reads per exposure, not CPU-bound.

using Distributed, ArgParse, Printf

function parse_args_local()
    s = ArgParseSettings(description = "Post-hoc exposure-type classifier sweep over a completed run")
    @add_arg_table s begin
        "--outdir"
        required = true
        help = "run root (the directory containing apred/ and almanac/)"
        "--model"
        required = true
        help = "trained classifier artifact (.jld2)"
        "--output"
        required = true
        help = "output TSV path"
        "--almanac"
        default = ""
        help = "almanac .h5 (default: auto-detect in <outdir>/almanac/)"
        "--nworkers"
        arg_type = Int
        default = 10
        help = "local worker processes (keep <= 10 on a shared workstation)"
        "--tele"
        default = ""
        help = "restrict to one telescope"
        "--mjd"
        default = ""
        help = "restrict to a comma-separated MJD list"
        "--limit"
        arg_type = Int
        default = 0
        help = "process at most N exposures (smoke test)"
        "--append"
        action = :store_true
        help = "append to --output instead of overwriting (resume)"
    end
    parse_args(s)
end

parg = parse_args_local()

const OUTDIR = abspath(parg["outdir"])
const APRED = joinpath(OUTDIR, "apred")
isdir(APRED) || error("no apred/ under $OUTDIR")

# ---- almanac -----------------------------------------------------------------
function find_almanac(outdir)
    d = joinpath(outdir, "almanac")
    isdir(d) || error("no almanac/ under $outdir")
    cands = filter(readdir(d)) do f
        endswith(f, ".h5") && !startswith(f, "runlist_") && !startswith(f, "valid_")
    end
    length(cands) == 1 ||
        error("cannot auto-detect almanac in $d (found $(length(cands)): $cands); pass --almanac")
    joinpath(d, cands[1])
end
const ALMANAC = parg["almanac"] == "" ? find_almanac(OUTDIR) : abspath(parg["almanac"])
isfile(ALMANAC) || error("almanac not found: $ALMANAC")

println("outdir  : ", OUTDIR)
println("almanac : ", ALMANAC)
println("model   : ", parg["model"])
flush(stdout)

# ---- enumerate delivered exposures from the products ------------------------
# Source of truth is the products on disk, not the runlist: we are auditing what
# the run actually delivered.  ar2D (not ar2Dcal) is what the in-pipeline check
# consumes, so we consume the same thing.
const CHIPS = ["R", "G", "B"]

"""Parse ar2D_<tele>_<mjd>_<exp>_<chip>_<imtype>.h5 -> (tele, mjd, exp, chip, imtype)."""
function parse_ar2D(fname)
    b = basename(fname)
    startswith(b, "ar2D_") && endswith(b, ".h5") || return nothing
    p = split(b[1:(end - 3)], "_")
    length(p) < 6 && return nothing
    # imtype never contains '_', so the last six fields are positionally exact
    _, tele, mjds, exps, chip, imtype = p[(end - 5):end]
    (tele in ("apo", "lco") && chip in CHIPS) || return nothing
    mjd = tryparse(Int, mjds)
    ex = tryparse(Int, exps)
    (mjd === nothing || ex === nothing) && return nothing
    (tele = tele, mjd = mjd, expnum = ex, chip = chip, imtype = imtype)
end

telefilter = parg["tele"]
mjdfilter = parg["mjd"] == "" ? Int[] :
            [parse(Int, strip(x)) for x in split(parg["mjd"], ",") if strip(x) != ""]

groups = Dict{NTuple{4, Any}, Dict{String, String}}()
for mjddir in sort(readdir(APRED))
    d = joinpath(APRED, mjddir)
    isdir(d) || continue
    isempty(mjdfilter) || parse(Int, mjddir) in mjdfilter || continue
    for f in readdir(d)
        startswith(f, "ar2D_") || continue
        p = parse_ar2D(f)
        p === nothing && continue
        telefilter == "" || p.tele == telefilter || continue
        key = (p.tele, p.mjd, p.expnum, p.imtype)
        get!(groups, key, Dict{String, String}())[p.chip] = joinpath(d, f)
    end
end
println("exposures with >=1 ar2D chip: ", length(groups))

# ---- attach almanac metadata -------------------------------------------------
# Join on the `exposure` VALUE, not on row position.  pipeline.jl uses
# `df[expnum, :]`, which is only correct while every night's exposure column is
# exactly 1:N.  That holds in this run (checked below) but it is not guaranteed
# by anything, so this script joins properly and reports if the assumption ever
# breaks.
using HDF5, DataFrames

almcache = Dict{Tuple{String, Int}, Any}()
function almrows(tele, mjd)
    get!(almcache, (tele, mjd)) do
        h5open(ALMANAC, "r") do f
            k = "raw/$tele/$mjd/exposures"
            haskey(f, k) ? DataFrame(read(f[k])) : nothing
        end
    end
end

getcol(df, name, i, default) = hasproperty(df, name) ? df[i, name] : default

tasks = Any[]
nmissing_alm = Ref(0)
noncontig = Set{Tuple{String, Int}}()
for (key, chipmap) in groups
    tele, mjd, expnum, imtype = key
    df = almrows(tele, mjd)
    if df === nothing
        nmissing_alm[] += 1
        continue
    end
    if hasproperty(df, :exposure) && df.exposure != collect(1:nrow(df))
        push!(noncontig, (tele, mjd))
    end
    idx = findfirst(==(expnum), df.exposure)
    if idx === nothing
        nmissing_alm[] += 1
        continue
    end
    prevtype = idx > 1 ? lowercase(string(strip(String(df[idx - 1, :image_type])))) : ""
    push!(tasks,
        (tele = tele, mjd = mjd, expnum = expnum, imtype = imtype,
            fnames = chipmap,
            alm_image_type = lowercase(string(strip(String(df[idx, :image_type])))),
            lamp_quartz = getcol(df, :lamp_quartz, idx, -1),
            lamp_thar = getcol(df, :lamp_thar, idx, -1),
            lamp_une = getcol(df, :lamp_une, idx, -1),
            flagged_bad = getcol(df, :flagged_bad, idx, -1),
            n_read = getcol(df, :n_read, idx, -1),
            prevtype = prevtype))
end
sort!(tasks, by = t -> (t.tele, t.mjd, t.expnum))
nmissing_alm[] > 0 && @warn "$(nmissing_alm[]) exposures had no almanac row (skipped)"
if !isempty(noncontig)
    @warn "almanac `exposure` is not 1:N for $(length(noncontig)) night(s); pipeline.jl's positional `df[expnum, :]` would be WRONG there" nights=sort(collect(noncontig))
end

# resume support: skip rows already in the output
done = Set{Tuple{String, Int, Int}}()
if parg["append"] && isfile(parg["output"])
    for ln in eachline(parg["output"])
        startswith(ln, "tele\t") && continue
        p = split(ln, '\t')
        length(p) >= 3 && push!(done, (p[1], parse(Int, p[2]), parse(Int, p[3])))
    end
    filter!(t -> !((t.tele, t.mjd, t.expnum) in done), tasks)
    println("resuming: ", length(done), " already done, ", length(tasks), " remaining")
end
parg["limit"] > 0 && (tasks = tasks[1:min(parg["limit"], length(tasks))])
println("exposures to classify: ", length(tasks))
flush(stdout)
isempty(tasks) && (println("nothing to do"); exit(0))

# ---- workers -----------------------------------------------------------------
nw = max(parg["nworkers"], 1)
proj = dirname(Base.active_project())
nw > 1 && addprocs(nw, exeflags = ["--project=$proj"])
println("workers: ", nworkers())
flush(stdout)

@everywhere begin
    using ApogeeReduction, JLD2
    using ApogeeReduction: exposure_class_features, load_exposure_classifier,
                           classify_exposure_type, exposure_class_label,
                           exposure_check_category, exposure_predicted_bad,
                           CHIP_LIST, CLASSIFIER_QUANTS, CLASSIFIER_PERSIST_SOURCES
    const _CLF = Ref{Any}(nothing)
    const _MODELPATH = Ref{String}("")
    function get_clf()
        _CLF[] === nothing && (_CLF[] = load_exposure_classifier(_MODELPATH[]))
        _CLF[]
    end

    function classify_one(t)
        # Mirrors pipeline.jl's exp_type_check_one exactly, including the
        # persistence_prior post-step, so that this post-hoc pass and the
        # in-pipeline check cannot silently diverge.
        labeled = exposure_class_label(t.alm_image_type, t.lamp_quartz, t.lamp_thar,
            t.lamp_une)
        if length(t.fnames) != length(CHIP_LIST)
            return (t..., labeled = labeled, pred = "nofiles", prob = NaN,
                status = "nofiles", predicted_bad = false, q90max = NaN,
                nchips = length(t.fnames))
        end
        try
            clf = get_clf()
            cf = Dict(c => exposure_class_features(JLD2.load(t.fnames[c], "dimage"))
            for c in CHIP_LIST)
            res = classify_exposure_type(clf, cf, t.tele)
            status = exposure_check_category(labeled, res, clf.flag_tau)
            if status == "ok" && res.pred == "dark_q0t0u0" && labeled == "dark_q0t0u0" &&
               t.prevtype in CLASSIFIER_PERSIST_SOURCES
                status = "persistence_prior"
            end
            q90i = findfirst(==(0.9), CLASSIFIER_QUANTS)
            q90max = maximum(cf[c][q90i] for c in CHIP_LIST)
            (t..., labeled = labeled, pred = res.pred, prob = res.prob,
                status = status,
                predicted_bad = exposure_predicted_bad(labeled, res.pred, status),
                q90max = q90max, nchips = length(t.fnames))
        catch e
            (t..., labeled = labeled, pred = "checkfail", prob = NaN,
                status = "checkfail", predicted_bad = false, q90max = NaN,
                nchips = length(t.fnames))
        end
    end
end
@everywhere _MODELPATH[] = $(abspath(parg["model"]))

const COLS = ["tele", "mjd", "exposure", "image_type", "labeled", "pred", "prob",
    "status", "predicted_bad", "flagged_bad", "lamp_quartz", "lamp_thar",
    "lamp_une", "n_read", "q90max", "nchips"]

fmt(x::AbstractFloat) = isfinite(x) ? @sprintf("%.4g", x) : "NaN"
fmt(x::Bool) = x ? "1" : "0"
fmt(x) = string(x)

newfile = !(parg["append"] && isfile(parg["output"]))
mkpath(dirname(abspath(parg["output"])))
# Chunked so the TSV is flushed as we go: a killed job can be resumed with
# --append rather than restarting an hour of I/O.
const CHUNK = 500
open(parg["output"], newfile ? "w" : "a") do io
    newfile && println(io, join(COLS, '\t'))
    t0 = time()
    ndone = 0
    for lo in 1:CHUNK:length(tasks)
        hi = min(lo + CHUNK - 1, length(tasks))
        for r in pmap(classify_one, tasks[lo:hi]; batch_size = 4)
            println(io,
                join([fmt(r.tele), fmt(r.mjd), fmt(r.expnum), fmt(r.alm_image_type),
                        fmt(r.labeled), fmt(r.pred), fmt(r.prob), fmt(r.status),
                        fmt(r.predicted_bad), fmt(r.flagged_bad), fmt(r.lamp_quartz),
                        fmt(r.lamp_thar), fmt(r.lamp_une), fmt(r.n_read), fmt(r.q90max),
                        fmt(r.nchips)], '\t'))
            ndone += 1
        end
        flush(io)
        el = (time() - t0) / 60
        @printf("  %6d / %6d  (%.1f min elapsed, ETA %.1f min)\n", ndone, length(tasks),
            el, el * (length(tasks) - ndone) / max(ndone, 1))
        flush(stdout)
    end
    @printf("classified %d exposures in %.1f min\n", ndone, (time() - t0) / 60)
end
println("wrote ", parg["output"])
