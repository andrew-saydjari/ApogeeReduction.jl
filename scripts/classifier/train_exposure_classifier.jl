# Train the exposure-type classifier (v7).
#
# Same architecture and feature set as the deployed v6 model (random forest,
# 300 trees, 0.7 feature subsampling; 21 features x 3 chips + 3 chip-color
# ratios + telescope flag; content-class taxonomy of train_classifier_v6.jl in
# the 2026_07_14 scratch dir). What changed is the TRAINING DATA, using truth
# sources the DR21 pass-1 sweep surfaced:
#   1. Twilight-census sun-altitude ground truth: census night-sky flats train
#      as content twilightflat_faint and true solar twilights as
#      twilightflat_q0t0u0/_faint by measured brightness — regardless of
#      commanded image_type. This removes the label noise that had hundreds of
#      sky flats typed "object"/"dark" training the object and dark classes
#      (LCO true twilights are ALL typed object and were invisible before).
#   2. Census daytime dome tests, engineering-carton exposures, FPS-era
#      config_id=-1 object frames, and the lco 57800 stale frames are
#      excluded from training (role=eval_only in the manifest).
#   3. Features come from the pass-1 reduction via the deployed
#      ApogeeReduction.exposure_class_features (bit-identical to inference).
# Everything else (content-relabel thresholds, dark_persist construction,
# unknown/flag taus) is v6 logic unchanged.
#
# Inputs: manifest CSV (build_training_manifest.jl), features HDF5
# (extract_training_features.jl), the almanac (for previous-exposure types
# feeding the dark_persist sequence prior), and optionally the round-1 human
# labels CSV. Outputs: results HDF5 (cross-validated probabilities and labels),
# the deployable model artifact (JLD2), and the final manifest CSV with
# content_label / v7 predictions filled in.
using Pkg;
Pkg.instantiate();
using HDF5, ArgParse, DataFrames, Random, Statistics, StatsBase, Printf, CSV, JLD2
using LinearAlgebra: diag
using DecisionTree
# run with --project=<this repo>: ApogeeReduction is needed for safe_jldsave, so
# the model artifact records the git provenance of the code that made it
using ApogeeReduction

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--manifest"
        required = true
        arg_type = String
        "--features"
        required = true
        arg_type = String
        "--almanac_file"
        required = true
        help = "almanac h5 (raw/ groups; previous-exposure types for dark_persist labels)"
        arg_type = String
        "--human_labels"
        default = ""
        help = "optional CSV of human verdicts (obs,mjd,expnum,given_class,human_verdict) to hold out"
        arg_type = String
        "--out_results"
        required = true
        arg_type = String
        "--out_confusion"
        required = true
        help = "cross-validated confusion matrix CSV (rows = content label, cols = predicted)"
        arg_type = String
        "--out_model"
        required = true
        arg_type = String
        "--out_manifest"
        required = true
        help = "final manifest CSV with content_label and v7 predictions"
        arg_type = String
    end
    return parse_args(s)
end
parg = parse_commandline()

const CONTENT_CLASSES = ["object_q0t0u0", "dark_q0t0u0", "dark_persist",
    "arclamp_q0t1u0", "arclamp_q0t0u1", "arclamp_q0t0u0", "quartzflat_q1t0u0",
    "domeflat_q0t0u0", "internalflat_q0t0u0", "twilightflat_q0t0u0",
    "twilightflat_faint"]
const CONT_ILLUM = ["quartzflat_q1t0u0", "domeflat_q0t0u0", "internalflat_q0t0u0",
    "arclamp_q0t0u0"]
const ARC_ILLUM = ["arclamp_q0t1u0", "arclamp_q0t0u1"]
const TWIFLAT_SPLIT = 50.0
const UNKNOWN_TAU = 0.5
const FLAG_TAU = 0.7
const NFOLD = 5
const PERSIST_SOURCES = ["internalflat", "quartzflat", "domeflat", "arclamp"]

man = CSV.read(parg["manifest"], DataFrame)
F, ok, featnames = h5open(parg["features"], "r") do f
    read(f["features"]), Bool.(read(f["ok"])), read(f["feature_names"])
end
nexp = nrow(man)
@assert size(F, 3) == nexp "manifest / features row mismatch"
obs, mjd, expnum = String.(man.obs), man.mjd, man.expnum
orig_class = String.(man.labeled)

q90i = findfirst(==("q0.9"), featnames)
smaxi = findfirst(==("spec_maxrat"), featnames)
q90max = [maximum(F[q90i, :, i]) for i in 1:nexp]
arcline = [maximum(F[smaxi, :, i]) for i in 1:nexp]

nchipfeat = size(F, 1)
X = zeros(nexp, nchipfeat * 3 + 3 + 1)
for i in 1:nexp
    x = vec(F[:, :, i])
    r, g, b = F[q90i, 1, i], F[q90i, 2, i], F[q90i, 3, i]
    colors = [(r - g) / (abs(r) + abs(g) + 1e-3), (b - g) / (abs(b) + abs(g) + 1e-3),
        (r - b) / (abs(r) + abs(b) + 1e-3)]
    X[i, :] = vcat(x, colors, obs[i] == "apo" ? 1.0 : 0.0)
end
X[.!isfinite.(X)] .= -9999.0

# ---- content labels ----
content = copy(orig_class)

# (1) census truth first (overrides the commanded label entirely)
isnight = man.truth_source .== "twilight_census_nightsky"
istwi = man.truth_source .== "twilight_census_twilight"
content[isnight] .= "twilightflat_faint"
content[istwi] .= ifelse.(q90max[istwi] .>= TWIFLAT_SPLIT,
    "twilightflat_q0t0u0", "twilightflat_faint")
census_locked = isnight .| istwi
println("census twilight relabels: nightsky ", sum(isnight), ", twilight ", sum(istwi),
    " (bright ", sum(content[istwi] .== "twilightflat_q0t0u0"), ")")

# (2) v6 brightness split for remaining twilightflat-typed rows
twfaint = .!census_locked .& (orig_class .== "twilightflat_q0t0u0") .&
          (q90max .< TWIFLAT_SPLIT)
content[twfaint] .= "twilightflat_faint"

# (3) dark-level frames under illuminated labels -> content = dark
#
# CENSUS PRECEDENCE. Census-truth rows are exempt from this clause in BOTH
# directions, and that exemption is load-bearing — it was measured, not assumed:
#   - they are excluded from the reference pool `dm` that sets the thresholds.
#     Commanded-"dark" exposures the census identifies as night-sky flats carry
#     real sky signal (median census signal ~120 counts vs a clean dark's
#     floor). They are a few dozen rows out of ~1,200 commanded darks per
#     observatory, so they land in the top percentiles and it is THEY, not
#     darks, that set the 0.95/0.99 quantiles. Leaving them in inflates
#     thr_cont / thr_arc and, through the same `cleanmask`, thr_persist in
#     clause (4): the first v7 build got 18 dark_persist labels with them
#     excluded and 16 with them in, against 56 for v6.
#   - their content is not overwritten. With them included, 77 census night-sky
#     flats were relabeled dark_q0t0u0 — and those 77 had the SAME median
#     measured sky signal (119.8) as the 681 kept as twilightflat_faint, i.e.
#     the cut was slicing a homogeneous population rather than separating
#     signal-free from signal-bearing frames. It cost 7 points of night-sky
#     accuracy against census truth (96.6% for v6 -> 89.5%; 98.4% once exempt).
# This is what makes clause (1)'s "overrides the commanded label entirely" true:
# the override has to cover the thresholds derived from that label too.
for o in ("apo", "lco")
    dm = (obs .== o) .& (orig_class .== "dark_q0t0u0") .& ok .& .!census_locked
    thr_cont = quantile(q90max[dm], 0.95)
    thr_arc = quantile(arcline[dm], 0.95)
    om = (obs .== o) .& .!census_locked
    isdarklike = om .& ok .& (
        (in.(content, Ref(CONT_ILLUM)) .& (q90max .< thr_cont)) .|
        (in.(content, Ref(ARC_ILLUM)) .& (arcline .< thr_arc) .& (q90max .< thr_cont)) .|
        ((content .== "twilightflat_faint") .& (arcline .< thr_arc) .& (q90max .< thr_cont)))
    content[isdarklike] .= "dark_q0t0u0"
    println("$o: thr_cont=", round(thr_cont, digits = 2), " thr_arc=",
        round(thr_arc, digits = 1), "; relabeled dark-content: ", sum(isdarklike),
        " (census rows exempt: ", sum((obs .== o) .& census_locked), ")")
end

# (4) dark_persist from sequence prior + flux evidence
prevtype = Dict{Tuple{String, Int, Int}, String}()
h5open(parg["almanac_file"], "r") do f
    for tele in keys(f["raw"])
        for m in keys(f["raw/$tele"])
            haskey(f["raw/$tele/$m"], "exposures") || continue
            g = f["raw/$tele/$m/exposures"]
            te, tt = read(g["exposure"]), read(g["image_type"])
            mi = parse(Int, m)
            for j in eachindex(te)
                prevtype[(tele, mi, te[j])] = tt[j]
            end
        end
    end
end
prev_bright(o, m, e) = get(prevtype, (o, m, e - 1), "") in PERSIST_SOURCES
# census rows exempt, same reasoning as clause (3): a commanded-dark night-sky
# flat is neither a clean dark (it has sky signal, so it corrupts the reference
# quantile) nor a persistence-contaminated dark (its excess flux is sky, not
# afterglow, so calling it dark_persist teaches the class the wrong thing).
darkmask = (orig_class .== "dark_q0t0u0") .& ok .& .!census_locked
cleanmask = darkmask .& .!prev_bright.(obs, mjd, expnum)
thr_persist = Dict(o => quantile(q90max[cleanmask .& (obs .== o)], 0.99)
for o in ("apo", "lco"))
ispersist = darkmask .& prev_bright.(obs, mjd, expnum) .&
            (q90max .> [thr_persist[o] for o in obs])
content[ispersist] .= "dark_persist"
println("dark_persist training labels: ", sum(ispersist), " (thr: ", thr_persist, ")")

println("\ncontent counts (train_pool only):")
show(sort(combine(groupby(DataFrame(c = content[man.role .== "train_pool"]), :c),
        nrow => :count), :count, rev = true), allrows = true)
println()

# ---- training pool ----
excluded = fill("", nexp)
excluded[man.role .== "eval_only"] .= "eval_only"
rare = .!in.(content, Ref(CONTENT_CLASSES))
excluded[rare .& (excluded .== "")] .= "excluded_rare"
if parg["human_labels"] != ""
    hl = CSV.read(parg["human_labels"], DataFrame)
    for r in eachrow(hl)
        i = findfirst((obs .== r.obs) .& (mjd .== r.mjd) .& (expnum .== r.expnum))
        i === nothing && continue
        excluded[i] = "human_" * r.human_verdict
    end
end
intrain = ok .& (excluded .== "")
println("training pool: ", sum(intrain), " (excluded_rare ", sum(rare), ", eval_only ",
    sum(excluded .== "eval_only"), ", human ", sum(startswith.(excluded, "human_")),
    ", readfail ", sum(.!ok), ")")

rng = MersenneTwister(42)
trainidx = findall(intrain)
folds = mod1.(randperm(rng, length(trainidx)), NFOLD)
classes = sort(CONTENT_CLASSES)
proba = fill(NaN, nexp, length(classes))
for k in 1:NFOLD
    tr = trainidx[folds .!= k]
    te = trainidx[folds .== k]
    model = build_forest(content[tr], X[tr, :], -1, 300, 0.7, -1; rng = 10 + k)
    proba[te, :] = apply_forest_proba(model, X[te, :], classes)
end
model_all = build_forest(content[trainidx], X[trainidx, :], -1, 300, 0.7, -1; rng = 7)
restidx = findall(ok .& (excluded .!= ""))
isempty(restidx) ||
    (proba[restidx, :] = apply_forest_proba(model_all, X[restidx, :], classes))

pred = [any(isnan, proba[i, :]) ? "" : classes[argmax(proba[i, :])] for i in 1:nexp]
maxp = [any(isnan, proba[i, :]) ? NaN : maximum(proba[i, :]) for i in 1:nexp]

println("\ncross-validated confusion (rows=content, cols=predicted):")
cm = zeros(Int, length(classes), length(classes))
for i in trainidx
    cm[findfirst(==(content[i]), classes), findfirst(==(pred[i]), classes)] += 1
end
shortname(c) = replace(c, "arclamp_q0t1u0" => "thar", "arclamp_q0t0u1" => "une",
    "arclamp_q0t0u0" => "fpi", "twilightflat_faint" => "twi_faint", r"_q\d.*" => "")
cmdf = DataFrame(cm, Symbol.(shortname.(classes)))
cmdf.given = classes
show(cmdf[:, vcat("given", names(cmdf)[1:(end - 1)])], allrows = true, allcols = true)
println("\ncross-validated accuracy: ", round(sum(diag(cm)) / sum(cm), digits = 5))
CSV.write(parg["out_confusion"], cmdf[:, vcat("given", names(cmdf)[1:(end - 1)])])
println("wrote ", parg["out_confusion"])

h5open(parg["out_results"], "w") do f
    f["proba"] = proba
    f["classes"] = classes
    f["pred"] = pred
    f["maxp"] = collect(maxp)
    f["excluded"] = String.(excluded)
    f["fold"] = begin
        fv = zeros(Int, nexp)
        fv[trainidx] = folds
        fv
    end
    f["observatory"] = obs
    f["mjd"] = collect(mjd)
    f["exposure"] = collect(expnum)
    f["class_label"] = content
    f["orig_label"] = orig_class
end
println("wrote ", parg["out_results"])

# safe_jldsave, not JLD2.jldsave: the artifact is a pinned calibration input, so
# it should carry the same git provenance as every other product. The extra
# "metadata/*" entries are additional top-level keys under JLD2.load, so
# load_exposure_classifier (which reads model/classes/unknown_tau/flag_tau by
# name) is unaffected — verified against a round-trip before switching.
ApogeeReduction.safe_jldsave(parg["out_model"];
    no_metadata = true,
    model = model_all,
    classes = classes,
    feature_names = featnames,
    unknown_tau = UNKNOWN_TAU,
    flag_tau = FLAG_TAU,
    train_date = "2026-09-13",
    train_pool = "pass1_dr21 reduction (allobs_57618_61230); v7: v6 content taxonomy + " *
                 "twilight-census sun-altitude ground-truth relabels (2026_09_11 census); " *
                 "engineering/config_id=-1/daytime/lco57800 excluded from training; " *
                 "manifest: metadata/exposure_classifier_training_manifest_v7.csv")
println("wrote ", parg["out_model"], " (",
    round(filesize(parg["out_model"]) / 1e6, digits = 1), " MB)")

man.content_label = content
man.v7_pred = pred
man.v7_prob = maxp
man.excluded = excluded
man.feat_ok = ok
CSV.write(parg["out_manifest"], man)
println("wrote ", parg["out_manifest"])
