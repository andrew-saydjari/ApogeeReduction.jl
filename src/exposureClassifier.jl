# Exposure-type classifier: predicts the exposure type from the ar2D image
# alone and compares it to the (possibly wrong) commanded image_type label.
# Mislabeled calibration exposures (e.g. ThAr/UNe swaps, FPI labeled arclamp)
# can silently poison wavelength calibration; this check runs right after the
# 2D stage and raises warnings — it never mutates labels on its own.
#
# Design notes:
# - Features are computed from dimage ONLY. Commanded metadata (n_read, lamp
#   flags, exptime) is deliberately excluded: it is generated alongside the
#   label being checked, so it would leak the label.
# - Class labels are image_type + lamp flags (e.g. "arclamp_q0t1u0" = ThAr,
#   "arclamp_q0t0u1" = UNe, "arclamp_q0t0u0" = FPI).
# - Predictions with max class probability < unknown_tau are reported as
#   "unknown" (image doesn't resemble any trained class).
# - The model artifact (random forest) is trained offline; see
#   scripts under the 2026_07_14 scratch dir (train_classifier.jl et al.).

using JLD2, Statistics, StatsBase
using LinearAlgebra: dot
using DecisionTree: apply_forest_proba

const CLASSIFIER_QUANTS = [0.001, 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 0.999]
const CLASSIFIER_NCHIPFEAT = length(CLASSIFIER_QUANTS) + 12
# labels that imply the detector was illuminated; predicted dark content under
# one of these labels means the lamp/LED/sky never delivered light
const CLASSIFIER_ILLUMINATED = ["arclamp_q0t1u0", "arclamp_q0t0u1",
    "arclamp_q0t0u0", "quartzflat_q1t0u0", "domeflat_q0t0u0",
    "internalflat_q0t0u0", "twilightflat_q0t0u0"]
# exposure types bright enough to leave persistence in the following dark
const CLASSIFIER_PERSIST_SOURCES = ["internalflat", "quartzflat", "domeflat",
    "arclamp"]

# Default exposure-type classifier artifact. The post-2D check runs with this
# unless it is explicitly disabled (`pipeline.jl --exp_class_model ""`).
#
# THE VERSION IS PINNED ON PURPOSE. This must never become "whatever the newest
# file in that directory is": the class list, the feature layout and the
# decision thresholds are all properties of one particular trained forest, and a
# reduction has to be able to say which forest judged it. Retraining means
# editing this line deliberately, not dropping a new file beside the old one.
#
# Location caveat, recorded honestly: this is a dated ANALYSIS directory. That
# is the same convention the other cal-artifact defaults already follow
# (caldir_darks, caldir_flats, gain_read_cal_dir in pipeline.jl), but it is not
# a stable artifact store and it would vanish under a scratch cleanup. The path
# is defined ONCE, here, precisely so that promoting the artifact to a canonical
# home is a one-line change and never a hunt through the DAG scripts. pipeline.jl
# fails fast if it is missing, so a cleanup can never silently downgrade a run to
# "no classification". See the PR discussion for the proposed canonical home.
const DEFAULT_EXP_CLASS_MODEL = "/mnt/ceph/users/sdssv/work/asaydjari/2026_07_14/meta/exposure_classifier_rf_v6.jld2"

# Tri-state exposure-class verdict carried into the 1D data products. -1 is a
# first-class value, not a filler: see the block comment above
# `exposure_class_unknown_metadata` for why this is not a Bool and not a bitmask.
const EXP_CLASS_BAD_UNKNOWN = Int8(-1)
const EXP_CLASS_BAD_FALSE = Int8(0)
const EXP_CLASS_BAD_TRUE = Int8(1)
# status sentinel meaning "the exposure-type check did not run for this exposure"
const EXP_CLASS_STATUS_NOTRUN = "notrun"
# pred/labeled sentinel for the same case (never "", which reads as a real class)
const EXP_CLASS_UNKNOWN_STR = "unknown"

"""
Count strict local maxima of profile `p` above `thresh`.
"""
classifier_countpeaks(p, thresh) = sum(@views (p[2:(end - 1)] .> p[1:(end - 2)]) .&
                                               (p[2:(end - 1)] .>= p[3:end]) .&
                                               (p[2:(end - 1)] .> thresh))

"""
Max normalized autocorrelation of profile `p` over lags 5:300 (FPI comb detector).
"""
function classifier_autocorr_peak(p)
    x = p .- mean(p)
    v = sum(abs2, x)
    v <= 0 && return 0.0
    best = 0.0
    for lag in 5:300
        c = @views dot(x[1:(end - lag)], x[(1 + lag):end]) / v
        c > best && (best = c)
    end
    best
end

"""
    exposure_class_features(dimage)

Compute the per-chip feature vector for the exposure-type classifier from a
2D image (science region rows 1:N_XPIX; dim1 = spectral, dim2 = fiber axis).
"""
function exposure_class_features(dimage)
    sci = @view dimage[1:N_XPIX, :]
    nanmask = isnan.(sci)
    nanfrac = mean(nanmask)
    s = copy(sci)
    s[nanmask] .= 0.0

    sub = vec(@view s[1:2:end, 1:2:end])
    q = quantile(sub, CLASSIFIER_QUANTS)
    madv = mad(sub; normalize = true)

    ps = vec(mean(s, dims = 2))
    ps_med = median(ps)
    ps_mad = mad(ps; normalize = true) + 1e-12
    spec_maxrat = (maximum(ps) - ps_med) / ps_mad
    spec_npk5 = classifier_countpeaks(ps, ps_med + 5 * ps_mad)
    spec_npk20 = classifier_countpeaks(ps, ps_med + 20 * ps_mad)
    pss = sort(ps, rev = true)
    spec_conc20 = sum(pss) > 0 ? sum(@view pss[1:20]) / (sum(pss) + 1e-12) : 0.0
    spec_acpk = classifier_autocorr_peak(ps)

    pf = vec(mean(s, dims = 1))
    pf_q10, pf_q90 = quantile(pf, [0.1, 0.9])
    fib_contrast = (pf_q90 - pf_q10) / (abs(pf_q90) + abs(pf_q10) + 1e-12)
    pf_med = median(pf)
    pf_mad = mad(pf; normalize = true) + 1e-12
    fib_npk2 = classifier_countpeaks(pf, pf_med + 2 * pf_mad)
    fib_acpk = classifier_autocorr_peak(pf)

    Float64[q..., madv, nanfrac,
        ps_med, spec_maxrat, spec_npk5, spec_npk20, spec_conc20, spec_acpk,
        pf_med, fib_contrast, fib_npk2, fib_acpk]
end

"""
    exposure_class_label(image_type, lamp_quartz, lamp_thar, lamp_une)

Build the class label string used by the classifier from almanac metadata,
e.g. ("arclamp", false, true, false) -> "arclamp_q0t1u0" (= ThAr).
Lamp values may be Bool, 0/1, or "T"/"F" strings; anything else maps to "?".
"""
function exposure_class_label(image_type, lamp_quartz, lamp_thar, lamp_une)
    lampstr(x) = (x in (true, 1, "T", "true", "True")) ? "1" :
                 ((x in (false, 0, "F", "false", "False")) ? "0" : "?")
    lowercase(string(image_type)) * "_q" * lampstr(lamp_quartz) *
    "t" * lampstr(lamp_thar) * "u" * lampstr(lamp_une)
end

"""
    load_exposure_classifier(model_path)

Load the trained exposure-type classifier artifact (random forest, classes,
and decision thresholds) saved by the offline training scripts.
"""
function load_exposure_classifier(model_path)
    d = JLD2.load(model_path)
    (model = d["model"], classes = d["classes"],
        unknown_tau = d["unknown_tau"], flag_tau = d["flag_tau"])
end

"""
    classify_exposure_type(clf, chip_features, tele)

Given per-chip feature vectors (Dict chip => features from
`exposure_class_features`), predict the exposure class label.
Returns (pred, prob, decision) where decision is one of
"ok"/"unknown", and pred is the predicted class label string.
"""
function classify_exposure_type(clf, chip_features::AbstractDict, tele)
    x = reduce(vcat, [chip_features[chip] for chip in CHIP_LIST])
    q90i = findfirst(==(0.9), CLASSIFIER_QUANTS)
    r = chip_features[CHIP_LIST[1]][q90i]
    g = chip_features[CHIP_LIST[2]][q90i]
    b = chip_features[CHIP_LIST[3]][q90i]
    colors = [(r - g) / (abs(r) + abs(g) + 1e-3), (b - g) / (abs(b) + abs(g) + 1e-3),
        (r - b) / (abs(r) + abs(b) + 1e-3)]
    x = vcat(x, colors, lowercase(tele) == "apo" ? 1.0 : 0.0)
    x[.!isfinite.(x)] .= -9999.0

    P = apply_forest_proba(clf.model, reshape(x, 1, :), clf.classes)
    maxp, mi = findmax(vec(P))
    pred = clf.classes[mi]
    decision = maxp < clf.unknown_tau ? "unknown" : "ok"
    (pred = pred, prob = maxp, decision = decision, proba = vec(P))
end

"""
    exposure_check_category(labeled_class, res)

Interpret a classification result against the commanded label (v5 content
taxonomy: the classifier predicts what the image IS; causes are expressed
here). Categories:
- "unknown": image doesn't resemble any trained class
- "faint_twilight": labeled twilightflat but content is faint line-sky
  (counts too low to flat-field with)
- "lamp_off_candidate": labeled as an illuminated exposure but content is a
  dark (lamp/LED off — equivalent to having taken a dark)
- "mislabel_candidate": confident content prediction that contradicts the
  label (e.g. ThAr/UNe swaps)
- "ok"
"""
function exposure_check_category(labeled_class, res, flag_tau = 0.7)
    if res.decision == "unknown"
        "unknown"
    elseif labeled_class == "twilightflat_q0t0u0" &&
           res.pred == "twilightflat_faint" && res.prob > flag_tau
        "faint_twilight"
    elseif labeled_class in CLASSIFIER_ILLUMINATED && res.pred == "dark_q0t0u0" &&
           res.prob > flag_tau
        "lamp_off_candidate"
    elseif res.pred != labeled_class && res.prob > flag_tau
        "mislabel_candidate"
    else
        "ok"
    end
end

"""
    exposure_predicted_bad(labeled_class, pred, status)

Masking policy for downstream consumers (cal runlists, wavecal arc/FPI
selection): should this exposure be excluded based on the classifier verdict?

- object_q0t0u0: never masked (science frames are handled downstream).
  n.b. this exemption is keyed on the EXACT label string, not on
  `image_type == "object"`. An object frame with anomalous lamp flags is
  labeled e.g. "object_q0t0u1" and does NOT take this branch, so it can come
  back masked (4 such frames exist in DR21, all on lco 57802). That is
  harmless today because the only consumer of the mask is the fiber-flat
  runlist builder, which selects on `image_type == "<flat_type>flat"` first and
  so can never see an object frame. If a future consumer masks science
  exposures with this, revisit the exemption before doing so.
- dark_q0t0u0: masked when the content prediction is anything but a clean
  dark (dark_persist, illuminated content) or the prediction is unknown.
  Sequence-only persistence risks whose image still classifies as a clean
  dark are kept.
- internalflat_q0t0u0 and arclamp_q0t0u0 (FPI): masked on any
  predicted-vs-labeled mismatch or unknown
- all other classes: masked whenever the check status is not "ok"
  (mislabel / lamp-off / faint-twilight / unknown / rare label)
- exposures without reduced 2D data ("nofiles"/"unclassified") are not
  masked here; missing products already exclude them downstream

This returns a Bool: it answers "did the classifier judge this bad?", and it
presumes the classifier actually ran. It cannot express "we do not know".
Consumers that must distinguish "judged good" from "never judged" want the
tri-state `exposure_class_metadata` / `EXP_CLASS_BAD_*` encoding instead —
see the note there on why a Bool is the wrong type at the file-metadata layer.
"""
function exposure_predicted_bad(labeled_class, pred, status)
    if status in ("nofiles", "unclassified", EXP_CLASS_STATUS_NOTRUN, "checkfail")
        false
    elseif labeled_class == "object_q0t0u0"
        false
    elseif labeled_class in ("dark_q0t0u0", "internalflat_q0t0u0", "arclamp_q0t0u0")
        (pred != labeled_class) || (status == "unknown")
    else
        status != "ok"
    end
end

##### Exposure-class provenance carried into the 1D data products #####
#
# Representation notes (deliberate, please read before extending):
#
# - `exp_class_predicted_bad` is a TRI-STATE Int8, not a Bool and not a bitmask:
#       1 = the classifier ran and judged this exposure bad
#       0 = the classifier ran and judged this exposure fine
#      -1 = UNKNOWN. The classifier never ran (it is off by default), its
#           artifact was absent, or the check errored on this exposure.
#   A missing verdict must never read as 0/"fine" to a downstream consumer.
#   That is the entire reason this is not a Bool.
#
# - It is a scalar owned by ONE producer (the exposure-type classifier), not a
#   shared bitmask. Other advisory per-exposure flags — e.g. the carton-derived
#   ENGINEERING bit — get their own namespaced scalar rather than a bit in this
#   one, so that two producers can never race over the same integer.
#
# - Field names are prefixed `exp_class_` so the exposure-level namespace stays
#   visibly separate from the per-pixel bitmasks documented in the README.

"""
    exposure_class_unknown_metadata()

The `exp_class_*` metadata block for an exposure the classifier never judged.
Every field is an explicit unknown sentinel: `predicted_bad = -1`,
`status = "notrun"`, `pred`/`labeled = "unknown"`, `prob = NaN`.

This is what a consumer sees when `--exp_class_model` was empty (the default),
when the per-MJD `exposureTypeCheck_*.h5` is absent, or when reading a 1D file
written before this field block existed and thus carrying none of it.
"""
exposure_class_unknown_metadata() = Dict{String, Any}(
    "exp_class_predicted_bad" => EXP_CLASS_BAD_UNKNOWN,
    "exp_class_status" => EXP_CLASS_STATUS_NOTRUN,
    "exp_class_pred" => EXP_CLASS_UNKNOWN_STR,
    "exp_class_labeled" => EXP_CLASS_UNKNOWN_STR,
    "exp_class_prob" => NaN)

"""
    exposure_class_metadata(labeled, pred, prob, status)

Build the `exp_class_*` metadata block for one exposure from a classifier
verdict. `status` is the category from `exposure_check_category` (plus
"persistence_prior", which `pipeline.jl` adds, and "checkfail").

A "checkfail" verdict maps to the UNKNOWN block, not to bad: an exception while
evaluating the check is a failure to form an opinion, and must not be recorded
as an adverse opinion about the data.
"""
function exposure_class_metadata(labeled, pred, prob, status)
    (status == "checkfail" || status == EXP_CLASS_STATUS_NOTRUN) &&
        return exposure_class_unknown_metadata()
    bad = exposure_predicted_bad(labeled, pred, status) ?
          EXP_CLASS_BAD_TRUE : EXP_CLASS_BAD_FALSE
    Dict{String, Any}(
        "exp_class_predicted_bad" => bad,
        "exp_class_status" => String(status),
        "exp_class_pred" => String(pred),
        "exp_class_labeled" => String(labeled),
        "exp_class_prob" => Float64(prob))
end

"""
    exposure_type_check_path(outdir, tele, mjd)

Path of the per-MJD exposure-type check table written by `pipeline.jl` between
the 2D and 1D stages.
"""
exposure_type_check_path(outdir, tele, mjd) = joinpath(
    outdir, "apred", string(mjd), "exposureTypeCheck_$(tele)_$(mjd).h5")

"""
    read_exposure_type_check(path) -> Dict{Int, Dict{String, Any}}

Read a per-MJD `exposureTypeCheck_*.h5` into expnum => `exp_class_*` metadata
block. Returns an empty Dict if the file is absent or unreadable, so that every
caller degrades to the explicit-unknown block rather than failing: this table is
advisory metadata and must never be able to break a reduction.

Tables written before `predicted_bad` was stored are handled by recomputing it
from (labeled, pred, flag) with `exposure_predicted_bad`.
"""
function read_exposure_type_check(path)
    out = Dict{Int, Dict{String, Any}}()
    isfile(path) || return out
    try
        d = load(path)
        expnum = d["expnum"]
        labeled, pred = d["labeled"], d["pred"]
        prob, flag = d["prob"], d["flag"]
        for i in eachindex(expnum)
            out[Int(expnum[i])] = exposure_class_metadata(
                labeled[i], pred[i], prob[i], flag[i])
        end
    catch e
        @warn "Could not read exposure-type check table $path; treating as unknown" exception=e
        return Dict{Int, Dict{String, Any}}()
    end
    out
end

"""
    check_exposure_type!(warnings, clf, fnames, tele, mjd, expnum, image_type, lamp_flags)

Post-2D hook: compute features for the three chip ar2D files `fnames`
(Dict chip => path), classify, and compare to the labeled class
(image_type + lamp flags, e.g. "arclamp_q0t1u0"). Pushes a NamedTuple row
into `warnings` when the prediction disagrees confidently or is unknown.
Returns the classification result.
"""
function check_exposure_type!(warnings, clf, fnames::AbstractDict, tele, mjd, expnum,
        labeled_class)
    chip_features = Dict(chip => exposure_class_features(load(fnames[chip], "dimage"))
    for chip in CHIP_LIST)
    res = classify_exposure_type(clf, chip_features, tele)
    flagged = exposure_check_category(labeled_class, res, clf.flag_tau)
    if flagged != "ok"
        @warn "Exposure-type check: $tele $mjd exp $expnum labeled '$labeled_class' but classified '$(res.pred)' (p=$(round(res.prob, digits = 3)), $flagged)"
        push!(warnings,
            (tele = tele, mjd = mjd, expnum = expnum, labeled = labeled_class,
                pred = res.pred, prob = res.prob, flag = flagged))
    end
    res
end
