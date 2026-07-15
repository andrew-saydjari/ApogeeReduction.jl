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
    flagged = if res.decision == "unknown"
        "unknown"
    elseif res.pred != labeled_class && res.prob > clf.flag_tau
        "mislabel_candidate"
    else
        "ok"
    end
    if flagged != "ok"
        @warn "Exposure-type check: $tele $mjd exp $expnum labeled '$labeled_class' but classified '$(res.pred)' (p=$(round(res.prob, digits = 3)), $flagged)"
        push!(warnings,
            (tele = tele, mjd = mjd, expnum = expnum, labeled = labeled_class,
                pred = res.pred, prob = res.prob, flag = flagged))
    end
    res
end
