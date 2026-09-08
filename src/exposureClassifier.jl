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

# The classifier artifact is a CALIBRATION INPUT, configured by path exactly
# like caldir_darks / caldir_flats / gain_read_cal_dir: pipeline.jl's arg table
# holds the default and airflow/dags/ar_common.py sets it for production. There
# is deliberately no model path constant in src/ — a calibration input does not
# belong baked into library source.

##### Exposure-level flag bitmask #####
#
# SHARED SCHEMA. These four definitions are byte-identical to the ones on
# feature/engineering-carton-flag (PR #397) and must stay that way: two
# producers write into one `exposure_flags` byte, and a disagreement about the
# bit numbering would be silent and unrecoverable. On merge, keep one copy.
#
# 2^0: predicted_bad — the exposure-type classifier's verdict
#      (policy: `exposure_predicted_bad`)
const EXPFLAG_PREDICTED_BAD = 0x01
# 2^1: engineering exposure — the configuration's science fibers were assigned
#      an engineering carton, so the frame exists to exercise the hardware, not
#      to do science (policy: `exposure_is_engineering`, PR #397)
const EXPFLAG_ENGINEERING = 0x02
# 2^2: NO VERDICT WAS FORMED for this exposure by the exposure-type classifier
#      — the check was deliberately disabled, its per-MJD table was missing, or
#      it threw while evaluating this frame.
#
#      This bit is why `exposure_flags == 0` is meaningful: zero means JUDGED
#      AND FINE, not "we never looked". "Never looked" is 2^2, a different
#      value, distinguishable from the byte alone.
#
#      MUTUALLY EXCLUSIVE with 2^0: if no verdict was formed there is no verdict
#      to be adverse. `exposure_flag_bits` asserts this; a byte carrying both is
#      a bug, not a state.
const EXPFLAG_NOTRUN = 0x04
# bits that mean "do not do science with this exposure"; prior builds and any
# other science-sample assembly must exclude these. Reduction must NOT.
#
# EXPFLAG_NOTRUN is deliberately NOT included: an unjudged exposure is not a
# known-bad one, and silently excluding everything we failed to look at would
# turn a monitoring gap into invisible data loss. A caller that wants "only
# frames positively cleared" must test the notrun bit itself, explicitly.
const EXPFLAG_NO_SCIENCE = EXPFLAG_PREDICTED_BAD | EXPFLAG_ENGINEERING

"""
    exposure_flag_bits(predicted_bad, engineering; notrun = false)

Pack the exposure-level verdicts into the `exposure_flags` bitmask
(`EXPFLAG_PREDICTED_BAD`, `EXPFLAG_ENGINEERING`, `EXPFLAG_NOTRUN`).

The positional signature is shared with PR #397; `notrun` is a keyword so that
producers which only know the science verdicts call it unchanged.

Throws if `notrun` and `predicted_bad` are both set: "no verdict was formed" and
"the verdict was adverse" cannot both be true, and a byte asserting both would
be silently misread by every consumer downstream.
"""
function exposure_flag_bits(predicted_bad::Bool, engineering::Bool; notrun::Bool = false)
    notrun && predicted_bad &&
        throw(ArgumentError("exposure_flags: EXPFLAG_NOTRUN and EXPFLAG_PREDICTED_BAD " *
                            "are mutually exclusive — no verdict cannot also be an adverse verdict"))
    b = 0x00
    predicted_bad && (b |= EXPFLAG_PREDICTED_BAD)
    engineering && (b |= EXPFLAG_ENGINEERING)
    notrun && (b |= EXPFLAG_NOTRUN)
    b
end

"""
    exposure_ok_for_science(flags)

True when none of the `EXPFLAG_NO_SCIENCE` bits are set. This is the single
predicate every science-sample assembler (prior builds, catalog construction)
should use. It is deliberately NOT consulted anywhere in the 2D/1D reduction:
engineering and predicted-bad exposures are still reduced to 1D.

n.b. this answers "is it flagged?", NOT "was it judged?" — `EXPFLAG_NOTRUN` is
deliberately not one of these bits, so an unjudged exposure reads as ok here. A
caller that needs to distinguish "judged fine" from "never judged" must use
`exposure_class_verdict`.
"""
exposure_ok_for_science(flags::Integer) = (UInt8(flags) & EXPFLAG_NO_SCIENCE) == 0x00

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
presumes the classifier actually ran. It cannot express "we do not know" — that
information lives one level up, in whether an `exposure_flags` field was written
at all. See `exposure_class_metadata` and `exposure_class_verdict`.
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
# REPRESENTATION (deliberate; read before extending).
#
# The stored verdict is the shared `exposure_flags` UInt8 bitmask, so that the
# exposure-type classifier and the engineering-carton check (PR #397) write into
# one byte with one agreed bit numbering, instead of each owning a private
# scalar that a consumer would have to know to combine.
#
# "No verdict was formed" is a STATE OF THE MASK, not the absence of the field:
# it is `EXPFLAG_NOTRUN` (2^2). The consequences a reader must internalise:
#
#   exposure_flags == 0            -> JUDGED, and nothing wrong. Not "unknown".
#   exposure_flags & NOTRUN != 0   -> never judged; bit 2^0 is guaranteed clear
#   exposure_flags & PREDICTED_BAD -> judged, and adverse
#
# so "judged fine" and "never judged" are different VALUES, distinguishable from
# the byte alone. That is the whole point of spending a bit on it.
#
# The companion `exp_class_status` string still carries WHY there is no verdict
# ("notrun" = disabled or table missing, "checkfail" = the check threw). Those
# share one bit deliberately: no consumer would act differently on them, the
# distinction is diagnostic rather than actionable, and bits in a byte shared
# between two producers are scarce enough not to spend on diagnostics.
#
# BACKWARD COMPATIBILITY: a product written before this PR has no
# `exposure_flags` field at all. `exposure_class_verdict` treats an absent field
# as unknown too, so both routes converge. That is a compatibility shim, not the
# design — new writers always emit the field.
#
# Field names are prefixed `exp_class_` (except the shared `exposure_flags`
# itself) so the exposure-level namespace stays visibly separate from the
# per-pixel bitmasks documented in the README.

"""
    exposure_class_unknown_metadata(status = EXP_CLASS_STATUS_NOTRUN)

The metadata block for an exposure the classifier never judged: `exposure_flags`
carries `EXPFLAG_NOTRUN` (and, by the mutual-exclusion invariant, never
`EXPFLAG_PREDICTED_BAD`).

`status` records WHY no verdict exists — "notrun" when the check was disabled or
its table was missing, "checkfail" when the check threw on this exposure. Both
set the same bit; only this string tells them apart.

This is what a consumer sees when the check was deliberately disabled
(`--exp_class_model ""`), when the per-MJD `exposureTypeCheck_*.h5` is absent,
or when the check errored. A 1D file written before these fields existed carries
no `exposure_flags` at all, and reads as unknown by the compatibility path in
`exposure_class_verdict`.
"""
exposure_class_unknown_metadata(status = EXP_CLASS_STATUS_NOTRUN) = Dict{String, Any}(
    "exposure_flags" => exposure_flag_bits(false, false; notrun = true),
    "exp_class_status" => String(status),
    "exp_class_pred" => EXP_CLASS_UNKNOWN_STR,
    "exp_class_labeled" => EXP_CLASS_UNKNOWN_STR,
    "exp_class_prob" => NaN)

"""
    exposure_class_metadata(labeled, pred, prob, status; engineering = false)

Build the metadata block for one exposure from a classifier verdict, including
the `exposure_flags` bitmask. `status` is the category from
`exposure_check_category` (plus "persistence_prior", which `pipeline.jl` adds,
and "checkfail").

A "checkfail" verdict returns the UNKNOWN block — `EXPFLAG_NOTRUN` set, never
`EXPFLAG_PREDICTED_BAD` — and keeps "checkfail" as its status string. An
exception while evaluating the check is a failure to form an opinion, and must
never be recorded as an adverse opinion, nor as a favourable one.

`engineering` is the hand-off point for PR #397's carton-derived bit. This
producer only knows the classifier's verdict, so it passes `false` by default;
the engineering check ORs its bit in here rather than writing a competing field.
"""
function exposure_class_metadata(labeled, pred, prob, status; engineering::Bool = false)
    (status == "checkfail" || status == EXP_CLASS_STATUS_NOTRUN) &&
        return exposure_class_unknown_metadata(status)
    bad = exposure_predicted_bad(labeled, pred, status)
    Dict{String, Any}(
        "exposure_flags" => exposure_flag_bits(bad, engineering),
        "exp_class_status" => String(status),
        "exp_class_pred" => String(pred),
        "exp_class_labeled" => String(labeled),
        "exp_class_prob" => Float64(prob))
end

"""
    exposure_class_verdict(md) -> Symbol

Three-way read of an exposure-class metadata block (as returned by
`exposure_class_metadata`, or read back from a 1D product's `metadata` group):

- `:unknown` — no verdict was formed: `EXPFLAG_NOTRUN` is set, or (compatibility
  path) the product predates `exposure_flags` and has no such field.
- `:bad` — a verdict was formed and `EXPFLAG_PREDICTED_BAD` is set.
- `:fine` — a verdict was formed and that bit is clear.

This is the ONLY correct way to ask "is this exposure fine?" — it is what makes
the notrun bit and the legacy no-field case give the same answer.
"""
function exposure_class_verdict(md::AbstractDict)
    haskey(md, "exposure_flags") || return :unknown   # pre-PR product
    f = UInt8(md["exposure_flags"])
    (f & EXPFLAG_NOTRUN) != 0x00 && return :unknown
    (f & EXPFLAG_PREDICTED_BAD) != 0x00 ? :bad : :fine
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

Read a per-MJD `exposureTypeCheck_*.h5` into expnum => metadata block. Returns
an empty Dict if the file is absent or unreadable, so that every caller degrades
to the explicit-unknown block rather than failing: this table is advisory
metadata and must never be able to break a reduction.

The block is rebuilt from `(labeled, pred, prob, flag)` through
`exposure_class_metadata`, so a table written by any version of the pipeline
yields the current encoding and the masking policy lives in exactly one place.
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
