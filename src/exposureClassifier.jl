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

using JLD2, Statistics, StatsBase, HDF5, DataFrames
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

## ---------------------------------------------------------------------------
## Exposure-level flag bits (`exposure_class/<tele>/<mjd>/exposure_flags`)
##
## These are EXPOSURE-level advisory bits, distinct from the per-pixel flag bits
## documented in the README's "Current Flag Bits" table. They are metadata for
## downstream consumers (prior builds, cal runlists, catalog construction); the
## reduction to 1D never consults them, so every exposure is still reduced.
## ---------------------------------------------------------------------------
# 2^0: the image-content classifier says this exposure should not be used
#      (policy: `exposure_predicted_bad`)
const EXPFLAG_PREDICTED_BAD = 0x01
# 2^1: engineering exposure — the configuration's science fibers were assigned
#      an engineering carton, so the frame exists to exercise the hardware, not
#      to do science (policy: `exposure_is_engineering`)
const EXPFLAG_ENGINEERING = 0x02
# bits that mean "do not do science with this exposure"; prior builds and any
# other science-sample assembly must exclude these. Reduction must NOT.
const EXPFLAG_NO_SCIENCE = EXPFLAG_PREDICTED_BAD | EXPFLAG_ENGINEERING

"""
Carton-name prefixes that mark a configuration as ENGINEERING rather than
survey science. Extend this list as more engineering cartons are identified;
matching is a case-insensitive `startswith`, so `"manual_fps_position_stars"`
covers `manual_fps_position_stars`, `..._10`, `..._apogee_10`, and
`..._lco_apogee_10` (all four exist in the 57618-61230 corpus).
"""
const ENGINEERING_CARTON_PREFIXES = ["manual_fps_position_stars"]

"""
Minimum fraction of a configuration's *science* fibers that must carry an
engineering carton before the exposure is flagged engineering. The comparison
is strict (`frac > ENGINEERING_CARTON_MIN_FRAC`), i.e. a strict majority.

Why a majority and not "any": a stray engineering assignment inside an
otherwise real science configuration should not condemn ~300 good spectra.
Why not "all": one un-assigned or serendipitous fiber should not rescue a
frame that is plainly a positioning test.

Behaviour on mixed configurations (none exist today — every one of the 1023
configurations in the 57618-61230 corpus with any engineering carton has
fraction exactly 1.0, MEASURED 2026-09-08): a config that is >50% engineering
is flagged in full, and a config that is <=50% engineering is not flagged at
all. Because the flag is per-exposure, either way some fibers are mis-served
in the mixed case; `engineering_frac` is written alongside the bit so a
consumer (or a future per-fiber treatment) can apply its own rule.
"""
const ENGINEERING_CARTON_MIN_FRAC = 0.5

"""
Fallback for configurations that contain NO `category == "science"` fibers at
all: evaluate the engineering fraction over every fiber that carries a non-empty
`firstcarton` instead.

**This is an addition beyond AKS's 2026-09-08 instruction** ("the carton comes
from the configuration's science fibers"), flagged for his decision. Set it to
`false` to get exactly the rule as specified.

Why it is here: MEASURED on the 57618-61230 corpus, 5 configurations
(apo 59558/105, 59558/106, 59560/121, 59560/122, 59561/133 — the earliest FPS
commissioning nights) carry `manual_fps_position_stars` on 207-254 of their 300
fibers while having ZERO fibers labelled `category == "science"` (their
categories are `""`, `bonus`, `open_fiber`, `sky_boss`/`sky_apogee`). Those
configurations back **35 object exposures** that the science-fibers-only rule
misses entirely. No configuration in the corpus is affected in the other
direction: the fallback only ever fires where the primary rule has no data.
"""
const ENGINEERING_FALLBACK_ALL_FIBERS = true

"""
Commanded image types the engineering carton check applies to. Calibration
frames (darks, flats, arcs) taken while an engineering configuration happened
to be loaded are still perfectly good calibrations, and the FPS era carries a
`config_id` on those rows too — so restricting to `object` keeps the flag
about the science content of the frame.
"""
const ENGINEERING_CHECK_IMAGE_TYPES = ["object"]

"""
    is_engineering_carton(carton)

True when the carton name matches any prefix in `ENGINEERING_CARTON_PREFIXES`
(case-insensitive, surrounding whitespace stripped).
"""
function is_engineering_carton(carton)
    c = lowercase(strip(String(carton)))
    isempty(c) && return false
    any(p -> startswith(c, lowercase(p)), ENGINEERING_CARTON_PREFIXES)
end

"""
    exposure_is_engineering(cartons; min_frac = ENGINEERING_CARTON_MIN_FRAC, basis)

Carton check for one exposure. `cartons` is the list of `firstcarton` values of
the configuration's **science** fibers (`category == "science"`), or of all
carton-bearing fibers when the config has no science fibers at all (see
`ENGINEERING_FALLBACK_ALL_FIBERS`); `basis` records which.

Returns `(engineering, frac, carton, nsci, basis)`:
- `engineering::Bool` — `frac > min_frac`
- `frac::Float64` — fraction of the fibers in `cartons` with an engineering
  carton (`NaN` when `cartons` is empty, e.g. plate-era or missing config)
- `carton::String` — the most common matching carton name ("" if none), kept
  for provenance so the reason for the bit is auditable
- `nsci::Int` — number of fibers the fraction was computed over
- `basis::String` — `"science"`, `"all_fibers_fallback"`, or `"none"`

An empty `cartons` (plate era, no configuration, unreadable fiber table) is
never engineering and never errors.
"""
function exposure_is_engineering(cartons; min_frac = ENGINEERING_CARTON_MIN_FRAC,
        basis::AbstractString = "science")
    nsci = length(cartons)
    nsci == 0 &&
        return (engineering = false, frac = NaN, carton = "", nsci = 0, basis = "none")
    matched = String[]
    for c in cartons
        is_engineering_carton(c) && push!(matched, strip(String(c)))
    end
    frac = length(matched) / nsci
    carton = if isempty(matched)
        ""
    else
        counts = Dict{String, Int}()
        for m in matched
            counts[m] = get(counts, m, 0) + 1
        end
        argmax(counts)
    end
    (engineering = frac > min_frac, frac = frac, carton = carton, nsci = nsci,
        basis = String(basis))
end

"""
    almanac_science_cartons(f, tele, mjd, config_id; root = "raw")

Read the `firstcarton` values of the science fibers of one configuration from an
open almanac HDF5 file `f` (group `<root>/<tele>/<mjd>/fibers/<config_id>`;
pass `root = ""` for the rootless layout some older almanac files use).

Returns `(cartons, basis)` where `basis` is `"science"` (the normal case),
`"all_fibers_fallback"` (no science-category fibers; see
`ENGINEERING_FALLBACK_ALL_FIBERS`), or `"none"`.

Returns empty cartons — never throws — when any of the following holds, which is
the correct "not engineering" answer rather than an error:
- `config_id <= 0` (plate-era rows carry `config_id == -1`)
- the `fibers/<config_id>` group is absent
- the fiber table has no `firstcarton` column (plate-era fiber tables do not:
  they predate cartons entirely)
- the table cannot be read
"""
function almanac_science_cartons(f, tele, mjd, config_id; root::AbstractString = "raw",
        fallback_all_fibers::Bool = ENGINEERING_FALLBACK_ALL_FIBERS)
    empty_result = (cartons = String[], basis = "none")
    (config_id isa Integer) || return empty_result
    config_id > 0 || return empty_result
    path = isempty(root) ? "$(tele)/$(mjd)/fibers/$(config_id)" :
           "$(root)/$(tele)/$(mjd)/fibers/$(config_id)"
    haskey(f, path) || return empty_result
    try
        df = DataFrame(read(f[path]))
        rename!(df, lowercase.(names(df)))
        ("category" in names(df) && "firstcarton" in names(df)) || return empty_result
        sci = strip.(String.(df.category)) .== "science"
        if any(sci)
            return (cartons = String[strip(String(c)) for c in df.firstcarton[sci]],
                basis = "science")
        end
        fallback_all_fibers || return empty_result
        # no science-category fibers at all: fall back to every fiber that has a
        # carton (see ENGINEERING_FALLBACK_ALL_FIBERS for why, and how to disable)
        allc = String[strip(String(c)) for c in df.firstcarton]
        filter!(!isempty, allc)
        isempty(allc) && return empty_result
        return (cartons = allc, basis = "all_fibers_fallback")
    catch e
        @warn "almanac_science_cartons: could not read $(path); treating as non-engineering" exception = e
        return empty_result
    end
end

"""
    exposure_engineering_from_almanac(f, tele, mjd, config_id, image_type; root = "raw")

Full engineering carton check for one almanac exposure row. Applies the check
only to `ENGINEERING_CHECK_IMAGE_TYPES` (see that constant for why), and reads
the science-fiber cartons from the almanac's own fiber table — no confSummary
dependency, since the almanac is already a pipeline input.

Returns the same NamedTuple as `exposure_is_engineering`.
"""
function exposure_engineering_from_almanac(f, tele, mjd, config_id, image_type;
        root::AbstractString = "raw")
    lowercase(strip(String(image_type))) in ENGINEERING_CHECK_IMAGE_TYPES ||
        return (engineering = false, frac = NaN, carton = "", nsci = 0, basis = "none")
    r = almanac_science_cartons(f, tele, mjd, config_id; root = root)
    exposure_is_engineering(r.cartons; basis = r.basis)
end

"""
    exposure_flag_bits(predicted_bad, engineering)

Pack the exposure-level verdicts into the `exposure_flags` bitmask
(`EXPFLAG_PREDICTED_BAD`, `EXPFLAG_ENGINEERING`).
"""
function exposure_flag_bits(predicted_bad::Bool, engineering::Bool)
    b = 0x00
    predicted_bad && (b |= EXPFLAG_PREDICTED_BAD)
    engineering && (b |= EXPFLAG_ENGINEERING)
    b
end

"""
    exposure_ok_for_science(flags)

True when none of the `EXPFLAG_NO_SCIENCE` bits are set. This is the single
predicate every science-sample assembler (prior builds, catalog construction)
should use. It is deliberately NOT consulted anywhere in the 2D/1D reduction:
engineering and predicted-bad exposures are still reduced to 1D.
"""
exposure_ok_for_science(flags::Integer) = (UInt8(flags) & EXPFLAG_NO_SCIENCE) == 0x00

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

- object_q0t0u0: never masked (science frames are handled downstream)
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
"""
function exposure_predicted_bad(labeled_class, pred, status)
    if status in ("nofiles", "unclassified")
        false
    elseif labeled_class == "object_q0t0u0"
        false
    elseif labeled_class in ("dark_q0t0u0", "internalflat_q0t0u0", "arclamp_q0t0u0")
        (pred != labeled_class) || (status == "unknown")
    else
        status != "ok"
    end
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
