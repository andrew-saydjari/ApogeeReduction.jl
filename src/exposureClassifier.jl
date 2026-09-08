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
# 2^1: engineering exposure — the configuration was assigned an engineering
#      carton, so the frame exists to exercise the hardware, not to do science
#      (policy: `engineering_verdict`, two clauses)
#
#      This is an EXPOSURE-level summary of a quantity that is fundamentally
#      per-fiber. Under clause 1 (purity) every science fiber agrees with the
#      exposure, so the summary is exact. Under clause 2 (science-less
#      positioning configurations) it is a deliberate whole-exposure call on a
#      configuration with NO science fibers to serve — 46-93 fibers per such
#      config are neither position-stars nor science and remain UNKNOWN.
#      A per-fiber notion can be added later WITHOUT a schema change: the bit
#      numbering is shared, so a future per-fiber `fiber_flags` array simply
#      reuses 2^1 and the exposure-level bit stays the summary. Do not renumber.
#      `engineering_frac` and `engineering_basis` are written alongside precisely
#      so the bit is never the only record of how the verdict was reached.
const EXPFLAG_ENGINEERING = 0x02
# 2^2 onward: UNASSIGNED here, reserved for the concurrent classifier-propagation
#      work (a NOT-RUN state as a bit rather than by field absence). Pass such a
#      bit through `exposure_flag_bits(...; extra = bit)` rather than renumbering.
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
"Dominated by" is defined as **purity**: ALL of the configuration's *science*
fibers must carry an engineering carton (`frac >= ENGINEERING_CARTON_PURITY`,
i.e. exactly 1.0) before the exposure is flagged engineering.

AKS 2026-09-08, after a full-corpus carton census. Purity costs nothing here —
MEASURED, all four `manual_fps_position_stars*` variants are always 100% of the
science fibers of every configuration they appear in, and they never appear as a
minority in someone else's configuration. Purity is what makes that free
property load-bearing rather than incidental:

- it automatically excludes the 27 other `manual_*` cartons, none of which owns
  a single 100%-pure configuration (median share 0.4-4.2%);
- it excludes `manual_mwm_crosscalib_apogee`, which owns exactly one pure
  configuration but which AKS decided is NOT engineering;
- a majority rule would additionally pull in configurations where crosscalib and
  validation_cool run 52-71%, which AKS does not want.

Zero `ops_*` cartons appear on science fibers anywhere in DR21, so the "except
standards and sky" exemption the ops team suggested is moot.

Behaviour on mixed configurations: a configuration that is *mostly but not
purely* an engineering carton is NOT flagged, and raises a loud warning (see
`ENGINEERING_CARTON_WARN_FRAC`). That has never happened in DR21, so if the
warning ever fires it is a real signal, not noise.
"""
const ENGINEERING_CARTON_PURITY = 1.0

"""
Fraction above which a non-pure configuration raises a loud warning. A
configuration whose engineering-carton share is in
`(ENGINEERING_CARTON_WARN_FRAC, ENGINEERING_CARTON_PURITY)` is NOT flagged
engineering — it fails the purity rule — but it is close enough to the boundary
that somebody should look. No configuration in the DR21 corpus lands here.
"""
const ENGINEERING_CARTON_WARN_FRAC = 0.5

"""
Second clause of the engineering rule, AKS 2026-09-08:

> "if any fibers are `manual_fps_position_stars*` and no science fibers, then
> reject."

i.e. `engineering := CLAUSE 1 OR CLAUSE 2` where

- **clause 1** — the configuration HAS science fibers and *all* of them carry an
  engineering carton (`ENGINEERING_CARTON_PURITY`);
- **clause 2** — the configuration has ZERO `category == "science"` fibers AND
  *any* fiber, of any category, carries an engineering carton.

This is deliberately NOT a general "when science is empty, fall back to all
fibers" rule. The distinction matters: a general fallback would evaluate purity
over whatever cartons a science-less configuration happened to carry and could
flag one that has nothing to do with positioning. Clause 2 is keyed on the
engineering carton itself, so a science-less configuration with different
cartons is untouched.

MEASURED over the full DR21 almanac before this was enabled: clause 2 matches
exactly 5 configurations / 35 object exposures, with ZERO overlap with clause 1,
so it cannot perturb the 2,448 that clause 1 produces (new total 2,483):

| config | position-stars fibers | object exposures |
| --- | --- | --- |
| apo 59558 cfg 105 | 245/300 | 5 |
| apo 59558 cfg 106 | 245/300 | 3 |
| apo 59560 cfg 121 | 207/300 | 1 |
| apo 59560 cfg 122 | 207/300 | 1 |
| apo 59561 cfg 133 | 254/300 | 25 |

All APO, MJD 59558-59561 — a three-night window at the very start of APO FPS
operations (the APO FPS divide is MJD 59423), three nights before
`manual_fps_position_stars` proper begins at 59564. Consistent with the earliest
positioning configurations predating the science-category convention.

Stated plainly rather than papered over: 46-93 fibers per configuration are
neither position-stars nor science, and what they are is UNKNOWN. The rule
deliberately does not care — "any position-stars fiber, no science fibers" — and
that indifference is exactly what makes it safe to apply without understanding
the rest of the configuration.
"""
const ENGINEERING_FLAG_SCIENCELESS_POSITION_STARS = true

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
    exposure_is_engineering(cartons; purity, warn_frac, basis, label)

Carton check for one exposure. `cartons` is the list of `firstcarton` values of
the configuration's **science** fibers (`category == "science"`), or of all
carton-bearing fibers when the config has no science fibers at all and the
configuration has none (clause 2); `basis` records which.

Returns `(engineering, frac, carton, nsci, basis)`:
- `engineering::Bool` — `frac >= purity` (purity rule: ALL science fibers)
- `frac::Float64` — fraction of the fibers in `cartons` with an engineering
  carton (`NaN` when `cartons` is empty, e.g. plate-era or missing config)
- `carton::String` — the most common matching carton name ("" if none), kept
  for provenance so the reason for the bit is auditable
- `nsci::Int` — number of fibers the fraction was computed over
- `basis::String` — `"science"`, `"scienceless_position_stars"`, or `"none"`

Raises a loud `@warn` for a configuration that is MOSTLY but not PURELY an
engineering carton (`warn_frac < frac < purity`). Such a configuration is NOT
flagged. No configuration in the DR21 corpus lands there, so the warning firing
is a real signal. `label` is used only to name the configuration in that
warning.

An empty `cartons` (plate era, no configuration, unreadable fiber table) is
never engineering and never errors.
"""
function exposure_is_engineering(cartons; purity = ENGINEERING_CARTON_PURITY,
        warn_frac = ENGINEERING_CARTON_WARN_FRAC,
        basis::AbstractString = "science", label::AbstractString = "")
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
    engineering = frac >= purity
    if !engineering && frac > warn_frac
        @warn "Engineering carton check: configuration $(label) is MOSTLY but not " *
              "PURELY an engineering carton ($(length(matched))/$(nsci) = " *
              "$(round(frac, digits = 4)) of $(basis) fibers are '$(carton)'). " *
              "It is NOT flagged engineering (the rule requires purity). This has " *
              "never happened in the DR21 corpus — treat it as a real signal and " *
              "decide whether the carton list or the purity rule needs to change."
    end
    (engineering = engineering, frac = frac, carton = carton, nsci = nsci,
        basis = String(basis))
end

"""
    engineering_verdict(sci_cartons, all_cartons; label, flag_scienceless)

Apply BOTH clauses of the engineering rule to one configuration and return the
verdict NamedTuple `(engineering, frac, carton, nsci, basis)`.

- `sci_cartons` — `firstcarton` of the `category == "science"` fibers
- `all_cartons` — `firstcarton` of EVERY fiber in the configuration

Clause 1 (science fibers exist): purity over `sci_cartons`, `basis = "science"`.
Clause 2 (`isempty(sci_cartons)` and any fiber carries an engineering carton):
flagged, `basis = "scienceless_position_stars"`. Clause 2 is keyed on the
engineering carton itself — a science-less configuration whose fibers carry some
OTHER carton is not flagged. See `ENGINEERING_FLAG_SCIENCELESS_POSITION_STARS`.

The two clauses are mutually exclusive by construction (one requires science
fibers, the other requires none), so enabling clause 2 cannot perturb clause 1's
verdicts.

For clause 2, `nsci == 0` (that is what triggered it) and `frac` is the share of
the configuration's CARTON-BEARING fibers that carry an engineering carton
(blank cartons are dropped upstream) — informational only, since the rule there
is "any", not a threshold, and no warning is raised on this path.
"""
function engineering_verdict(sci_cartons, all_cartons; label::AbstractString = "",
        flag_scienceless::Bool = ENGINEERING_FLAG_SCIENCELESS_POSITION_STARS,
        purity = ENGINEERING_CARTON_PURITY,
        warn_frac = ENGINEERING_CARTON_WARN_FRAC)
    none = (engineering = false, frac = NaN, carton = "", nsci = 0, basis = "none")
    # CLAUSE 1: the configuration has science fibers -> purity over those
    if !isempty(sci_cartons)
        return exposure_is_engineering(sci_cartons; purity = purity,
            warn_frac = warn_frac, basis = "science", label = label)
    end
    # CLAUSE 2: zero science fibers AND any fiber carries an engineering carton
    flag_scienceless || return none
    isempty(all_cartons) && return none
    matched = String[strip(String(c)) for c in all_cartons if is_engineering_carton(c)]
    isempty(matched) && return none
    counts = Dict{String, Int}()
    for m in matched
        counts[m] = get(counts, m, 0) + 1
    end
    (engineering = true, frac = length(matched) / length(all_cartons),
        carton = argmax(counts), nsci = 0, basis = "scienceless_position_stars")
end

"""
    almanac_config_cartons(f, tele, mjd, config_id; root = "raw")

Read one configuration's `firstcarton` column from an open almanac HDF5 file `f`
(group `<root>/<tele>/<mjd>/fibers/<config_id>`; pass `root = ""` for the
rootless layout some older almanac files use).

Returns `(science, all)`: the cartons of the `category == "science"` fibers, and
the cartons of every fiber (empty carton strings dropped from `all`, since a
blank carton carries no information for either clause).

Returns both empty — never throws — when any of the following holds, which is
the correct "not engineering" answer rather than an error:
- `config_id <= 0` (plate-era rows carry `config_id == -1`)
- the `fibers/<config_id>` group is absent
- the fiber table has no `firstcarton` column (plate-era fiber tables do not:
  they predate cartons entirely)
- the table cannot be read
"""
function almanac_config_cartons(f, tele, mjd, config_id; root::AbstractString = "raw")
    empty_result = (science = String[], all = String[])
    (config_id isa Integer) || return empty_result
    config_id > 0 || return empty_result
    path = isempty(root) ? "$(tele)/$(mjd)/fibers/$(config_id)" :
           "$(root)/$(tele)/$(mjd)/fibers/$(config_id)"
    haskey(f, path) || return empty_result
    try
        df = DataFrame(read(f[path]))
        rename!(df, lowercase.(names(df)))
        ("category" in names(df) && "firstcarton" in names(df)) || return empty_result
        cart = String[strip(String(c)) for c in df.firstcarton]
        sci = strip.(String.(df.category)) .== "science"
        return (science = cart[sci], all = filter(!isempty, cart))
    catch e
        @warn "almanac_config_cartons: could not read $(path); treating as non-engineering" exception = e
        return empty_result
    end
end

"""
    exposure_engineering_from_almanac(f, tele, mjd, config_id, image_type; root = "raw")

Full engineering carton check for one almanac exposure row. Applies the check
only to `ENGINEERING_CHECK_IMAGE_TYPES` (see that constant for why), and reads
the cartons from the almanac's own fiber table — no confSummary dependency,
since the almanac is already a pipeline input.

Returns the `engineering_verdict` NamedTuple.
"""
function exposure_engineering_from_almanac(f, tele, mjd, config_id, image_type;
        root::AbstractString = "raw")
    lowercase(strip(String(image_type))) in ENGINEERING_CHECK_IMAGE_TYPES ||
        return (engineering = false, frac = NaN, carton = "", nsci = 0, basis = "none")
    r = almanac_config_cartons(f, tele, mjd, config_id; root = root)
    engineering_verdict(r.science, r.all; label = "$(tele)/$(mjd)/config $(config_id)")
end

"""
    exposure_flag_bits(predicted_bad, engineering; extra = 0x00)

Pack the exposure-level verdicts into the `exposure_flags` bitmask
(`EXPFLAG_PREDICTED_BAD`, `EXPFLAG_ENGINEERING`).

This is the SHARED packer: every producer of `exposure_flags` should go through
it so the bit numbering has exactly one definition. Bit 2 onward are unassigned
here and reserved for the concurrent classifier-propagation work (which adds a
NOT-RUN state as a bit rather than by field absence); `extra` lets a caller OR in
such a bit without this function needing to know about it, and without either
side renumbering. Bits 0 and 1 are fixed — do not renumber them.
"""
function exposure_flag_bits(predicted_bad::Bool, engineering::Bool;
        extra::Integer = 0x00)
    b = UInt8(extra)
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
