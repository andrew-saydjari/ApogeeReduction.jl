# Exposure-type classifier: predicts the CONTENT of an exposure from its ar2D
# images and compares it to the (possibly wrong) commanded image_type label.
# Mislabeled calibration exposures (e.g. ThAr/UNe swaps, FPI labeled arclamp)
# can silently poison wavelength calibration, and commanded-object frames with
# no sky light or non-sky content are not science; this check runs right after
# the 2D stage. It never mutates labels and never drops an exposure from the
# reduction: its verdict is the `exposure_flags` metadata below.
#
# Design notes:
# - Since AKS 2026-11-22 the model is the 2D content forest
#   (src/exposureClassifier2D.jl): ten content classes (`CLASSIFIER2D_CLASSES`:
#   sky, twilight, dark, dark_persist, ledflat, quartzflat, domeflat, thar, une,
#   fpi), 173 features. It replaced the v6 forest (11 commanded-label classes,
#   67 features); v6 verdict strings are still readable via
#   `classifier_content_class`.
# - Features come from the ar2D image products only (v2, AKS 2026-11-22: the read
#   counts n_read / n_ndiff of the approved prototype were removed). The one
#   non-image input is the telescope flag, as in v6; see exposureClassifier2D.jl.
# - Commanded labels are image_type + lamp flags (e.g. "arclamp_q0t1u0" = ThAr,
#   "arclamp_q0t0u1" = UNe, "arclamp_q0t0u0" = FPI); `expected_content_class`
#   maps a label to the content the forest should see.
# - Predictions with max class probability < unknown_tau are reported as
#   "unknown" (image doesn't resemble any trained class).
# - The model artifact is trained offline by scripts/classifier/ (see the
#   README there); it is a calibration input configured by path.

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

# The classifier artifact is a CALIBRATION INPUT, configured by path exactly
# like caldir_darks / caldir_flats / gain_read_cal_dir: pipeline.jl's arg table
# holds the default and airflow/dags/ar_common.py sets it for production. There
# is deliberately no model path constant in src/ — a calibration input does not
# belong baked into library source.

## ---------------------------------------------------------------------------
## Exposure-level flag bits (`exposure_class/<tele>/<mjd>/exposure_flags`)
##
## These are EXPOSURE-level advisory bits, distinct from the per-pixel flag bits
## documented in the README's "Current Flag Bits" table. They are metadata for
## downstream consumers (prior builds, cal runlists, catalog construction); the
## reduction to 1D never consults them, so every exposure is still reduced.
##
## SHARED SCHEMA (PR #397 + PR #398). Two producers write into this one byte, so
## the bit numbering has exactly one definition and `exposure_flag_bits` is the
## only packer. Bits 2^0 and 2^1 are fixed; 2^2 is EXPFLAG_NOTRUN. Do not
## renumber any of them — a disagreement here would be silent and unrecoverable.
## A new producer ORs its bit in through `exposure_flag_bits(...; extra = bit)`.
## ---------------------------------------------------------------------------
# 2^0: the image-content classifier says this exposure should not be used
#      (policy: `exposure_predicted_bad`, PR #398)
const EXPFLAG_PREDICTED_BAD = 0x01
# 2^1: engineering exposure — the configuration was assigned an engineering
#      carton, so the frame exists to exercise the hardware, not to do science
#      (policy: `engineering_verdict`, three clauses, PR #397)
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
# 2^2: NO VERDICT WAS FORMED for this exposure by the exposure-type classifier
#      — the check was deliberately disabled, its per-MJD table was missing, or
#      it threw while evaluating this frame (PR #398).
#
#      This bit is why `exposure_flags == 0` is meaningful: zero means JUDGED
#      AND FINE, not "we never looked". "Never looked" is 2^2, a different
#      value, distinguishable from the byte alone.
#
#      MUTUALLY EXCLUSIVE with 2^0: if no verdict was formed there is no verdict
#      to be adverse. `exposure_flag_bits` asserts this; a byte carrying both is
#      a bug, not a state.
const EXPFLAG_NOTRUN = 0x04
# 2^3 onward: UNASSIGNED. Claim one by adding a const here and ORing it in via
#      `exposure_flag_bits(...; extra = bit)` rather than renumbering.
# bits that mean "do not do science with this exposure"; prior builds and any
# other science-sample assembly must exclude these. Reduction must NOT.
#
# EXPFLAG_NOTRUN is deliberately NOT included: an unjudged exposure is not a
# known-bad one, and silently excluding everything we failed to look at would
# turn a monitoring gap into invisible data loss. A caller that wants "only
# frames positively cleared" must test the notrun bit itself, explicitly.
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
almanac fiber categories that count as SCIENCE fibers for the engineering carton
check (`almanac_config_cartons`). `veto_location_boss` / `veto_location_apogee`
are real catalog stars placed at the veto locations of manually built FPS designs
and are reduced as `"sci"` (see `ALMANAC_CATEGORY_TO_FIBER_TYPE`). A design made
almost entirely of them (15 standards, ~35 sky, ~240 veto_location, zero
`science`) is therefore judged by clause 1 over their cartons
(`ops_gaia_brightneighbors`, `ops_2mass_psc_brightneighbors`), not treated as a
science-less configuration.
"""
const ALMANAC_SCIENCE_CATEGORIES = ("science", "veto_location_boss", "veto_location_apogee")

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
- **clause 2** — the configuration has ZERO science fibers (`ALMANAC_SCIENCE_CATEGORIES`) AND
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
The `design_id` an FPS configuration carries when it was built by jaeger WITHOUT
a robostrategy design. This is the SDSS-V convention, not an almanac sentinel —
established 2026-09-09 and worth restating, because the two are easy to confuse
and the conclusion flips depending on which it is:

- almanac's own missing-value sentinel for `design_id` is `-1`, set at
  `almanac/src/almanac/data_models/exposure.py:46`
  (`Field(default=-1, alias="designid")`) and coerced through
  `empty_string_to_int(v, -1)`. A whole-package grep for the literal `999` in
  `almanac/src` returns exactly one hit, unrelated (a `plate = -999` fallback in
  the bad-exposure CSV loader).
- `-999` is present in the RAW FITS HEADERS and in the confSummary files:
  `apo/59765/apR-a-42030010.apz` carries `DESIGNID = -999`, `FIELDID = '-999'`,
  `CONFIGID = 6106`, `OBSCMNT = 'FPI test'`, and `confSummary-6106.par` has
  `design_id -999`, `field_id -999`, `robostrategy_run NA`, `raCen/decCen -999.0`.

So there is no parser bug to fix here: `-999` is genuine observatory output
meaning "no design was ever generated for this configuration".

`design_id == -999` and `field_id == -999` co-occur perfectly over the DR21
corpus (17,311 rows each), so this is one signal, not two, and `design_id == 0`
never occurs.
"""
const ENGINEERING_DESIGNLESS_DESIGN_ID = -999

"""
Third clause of the engineering rule, AKS 2026-09-09 ("implement the -999
clause"), after a full-corpus census:

> **clause 3** — an `object` exposure whose configuration was built without a
> robostrategy design (`design_id == ENGINEERING_DESIGNLESS_DESIGN_ID`).

It closes a real hole rather than duplicating clauses 1 and 2: those catch
engineering configurations built *from a position-stars design*, and miss
configurations built from *no design at all*.

MEASURED over the full DR21 almanac (`allobs_57618_61230.h5`, 2026-09-09):

- 17,311 rows carry `design_id == -999`, but **17,296 of them (99.91%) are
  CALIBRATION frames** — 7,932 arclamp, 3,540 dark, 3,082 domeflat, 1,767
  quartzflat, 956 internalflat, 19 twilightflat — spread over 1,323 nights,
  continuous from MJD 59637 to 61227.
- exactly **15 are `object`**: apo 59697 (2), apo 59765 (11), apo 60212 (2).

The `ENGINEERING_CHECK_IMAGE_TYPES` restriction is therefore NOT cosmetic and
NOT optional. Implementing this clause anywhere that does not inherit that gate
— a standalone almanac-row filter, say — would destroy 17,296 perfectly good
calibration frames, including most FPS-era APO arclamps and darks. Keep the
clause inside `exposure_engineering_from_almanac`, which gates on image type at
its first line.

All 15 have a POSITIVE `config_id` (5055, 6106, 10680) and identical fiber
tables: 300 fibers, 0 science, 0 sky, 0 standard, 0 non-empty cartons, 298
blank category + 2 `bonus`, `assigned = 0`, `on_target = 0`. A control
configuration on the same night (10681) has 175 science / 91 sky_apogee / 15
standard_apogee, so the emptiness is real and not an ingest failure.

NET NEW: none of the 15 is reachable by clause 1 (needs science fibers) or
clause 2 (needs a `manual_fps_position_stars*` fiber), and they have neither.
The engineering total therefore goes 2,483 -> 2,498 with no reclassification of
anything already flagged.

Why it matters despite being 0.015% of 100,007 FPS-era object exposures: these
are the ONLY `-999` exposures that reach arMADGICS. arM ingests on
`image_type == "object"` (`arMADGICS.jl/src/ingest.jl`), so the 17,296
calibration frames are already invisible to it, while these 15 would be solved
as science against a fiber table with no targets. apo 60212 exp 12/13 are the
clearest case: they sit inside a six-rung twilight ladder (nread 15, 15, **20**,
**25**, 30, 40) whose other four rungs are correctly typed `twilightflat` with
the comment "twilight flats", while these two carry `IMAGETYP = object` and a
blank comment. Without this clause, arM would fit twilight sky continuum with
stellar and sky models. The other 13 announce themselves in
`observer_comment`: 'FPI test', 'test', 'testok', and
'test w/ BOSS FF lamps on for 30s during exp'.

See `metadata/special_cal_obs.txt` for the twilight sequence itself, which is
recorded there as a calibration observation independent of this flag.
"""
const ENGINEERING_FLAG_DESIGNLESS = true

"""
    is_designless_design(design_id)

True when `design_id` marks a configuration built without a robostrategy design
(`ENGINEERING_DESIGNLESS_DESIGN_ID`), and clause 3 is enabled.

Deliberately total: a `nothing`, a missing column, or a non-integer gives
`false` — "we do not know the design" is not the same claim as "there was no
design", and only the latter is evidence of an engineering frame.
"""
function is_designless_design(design_id)
    ENGINEERING_FLAG_DESIGNLESS || return false
    (design_id isa Integer) || return false
    Int(design_id) == ENGINEERING_DESIGNLESS_DESIGN_ID
end

"""
    apply_designless_clause(verdict, design_id)

OR clause 3 into an `engineering_verdict` NamedTuple.

Separate from `engineering_verdict` on purpose. The carton clauses need the
configuration's fiber table, which callers cache per `config_id`; `design_id`
lives on the EXPOSURE row. Folding clause 3 into the cached call would let one
exposure's design leak onto every other exposure sharing its configuration.
Applying it here, per row, after the cache lookup, makes that impossible.

A verdict already flagged by clause 1 or 2 keeps its more specific `basis`; only
an unflagged verdict can become `basis = "designless"`.
"""
function apply_designless_clause(verdict, design_id)
    (verdict.engineering || !is_designless_design(design_id)) && return verdict
    (engineering = true, frac = verdict.frac, carton = verdict.carton,
        nsci = verdict.nsci, basis = "designless")
end

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
the configuration's **science** fibers (category in `ALMANAC_SCIENCE_CATEGORIES`), or of all
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

- `sci_cartons` — `firstcarton` of the science fibers (category in `ALMANAC_SCIENCE_CATEGORIES`)
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

Returns `(science, all)`: the cartons of the science fibers (category in
`ALMANAC_SCIENCE_CATEGORIES`), and the cartons of every fiber (empty carton
strings dropped from `all`, since a blank carton carries no information for
either clause).

A fiber that almanac relabelled `"offtarget"` (assigned but not on target) is
classified by its `original_category` when that column is present: this check
is about what the configuration was DESIGNED to observe, which a positioner
failing to reach its target does not change.

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
        cat = strip.(String.(df.category))
        if "original_category" in names(df)
            orig = strip.(String.(df.original_category))
            cat = ifelse.(cat .== "offtarget", orig, cat)
        end
        sci = in.(cat, Ref(ALMANAC_SCIENCE_CATEGORIES))
        return (science = cart[sci], all = filter(!isempty, cart))
    catch e
        @warn "almanac_config_cartons: could not read $(path); treating as non-engineering" exception = e
        return empty_result
    end
end

"""
    exposure_engineering_from_almanac(f, tele, mjd, config_id, image_type; root = "raw")

Full engineering check for one almanac exposure row. Applies the check only to
`ENGINEERING_CHECK_IMAGE_TYPES` (see that constant for why), and reads the
cartons from the almanac's own fiber table — no confSummary dependency, since
the almanac is already a pipeline input.

Pass `design_id` to include clause 3 (`ENGINEERING_FLAG_DESIGNLESS`). Omitting it
evaluates the carton clauses alone, which is the correct behaviour for a caller
that has no design column rather than a silent "not engineering".

Callers that CACHE this result per `config_id` must leave `design_id` unset here
and apply `apply_designless_clause` per exposure row afterwards — `design_id`
lives on the exposure, not the configuration, so caching it would leak one
exposure's design onto its config-mates.

Returns the `engineering_verdict` NamedTuple.
"""
function exposure_engineering_from_almanac(f, tele, mjd, config_id, image_type;
        design_id = nothing, root::AbstractString = "raw")
    lowercase(strip(String(image_type))) in ENGINEERING_CHECK_IMAGE_TYPES ||
        return (engineering = false, frac = NaN, carton = "", nsci = 0, basis = "none")
    r = almanac_config_cartons(f, tele, mjd, config_id; root = root)
    v = engineering_verdict(r.science, r.all; label = "$(tele)/$(mjd)/config $(config_id)")
    apply_designless_clause(v, design_id)
end

"""
    exposure_flag_bits(predicted_bad, engineering; notrun = false, extra = 0x00)

Pack the exposure-level verdicts into the `exposure_flags` bitmask
(`EXPFLAG_PREDICTED_BAD`, `EXPFLAG_ENGINEERING`, `EXPFLAG_NOTRUN`).

This is the SHARED packer: every producer of `exposure_flags` should go through
it so the bit numbering has exactly one definition. Bits 0 and 1 are fixed and
bit 2 is `EXPFLAG_NOTRUN` — do not renumber any of them.

The positional signature is shared by PR #397 and PR #398; both extra states are
keywords, so a producer that only knows the science verdicts calls it unchanged.
`extra` remains the escape hatch for a FUTURE producer to OR in a bit (2^3
onward) without this function needing to know about it, and without any side
renumbering.

Throws if the composed byte would carry both `EXPFLAG_NOTRUN` and
`EXPFLAG_PREDICTED_BAD`: "no verdict was formed" and "the verdict was adverse"
cannot both be true, and a byte asserting both would be silently misread by
every consumer downstream. The check is on the COMPOSED BYTE, not just on the
`notrun` keyword, so routing the bit in through `extra` cannot evade it.
"""
function exposure_flag_bits(predicted_bad::Bool, engineering::Bool;
        notrun::Bool = false, extra::Integer = 0x00)
    b = UInt8(extra)
    predicted_bad && (b |= EXPFLAG_PREDICTED_BAD)
    engineering && (b |= EXPFLAG_ENGINEERING)
    notrun && (b |= EXPFLAG_NOTRUN)
    (b & EXPFLAG_NOTRUN) != 0x00 && (b & EXPFLAG_PREDICTED_BAD) != 0x00 &&
        throw(ArgumentError("exposure_flags: EXPFLAG_NOTRUN and EXPFLAG_PREDICTED_BAD " *
                            "are mutually exclusive — no verdict cannot also be an adverse verdict"))
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

Build the COMMANDED class label string from almanac metadata, e.g.
("arclamp", false, true, false) -> "arclamp_q0t1u0" (= ThAr). This is the label
being checked; the forest itself predicts CONTENT classes
(`CLASSIFIER2D_CLASSES`), and `expected_content_class` maps one to the other.
Lamp values may be Bool, 0/1, or "T"/"F" strings; anything else maps to "?".
"""
function exposure_class_label(image_type, lamp_quartz, lamp_thar, lamp_une)
    lampstr(x) = (x in (true, 1, "T", "true", "True")) ? "1" :
                 ((x in (false, 0, "F", "false", "False")) ? "0" : "?")
    lowercase(string(image_type)) * "_q" * lampstr(lamp_quartz) *
    "t" * lampstr(lamp_thar) * "u" * lampstr(lamp_une)
end

"""
Commanded label -> the content class the 2D forest should predict for it.
The commanded labels are unchanged; only the forest's vocabulary changed when
the 2D content forest replaced v6 (AKS 2026-11-22). A commanded label not listed
here (rare lamp combinations) has no expected content, so any confident
prediction reads as a mislabel, exactly as v6 never predicted those labels.
"""
const CLASSIFIER_LABEL_CONTENT = Dict(
    "object_q0t0u0" => "sky",
    "dark_q0t0u0" => "dark",
    "arclamp_q0t1u0" => "thar",
    "arclamp_q0t0u1" => "une",
    "arclamp_q0t0u0" => "fpi",
    "quartzflat_q1t0u0" => "quartzflat",
    "domeflat_q0t0u0" => "domeflat",
    "internalflat_q0t0u0" => "ledflat",
    "twilightflat_q0t0u0" => "twilight")

"""
v6 (11 label-string classes) -> 2D content class. Used ONLY to read verdicts
written by the v6 model (old per-MJD exposureTypeCheck tables, old almanac
decorations) under the current policy; the deployed forest predicts the new
names directly. `twilightflat_faint` (night-sky "flats") is sky content.
"""
const CLASSIFIER_LEGACY_PRED = Dict(
    "object_q0t0u0" => "sky", "twilightflat_faint" => "sky",
    "twilightflat_q0t0u0" => "twilight", "dark_q0t0u0" => "dark",
    "dark_persist" => "dark_persist", "internalflat_q0t0u0" => "ledflat",
    "quartzflat_q1t0u0" => "quartzflat", "domeflat_q0t0u0" => "domeflat",
    "arclamp_q0t1u0" => "thar", "arclamp_q0t0u1" => "une", "arclamp_q0t0u0" => "fpi")

"Content class of a prediction string (identity for 2D-forest names; v6 strings translated)."
classifier_content_class(pred) = get(CLASSIFIER_LEGACY_PRED, String(pred), String(pred))
"Expected content class of a commanded label, or `nothing` when the label has none."
expected_content_class(labeled) = get(CLASSIFIER_LABEL_CONTENT, String(labeled), nothing)
"True for a COMMANDED object exposure (any lamp flags): the object-frame rule applies."
is_commanded_object(labeled) = startswith(lowercase(String(labeled)), "object_")

# decision thresholds carried over from v6 unchanged (AKS 2026-11-22: calibration
# masking intent unchanged; no new calibration cuts). The training script writes
# them into the artifact; these are the defaults for an artifact without them.
"max class probability below this -> prediction \"unknown\" (v6 value)"
const CLASSIFIER_UNKNOWN_TAU = 0.5
"confidence above which a content/label disagreement is flagged (v6 value)"
const CLASSIFIER_FLAG_TAU = 0.7

"""
    load_exposure_classifier(model_path, templates_path, light_reference_path)

Load THE exposure classifier: the 2D content forest (`load_exposure_classifier_2d`)
plus the v6-era decision thresholds (`unknown_tau`, `flag_tau`, from the artifact
or `CLASSIFIER_UNKNOWN_TAU` / `CLASSIFIER_FLAG_TAU`). Returns `nothing` (logged:
check DISABLED) when the forest or either stored reference is absent.

A v6 artifact (no `features` list) is rejected with an error: its 11
label-string classes and 67-feature design row are no longer supported.
"""
function load_exposure_classifier(model_path, templates_path, light_reference_path)
    if model_path !== nothing && !isempty(model_path) && isfile(model_path)
        keys_ = JLD2.jldopen(model_path, "r") do f
            collect(keys(f))
        end
        "features" in keys_ ||
            error("$model_path looks like a v6 exposure-classifier artifact (no `features` " *
                  "list). The 2D content forest replaced v6 (AKS 2026-11-22); point " *
                  "--exp_class_model at a classifier_2d artifact " *
                  "(scripts/classifier/train_classifier_2d.jl).")
    end
    c = load_exposure_classifier_2d(model_path, templates_path, light_reference_path)
    c === nothing && return nothing
    d = JLD2.jldopen(model_path, "r") do f
        (haskey(f, "unknown_tau") ? f["unknown_tau"] : CLASSIFIER_UNKNOWN_TAU,
            haskey(f, "flag_tau") ? f["flag_tau"] : CLASSIFIER_FLAG_TAU)
    end
    merge(c, (unknown_tau = Float64(d[1]), flag_tau = Float64(d[2])))
end

"""
    classify_exposure_type(clf, features::AbstractDict)

Score one exposure's 2D feature Dict (`exposure_class_features_2d`).
Returns `(pred, prob, decision, proba, p_sky)` where `pred` is a content class,
`decision` is "unknown" when `prob < clf.unknown_tau`, else "ok".
"""
function classify_exposure_type(clf, features::AbstractDict)
    s = classify_exposure_2d(clf, features)
    decision = s.prob < clf.unknown_tau ? "unknown" : "ok"
    (pred = s.pred, prob = s.prob, decision = decision, proba = s.proba, p_sky = s.p_sky)
end

"""
    exposure_check_category(labeled_class, res, flag_tau = CLASSIFIER_FLAG_TAU)

Interpret a content prediction against the commanded label (the v6 categories,
translated to the 2D content classes; intent unchanged, AKS 2026-11-22):
- "unknown": max probability below `unknown_tau`
- "faint_twilight": labeled twilightflat but the content is night sky
  (v6 predicted `twilightflat_faint`, which is `sky` content now)
- "lamp_off_candidate": labeled as an illuminated exposure but the content is
  `dark` (lamp/LED off: equivalent to having taken a dark)
- "mislabel_candidate": confident content that is not the label's expected
  content (`expected_content_class`), e.g. ThAr/UNe swaps, `dark_persist` on a
  commanded dark
- "ok"
For commanded objects this category is informational only; their mask comes from
`object_frame_rule` (see `exposure_predicted_bad`).
"""
function exposure_check_category(labeled_class, res, flag_tau = CLASSIFIER_FLAG_TAU)
    pred = classifier_content_class(res.pred)
    if res.decision == "unknown"
        "unknown"
    elseif labeled_class == "twilightflat_q0t0u0" && pred == "sky" && res.prob > flag_tau
        "faint_twilight"
    elseif labeled_class in CLASSIFIER_ILLUMINATED && pred == "dark" && res.prob > flag_tau
        "lamp_off_candidate"
    elseif pred != expected_content_class(labeled_class) && res.prob > flag_tau
        "mislabel_candidate"
    else
        "ok"
    end
end

"""
    exposure_predicted_bad(labeled_class, pred, status; p_sky, oh_rel, line_z, tr_dyn)

Masking policy behind `EXPFLAG_PREDICTED_BAD`: should this exposure be excluded?

- COMMANDED OBJECTS (`is_commanded_object`, any lamp flags): bad exactly when the
  approved object-frame rule acts (AKS 2026-11-22, thresholds approved
  2026-11-21/22): `object_frame_rule(1 - p_sky, oh_rel, line_z, tr_dyn).act`, i.e.
  (P(non-sky) >= 0.8 OR (oh_rel < -0.6 AND line_z < 2.2)) AND NOT tr_dyn >= 0.9.
  There is NO object exemption any more. Missing evidence (NaN, e.g. a table
  written by the v6 model) means the rule cannot fire -> not bad.
- commanded dark, internal (LED) flat and FPI (`arclamp_q0t0u0`): masked when the
  content prediction is not the label's expected content (`dark`, `ledflat`,
  `fpi`; so `dark_persist` masks a dark) or the prediction is unknown.
- every other commanded class: masked whenever the check status is not "ok"
  (mislabel / lamp-off / faint-twilight / unknown / rare label).
- no verdict ("nofiles", "unclassified", "notrun", "checkfail"): never masked.

`pred` may be a 2D content class or a legacy v6 string (translated by
`classifier_content_class`). Returns a Bool; "we do not know" is carried by
`EXPFLAG_NOTRUN`, see `exposure_class_metadata` and `exposure_class_verdict`.
"""
function exposure_predicted_bad(labeled_class, pred, status; p_sky = NaN, oh_rel = NaN,
        line_z = NaN, tr_dyn = NaN)
    if status in ("nofiles", "unclassified", EXP_CLASS_STATUS_NOTRUN, "checkfail")
        false
    elseif is_commanded_object(labeled_class)
        object_frame_rule(1 - _rule_float(p_sky), oh_rel, line_z, tr_dyn).act
    elseif labeled_class in ("dark_q0t0u0", "internalflat_q0t0u0", "arclamp_q0t0u0")
        (classifier_content_class(pred) != expected_content_class(labeled_class)) ||
            (status == "unknown")
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
exposure_class_unknown_metadata(status = EXP_CLASS_STATUS_NOTRUN) = merge(
    Dict{String, Any}(
        "exposure_flags" => exposure_flag_bits(false, false; notrun = true),
        "exp_class_status" => String(status),
        "exp_class_pred" => EXP_CLASS_UNKNOWN_STR,
        "exp_class_labeled" => EXP_CLASS_UNKNOWN_STR,
        "exp_class_prob" => NaN),
    exposure_object_rule_metadata(false, NaN, NaN, NaN, NaN))

"""
    exposure_object_rule_metadata(is_object, p_sky, oh_rel, line_z, tr_dyn)

The object-frame rule's provenance fields (written for EVERY exposure, so the
columns exist uniformly in the 1D products and the almanac decoration):
- `exp_class_p_sky`, `exp_class_oh_rel`, `exp_class_line_z`, `exp_class_tr_dyn`:
  the measured inputs (NaN when not measured);
- `exp_class_object_rule`: `object_frame_rule(...).reason` for a commanded object
  ("forest", "no_sky_light", "forest+no_sky_light", "guard_exempt", "none"),
  "n/a" for any other commanded type;
- `exp_class_rule_forest`, `exp_class_rule_no_sky_light`, `exp_class_rule_guard`:
  the three tests individually (false for non-objects).
"""
function exposure_object_rule_metadata(is_object::Bool, p_sky, oh_rel, line_z, tr_dyn)
    p, o, z, t = _rule_float(p_sky), _rule_float(oh_rel), _rule_float(line_z),
    _rule_float(tr_dyn)
    r = is_object ? object_frame_rule(1 - p, o, z, t) : nothing
    Dict{String, Any}(
        "exp_class_p_sky" => p, "exp_class_oh_rel" => o,
        "exp_class_line_z" => z, "exp_class_tr_dyn" => t,
        "exp_class_object_rule" => r === nothing ? "n/a" : r.reason,
        "exp_class_rule_forest" => r === nothing ? false : r.forest,
        "exp_class_rule_no_sky_light" => r === nothing ? false : r.no_sky_light,
        "exp_class_rule_guard" => r === nothing ? false : r.guard)
end

"""
    exposure_class_metadata(labeled, pred, prob, status; engineering = false,
                            p_sky = NaN, oh_rel = NaN, line_z = NaN, tr_dyn = NaN)

Build the metadata block for one exposure from a classifier verdict, including
the `exposure_flags` bitmask. `status` is the category from
`exposure_check_category` (plus "persistence_prior", which `pipeline.jl` adds,
and "checkfail"). `p_sky`, `oh_rel`, `line_z`, `tr_dyn` are the object-frame
rule's inputs; for a commanded object they decide `EXPFLAG_PREDICTED_BAD`
(see `exposure_predicted_bad`) and are recorded with the test that fired
(`exposure_object_rule_metadata`).

A "checkfail" verdict returns the UNKNOWN block — `EXPFLAG_NOTRUN` set, never
`EXPFLAG_PREDICTED_BAD` — and keeps "checkfail" as its status string. An
exception while evaluating the check is a failure to form an opinion, and must
never be recorded as an adverse opinion, nor as a favourable one.

`engineering` is the hand-off point for PR #397's carton-derived bit. This
producer only knows the classifier's verdict, so it passes `false` by default;
the engineering check ORs its bit in here rather than writing a competing field.
"""
function exposure_class_metadata(labeled, pred, prob, status; engineering::Bool = false,
        p_sky = NaN, oh_rel = NaN, line_z = NaN, tr_dyn = NaN)
    (status == "checkfail" || status == EXP_CLASS_STATUS_NOTRUN) &&
        return exposure_class_unknown_metadata(status)
    bad = exposure_predicted_bad(labeled, pred, status; p_sky = p_sky, oh_rel = oh_rel,
        line_z = line_z, tr_dyn = tr_dyn)
    merge(
        Dict{String, Any}(
            "exposure_flags" => exposure_flag_bits(bad, engineering),
            "exp_class_status" => String(status),
            "exp_class_pred" => String(pred),
            "exp_class_labeled" => String(labeled),
            "exp_class_prob" => Float64(prob)),
        exposure_object_rule_metadata(is_commanded_object(labeled), p_sky, oh_rel, line_z,
            tr_dyn))
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

The block is rebuilt from `(labeled, pred, prob, flag)` and, when present, the
object-rule inputs `(p_sky, oh_rel, line_z, tr_dyn)` through
`exposure_class_metadata`, so a table written by any version of the pipeline
yields the current encoding and the masking policy lives in exactly one place.
A table written by the v6 model has no rule inputs: its commanded objects read
back with NaN evidence, which cannot fire the rule (not bad).
"""
function read_exposure_type_check(path)
    out = Dict{Int, Dict{String, Any}}()
    isfile(path) || return out
    try
        d = load(path)
        expnum = d["expnum"]
        labeled, pred = d["labeled"], d["pred"]
        prob, flag = d["prob"], d["flag"]
        n = length(expnum)
        col(k) = haskey(d, k) ? d[k] : fill(NaN, n)
        p_sky, oh_rel, line_z, tr_dyn = col("p_sky"), col("oh_rel"), col("line_z"),
        col("tr_dyn")
        for i in eachindex(expnum)
            out[Int(expnum[i])] = exposure_class_metadata(
                labeled[i], pred[i], prob[i], flag[i]; p_sky = p_sky[i],
                oh_rel = oh_rel[i], line_z = line_z[i], tr_dyn = tr_dyn[i])
        end
    catch e
        @warn "Could not read exposure-type check table $path; treating as unknown" exception=e
        return Dict{Int, Dict{String, Any}}()
    end
    out
end
