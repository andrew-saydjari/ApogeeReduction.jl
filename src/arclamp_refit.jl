# Arc-lamp + FPI wavelength REFIT, design stages 3-6 (ARCLAMP_PRODUCTION_DESIGN.md, revision 2,
# with the PI decisions of 2026-11-17 that override it):
#
#   stage 3  arc-anchored comb per night: the comb's integer orders are assigned from the arc
#            solution at the first comb exposure's MEASURED position and its cavity (L, eps) is
#            re-initialised from the arc per night ("ok, let's keep the cavity refit per night"),
#            cavity only, no per-chip zero points; the per-chip arc - comb residual is QA.
#   stage 4  FPI-fiber per-exposure transform (FPS era): the two FPI fibers of every exposure,
#            measured against the night comb, give the exposure's position and drift.
#   stage 5  pooled sky per-exposure transform (no FPI fibers): a robust, line-balanced pooled
#            OH-line offset per chip, anchored to the arcs (never an absolute sky level).
#   plate    era keyed on the plate/FPS MJD boundary (`refit_era`, `get_fps_plate_divide`), never on
#            FPI light; plate-era shape = the arc shape at BOTH telescopes (PI 2026-11-18,
#            `plate_shape`). The comb SHAPE template (file + rebuild script,
#            `scripts/refit/build_comb_shape_template.jl`) serves FPS-era exposures without a comb
#            and stays a switch for the plate era; arc shape is flagged WAVE_SHAPE_ARC.
#   anchor   the comb's arc anchor block is chosen by time to the nearest science exposure,
#            preferring blocks at the science exposures' measured position or closing an object
#            sequence, with gate quality first (`anchor_block = "nearest_science"`; clock-free).
#   stage 6  assembly per exposure: absolute + relative + per-exposure terms, `wavecal_type =
#            <absolute>_<relative>`, `bitmsk_wavecal` (design section 2.4), QA (per-chip arc - comb
#            residual, cross-checks against the sky and against pass 1).
#
# Everything undecided is a field of `RefitConfig` with a documented default, so a variant is a
# rerun of the affected stage only. Inputs are opened read-only. Recorded dither positions are
# never used to place a frame: every position is measured (arc sequences, FPI fibers, sky lines).
# Nothing in the production pipeline calls this file.

using HDF5
using CSV
using DataFrames
using LinearAlgebra
using Statistics: median, mean, quantile

const REFIT_C_MPS = 299792458.0
const REFIT_PXG = 65:20:1985                       # 97-point flat pixel grid (HOLDOUT_RV / AV2)
const REFIT_FORMAT_VERSION = 1

# ------------------------------------------------------------------------------------------
# per-exposure wavelength bitmask (design section 2.4; bits 2^0-2^17 as in the design table)
# ------------------------------------------------------------------------------------------
const WAVE_SKY = 1 << 0
const WAVE_NIGHT_AVE = 1 << 1
const WAVE_ROUGH = 1 << 2
const WAVE_FPI_GATED_OFF = 1 << 3
const WAVE_NIGHTAVE_POS_MISMATCH = 1 << 4
const WAVE_ABS_ARC_WINDOW = 1 << 5
const WAVE_REL_ARC_TIME = 1 << 6
const WAVE_REL_ARC_NEAREST = 1 << 7
const WAVE_FPIFIB_UNAVAILABLE = 1 << 8
const WAVE_ARC_SEQ_GATED = 1 << 9       # an arc sequence was rejected (#417 arc-only rejection)
const WAVE_ARC_ALL_GATED = 1 << 10      # final-fit sequences fail the absolute ceilings (ARCNIGHT_ABS_FAIL)
const WAVE_ARC_UNGATED = 1 << 11        # too few sequences to compare: rejection not possible (ARCNIGHT_FEW_SEQ)
const WAVE_ARC_MULTIPOS = 1 << 12
const WAVE_COMB_MODE_SLIP = 1 << 13
const WAVE_LEVEL_INCONSISTENT = 1 << 14
const WAVE_HIGHORDER_EVENT = 1 << 15       # needs the cross-night pass; not set per night
const WAVE_CHIPSCALE_FREE = 1 << 16
const WAVE_FPIFIB_ROTATION = 1 << 17
# refit-only bits (proposed here; not yet in the design table)
const WAVE_PASS1_COPY = 1 << 18            # no new solution possible: the pass-1 adopted solution is copied
const WAVE_SHAPE_ARC = 1 << 19             # no comb shape (no comb, no template file): arc shape used
const WAVE_SKY_ANCHOR_NIGHT = 1 << 20      # pooled-sky tweak anchored on the whole night (no exposure near an arc)
const WAVE_POS_UNMATCHED = 1 << 21         # measured position matched no arc position cluster
const WAVE_FPIFIB_ONE_FIBER = 1 << 22      # only one FPI fiber usable (no across-slit term)
const WAVE_COMB_ARC_INCONSISTENT = 1 << 23 # |arc - comb| per chip above the QA limit (no slip)
const WAVE_ARC_QUARTZ_PATTERN = 1 << 24    # an arc sequence right after a same-dither quartz flat (no
                                           # commanded move) failed the 0.02 px gate (Q7 pattern)
const WAVE_ANCHOR_FALLBACK = 1 << 25       # comb anchor is not the anchor rule's first-ranked arc block (it was bad,
                                           # unverified or had no matching comb exposure); "morning": no gated morning block
const WAVE_ANCHOR_NOT_SAME_BLOCK = 1 << 26 # no comb exposure of the anchoring arc block at its measured position
const WAVE_ANCHOR_NOT_PREFERRED = 1 << 27  # anchor block of a non-preferred kind: "nearest_science": neither at the nearest
                                           # science exposure's measured position nor closing an object sequence (or no
                                           # science exposures); "morning": an evening block
const WAVE_ANCHOR_UNVERIFIED = 1 << 28     # anchoring arc sequence not verified by the 0.02 px gate (ungated, few lines,
                                           # lamp mismatch, or only gate-failed sequences left)
const WAVE_FPS_NO_COMB = 1 << 29           # FPS era (plate/FPS MJD boundary) but no night comb: arc-route fallback
const WAVE_ANCHOR_TREND_OUTLIER = 1 << 30  # verified_pair: the kept anchor's fixed-mode comb velocity is more than
                                           # anchor_trend_ms from the neighbouring nights' strict-pair trend
# in-order nightly pipeline (`--wavecal_order arc_first`, src/wavecal_arcfirst.jl)
const WAVE_FIBER_FILLED = 1 << 31          # >= 1 fiber has no own fit at its route (arc fiber or template fit
                                           # missing): filled across the slit (interpolate_wave_params) or by
                                           # the arc shape / nearest finite fiber
const WAVE_NOT_ARCFIRST = 1 << 32          # arc_first run, but this exposure has no arc-first solution (no usable
                                           # arcs on the night, or the night's arc-first chain failed): the legacy
                                           # chain's product (or its night-average / rough fallback) was resampled
const WAVE_REANCHORED = 1 << 33            # arc_first: the end-of-night check rejected the first comb anchor;
                                           # the comb was re-anchored on the next candidate (conditional 2nd pass)
const WAVE_ANCHOR_CHECK_FAILED = 1 << 34   # arc_first: no candidate passed the end-of-night check (cap reached);
                                           # the candidate with the smallest check metric is kept
const WAVE_ARC_SKY_TIE_REPLACED = 1 << 31  # plate era: the night's arc-sky tie was > plate_tie_ms from the telescope's
                                           # arc-sky ensemble offset (arcs not at the science position); replaced by it

"""Plate-era arc - sky tie reference (pooled sky - arc model zero point, R/G/B m/s) per telescope: the median
over run1's plate nights (2026-11-18, `anchor_verified/thA_dist.jl`, robust sigma 29-39 m/s per chip).
Rebuild from any run with that script."""
const REFIT_PLATE_TIE_REF = Dict("apo" => [-6.0, 156.0, -10.0], "lco" => [10.0, 58.0, 77.0])

const WAVE_BIT_NAMES = [
    "WAVE_SKY", "WAVE_NIGHT_AVE", "WAVE_ROUGH", "WAVE_FPI_GATED_OFF", "WAVE_NIGHTAVE_POS_MISMATCH",
    "WAVE_ABS_ARC_WINDOW", "WAVE_REL_ARC_TIME", "WAVE_REL_ARC_NEAREST", "WAVE_FPIFIB_UNAVAILABLE",
    "WAVE_ARC_SEQ_GATED", "WAVE_ARC_ALL_GATED", "WAVE_ARC_UNGATED", "WAVE_ARC_MULTIPOS",
    "WAVE_COMB_MODE_SLIP", "WAVE_LEVEL_INCONSISTENT", "WAVE_HIGHORDER_EVENT", "WAVE_CHIPSCALE_FREE",
    "WAVE_FPIFIB_ROTATION", "WAVE_PASS1_COPY", "WAVE_SHAPE_ARC", "WAVE_SKY_ANCHOR_NIGHT",
    "WAVE_POS_UNMATCHED", "WAVE_FPIFIB_ONE_FIBER", "WAVE_COMB_ARC_INCONSISTENT",
    "WAVE_ARC_QUARTZ_PATTERN", "WAVE_ANCHOR_FALLBACK", "WAVE_ANCHOR_NOT_SAME_BLOCK",
    "WAVE_ANCHOR_NOT_PREFERRED", "WAVE_ANCHOR_UNVERIFIED", "WAVE_FPS_NO_COMB", "WAVE_ANCHOR_TREND_OUTLIER",
    "WAVE_FIBER_FILLED", "WAVE_NOT_ARCFIRST", "WAVE_REANCHORED", "WAVE_ANCHOR_CHECK_FAILED", "WAVE_ARC_SKY_TIE_REPLACED"]

"""Era of a telescope-night, keyed ONLY on the canonical plate/FPS MJD boundary
(`get_fps_plate_divide`: FPS is mjd > divide), never on whether FPI light exists that night."""
refit_era(tele::AbstractString, mjd::Integer) = mjd > get_fps_plate_divide(String(tele)) ? "fps" : "plate"

"Night-level arc flags of the stage-1 product mapped onto the per-exposure bits 2^9-2^12, 2^16."
function refit_arc_night_bits(night_flags::Integer, seq_flags::AbstractVector{<:Integer})
    b = 0
    # #417 arc-only rejection (gate B removed): a rejected sequence -> 2^9; ARCNIGHT_ABS_FAIL (final-fit
    # sequences fail the absolute ceilings, the former all-gated screen) -> 2^10; ARCNIGHT_FEW_SEQ
    # (< min_seq_compare sequences, rejection not possible) -> 2^11
    any(s -> (s & ARCSEQ_REJECTED) != 0, seq_flags) && (b |= WAVE_ARC_SEQ_GATED)
    (night_flags & ARCNIGHT_ABS_FAIL) != 0 && (b |= WAVE_ARC_ALL_GATED)
    (night_flags & ARCNIGHT_FEW_SEQ) != 0 && (b |= WAVE_ARC_UNGATED)
    (night_flags & ARCNIGHT_MULTIPOS) != 0 && (b |= WAVE_ARC_MULTIPOS)
    (night_flags & ARCNIGHT_CHIPSCALE_FREE) != 0 && (b |= WAVE_CHIPSCALE_FREE)
    return b
end

"Pass-1 `wavecal_type` -> (new-scheme name, bits) for exposures that keep the pass-1 solution."
function refit_pass1_type(wt::AbstractString)
    wt == "waveCalNightfpiDither" && return ("sky_skydith", WAVE_SKY)
    wt == "waveCalNightskyDither" && return ("sky_skydith", WAVE_SKY)
    startswith(wt, "error_night_ave") &&
        return ("sky_nightpos", WAVE_SKY | WAVE_NIGHT_AVE | WAVE_NIGHTAVE_POS_MISMATCH)
    wt == "error_fixed_fallback" && return ("none_rough", WAVE_ROUGH)
    return ("none_rough", WAVE_ROUGH)
end

# ------------------------------------------------------------------------------------------
# configuration: every undecided choice is a switch with a documented default
# ------------------------------------------------------------------------------------------
"""
    RefitConfig

Switches of the refit (defaults = the design's recommendation where it has one):

  * `neon` ("none"): arc line list neon handling, design L4 (PI decision pending). "none" (N0) =
    shape-clean Th/U/Ar (`use_default`, the stage-1 default); "lab" (N1) = plus lab neon
    (`use_with_neon`); "lab_offset" (N3) = lab neon plus a measured per-lamp, per-telescope offset
    from `neon_offset_file` (learned set `neon_offset_set`); "lab_line_offset" (N3L) = lab neon
    plus a learned per-line, per-telescope offset from `neon_line_offset_file` (CSV: tele, set,
    lamp, lab_wave_A, offset_A); a missing file reduces to "lab" (recorded).
  * `plate_shape` ("arc"): relative shape of every plate-era exposure, the SAME at APO and LCO (PI
    2026-11-18: "I don't want to do something different at APO versus LCO ... not use the comb for
    the plate era shape"); "template" (the comb shape template, needs `shape_template_file`) stays
    a non-default switch. The era is keyed on the plate/FPS MJD boundary (`refit_era`), never on
    whether the night has FPI light.
  * `fps_shape_apo` ("template"), `fps_shape_lco` ("template"): relative shape for FPS-era
    exposures that fall back to the arc route (FPS nights without a comb, exposures without FPI
    light); FPS nights without a comb are flagged WAVE_FPS_NO_COMB (ARCFPI_VALIDATION Q6: the comb
    template helps LCO FPS, -26 m/s held-out on R).
  * `anchor_block` ("verified_pair", see below): which arc block anchors the comb cavity. Never
    clock-based (PI 2026-11-18). "nearest_science" (switch): the used arc sequences are ranked by (1) quality: gated
    and clean, then unverified (ungated / few lines / lamp mismatch), then gate-failed (present
    only when every gated sequence failed); (2) preference (`anchor_prefer`): sequences whose
    MEASURED state matches the nearest science exposure's measured position (pass-1 per-exposure
    solution, < `pos_match_px`) or that close an object sequence (the nearest science exposure
    precedes them); (3) time to the nearest science exposure. The first-ranked sequence with a
    same-block comb exposure (|dt| < `anchor_block_h`) anchors; if none has one, the first-ranked
    sequence with any matching comb exposure. Anything but the first-ranked sequence is flagged
    WAVE_ANCHOR_FALLBACK, a non-preferred kind WAVE_ANCHOR_NOT_PREFERRED, an unverified sequence
    WAVE_ANCHOR_UNVERIFIED; a night without science exposures falls back to "any" (FALLBACK +
    NOT_PREFERRED). "morning" (switch): gated sequences after local midnight, nearest comb first
    (WAVE_ANCHOR_FALLBACK + NOT_PREFERRED when there is none); "any" (switch): every used sequence.
  * `anchor_prefer` ("strict"): "strict" = preferred kinds (same measured position, or closing an
    object sequence) rank before all others of the same quality; "none" = time to the nearest
    science exposure only (quality tiers kept).
  * `anchor_block = "verified_pair"` (DEFAULT; PI 2026-11-18: the anchor must be a pair that did not move): the anchor is an (arc
    sequence, full-frame FPI exposure) PAIR between which the instrument did not move, from the
    exposure SEQUENCE (almanac; commanded dither used only as evidence of a commanded move):
    candidates are "strict" pairs (only darks between, no commanded dither change), then, with
    `anchor_adjacency = "strict_then_nomove"`, "nomove" pairs (flats allowed between, no object
    exposure, no commanded change); at most `anchor_max_candidates`, each costing one comb fit.
    Checks (in-order inputs only: raw arc-line and comb-peak positions, timestamps, sequence):
    the sequence's ThAr - UNe line-position offset |.| < `anchor_lamp_tol_px`, and the per-chip
    arc - comb residual pattern after the cavity fit expressed as a pixel shift b (z_c = a + b v_c,
    a rigid shift is absorbed by L, a motion is not) agreeing with the other candidates' within
    `anchor_pattern_tol_px`. Choice among verified pairs by check quality (smallest |b - median b|,
    strict before nomove), never by clock or science proximity. Fallbacks, each flagged: the first
    strict/nomove candidate unverified (WAVE_ANCHOR_UNVERIFIED); no adjacent pair at all: the
    "any" rule (WAVE_ANCHOR_FALLBACK | NOT_PREFERRED | UNVERIFIED).
  * `anchor_trend` ("both"), `anchor_trend_nights` (10), `anchor_trend_ms` (20), `anchor_trend_min_nights`
    (4): verified_pair cross-night continuity (phase 2, `refit_night_anchorcheck`): the fixed-mode
    comb velocity (`refit_comb_mode_velocity`) of the night's chosen pair against the median of the
    neighbouring nights' strict candidates within `anchor_trend_nights` ("both" sides for bulk
    runs; "trailing" = earlier nights only, for a daily run; "off"). More than `anchor_trend_ms` off
    (about 5 robust sigma of strict pairs, 4.1-4.2 m/s): the night's other candidate closest to the
    trend is used if it is within `anchor_trend_ms` (comb refit, kind "...:trend_swap",
    WAVE_ANCHOR_FALLBACK), else the anchor is kept and flagged WAVE_ANCHOR_TREND_OUTLIER.
  * `fps_route` ("fpifib"): per-exposure position of FPS exposures with FPI fibers: "fpifib" = the night comb
    transported by the exposure's two FPI fibers; "sky" (experiment switch) = ignore the FPI fibers and use the
    arc route (arc time model + pooled-sky per-exposure terms per `prefps_time`), i.e. position from the sky.
  * `plate_tie_check` ("replace"), `plate_tie_ms` (300): plate era. The arc - sky tie (pooled sky - arc model
    at the reference exposures near the arcs; chip G) is the arc-sky line-ensemble offset when arcs and science
    sit at the same position: run1 APO median +156 m/s, robust sigma 34. A tie more than `plate_tie_ms` from the
    telescope's reference (`REFIT_PLATE_TIE_REF`) means the arcs were not where the science was (dither mechanism
    landed differently); the reference replaces it (WAVE_ARC_SKY_TIE_REPLACED, 2^31). run1: |tie - median| > 300
    on 49 APO plate nights, 48 of the 51 tail nights (|night shift| > 250 m/s) and none of the normal ones;
    removing the tie term takes the APO plate night-shift rms from 134 to 33 m/s. "off" = run1 behaviour.
  * `anchor_samepos` ("pass1"): source of the "same measured position as the science" input:
    "pass1" = the science exposures' pass-1 per-exposure solutions (`refit_anchor_samepos`; a
    later-stage product, fine for a refit of an existing run, NOT available in an in-order
    pipeline); "none" = not used. The quality tiers (`refit_anchor_quality`) read the stage-1
    gate flags; both are the single swap points if gate B is replaced.
  * `anchor_pairing` ("block"): how the anchoring arc sequence is paired with a comb exposure:
    "block" = the time-adjacent comb exposure of the same calibration block (|dt| <
    `anchor_block_h`) whose MEASURED position matches the arc sequence's cluster, the arc
    sequence's own measured state seeding the comb; without one, the nearest-in-time matching comb
    exposure (flagged WAVE_ANCHOR_NOT_SAME_BLOCK); "position" = the first comb exposure and the arc
    time model of the cluster that matches it by measured position (the earlier behaviour).
  * `chip_transform_file`: comb chip-scale calibration (PR #415 format); "" = free chip scale (flagged).
  * `exp_terms` ("chipzp_stretch_tilt"): per-exposure position transform from the FPI fibers, design
    D2/D3 (OPEN-UNDER-TEST): "shift" (one pixel shift, all chips), "chipzp" (per-chip zero point),
    "chipzp_stretch" (+ per-chip stretch; design D2'), "chipzp_stretch_tilt" (+ linear in fiber
    through the two FPI fibers; design D3).
  * `sky_fiber_term` ("full"): per-fiber sky term, design D6 (OPEN-UNDER-TEST, ARCFPI_VALIDATION Q5
    pending): "full" (per-fiber per-exposure sky zero point + stretch, fiber median removed, i.e.
    fiber-relative only; the design's D6 recommendation, pass 1's per-fiber dporder-1 behaviour
    without its absolute level); "static" (per-fiber night-median sky offset, fiber median
    removed); "pooled" (fiber-common per-exposure sky tweak in every era, arc-anchored); "none".
    Smoke test (7 nights with stable-star rows, first-order RV): dropping the per-fiber term
    raises the row-level repeat scatter by +61 m/s at APO FPS ("none") vs +4 ("full") / +1
    ("static"); a small sample, so the choice stays a switch.
  * `prefps_time` ("regularised"): time dependence where no FPI fibers are usable (plate era, FPS
    exposures without FPI light), ARCFPI_VALIDATION Q4: "regularised" = option (c), the arc time
    model of the night's REFERENCE position cluster plus `reg_weight_<tele>` x the pooled-sky
    departure (APO 0.825, LCO 0.65, chosen on held-out nights); "sky" = option (a), a fixed arc
    reference state plus the full pooled-sky departure; "arc" = option (b), the arc time model of
    the exposure's own cluster, no sky. In (a) and (c) the step BETWEEN position clusters comes
    from the sky (median sky departure of the cluster's exposures), never from the arcs (the
    arc-measured step is biased on APO chip G: +59..+72 m/s zero point, -127..-133 m/s stretch).
  * `sky_exp_terms` ("chipzp_stretch_tilt"): terms of the pooled-sky per-exposure transform,
    per chip: "chipzp", "chipzp_stretch", "chipzp_stretch_tilt" (the APO G/B tilt moves in the
    plate era too: 93/45 m/s within a position against noise of 8).
  * `reg_weight_apo` (0.825), `reg_weight_lco` (0.65): weights of option (c).
  * `shape_template_file`, `shape_template_block` ("all"): comb shape template for arc-route
    exposures whose shape switch is "template"; "" = arc shape (flagged WAVE_SHAPE_ARC).
  * `arcwin_nights` (15): neighbour-night window for the comb cavity when a night has no arc anchor.
  * `nsig` (8): arc and FPI-fiber peak detection threshold (stage-1 finder).
  * `sky_anchor_window_h` (1.5): exposures within this time of a used arc sequence anchor the
    pooled sky tweak to the arcs.
  * `pos_match_px` (0.25): position-cluster match radius (stage-1 `group_px`).
  * `xcheck_threshold_ms` (60): `WAVE_LEVEL_INCONSISTENT` threshold on the sky cross-check.
  * `comb_arc_qa_ms` (150): `WAVE_COMB_ARC_INCONSISTENT` threshold on |arc - comb| per chip.
  * `coeff_degree` (9): Legendre degree of the stored per-fiber, per-chip solution.
"""
Base.@kwdef struct RefitConfig
    neon::String = "none"
    neon_offset_file::String = ""
    neon_offset_set::String = "all"
    neon_line_offset_file::String = ""
    chip_transform_file::String = ""
    exp_terms::String = "chipzp_stretch_tilt"
    sky_fiber_term::String = "full"
    prefps_time::String = "regularised"
    sky_exp_terms::String = "chipzp_stretch_tilt"
    reg_weight_apo::Float64 = 0.825
    reg_weight_lco::Float64 = 0.65
    shape_template_file::String = ""
    shape_template_block::String = "all"
    plate_shape::String = "arc"
    fps_shape_apo::String = "template"
    fps_shape_lco::String = "template"
    anchor_block::String = "verified_pair"
    anchor_prefer::String = "strict"
    anchor_samepos::String = "pass1"
    anchor_adjacency::String = "strict_then_nomove"
    anchor_max_candidates::Int = 3
    anchor_pattern_tol_px::Float64 = 0.02
    anchor_lamp_tol_px::Float64 = 0.1
    anchor_trend::String = "both"
    anchor_trend_nights::Int = 10
    anchor_trend_ms::Float64 = 20.0
    anchor_trend_min_nights::Int = 4
    fps_route::String = "fpifib"
    plate_tie_check::String = "replace"
    plate_tie_ms::Float64 = 300.0
    anchor_pairing::String = "block"
    anchor_block_h::Float64 = 1.5
    arcwin_nights::Int = 15
    nsig::Float64 = 8.0
    sky_anchor_window_h::Float64 = 1.5
    pos_match_px::Float64 = 0.25
    xcheck_threshold_ms::Float64 = 60.0
    comb_arc_qa_ms::Float64 = 150.0
    coeff_degree::Int = 9
end

const REFIT_SWITCH_VALUES = Dict(
    :neon => ("none", "lab", "lab_offset", "lab_line_offset"),
    :plate_shape => ("arc", "template"),
    :fps_shape_apo => ("template", "arc"), :fps_shape_lco => ("template", "arc"),
    :anchor_block => ("nearest_science", "morning", "any", "verified_pair"),
    :anchor_adjacency => ("strict_then_nomove", "strict"), :anchor_trend => ("both", "trailing", "off"),
    :fps_route => ("fpifib", "sky"), :plate_tie_check => ("replace", "off"), :anchor_prefer => ("strict", "none"),
    :anchor_samepos => ("pass1", "none"),
    :anchor_pairing => ("block", "position"),
    :exp_terms => ("shift", "chipzp", "chipzp_stretch", "chipzp_stretch_tilt"),
    :sky_fiber_term => ("none", "pooled", "static", "full"),
    :prefps_time => ("regularised", "sky", "arc"),
    :sky_exp_terms => ("chipzp", "chipzp_stretch", "chipzp_stretch_tilt"))

"Check the switch values; throws an ArgumentError naming the bad field."
function refit_check_config(cfg::RefitConfig)
    for (k, vals) in REFIT_SWITCH_VALUES
        v = getfield(cfg, k)
        v in vals || throw(ArgumentError("RefitConfig.$(k) = \"$(v)\" not in $(vals)"))
    end
    cfg.coeff_degree >= 4 || throw(ArgumentError("coeff_degree must be >= 4"))
    return true
end

"All switch values as strings (recorded in every product)."
refit_config_dict(cfg::RefitConfig) =
    Dict{String, String}(string(k) => string(getfield(cfg, k)) for k in fieldnames(RefitConfig))

"Which stages a changed switch invalidates (for the plan; the driver uses the same table)."
const REFIT_SWITCH_STAGES = Dict(
    "neon" => "arc,comb,fpifib,assemble", "neon_offset_file" => "arc,comb,fpifib,assemble",
    "neon_offset_set" => "arc,comb,fpifib,assemble", "chip_transform_file" => "arc,comb,fpifib,assemble",
    "neon_line_offset_file" => "arc,comb,fpifib,assemble", "anchor_block" => "comb,fpifib,assemble",
    "anchor_prefer" => "comb,fpifib,assemble", "anchor_samepos" => "comb,fpifib,assemble",
    "anchor_adjacency" => "comb,fpifib,assemble", "anchor_max_candidates" => "comb,fpifib,assemble",
    "anchor_pattern_tol_px" => "comb,fpifib,assemble", "anchor_lamp_tol_px" => "comb,fpifib,assemble",
    "anchor_trend" => "anchorcheck,fpifib,assemble", "anchor_trend_nights" => "anchorcheck,fpifib,assemble",
    "anchor_trend_ms" => "anchorcheck,fpifib,assemble", "anchor_trend_min_nights" => "anchorcheck,fpifib,assemble",
    "fps_route" => "assemble", "plate_tie_check" => "assemble", "plate_tie_ms" => "assemble",
    "anchor_pairing" => "comb,fpifib,assemble", "anchor_block_h" => "comb,fpifib,assemble",
    "plate_shape" => "assemble", "fps_shape_apo" => "assemble", "fps_shape_lco" => "assemble",
    "nsig" => "arc,comb,fpifib,assemble", "arcwin_nights" => "comb,fpifib,assemble",
    "pos_match_px" => "comb,fpifib,assemble", "comb_arc_qa_ms" => "comb,fpifib,assemble",
    "exp_terms" => "assemble", "sky_fiber_term" => "assemble", "prefps_time" => "assemble",
    "sky_exp_terms" => "assemble", "reg_weight_apo" => "assemble", "reg_weight_lco" => "assemble",
    "shape_template_file" => "assemble", "shape_template_block" => "assemble",
    "sky_anchor_window_h" => "assemble", "xcheck_threshold_ms" => "assemble",
    "coeff_degree" => "assemble")

"""Value of each switch added after the first archive run that reproduces the behaviour before it existed
(used by `--reuse_from` when the base run's config.txt lacks the switch; a switch not listed here, or set to
another value, invalidates its stages)."""
const REFIT_SWITCH_LEGACY = Dict("anchor_prefer" => "strict", "anchor_samepos" => "pass1",
    "anchor_adjacency" => "strict_then_nomove", "anchor_max_candidates" => "3", "anchor_pattern_tol_px" => "0.02",
    "anchor_lamp_tol_px" => "0.1", "anchor_trend" => "both", "anchor_trend_nights" => "10", "anchor_trend_ms" => "20.0",
    "anchor_trend_min_nights" => "4", "fps_route" => "fpifib", "plate_tie_check" => "off", "plate_tie_ms" => "300.0")

"FNV-1a fingerprint of a file, or \"none\" for an empty path / missing file."
refit_fingerprint(path::AbstractString) =
    (isempty(path) || !isfile(path)) ? "none" : arclamp_file_fingerprint(path)

# ------------------------------------------------------------------------------------------
# calibration inputs (all optional; the code runs without them and records/flags that)
# ------------------------------------------------------------------------------------------
"""Comb chip scales (N_FIBERS x 2: b_R, b_B) for (tele, mjd) from the PR #415 chip-transform
calibration (`load_chip_transform_cal` + `chip_scales`), used in BOTH eras (PI 2026-11-18: the plate
era keeps the comb chip scale; only the higher-order shape template is not used there). `""` =
free chip scale (`nothing`, flagged downstream); a configured file that is missing or unreadable is
an error (as in #415), and an MJD without a block gives `nothing`."""
refit_chip_scales(path::AbstractString, tele::AbstractString, mjd::Integer) =
    isempty(path) ? nothing : chip_scales(load_chip_transform_cal(path), String(tele), Int(mjd))

"""Per-lamp neon offsets (Dict "ThAr"/"UNe" => Å) for a telescope from the neon calibration file
(`learned/<tele>/<set>/lamp_offset_A`, ARCFPI_VALIDATION build_neon_offsets.jl format), or `nothing`."""
function refit_neon_offsets(path::AbstractString, tele::AbstractString, set::AbstractString)
    (isempty(path) || !isfile(path)) && return nothing
    h5open(path, "r") do f
        k = "learned/$(tele)/$(set)/lamp_offset_A"
        haskey(f, k) || return nothing
        v = read(f[k])
        return Dict("ThAr" => Float64(v[1]), "UNe" => Float64(v[2]))
    end
end

"""Per-line neon offsets (Dict (lamp, lab_wave) => Å) for a telescope and learned set from the
N3L table (`neon_line_offsets_v1.csv`), or `nothing`."""
function refit_neon_line_offsets(path::AbstractString, tele::AbstractString, set::AbstractString)
    (isempty(path) || !isfile(path)) && return nothing
    T = CSV.read(path, DataFrame)
    T = T[(T.tele .== tele) .& (T.set .== set) .& isfinite.(T.offset_A), :]
    nrow(T) == 0 && return nothing
    return Dict((String(r.lamp), Float64(r.lab_wave_A)) => Float64(r.offset_A) for r in eachrow(T))
end

"The arc line list for a configuration and telescope; returns (list, neon_mode_used)."
function refit_linelist(cfg::RefitConfig, tele::AbstractString)
    cfg.neon == "none" && return load_arclamp_linelist(selection = "use_default"), "none"
    if cfg.neon == "lab_offset"
        off = refit_neon_offsets(cfg.neon_offset_file, tele, cfg.neon_offset_set)
        off === nothing && return load_arclamp_linelist(selection = "use_with_neon"), "lab(no offset file)"
        return load_arclamp_linelist(selection = "use_with_neon", neon_offset_A = off), "lab_offset"
    end
    ll = load_arclamp_linelist(selection = "use_with_neon")
    if cfg.neon == "lab_line_offset"
        off = refit_neon_line_offsets(cfg.neon_line_offset_file, tele, cfg.neon_offset_set)
        off === nothing && return ll, "lab(no line offset file)"
        lw = copy(ll.labwave)
        nmatch = 0
        df = ll.df
        for i in 1:nrow(df)
            (Bool(df.neon[i]) && df.lab_id[i] > 0) || continue
            lam = Float64(df.lab_wave[i])
            k = findfirst(key -> key[1] == String(df.lamp[i]) && abs(key[2] - lam) < 0.005, collect(keys(off)))
            k === nothing && continue
            lw[df.lab_id[i]] = lam + off[collect(keys(off))[k]]
            nmatch += 1
        end
        sel = ll.selection * "+neon_line_offset($(tele),$(cfg.neon_offset_set),n=$(nmatch))"
        return ArcLineList(ll.path, ll.version, ll.fingerprint, sel, ll.df, ll.slices, ll.seed, lw,
            ll.lamp_of_id, ll.fit_ids), "lab_line_offset(n=$(nmatch))"
    end
    return ll, "lab"
end

# ------------------------------------------------------------------------------------------
# stage-1 arc night product: loader and arc states
# ------------------------------------------------------------------------------------------
"""Stored stage-1 arc night fit, with the final-fit layout rebuilt. Sequence arrays with suffix
`_u` are over the sequences used in the final fit (in almanac order)."""
struct RefitArc
    fname::String
    TH::Matrix{Float64}
    okf::Vector{Bool}
    g::Vector{Float64}
    lay::ArcGlobalLayout
    used::Vector{Int}
    grp_u::Vector{Int}
    tau_u::Vector{Float64}
    tseq_u::Vector{Float64}
    clus_u::Vector{Int}
    tmean::Float64
    night_flags::Int
    seq_flags::Vector{Int}
    chip_scale_source::String
    lines::NamedTuple
end

function load_refit_arc(fname::AbstractString)
    h5open(fname, "r") do f
        TH = read(f["fiber_params"])
        okf = Bool.(read(f["fiber_ok"]))
        g = Float64.(read(f["global_params"]))
        used = findall(Bool.(read(f["seq_used"])))
        has2 = Bool(read(f["has_pos2"]))
        grp = Int.(read(f["seq_pos_group"]))
        tau = Float64.(read(f["seq_tau"]))
        grp_u = has2 ? grp[used] : ones(Int, length(used))
        lay = ArcGlobalLayout(length(used), grp_u, tau[used])
        lay.n == length(g) || error("$(fname): stored global terms do not match the rebuilt layout")
        css = haskey(f, "metadata/chip_scale_source") ? String(read(f["metadata/chip_scale_source"])) : "unknown"
        lines = (seq = Int.(read(f["lines_seq"])), fiber = Int.(read(f["lines_fiber"])),
            chip = Int.(read(f["lines_chip"])), id = Int.(read(f["lines_id"])),
            pix = Float64.(read(f["lines_pix"])), pixvar = Float64.(read(f["lines_pixvar"])),
            keep = Bool.(read(f["lines_keep"])))
        RefitArc(String(fname), TH, okf, g, lay, used, grp_u, tau[used],
            Float64.(read(f["seq_mjd_mid"]))[used], Int.(read(f["seq_pos_cluster"]))[used],
            Float64(read(f["tau_mjd_mean"])), Int(read(f["night_flags"])),
            Int.(read(f["seq_flags"])), css, lines)
    end
end

"""
    refit_arc_fill_fibers(A) -> (A_filled, n_filled)

Fibers without an arc fit (`okf` false or non-finite parameters) get their 9 per-fiber parameters
(a_R, b_R, a_B, b_B, p0..p4) from the existing across-slit regulariser `interpolate_wave_params`
(quadratic in fiber with per-MTP-block steps; the constant term per MTP block), evaluated on the
fitted fibers. Fitted fibers are NOT changed ("fill_only"). Used by the in-order nightly path
(`fiber_fill = "fill_only"`) so every fiber has a solution; the refit keeps NaN (`"none"`).
"""
function refit_arc_fill_fibers(A::RefitArc)
    bad = [!(A.okf[f] && all(isfinite, view(A.TH, f, :))) for f in 1:N_FIBERS]
    nb = count(bad)
    (nb == 0 || nb == N_FIBERS) && return A, 0
    TH = copy(A.TH)
    TH[bad, :] .= NaN
    inl, ilin = interpolate_wave_params(collect(1:N_FIBERS), TH[:, 1:4], TH[:, 5:9]; block_steps = true)
    for f in findall(bad)
        TH[f, 1:4] .= inl[f, :]
        TH[f, 5:9] .= ilin[f, :]
    end
    # an MTP block without any fitted fiber can leave NaN: nearest filled-or-fitted fiber
    good = [all(isfinite, view(TH, f, :)) for f in 1:N_FIBERS]
    for f in findall(.!good)
        g = findall(good)
        TH[f, :] .= TH[g[argmin(abs.(g .- f))], :]
    end
    okf = trues(N_FIBERS)
    return RefitArc(A.fname, TH, Vector{Bool}(okf), A.g, A.lay, A.used, A.grp_u, A.tau_u, A.tseq_u,
        A.clus_u, A.tmean, A.night_flags, A.seq_flags, A.chip_scale_source, A.lines), nb
end

"An arc night fit is usable: at least one sequence in the final fit and >= 100 fitted fibers."
refit_arc_usable(A::RefitArc) = !isempty(A.used) && count(A.okf) >= 100

"Position clusters present among the used sequences."
refit_arc_clusters(A::RefitArc) = sort(unique(A.clus_u))

"Linear interpolation in time of (d, t) between sequences (constant outside their time range)."
function refit_interp_dt(ts::AbstractVector, D::AbstractVector, T::AbstractVector, t::Real)
    o = sortperm(ts)
    ts = ts[o]
    D = D[o]
    T = T[o]
    n = length(ts)
    (n == 1 || t <= ts[1]) && return D[1], T[1]
    t >= ts[n] && return D[n], T[n]
    k = searchsortedlast(ts, t)
    w = (t - ts[k]) / (ts[k + 1] - ts[k])
    return (1 - w) * D[k] + w * D[k + 1], (1 - w) * T[k] + w * T[k + 1]
end

"""Arc state (grpv, d, t, tau) at time `t` (MJD) for position cluster `cl`.
`mode = :time`: (d, t) interpolated in time between the cluster's sequences, time terms at the
exposure time (clamped to the sequence range) = the arc time model. `mode = :nearest`: the state
of the cluster's sequence nearest in time. Returns `nothing` if the cluster has no sequence."""
function refit_arc_state(A::RefitArc, t::Real, cl::Integer; mode::Symbol = :time)
    js = [j for j in eachindex(A.used) if A.clus_u[j] == cl]
    isempty(js) && return nothing
    if mode == :nearest || !isfinite(t)
        j = isfinite(t) ? js[argmin(abs.(A.tseq_u[js] .- t))] : js[1]
        d, tt = arcfit_seq_state(A.g, A.lay, j)
        return (A.grp_u[j], d, tt, A.tau_u[j])
    end
    Ds = [arcfit_seq_state(A.g, A.lay, j)[1] for j in js]
    Ts = [arcfit_seq_state(A.g, A.lay, j)[2] for j in js]
    d, tt = refit_interp_dt(A.tseq_u[js], Ds, Ts, t)
    lo, hi = extrema(A.tau_u)
    return (A.grp_u[js[1]], d, tt, clamp((t - A.tmean) * 2.0, lo, hi))
end

"Shift of the arc state on chip c of fiber f as A + B x (x = (pix-1024)/2048), in z units."
function refit_arc_shift_AB(A::RefitArc, f::Integer, c::Integer, st)
    grpv, d, t, tau = st
    lay = A.lay
    a = d + t * arcfit_phi(f)
    b = 0.0
    if lay.has2 && grpv == 2
        ia = lay.idx_alpha[c]
        ib = lay.idx_beta[c]
        ia > 0 && (a += A.g[ia])
        ib > 0 && (b += A.g[ib])
    end
    if lay.timeok
        iu = lay.idx_u[c]
        iv = lay.idx_v[c]
        iu > 0 && (a += A.g[iu] * tau)
        iv > 0 && (b += A.g[iv] * tau)
    end
    return a, b
end

"Arc model wavelengths (px, fibers, chips) at a state on pixel vector `px`."
function refit_arc_grid(A::RefitArc, st, px::AbstractVector{<:Real})
    grpv, d, t, tau = st
    G = fill(NaN, length(px), N_FIBERS, N_CHIPS)
    pxv = collect(Float64, px)
    for f in 1:N_FIBERS
        A.okf[f] || continue
        th = view(A.TH, f, :)
        for c in 1:N_CHIPS
            G[:, f, c] .= arcfit_eval(th, A.g, A.lay, f, c, pxv; grpv = grpv, d = d, t = t, tau = tau)
        end
    end
    return G
end

"""
    refit_arc_lin_nl(A, st) -> (linParams N_FIBERS x 5, nlParams N_FIBERS x 4, n_filled)

The arc solution at state `st` rewritten EXACTLY in the production comb form (quartic in the
green-chip coordinate, chip G the identity, chips R/B linear: `transform_x_chips` with cporder
1), so it can seed `comb_fpi_fit_core`. With xs_c = O_c + S_c x the arc's argument on chip c and
lambda = P(xs), the production form is Q(u) = P(O_G + S_G u) with o_c = (O_c - O_G)/S_G and
s_c = S_c/S_G. Fibers without an arc fit take the nearest good fiber's parameters (seed only).
"""
function refit_arc_lin_nl(A::RefitArc, st)
    lin = fill(NaN, N_FIBERS, 5)
    nl = fill(NaN, N_FIBERS, 4)
    for f in 1:N_FIBERS
        A.okf[f] || continue
        th = A.TH[f, :]
        all(isfinite, th) || continue
        aR, bR = refit_arc_shift_AB(A, f, 1, st)
        aG, bG = refit_arc_shift_AB(A, f, 2, st)
        aB, bB = refit_arc_shift_AB(A, f, 3, st)
        O = (th[1] + aR, aG, th[3] + aB)
        S = (th[2] + bR, 1.0 + bG, th[4] + bB)
        q = Polynomial(th[5:9])(Polynomial([O[2], S[2]])).coeffs
        lin[f, :] .= 0.0
        lin[f, 1:min(5, length(q))] .= q[1:min(5, length(q))]
        nl[f, :] .= ((O[1] - O[2]) / S[2], S[1] / S[2], (O[3] - O[2]) / S[2], S[3] / S[2])
    end
    good = findall(f -> all(isfinite, view(lin, f, :)), 1:N_FIBERS)
    isempty(good) && return lin, nl, N_FIBERS
    nfill = 0
    for f in 1:N_FIBERS
        f in good && continue
        j = good[argmin(abs.(good .- f))]
        lin[f, :] .= lin[j, :]
        nl[f, :] .= nl[j, :]
        nfill += 1
    end
    return lin, nl, nfill
end

"Production-form (cporder 1 or 0) solution of fiber f, chip c at pixels px."
function refit_prod_wave(lin::AbstractMatrix, nl::AbstractMatrix, f::Integer, c::Integer, px)
    x = (px .- (N_XPIX ÷ 2)) ./ N_XPIX
    cp = size(nl, 2) ÷ 2 - 1
    if c == 2
        xt = x
    else
        o = c == 1 ? 0 : cp + 1
        xt = cp == 0 ? nl[f, o + 1] .+ x : nl[f, o + 1] .+ nl[f, o + 2] .* x
    end
    return positional_poly_mat(xt, porder = size(lin, 2) - 1) * lin[f, :]
end

# ------------------------------------------------------------------------------------------
# stage 3: arc-anchored comb
# ------------------------------------------------------------------------------------------
const COMB_ANCHOR_ARC = "arc"
const COMB_ANCHOR_ARCWIN = "arcwin"

"""Per fiber-chip comb - arc difference in units of the local FSR (λ²/2L), median over the pixel
grid, and the flat-mean velocity difference arc - comb (m/s). Returns (dfsr, dv) each N_FIBERS x 3."""
function refit_comb_vs_arc(lin, nl, arcG::AbstractArray, L::Real; px = REFIT_PXG)
    dfsr = fill(NaN, N_FIBERS, N_CHIPS)
    dv = fill(NaN, N_FIBERS, N_CHIPS)
    pxv = collect(Float64, px)
    for f in 1:N_FIBERS, c in 1:N_CHIPS
        a = view(arcG, :, f, c)
        all(isfinite, a) || continue
        lc = refit_prod_wave(lin, nl, f, c, pxv)
        all(isfinite, lc) || continue
        fsr = lc .^ 2 ./ (2L)
        dfsr[f, c] = median((lc .- a) ./ fsr)
        dv[f, c] = mean((a .- lc) ./ lc) * REFIT_C_MPS
    end
    return dfsr, dv
end

"""Mode-slip check C3 (i): fiber-chips whose comb - arc difference is within `tol` FSR of a
NONZERO integer. Returns (count, mask N_FIBERS x 3)."""
function refit_mode_slips(dfsr::AbstractMatrix; tol = 0.2)
    m = falses(size(dfsr))
    for i in eachindex(dfsr)
        x = dfsr[i]
        isfinite(x) || continue
        k = round(x)
        (k != 0 && abs(x - k) < tol) && (m[i] = true)
    end
    return count(m), m
end

"""Integer-order change per fiber-chip between two peak-integer tables of the SAME peak ingest
(median over the fiber-chip's peaks). Returns an Int16 matrix N_FIBERS x 3 (0 = unchanged)."""
function refit_order_change(pk_new::AbstractMatrix, pk_old::AbstractMatrix, chipInt::AbstractMatrix,
        ok::AbstractMatrix)
    out = zeros(Int16, N_FIBERS, N_CHIPS)
    size(pk_new) == size(pk_old) || return fill(Int16(-999), N_FIBERS, N_CHIPS)
    for f in 1:N_FIBERS, c in 1:N_CHIPS
        m = [i for i in axes(pk_new, 1) if chipInt[i, f] == c && ok[i, f]]
        isempty(m) && continue
        out[f, c] = Int16(round(median(pk_new[m, f] .- pk_old[m, f])))
    end
    return out
end

"Local midnight as a fraction of the UT day (APO UT-7, LCO about UT-4): splits evening from morning arcs."
const REFIT_LOCAL_MIDNIGHT_FRAC = Dict("apo" => 0.29, "lco" => 0.18)
"Arc taken after local midnight (the UT day wraps during the night, so compare in local time)."
refit_is_morning(tele, t) = isfinite(t) && (mod(t - REFIT_LOCAL_MIDNIGHT_FRAC[String(tele)] + 0.5, 1.0) - 0.5) > 0

# ---- POSITION REFERENCES (re-homed from #417, which removed them with gate B) ----------------
"""A (pixel window x fibers 100:200) slice of an exposure's chip-G wavelength solution, with its
mid-exposure MJD: the measured-position reference of a comb (or science) exposure for the arc
pairing. The refit builds it from a previous run's per-exposure solution, the in-order pipeline from
its provisional comb fit (`arcfirst_comb_refs_provisional`)."""
struct ArcObjRef
    expnum::Int
    mjd_mid::Float64
    kind::String                     # "fpi" or "sky" (which per-exposure solution)
    px0::Int                         # first pixel of the slice
    fibers::UnitRange{Int}
    wave::Matrix{Float64}            # (pixels, fibers) chip G
end

const ARCGATE_FIBERS = 100:200
const ARCGATE_PX = collect(900.0:25.0:1150.0)

"""Measured offset (px) of arc state (`th`, `g`, ...) from an exposure's chip-G solution slice:
median over fibers 100-200 and pixels 900-1150 of (λ_arc − λ_ref)/(dλ_ref/dpix)."""
function arclamp_offset_to_ref(TH, okf, g, lay, grpv, d, t, tau, ref::ArcObjRef)
    v = Float64[]
    for (jf, f) in enumerate(ref.fibers)
        okf[f] || continue
        la = arcfit_eval(view(TH, f, :), g, lay, f, 2, ARCGATE_PX; grpv = grpv, d = d, t = t,
            tau = tau)
        for (k, p) in enumerate(ARCGATE_PX)
            i = Int(p) - ref.px0 + 1
            (2 <= i <= size(ref.wave, 1) - 1) || continue
            lr = ref.wave[i, jf]
            dl = (ref.wave[i + 1, jf] - ref.wave[i - 1, jf]) / 2
            isfinite(lr) && isfinite(dl) && dl != 0 && push!(v, (la[k] - lr) / dl)
        end
    end
    return isempty(v) ? NaN : median(v)
end

# ---- SEQUENCE-QUALITY ADAPTER ---------------------------------------------------------------
# The only places that interpret the stage-1 per-sequence / per-night arc flags (ARCSEQ_*,
# ARCNIGHT_*): `refit_anchor_quality`, `refit_seq_gated` (here) and `refit_arc_night_bits` (top of
# the file).
"""Sequence verified by the stage-1 arc-only check (#417) for the legacy anchor rules' "verified"
input: clean, i.e. anchor quality tier 0 (not rejected, not a position outlier, compared, enough
lines, lamps agree, within the absolute ceilings)."""
refit_seq_gated(sf::Integer) = refit_anchor_quality(sf) == 0

"""Quality tier of an arc sequence for anchoring, from its stage-1 flags (#417 arc-only rejection):
0 = clean; 1 = few lines / lamp mismatch / not compared (night too small to judge) / absolute ceiling;
2 = rejected (ARCSEQ_REJECTED) or position outlier (never an anchor candidate of verified_pair)."""
function refit_anchor_quality(sf::Integer)
    (sf & (ARCSEQ_REJECTED | ARCSEQ_POS_OUTLIER)) != 0 && return 2
    (sf & (ARCSEQ_FEW_LINES | ARCSEQ_LAMP_MISMATCH | ARCSEQ_NOT_COMPARED | ARCSEQ_ABS_FAIL)) != 0 && return 1
    return 0
end

"""
    refit_anchor_rank(tseq, quality, sci_t, samepos; prefer = "strict") -> NamedTuple

Clock-free ranking of arc sequences for the comb anchor ("nearest_science"). Per sequence j:
`dt_sci_h[j]` = time to the nearest science exposure (h); `trailing[j]` = that exposure precedes
the sequence (the block closes an object sequence); `pref[j] = samepos[j] | trailing[j]` (samepos:
the sequence's measured state matches the nearest science exposure's measured position, supplied
by the caller). `order` sorts by (quality, !pref if `prefer == "strict"`, dt_sci_h, j). Without
science times every dt is Inf and `pref` is false.
"""
function refit_anchor_rank(tseq::AbstractVector{<:Real}, quality::AbstractVector{<:Integer},
        sci_t::AbstractVector{<:Real}, samepos::AbstractVector{Bool}; prefer::AbstractString = "strict")
    n = length(tseq)
    st = filter(isfinite, collect(Float64, sci_t))
    dts = fill(Inf, n)
    trailing = falses(n)
    for j in 1:n
        (isempty(st) || !isfinite(tseq[j])) && continue
        d = st .- tseq[j]
        i = argmin(abs.(d))
        dts[j] = abs(d[i]) * 24
        trailing[j] = d[i] < 0
    end
    pref = (samepos .| trailing) .& isfinite.(dts)
    key(j) = (quality[j], prefer == "strict" ? !pref[j] : false, dts[j], j)
    return (; order = sort(collect(1:n), by = key), dt_sci_h = dts, trailing, pref)
end

"""
    refit_anchor_select(anchor_block, tseq, quality, comb_dt_h, comb_k, sci_t, samepos, gated, morning;
        prefer = "strict", anchor_block_h = 1.5) -> NamedTuple

Pure selection kernel of `refit_choose_anchor` (also used by the rule comparison scripts). Per
used arc sequence j: `quality[j]` (`refit_anchor_quality`), its nearest comb exposure at a matching
MEASURED position (`comb_dt_h[j]` h, `comb_k[j]`, 0 = none), whether it matches the nearest science
exposure's measured position (`samepos[j]`), gate-verified (`gated[j]`) and after local midnight
(`morning[j]`, used by the "morning" switch only). Returns (ok, j, bits, kind, rank, dt_sci_h)
with all anchor bits (FALLBACK, NOT_SAME_BLOCK, NOT_PREFERRED, UNVERIFIED).
"""
function refit_anchor_select(anchor_block::AbstractString, tseq::AbstractVector{<:Real},
        quality::AbstractVector{<:Integer}, cdt::AbstractVector{<:Real}, ck::AbstractVector{<:Integer},
        sci_t::AbstractVector{<:Real}, samepos::AbstractVector{Bool}, gated::AbstractVector{Bool},
        morning::AbstractVector{Bool}; prefer::AbstractString = "strict", anchor_block_h::Real = 1.5)
    nU = length(tseq)
    bits = 0
    kind = String(anchor_block)
    rk = refit_anchor_rank(tseq, quality, sci_t, samepos; prefer)
    fail = (; ok = false, j = 0, bits = 0, kind = "none", rank = 0, dt_sci_h = NaN)
    if anchor_block == "nearest_science" && any(isfinite, sci_t)
        order = rk.order
        r = something(findfirst(j -> ck[j] > 0 && cdt[j] <= anchor_block_h, order),
            findfirst(j -> ck[j] > 0, order), 0)
        r == 0 && return fail
        jj = order[r]
        r > 1 && (bits |= WAVE_ANCHOR_FALLBACK)
        rk.pref[jj] || (bits |= WAVE_ANCHOR_NOT_PREFERRED)
        kind = "nearest_science:" * (samepos[jj] ? "same_position" : rk.trailing[jj] ? "end_of_sequence" : "other")
        rank = r
    else
        # "morning", "any", or "nearest_science" without science exposures: smallest arc - comb |dt|
        cand = collect(1:nU)
        if anchor_block == "nearest_science"
            bits |= WAVE_ANCHOR_FALLBACK | WAVE_ANCHOR_NOT_PREFERRED
            kind = "nearest_science:no_science"
        elseif anchor_block == "morning"
            mo = [j for j in 1:nU if morning[j] && gated[j]]
            isempty(mo) ? (bits |= WAVE_ANCHOR_FALLBACK | WAVE_ANCHOR_NOT_PREFERRED) : (cand = mo)
        end
        pick(js) = isempty(js) ? 0 : js[argmin([(cdt[j], j) for j in js])]
        jj = pick([j for j in cand if ck[j] > 0])
        if jj == 0 && length(cand) < nU
            bits |= WAVE_ANCHOR_FALLBACK | WAVE_ANCHOR_NOT_PREFERRED
            jj = pick([j for j in 1:nU if ck[j] > 0])
        end
        jj == 0 && return fail
        anchor_block == "morning" && (kind = morning[jj] ? "morning" : "evening")
        rank = 1 + count(j -> ck[j] > 0 && (cdt[j], j) < (cdt[jj], jj), 1:nU)
    end
    quality[jj] > 0 && (bits |= WAVE_ANCHOR_UNVERIFIED)
    cdt[jj] > anchor_block_h && (bits |= WAVE_ANCHOR_NOT_SAME_BLOCK)
    return (; ok = true, j = jj, bits, kind, rank, dt_sci_h = rk.dt_sci_h[jj])
end

"""
    refit_anchor_samepos(A, scis, pos_match_px) -> Vector{Bool}

"Same dither as the science" input of the nearest_science rule, isolated here so it can be swapped:
per used arc sequence, whether its own measured state matches the MEASURED position of the
nearest-in-time science exposure (`scis[i].ref`, a slice of that exposure's pass-1 per-exposure
solution) within `pos_match_px`. This reads a LATER-stage product (per-exposure sky/FPI solutions),
available to a refit of an existing run but not to an in-order pipeline (arcs -> comb -> per
exposure); `anchor_samepos = "none"` switches it off (then only "closes an object sequence"
counts as preferred). An in-order replacement would measure positions from raw sky-line / FPI peak
positions.
"""
function refit_anchor_samepos(A::RefitArc, scis, pos_match_px::Real)
    nU = length(A.used)
    samepos = falses(nU)
    st_sci = Float64[s.t for s in scis]
    any(isfinite, st_sci) || return samepos
    for j in 1:nU
        isfinite(A.tseq_u[j]) || continue
        i = argmin([isfinite(t) ? abs(t - A.tseq_u[j]) : Inf for t in st_sci])
        scis[i].ref === nothing && continue
        d, t = arcfit_seq_state(A.g, A.lay, j)
        o = arclamp_offset_to_ref(A.TH, A.okf, A.g, A.lay, A.grp_u[j], d, t, A.tau_u[j], scis[i].ref)
        samepos[j] = isfinite(o) && abs(o) < pos_match_px
    end
    return samepos
end

"""
    refit_choose_anchor(A, seq_flags, combs, tele; scis, anchor_block, anchor_prefer, anchor_pairing,
        anchor_block_h, pos_match_px)

Pick the arc sequence that anchors the comb and the comb exposure it is paired with.
`combs`: per comb exposure (order of the FPI file list) a NamedTuple (ref::ArcObjRef or nothing,
t); `scis`: the same for the night's science (object) exposures with a pass-1 per-exposure
solution (`refit_object_refs_pass1`). Positions are always MEASURED: a pair is allowed only if the arc
sequence's own state matches the comb exposure's measured solution within `pos_match_px`, and
"same position as the science" compares the arc state with the science exposure's measured
solution (never the recorded dither).
  * `anchor_block = "nearest_science"` (default): `refit_anchor_rank`; the first-ranked sequence
    with a same-block comb exposure (|dt| < `anchor_block_h`), else the first-ranked sequence with
    any matching comb exposure (WAVE_ANCHOR_NOT_SAME_BLOCK). Not the first-ranked sequence:
    WAVE_ANCHOR_FALLBACK. No science exposures: "any" + FALLBACK + NOT_PREFERRED.
  * `"morning"`: gated sequences after local midnight if any (else every used sequence,
    WAVE_ANCHOR_FALLBACK + NOT_PREFERRED), the matching pair with the smallest |dt|.
  * `"any"`: every used sequence, the matching pair with the smallest |dt|.
  * `anchor_pairing = "position"`: the first comb exposure and the arc TIME MODEL of the cluster
    that matches it (the earlier behaviour; no anchor block).
Every rule sets WAVE_ANCHOR_UNVERIFIED when the chosen sequence's quality tier is > 0.
Returns (ok, state, k_comb, j_seq, bits, dt_h, offset_px, reason, offsets, kind, rank, dt_sci_h,
quality); `offsets` holds the measured arc - comb offset (px) of every (used sequence, comb
exposure) pair (QA); `kind` names the anchor kind, `rank` the chosen sequence's rank (1 = first).
"""
function refit_choose_anchor(A::RefitArc, seq_flags::AbstractVector{<:Integer}, combs, tele;
        scis = NamedTuple[], anchor_block = "nearest_science", anchor_prefer = "strict",
        anchor_samepos = "pass1", anchor_pairing = "block", anchor_block_h = 1.5, pos_match_px = 0.25,
        exclude_seqs = Int[])
    nU = length(A.used)
    nC = length(combs)
    offs = fill(NaN, nU, nC)
    for j in 1:nU, k in 1:nC
        combs[k].ref === nothing && continue
        d, t = arcfit_seq_state(A.g, A.lay, j)
        offs[j, k] = arclamp_offset_to_ref(A.TH, A.okf, A.g, A.lay, A.grp_u[j], d, t, A.tau_u[j], combs[k].ref)
    end
    if anchor_pairing == "position"
        cls = refit_arc_clusters(A)
        best = (Inf, 0, nothing)
        for cl in cls
            st = refit_arc_state(A, combs[1].t, cl; mode = :time)
            st === nothing && continue
            o = combs[1].ref === nothing ? (length(cls) == 1 ? 0.0 : NaN) :
                arclamp_offset_to_ref(A.TH, A.okf, A.g, A.lay, st[1], st[2], st[3], st[4], combs[1].ref)
            (isfinite(o) && abs(o) < abs(best[1])) && (best = (o, cl, st))
        end
        (best[3] === nothing || abs(best[1]) >= pos_match_px) &&
            return (; ok = false, reason = "no arc position cluster at the frame exposure", offsets = offs)
        return (; ok = true, state = best[3], k_comb = 1, j_seq = 0, bits = 0, dt_h = NaN,
            offset_px = best[1], reason = "", offsets = offs, kind = "position", rank = 0,
            dt_sci_h = NaN, quality = -1)
    end
    q = [refit_anchor_quality(seq_flags[A.used[j]]) for j in 1:nU]
    # nearest matching comb exposure (|dt| h, k) of sequence j, (Inf, 0) if none matches
    function bestcomb(j)
        b = (Inf, 0)
        for k in 1:nC
            o = offs[j, k]
            (isfinite(o) && abs(o) < pos_match_px) || continue
            dt = abs(combs[k].t - A.tseq_u[j]) * 24
            isfinite(dt) || (dt = 1e3)
            dt < b[1] && (b = (dt, k))
        end
        return b
    end
    # `exclude_seqs` (used-sequence indices j): candidates rejected by a later verification (the
    # in-order chain's end-of-night check); they take no comb, so the rule falls to the next one
    bc = [j in exclude_seqs ? (Inf, 0) : bestcomb(j) for j in 1:nU]
    # science exposures: times and, for the nearest one, the measured arc - science offset
    st_sci = Float64[s.t for s in scis]
    dts = fill(Inf, nU)
    samepos = anchor_samepos == "pass1" ? refit_anchor_samepos(A, scis, pos_match_px) : falses(nU)
    S = refit_anchor_select(anchor_block, A.tseq_u, q, [x[1] for x in bc], [x[2] for x in bc], st_sci,
        Vector{Bool}(samepos), [refit_seq_gated(seq_flags[A.used[j]]) for j in 1:nU],
        [refit_is_morning(tele, t) for t in A.tseq_u]; prefer = anchor_prefer, anchor_block_h)
    S.ok || return (; ok = false,
        reason = "no arc sequence matches any comb exposure by measured position", offsets = offs)
    jj, bits, kind, rank = S.j, S.bits, S.kind, S.rank
    dt, k = bc[jj]
    d, t = arcfit_seq_state(A.g, A.lay, jj)
    return (; ok = true, state = (A.grp_u[jj], d, t, A.tau_u[jj]), k_comb = k, j_seq = jj, bits,
        dt_h = dt, offset_px = offs[jj, k], reason = "", offsets = offs, kind, rank,
        dt_sci_h = S.dt_sci_h, quality = q[jj])
end

"""
    refit_anchor_candidates(A, seq_flags, combs, tele; max_n = 10, kw...) -> Vector{NamedTuple}

The ordered list of anchor candidates (arc sequence, paired comb exposure) that the configured rule
(`anchor_block`, `anchor_prefer`, `anchor_pairing`; dispatched in `refit_choose_anchor` /
`refit_anchor_select`, where a new rule such as a verified-pair rule drops in) would take in turn
if each earlier one were rejected. Each element: (j_seq (used-sequence index), seq (almanac
sequence), k_comb, verified = quality tier 0 and same calibration block, bits, kind, rank,
offset_px, dt_h). The in-order chain re-anchors on the next element when its end-of-night check
fails.
"""
function refit_anchor_candidates(A::RefitArc, seq_flags::AbstractVector{<:Integer}, combs, tele;
        max_n::Integer = 10, kw...)
    out = NamedTuple[]
    ex = Int[]
    for _ in 1:max_n
        P = refit_choose_anchor(A, seq_flags, combs, tele; kw..., exclude_seqs = copy(ex))
        (P.ok && P.j_seq > 0) || break
        push!(out, (; j_seq = P.j_seq, seq = A.used[P.j_seq], k_comb = P.k_comb,
            verified = P.quality == 0 && (P.bits & WAVE_ANCHOR_NOT_SAME_BLOCK) == 0, bits = P.bits,
            kind = P.kind, rank = P.rank, offset_px = P.offset_px, dt_h = P.dt_h))
        push!(ex, P.j_seq)
    end
    return out
end

"""Integer-order change per fiber-chip between the new comb and the pass-1 comb solution (both in
production form), from their wavelength difference in units of the local FSR at the same pixels
(robust to a different frame exposure: a dither step is a small fraction of an FSR)."""
function refit_order_change_lambda(lin, nl, p1lin, p1nl, L; px = REFIT_PXG)
    out = zeros(Int16, N_FIBERS, N_CHIPS)
    pxv = collect(Float64, px)
    for f in 1:N_FIBERS, c in 1:N_CHIPS
        a = refit_prod_wave(lin, nl, f, c, pxv)
        b = refit_prod_wave(p1lin, p1nl, f, c, pxv)
        (all(isfinite, a) && all(isfinite, b)) || continue
        out[f, c] = Int16(round(median((a .- b) ./ (a .^ 2 ./ (2L)))))
    end
    return out
end

"""Exposure sequence of a night from the almanac (read-only): sorted exposure numbers, a one-letter
type code per exposure (T ThAr, U UNe, F full-frame FPI = arclamp with no lamp, O object, D dark,
Q quartz flat, I internal flat, d dome flat, x other) and the COMMANDED dither (used only to detect
a commanded move, never as a position). Empty vectors if the night is not in the almanac."""
function refit_exposure_sequence(almanac_file::AbstractString, tele::AbstractString, mjd::Integer)
    (isempty(almanac_file) || !isfile(almanac_file)) && return (Int[], String[], Float64[])
    h5open(almanac_file, "r") do f
        k = "raw/$(tele)/$(mjd)/exposures"
        haskey(f, k) || return (Int[], String[], Float64[])
        g = f[k]
        ex = Int.(read(g["exposure"])); o = sortperm(ex)
        it = String.(read(g["image_type"]))[o]; th = Int.(read(g["lamp_thar"]))[o]
        un = Int.(read(g["lamp_une"]))[o]; d = Float64.(read(g["dithered_pixels"]))[o]
        c = [it[i] == "arclamp" ? (th[i] == 1 ? "T" : un[i] == 1 ? "U" : "F") : it[i] == "object" ? "O" :
             it[i] == "dark" ? "D" : it[i] == "quartzflat" ? "Q" : it[i] == "internalflat" ? "I" :
             it[i] == "domeflat" ? "d" : "x" for i in eachindex(it)]
        (ex[o], c, d)
    end
end

"""Sequence adjacency of an arc sequence (its ThAr and UNe exposure numbers) and an FPI exposure:
"strict" (every exposure between them a dark or another arc-lamp / FPI frame, e.g. the second
ThAr/UNe of a `TTUU` sequence; no commanded dither change from the sequence through the FPI exposure), "nomove" (no object exposure between and no commanded change; flats allowed),
else "other"; "unknown" if an exposure is missing from the sequence."""
function refit_pair_adjacency(seqexp::NTuple{2, <:Integer}, fpiexp::Integer, ex, code, dith)
    ia = findfirst(==(minimum(seqexp)), ex); ib = findfirst(==(maximum(seqexp)), ex)
    ik = findfirst(==(fpiexp), ex)
    (ia === nothing || ib === nothing || ik === nothing) && return "unknown"
    lo, hi = ik > ib ? (ib, ik) : (ik, ia)
    nomove = all(abs(dith[j] - dith[ia]) < 0.01 for j in lo:hi)
    btw = (lo + 1):(hi - 1)
    nomove && all(code[j] in ("D", "T", "U", "F") for j in btw) && return "strict"
    nomove && !any(code[j] == "O" for j in btw) && return "nomove"
    return "other"
end

"""Per-chip arc - comb residual (m/s, after the cavity fit) -> (a, b): z_c = a + b v_c with v_c the
velocity of one pixel on chip c. A rigid velocity shift between arc and comb only changes a (and is
absorbed by the cavity L); a motion of the spectrum between the two exposures by b px gives b."""
function refit_chip_pattern_px(z::AbstractVector{<:Real}, tele::AbstractString)
    v = tele == "apo" ? [4231.0, 5253.0, 6328.0] : [4256.0, 5281.0, 6359.0]
    all(isfinite, z) || return (NaN, NaN)
    ab = hcat(ones(3), v) \ collect(Float64, z)
    return (ab[1], ab[2])
end

"""
    refit_arc_anchored_comb(fpi_fnames, A, combs, tele; seq_flags, pass1_comb = nothing, ...)

Stage 3 for a night with an arc anchor. `fpi_fnames`: the night's FPI peak files (R chip names,
read-only); `combs`: per file (ref, t) for the measured-position match (`refit_choose_anchor`).
The paired comb exposure becomes the comb's frame (first file); the anchoring arc sequence's own
measured state seeds the comb: integer orders from the arc, (L, eps) by least squares to the arc,
then the production comb loop (L, eps refit per night). Returns a NamedTuple (`ok = false` with a
`reason` if it cannot anchor).
"""
function refit_arc_anchored_comb(fpi_fnames::AbstractVector{<:AbstractString}, A::RefitArc,
        combs, tele; seq_flags, pass1_comb = nothing, pos_match_px = 0.25, comb_arc_qa_ms = 150.0,
        scis = NamedTuple[], anchor_block = "nearest_science", anchor_prefer = "strict",
        anchor_samepos = "pass1", anchor_pairing = "block", anchor_block_h = 1.5, n_sigma = 4, max_iter = 2,
        exclude_seqs = Int[],
        adjacency = nothing, lamp_off_px = nothing, anchor_max_candidates = 3, anchor_pattern_tol_px = 0.02,
        anchor_lamp_tol_px = 0.1, anchor_adjacency = "strict_then_nomove", force = nothing, seq_score = nothing)
    if anchor_block == "verified_pair" && adjacency !== nothing
        Q = refit_verified_pair_comb(fpi_fnames, A, tele; adjacency, lamp_off_px, pass1_comb, comb_arc_qa_ms,
            anchor_max_candidates, anchor_pattern_tol_px, anchor_lamp_tol_px, anchor_adjacency, n_sigma, max_iter,
            exclude_seqs, force, seq_flags = [seq_flags[s] for s in A.used],
            seq_score = seq_score === nothing ? nothing : seq_score[A.used])
        Q === nothing || return Q
        force === nothing || return (; ok = false, reason = "forced pair $(force) not refittable")
        # no sequence-adjacent pair at all: the "any" rule, flagged
        Qf = refit_arc_anchored_comb(fpi_fnames, A, combs, tele; seq_flags, pass1_comb, pos_match_px, comb_arc_qa_ms,
            scis, anchor_block = "any", anchor_pairing, anchor_block_h, n_sigma, max_iter, exclude_seqs)
        Qf.ok || return Qf
        return merge(Qf, (; anchor_bits = Qf.anchor_bits | WAVE_ANCHOR_FALLBACK | WAVE_ANCHOR_NOT_PREFERRED |
            WAVE_ANCHOR_UNVERIFIED, anchor_kind = "verified_pair:fallback_any"))
    end
    P = refit_choose_anchor(A, seq_flags, combs, tele; scis, anchor_block, anchor_prefer, anchor_samepos, anchor_pairing,
        anchor_block_h, pos_match_px, exclude_seqs)
    P.ok || return (; ok = false, reason = P.reason, cluster_offsets_px = vec(P.offsets))
    st = P.state
    order = vcat(P.k_comb, setdiff(1:length(fpi_fnames), P.k_comb))
    fn = fpi_fnames[order]
    lin0, nl0, nfill = refit_arc_lin_nl(A, st)
    R = comb_fpi_fit_core(fn, lin0, nl0; cporder = 1, wporder = 4, dporder = 2,
        n_sigma = n_sigma, max_ang_sigma = 0.2, max_iter = max_iter, verbose = false,
        seed_LE_from_initial = true)
    R === nothing && return (; ok = false, reason = "comb fit returned nothing", cluster_offsets_px = vec(P.offsets))
    arcG = refit_arc_grid(A, st, REFIT_PXG)
    return refit_comb_qa(R, arcG; anchor = COMB_ANCHOR_ARC, pass1_comb = pass1_comb,
        comb_arc_qa_ms = comb_arc_qa_ms, extra = (; arc_cluster = P.j_seq == 0 ? 0 : A.clus_u[P.j_seq],
            arc_cluster_offset_px = P.offset_px, cluster_offsets_px = vec(P.offsets),
            arc_state = collect(Float64, st), n_seed_fibers_filled = nfill,
            frame_offset_verified = combs[P.k_comb].ref !== nothing, anchor_bits = P.bits,
            anchor_seq = P.j_seq == 0 ? 0 : A.used[P.j_seq], anchor_dt_h = P.dt_h, fnames_ordered = fn,
            anchor_morning = P.j_seq == 0 ? false : refit_is_morning(tele, A.tseq_u[P.j_seq]),
            anchor_kind = P.kind, anchor_rank = P.rank, anchor_dt_sci_h = P.dt_sci_h,
            anchor_quality = P.quality))
end

"Reference cavity (L, eps) per telescope for the fixed-mode comb velocity (any fixed value works; it only sets the zero)."
const REFIT_MODE_REF = Dict("apo" => (3.7363e7, 0.23), "lco" => (3.73611e7, 0.22))
const REFIT_MODE_M0 = 4670

"""Fixed-mode comb velocity (m/s), independent of the integer-order convention: the continuous mode
number at the fixed wavelength lambda_ref (= the wavelength of mode REFIT_MODE_M0 of the reference cavity),
phi = 2L/lambda_ref - eps, compared with REFIT_MODE_M0, wrapped to +-1/2 order, times c/(M0 + eps0). A rigid
scaling of all comb wavelengths by (1 + v/c) gives v; an L - eps trade that keeps the modes near lambda_ref
fixed gives 0; a change of the integer convention by whole orders gives 0. (c dL/L alone is NOT the comb's
velocity: L and eps trade against each other in a fit.)"""
function refit_comb_mode_velocity(L::Real, eps::Real, tele::AbstractString)
    L0, e0 = REFIT_MODE_REF[String(tele)]
    lref = 2L0 / (REFIT_MODE_M0 + e0)
    d = 2L / lref - eps - REFIT_MODE_M0
    d -= round(d)
    return REFIT_C_MPS * d / (REFIT_MODE_M0 + e0)
end

"""
    refit_verified_pair_comb(fpi_fnames, A, tele; adjacency, lamp_off_px, ...) -> NamedTuple or nothing

`anchor_block = "verified_pair"`: candidates = sequence-adjacent (used arc sequence j, FPI file k)
pairs (`adjacency[j, k]` from `refit_pair_adjacency`), strict first, then nomove, at most
`anchor_max_candidates`; each gets the arc-seeded comb fit with its FPI exposure as the frame, then
the checks: lamp offset (ThAr - UNe, from the arc line positions) and the per-chip residual pattern
as a pixel shift b (`refit_chip_pattern_px`), compared across candidates. `nothing` if the night has
no adjacent pair. The chosen fit carries `anchor_kind` ("verified_pair:strict|nomove" or
"verified_pair:unverified_<class>") and per-candidate QA (`vp_*`).
"""
function refit_verified_pair_comb(fpi_fnames, A::RefitArc, tele; adjacency, lamp_off_px, pass1_comb = nothing,
        comb_arc_qa_ms = 150.0, anchor_max_candidates = 3, anchor_pattern_tol_px = 0.02, anchor_lamp_tol_px = 0.1,
        anchor_adjacency = "strict_then_nomove", n_sigma = 4, max_iter = 2, exclude_seqs = Int[], force = nothing,
        seq_flags = nothing, seq_score = nothing)
    nU, nC = size(adjacency)
    cls_ok = anchor_adjacency == "strict" ? ("strict",) : ("strict", "nomove")
    # candidates: never a rejected sequence or a position outlier, nor one rejected by the in-order chain's
    # end-of-night check (`exclude_seqs`, re-anchoring); within a class, cleaner sequences first
    # (#417 seq_reject_score, < 1 passes, lower is cleaner)
    sf = seq_flags === nothing ? zeros(Int, nU) : seq_flags
    sc = seq_score === nothing ? zeros(nU) : [isfinite(x) ? x : Inf for x in seq_score]
    cand = [(j, k, adjacency[j, k]) for c in cls_ok for j in sortperm(sc) for k in 1:nC
            if adjacency[j, k] == c && refit_anchor_quality(sf[j]) < 2 && !(j in exclude_seqs)]
    isempty(cand) && return nothing
    cand = cand[1:min(end, anchor_max_candidates)]
    # force = (stage-1 sequence index, FPI file basename): refit only that pair (continuity swap)
    force === nothing || (cand = [x for x in cand if A.used[x[1]] == force[1] && basename(fpi_fnames[x[2]]) == force[2]])
    isempty(cand) && return nothing
    fits = Any[]
    for (j, k, c) in cand
        d, t = arcfit_seq_state(A.g, A.lay, j)
        st = (A.grp_u[j], d, t, A.tau_u[j])
        order = vcat(k, setdiff(1:length(fpi_fnames), k))
        fn = fpi_fnames[order]
        lin0, nl0, nfill = refit_arc_lin_nl(A, st)
        R = comb_fpi_fit_core(fn, lin0, nl0; cporder = 1, wporder = 4, dporder = 2,
            n_sigma = n_sigma, max_ang_sigma = 0.2, max_iter = max_iter, verbose = false, seed_LE_from_initial = true)
        R === nothing && continue
        arcG = refit_arc_grid(A, st, REFIT_PXG)
        Q = refit_comb_qa(R, arcG; anchor = COMB_ANCHOR_ARC, pass1_comb = pass1_comb, comb_arc_qa_ms = comb_arc_qa_ms,
            extra = (; arc_cluster = A.clus_u[j], arc_cluster_offset_px = NaN, cluster_offsets_px = Float64[],
                arc_state = collect(Float64, st), n_seed_fibers_filled = nfill, frame_offset_verified = false,
                anchor_seq = A.used[j], anchor_dt_h = NaN, fnames_ordered = fn,
                anchor_morning = refit_is_morning(tele, A.tseq_u[j]), anchor_rank = 0, anchor_dt_sci_h = NaN,
                anchor_quality = -1))
        a, b = refit_chip_pattern_px(Q.chip_arc_minus_comb_ms, tele)
        lo = lamp_off_px === nothing ? NaN : lamp_off_px[j]
        push!(fits, (; j, k, c, Q, a, b, lamp_ok = !isfinite(lo) || abs(lo) < anchor_lamp_tol_px, lamp = lo))
    end
    isempty(fits) && return nothing
    bs = [x.b for x in fits if x.lamp_ok && isfinite(x.b)]
    bref = isempty(bs) ? NaN : median(bs)
    # pattern check: needs >= 2 candidates to compare (a single candidate cannot verify itself)
    verified = [x.lamp_ok && length(bs) >= 2 && isfinite(x.b) && abs(x.b - bref) < anchor_pattern_tol_px for x in fits]
    score(i) = (verified[i] ? 0 : 1, fits[i].c == "strict" ? 0 : 1, isfinite(fits[i].b) && isfinite(bref) ? abs(fits[i].b - bref) : Inf, i)
    ib = sort(collect(eachindex(fits)), by = score)[1]
    x = fits[ib]
    bits = verified[ib] ? 0 : WAVE_ANCHOR_UNVERIFIED
    x.c == "strict" || (bits |= WAVE_ANCHOR_NOT_PREFERRED)
    ib == 1 || (bits |= WAVE_ANCHOR_FALLBACK)
    kind = "verified_pair:" * (verified[ib] ? "" : "unverified_") * x.c
    return merge(x.Q, (; anchor_bits = bits, anchor_kind = kind, anchor_rank = ib,
        vp_seq = [A.used[f.j] for f in fits], vp_fpi = [basename(fpi_fnames[f.k]) for f in fits],
        vp_class = [f.c for f in fits], vp_b_px = [f.b for f in fits], vp_a_ms = [f.a for f in fits],
        vp_lamp_px = [f.lamp for f in fits], vp_verified = verified,
        vp_L = [f.Q.R.cavity_size for f in fits], vp_eps = [f.Q.R.m_offset for f in fits],
        vp_vmode_ms = [refit_comb_mode_velocity(f.Q.R.cavity_size, f.Q.R.m_offset, tele) for f in fits]))
end

"""Stage 3 for a night WITHOUT an arc anchor: integer orders from the pass-1 seed (`lin0`, `nl0`,
the night-average sky solution), cavity (L, eps) FIXED at a neighbouring night's arc-anchored
values (`arcwin`). `arcG` (optional) is an arc grid at the frame for QA only."""
function refit_arcwin_comb(fpi_fnames, lin0, nl0, L::Real, eps::Real; pass1_comb = nothing,
        comb_arc_qa_ms = 150.0, n_sigma = 4, max_iter = 2, neighbour_mjd = -1)
    R = comb_fpi_fit_core(fpi_fnames, lin0, nl0; cporder = 1, wporder = 4, dporder = 2,
        n_sigma = n_sigma, max_ang_sigma = 0.2, max_iter = max_iter, verbose = false,
        cavity_override = L, m_offset_override = eps, fix_LE = true)
    R === nothing && return (; ok = false, reason = "comb fit returned nothing")
    return refit_comb_qa(R, nothing; anchor = COMB_ANCHOR_ARCWIN, pass1_comb = pass1_comb,
        comb_arc_qa_ms = comb_arc_qa_ms, extra = (; arc_cluster = 0, arc_cluster_offset_px = NaN,
            cluster_offsets_px = Float64[], arc_state = Float64[], n_seed_fibers_filled = 0,
            frame_offset_verified = false, neighbour_mjd = neighbour_mjd, anchor_bits = 0,
            anchor_seq = 0, anchor_dt_h = NaN, fnames_ordered = fpi_fnames, anchor_morning = false,
            anchor_kind = "arcwin", anchor_rank = 0, anchor_dt_sci_h = NaN, anchor_quality = -1))
end

"QA of a comb fit result (C3 (i), arc - comb per chip, order changes vs pass 1)."
function refit_comb_qa(R, arcG; anchor, pass1_comb, comb_arc_qa_ms, extra)
    gp = R.good_peaks
    resid = filter(isfinite, R.resid_vec'[gp])
    rms_mA = isempty(resid) ? NaN : 1000 * sqrt(mean(resid .^ 2))
    if arcG !== nothing
        dfsr, dv = refit_comb_vs_arc(R.linParams, R.nlParams, arcG, R.cavity_size)
        nslip, slipmask = refit_mode_slips(dfsr)
        chip_dv = [_arc_nanmedian(dv[:, c]) for c in 1:N_CHIPS]
        chip_dv_rsig = [_arc_robstd(dv[:, c]) for c in 1:N_CHIPS]
        n_far = count(x -> isfinite(x) && abs(x) > 0.3 && abs(x - round(x)) >= 0.2, dfsr)
    else
        dfsr = fill(NaN, N_FIBERS, N_CHIPS)
        dv = fill(NaN, N_FIBERS, N_CHIPS)
        nslip, slipmask = 0, falses(N_FIBERS, N_CHIPS)
        chip_dv = fill(NaN, 3)
        chip_dv_rsig = fill(NaN, 3)
        n_far = 0
    end
    ordch = fill(Int16(-999), N_FIBERS, N_CHIPS)
    p1_L = NaN
    p1_eps = NaN
    if pass1_comb !== nothing
        ordch = refit_order_change_lambda(R.linParams, R.nlParams, pass1_comb.lin, pass1_comb.nl,
            R.cavity_size)
        p1_L = pass1_comb.L
        p1_eps = pass1_comb.eps
    end
    n_ord = count(x -> x != 0 && x != -999, ordch)
    inconsistent = any(x -> isfinite(x) && abs(x) > comb_arc_qa_ms, chip_dv)
    return merge((; ok = true, reason = "", anchor, R, rms_mA, dfsr, dv_arc_minus_comb = dv,
        chip_arc_minus_comb_ms = chip_dv, chip_arc_minus_comb_rsig_ms = chip_dv_rsig, n_slip = nslip,
        slipmask, n_far_fiberchips = n_far, order_change_vs_pass1 = ordch,
        n_fiberchip_order_changed = n_ord, pass1_L = p1_L, pass1_eps = p1_eps,
        comb_arc_inconsistent = inconsistent), extra)
end

"Pass-1 comb (waveCalFPI) integers, use flags and cavity, opened read-only; `nothing` if absent."
function refit_read_pass1_comb(fname::AbstractString)
    isfile(fname) || return nothing
    h5open(fname, "r") do f
        (; peak_ints = read(f["resid_peak_ints"]), used = Bool.(read(f["resid_used_in_fit"])),
            L = Float64(read(f["fpi_cavity_size"])), eps = Float64(read(f["fpi_m0_offset"])),
            lin = read(f["linParams"]), nl = read(f["nlParams"]))
    end
end

"Write the stage-3 night product `waveCalFPIArc_<tele>_<mjd>.h5` (atomic)."
function write_refit_comb(outname::AbstractString, Q, fpi_fnames; tele, mjd, cfg::RefitConfig,
        arc_fname = "", provenance = Dict{String, String}())
    R = Q.R
    meta = Dict{String, Any}("product" => "waveCalFPIArc", "stage" => 3, "tele" => String(tele),
        "mjd" => Int(mjd), "anchor" => Q.anchor, "arc_product" => String(arc_fname),
        "format_version" => REFIT_FORMAT_VERSION)
    for (k, v) in refit_config_dict(cfg)
        meta["cfg_" * k] = v
    end
    for (k, v) in provenance
        meta[k] = v
    end
    tmp = outname * ".tmp$(getpid())"
    safe_jldsave(tmp, meta; linParams = R.linParams, nlParams = R.nlParams,
        ditherParams = R.ditherParams, ditherPolyParams = R.ditherPolyParams,
        exposure_names = [basename(x) for x in Q.fnames_ordered], fpi_cavity_size = R.cavity_size,
        anchor_bits = Q.anchor_bits, anchor_seq = Q.anchor_seq, anchor_dt_h = Q.anchor_dt_h,
        anchor_morning = Q.anchor_morning, anchor_kind = String(Q.anchor_kind), anchor_rank = Q.anchor_rank,
        anchor_dt_sci_h = Q.anchor_dt_sci_h, anchor_quality = Q.anchor_quality,
        fpi_m0_offset = R.m_offset, fpi_m0 = R.med_m0, seed_cavity_size = R.seed_cavity_size,
        seed_m0_offset = R.seed_m_offset, init_cavity_size = R.init_cavity_size,
        init_m0_offset = R.init_m_offset, first_dither_offset = R.first_dither_offset,
        comb_resid_rms_mA = Q.rms_mA, n_good_peaks = count(R.good_peaks),
        dfsr_comb_minus_arc = Q.dfsr, dv_arc_minus_comb_ms = Q.dv_arc_minus_comb,
        chip_arc_minus_comb_ms = Q.chip_arc_minus_comb_ms,
        chip_arc_minus_comb_rsig_ms = Q.chip_arc_minus_comb_rsig_ms, n_mode_slip = Q.n_slip,
        mode_slip_mask = Matrix{Bool}(Q.slipmask), n_far_fiberchips = Q.n_far_fiberchips,
        order_change_vs_pass1 = Q.order_change_vs_pass1,
        n_fiberchip_order_changed = Q.n_fiberchip_order_changed, pass1_cavity_size = Q.pass1_L,
        pass1_m0_offset = Q.pass1_eps, comb_arc_inconsistent = Q.comb_arc_inconsistent,
        arc_cluster = Q.arc_cluster, arc_cluster_offset_px = Q.arc_cluster_offset_px,
        cluster_offsets_px = Float64.(Q.cluster_offsets_px), arc_state = Float64.(Q.arc_state),
        n_seed_fibers_filled = Q.n_seed_fibers_filled,
        frame_offset_verified = Q.frame_offset_verified,
        neighbour_mjd = hasproperty(Q, :neighbour_mjd) ? Q.neighbour_mjd : -1,
        vp_seq = hasproperty(Q, :vp_seq) ? Q.vp_seq : Int[], vp_fpi = hasproperty(Q, :vp_fpi) ? Q.vp_fpi : String[],
        vp_class = hasproperty(Q, :vp_class) ? Q.vp_class : String[],
        vp_b_px = hasproperty(Q, :vp_b_px) ? Q.vp_b_px : Float64[], vp_a_ms = hasproperty(Q, :vp_a_ms) ? Q.vp_a_ms : Float64[],
        vp_lamp_px = hasproperty(Q, :vp_lamp_px) ? Q.vp_lamp_px : Float64[],
        vp_verified = hasproperty(Q, :vp_verified) ? Vector{Bool}(Q.vp_verified) : Bool[],
        vp_L = hasproperty(Q, :vp_L) ? Q.vp_L : Float64[], vp_eps = hasproperty(Q, :vp_eps) ? Q.vp_eps : Float64[],
        vp_vmode_ms = hasproperty(Q, :vp_vmode_ms) ? Q.vp_vmode_ms : Float64[],
        vmode_ms = refit_comb_mode_velocity(R.cavity_size, R.m_offset, String(tele)))
    mv(tmp, outname; force = true)
    return outname
end

"Stored stage-3 comb: (lin, nl, L, eps, anchor, flags...)."
function load_refit_comb(fname::AbstractString)
    h5open(fname, "r") do f
        (; lin = read(f["linParams"]), nl = read(f["nlParams"]), L = Float64(read(f["fpi_cavity_size"])),
            eps = Float64(read(f["fpi_m0_offset"])), anchor = String(read(f["metadata/anchor"])),
            n_slip = Int(read(f["n_mode_slip"])), inconsistent = Bool(read(f["comb_arc_inconsistent"])),
            anchor_bits = haskey(f, "anchor_bits") ? Int(read(f["anchor_bits"])) : 0,
            anchor_kind = haskey(f, "anchor_kind") ? String(read(f["anchor_kind"])) : "unknown",
            anchor_rank = haskey(f, "anchor_rank") ? Int(read(f["anchor_rank"])) : 0,
            vp_class = haskey(f, "vp_class") ? String.(read(f["vp_class"])) : String[],
            vp_seq = haskey(f, "vp_seq") ? Int.(read(f["vp_seq"])) : Int[],
            vp_fpi = haskey(f, "vp_fpi") ? String.(read(f["vp_fpi"])) : String[],
            vp_vmode_ms = haskey(f, "vp_vmode_ms") ? Float64.(read(f["vp_vmode_ms"])) : Float64[],
            chip_arc_minus_comb_ms = Float64.(read(f["chip_arc_minus_comb_ms"])),
            exposure_names = String.(read(f["exposure_names"])))
    end
end

# ------------------------------------------------------------------------------------------
# stage 4: FPI fibers in every exposure (position and drift against the night comb)
# ------------------------------------------------------------------------------------------
"""Weighted line fit d = a + b t (t = (x - 1024)/1024) with iterative clipping (error floor
0.003 px), as in the FPI_FIBER_DITHER study. Returns (a, b, ea, eb, chi2r, n, rms)."""
function refit_fit_ab(x, d, err; clip = 4.0, floor = 0.003, iters = 5, minn = 8)
    t = (x .- 1024.0) ./ 1024.0
    e = sqrt.(err .^ 2 .+ floor^2)
    ok = isfinite.(d) .& isfinite.(e)
    count(ok) < minn && return (NaN, NaN, NaN, NaN, NaN, 0, NaN)
    med = median(d[ok])
    s = 1.4826 * median(abs.(d[ok] .- med))
    s = (isfinite(s) && s > 0) ? s : 0.01
    ok .&= abs.(d .- med) .< max(8 * s, 0.05)
    a = b = ea = eb = chi = NaN
    for _ in 1:iters
        count(ok) < minn && return (NaN, NaN, NaN, NaN, NaN, 0, NaN)
        w = 1 ./ e[ok] .^ 2
        X = hcat(ones(count(ok)), t[ok])
        M = X' * (w .* X)
        p = M \ (X' * (w .* d[ok]))
        C = inv(M)
        a, b = p
        ea, eb = sqrt(C[1, 1]), sqrt(C[2, 2])
        r = (d .- (a .+ b .* t)) ./ e
        chi = sum(abs2, r[ok]) / max(1, count(ok) - 2)
        sc = sqrt(max(chi, 1.0))
        newok = isfinite.(d) .& isfinite.(e) .& (abs.(r) .< clip * sc) .&
                (abs.(d .- med) .< max(8 * s, 0.05))
        newok == ok && break
        ok = newok
    end
    rms = sqrt(mean(abs2, (d .- (a .+ b .* t))[ok]))
    return (a, b, ea, eb, chi, count(ok), rms)
end

"""Per-peak displacement (px) of comb peaks at pixels `x` on fiber f, chip c against the comb
night solution (lin, nl, L, eps): delta = (λ(x) - λ_m)/(dλ/dx), λ_m = 2L/(m + eps), m the nearest
integer. Returns (delta, m_frac)."""
function refit_comb_delta(x::AbstractVector, f::Integer, c::Integer, lin, nl, L, eps)
    n = length(x)
    dl = fill(NaN, n)
    mf = fill(NaN, n)
    n == 0 && return dl, mf
    l0 = refit_prod_wave(lin, nl, f, c, x)
    lp = refit_prod_wave(lin, nl, f, c, x .+ 0.5)
    lm = refit_prod_wave(lin, nl, f, c, x .- 0.5)
    for i in 1:n
        isfinite(l0[i]) || continue
        mc = 2L / l0[i] - eps
        m = round(mc)
        mf[i] = mc - m
        dl[i] = (l0[i] - 2L / (m + eps)) / (lp[i] - lm[i])
    end
    return dl, mf
end

"""
    refit_fpifib_exposure(fnames_by_chip, fibers, comb; nsig = 8.0) -> NamedTuple

Stage 4 for one exposure: the local-noise finder on the FPI fibers (`fibers`, up to 2) of the
three ar1Dcal files (read-only), per-peak displacement against the night comb, and per (fiber,
chip) a + b t fits. Arrays are (slot, chip): a, b, ea, eb, n, rms, mfrac_max.
"""
function refit_fpifib_exposure(fnames_by_chip::AbstractVector{<:AbstractString},
        fibers::AbstractVector{<:Integer}, comb; nsig = 8.0)
    ns = 2
    a = fill(NaN, ns, N_CHIPS)
    b = fill(NaN, ns, N_CHIPS)
    ea = fill(NaN, ns, N_CHIPS)
    eb = fill(NaN, ns, N_CHIPS)
    n = zeros(Int, ns, N_CHIPS)
    rms = fill(NaN, ns, N_CHIPS)
    mfx = fill(NaN, ns, N_CHIPS)
    for (c, fn) in enumerate(fnames_by_chip)
        isfile(fn) || continue
        h5open(fn, "r") do fh
            for (s, fb) in enumerate(fibers[1:min(end, ns)])
                fl = fh["flux_1d"][:, fb]
                iv = fh["ivar_1d"][:, fb]
                mk = fh["mask_1d"][:, fb]
                _, prm, cov, _ = get_arclamp_peaks_local(fl, iv, mk; nsig = nsig)
                size(prm, 1) == 0 && continue
                x = prm[:, 2]
                er = sqrt.(max.(cov[:, 2, 2], 0.0))
                dl, mf = refit_comb_delta(x, fb, c, comb.lin, comb.nl, comb.L, comb.eps)
                r = refit_fit_ab(x, dl, er)
                a[s, c], b[s, c], ea[s, c], eb[s, c] = r[1], r[2], r[3], r[4]
                n[s, c] = r[6]
                rms[s, c] = r[7]
                fm = filter(isfinite, mf)
                mfx[s, c] = isempty(fm) ? NaN : maximum(abs.(fm))
            end
        end
    end
    return (; a, b, ea, eb, n, rms, mfrac_max = mfx, fibers = collect(Int, fibers[1:min(end, ns)]))
end

"""A (fiber-slot, chip) FPI-fiber fit is usable: it looks like a comb (>= 80 clipped peaks on the
chip; a lit comb gives 120-210) with a comb-like residual (rms < 0.1 px; arc-lamp lines landing on
an FPI fiber give 3-5 px) and a centroid-level zero-point error (< 0.005 px)."""
refit_fpifib_usable(P, s, c) = P.n[s, c] >= 80 && isfinite(P.a[s, c]) && P.ea[s, c] < 0.005 &&
                               (!hasproperty(P, :rms) || P.rms[s, c] < 0.1)

"""Per-exposure displacement model from the FPI-fiber fits for a term set (`exp_terms`).
Returns a function (f, c, p) -> δ px, plus (n_fibers_used, tilt_px per chip, flags)."""
function refit_fpifib_transform(P, exp_terms::AbstractString)
    use = [all(refit_fpifib_usable(P, s, c) for c in 1:N_CHIPS) for s in 1:length(P.fibers)]
    nf = count(use)
    nf == 0 && return nothing, 0, fill(NaN, 3), WAVE_FPIFIB_UNAVAILABLE
    flags = 0
    ss = findall(use)
    za = [mean(P.a[ss, c]) for c in 1:N_CHIPS]
    zb = [mean(P.b[ss, c]) for c in 1:N_CHIPS]
    tilt = fill(NaN, 3)
    if exp_terms == "shift"
        sh = mean(za)
        return ((f, c, p) -> sh), nf, tilt, flags
    elseif exp_terms == "chipzp"
        return ((f, c, p) -> za[c]), nf, tilt, flags
    elseif exp_terms == "chipzp_stretch" || nf < 2
        exp_terms == "chipzp_stretch_tilt" && (flags |= WAVE_FPIFIB_ONE_FIBER)
        return ((f, c, p) -> za[c] + zb[c] * (p - 1024.0) / 1024.0), nf, tilt, flags
    end
    f1, f2 = P.fibers[1], P.fibers[2]
    for c in 1:N_CHIPS
        tilt[c] = P.a[2, c] - P.a[1, c]      # px between the two FPI fibers
        er = sqrt(P.ea[1, c]^2 + P.ea[2, c]^2)
        abs(tilt[c]) > 3 * max(er, 0.001) && (flags |= WAVE_FPIFIB_ROTATION)
    end
    a1, b1, a2, b2 = P.a[1, :], P.b[1, :], P.a[2, :], P.b[2, :]
    fn = function (f, c, p)
        w = (f - f1) / (f2 - f1)
        u = (p - 1024.0) / 1024.0
        return (1 - w) * (a1[c] + b1[c] * u) + w * (a2[c] + b2[c] * u)
    end
    return fn, nf, tilt, flags
end

# ------------------------------------------------------------------------------------------
# stage 5: sky lines (pooled, robust to misassigned OH lines) and per-fiber sky terms
# ------------------------------------------------------------------------------------------
"""OH-line rows of an object exposure from the pass-1 `skyLinePeaks_*` files (read-only), with
pass-1's ingest mask (a line whose across-fiber IQR in x exceeds 0.002 is dropped). Returns a
NamedTuple of vectors (fiber, chip, line, pix, lam) or `nothing` if a chip file is missing."""
function refit_read_sky(mdir::AbstractString, tele::AbstractString, mjd::Integer, e::Integer)
    rows = (fiber = Int[], chip = Int[], line = Int[], pix = Float64[], lam = Float64[])
    for (c, chip) in enumerate(CHIP_LIST)
        fn = joinpath(mdir, "skyLinePeaks_$(tele)_$(mjd)_$(_arc_expstr(e))_$(chip)_object.h5")
        isfile(fn) || return nothing
        M = h5open(f -> read(f["sky_line_mat_clean"]), fn, "r")
        size(M, 3) == N_FIBERS || return nothing
        for l in 1:size(M, 1)
            xx = filter(isfinite, (M[l, 1, :] .- 1024) ./ 2048)
            xx = filter(!=(0), xx)
            length(xx) < 10 && continue
            q = quantile(xx, [0.25, 0.75])
            q[2] - q[1] > 0.002 && continue
            for f in 1:N_FIBERS
                p = M[l, 1, f]
                w = M[l, 2, f]
                (isfinite(p) && isfinite(w) && p > 1) || continue
                push!(rows.fiber, f)
                push!(rows.chip, c)
                push!(rows.line, round(Int, M[l, 9, f]))
                push!(rows.pix, p)
                push!(rows.lam, w)
            end
        end
    end
    return rows
end

"Velocity (m/s) of OH rows against a model: c (λ_OH - λ_model(pix))/λ_model, model = callback."
function refit_sky_velocities(rows, modelfun)
    v = fill(NaN, length(rows.pix))
    for i in eachindex(v)
        lm = modelfun(rows.fiber[i], rows.chip[i], rows.pix[i])
        isfinite(lm) && (v[i] = (rows.lam[i] - lm) / lm * REFIT_C_MPS)
    end
    return v
end

"""Pooled per-chip sky offset (m/s), robust to misassigned OH lines (the ARCFPI_VALIDATION
estimator): per line, rows farther than max(5 NMAD, 600 m/s) from the line median are dropped (a
wrong OH component sits >= 0.1 Å = 1.8 km/s away); a line needs `minfib` surviving fibers and
contributes its median (line-balanced); the chip value is a clipped mean over lines (2 rounds of
2.5 sigma about the median, sigma floor 5 m/s); a chip needs `minlines` lines.
Returns (zp[3], nlines[3], keep::Vector{Bool})."""
function refit_pooled_sky(rows, v; minfib = 30, minlines = 3)
    zp = fill(NaN, 3)
    nl = zeros(Int, 3)
    keep = falses(length(v))
    for c in 1:3
        ii = findall(==(c), rows.chip)
        isempty(ii) && continue
        lv = Float64[]
        lrows = Vector{Int}[]
        for l in unique(rows.line[ii])
            jj = ii[rows.line[ii] .== l]
            vv = v[jj]
            ok = isfinite.(vv)
            count(ok) < minfib && continue
            m = median(vv[ok])
            s = 1.4826 * median(abs.(vv[ok] .- m))
            gsel = jj[ok .& (abs.(vv .- m) .< max(5s, 600.0))]
            length(gsel) < minfib && continue
            push!(lrows, gsel)
            push!(lv, median(v[gsel]))
        end
        length(lv) < minlines && continue
        use = trues(length(lv))
        for _ in 1:2
            m = median(lv[use])
            s = max(1.4826 * median(abs.(lv[use] .- m)), 5.0)
            nu = abs.(lv .- m) .< 2.5s
            count(nu) == 0 && break
            use = nu
        end
        nl[c] = count(use)
        nl[c] < minlines && continue
        zp[c] = mean(lv[use])
        for k in findall(use)          # rows of whole-line misassignments are not kept either
            keep[lrows[k]] .= true
        end
    end
    return zp, nl, keep
end

"""Per-fiber sky terms (m/s) on the rows kept by `refit_pooled_sky`: per fiber and chip a robust
zero point and stretch (v = z + s u, u = (pix - 1024.5)/1024), fiber median of z removed (fiber-
relative only). Returns (z, s) each N_FIBERS x 3 (NaN where < 3 lines)."""
function refit_fiber_sky_terms(rows, v, keep)
    z = fill(NaN, N_FIBERS, N_CHIPS)
    s = fill(NaN, N_FIBERS, N_CHIPS)
    for c in 1:N_CHIPS
        for f in 1:N_FIBERS
            ii = [i for i in eachindex(v) if keep[i] && rows.chip[i] == c && rows.fiber[i] == f]
            length(ii) < 3 && continue
            u = (rows.pix[ii] .- 1024.5) ./ 1024
            y = v[ii]
            if length(ii) >= 5 && (maximum(u) - minimum(u)) > 0.3
                X = hcat(ones(length(y)), u)
                w = ones(length(y))
                p = zeros(2)
                for _ in 1:3
                    p = (X .* w) \ (y .* w)
                    r = y .- X * p
                    sg = max(1.4826 * median(abs.(r)), 5.0)
                    w = sqrt.(min.(1.0, 1.5sg ./ max.(abs.(r), 1e-9)))
                end
                z[f, c], s[f, c] = p[1], p[2]
            else
                z[f, c] = median(y)
                s[f, c] = 0.0
            end
        end
        m = _arc_nanmedian(z[:, c])
        z[:, c] .-= m
    end
    return z, s
end

"""Per-chip pooled-sky per-exposure transform (m/s) on the rows kept by `refit_pooled_sky` (the QA
that removes wrong-peak matches): theta[:, c] = (zp, stretch, tilt) of v = zp + st u + tl phi,
u = (pix - 1024.5)/1024, phi = (fiber - 150.5)/149.5. "chipzp": zp = the line-balanced pooled
value, st = tl = 0. Otherwise a robust (Huber IRLS) fit on the kept rows of the terms requested.
Chips with fewer than 3 lines are NaN."""
function refit_sky_terms(rows, v, keep, zp_pooled, nlines, terms::AbstractString; minlines = 3)
    th = fill(NaN, 3, 3)
    for c in 1:3
        nlines[c] >= minlines || continue
        if terms == "chipzp"
            th[:, c] .= (zp_pooled[c], 0.0, 0.0)
            continue
        end
        ii = [i for i in eachindex(v) if keep[i] && rows.chip[i] == c && isfinite(v[i])]
        length(ii) < 20 && continue
        u = (rows.pix[ii] .- 1024.5) ./ 1024
        ph = (rows.fiber[ii] .- 150.5) ./ 149.5
        y = v[ii]
        X = terms == "chipzp_stretch" ? hcat(ones(length(y)), u) : hcat(ones(length(y)), u, ph)
        w = ones(length(y))
        b = zeros(size(X, 2))
        for _ in 1:4
            b = (X .* w) \ (y .* w)
            r = y .- X * b
            sg = max(1.4826 * median(abs.(r)), 5.0)
            w = sqrt.(min.(1.0, 1.5sg ./ max.(abs.(r), 1e-9)))
        end
        th[1, c] = b[1]
        th[2, c] = b[2]
        th[3, c] = size(X, 2) == 3 ? b[3] : 0.0
    end
    return th
end

"Apply per-chip sky terms theta (3 x 3, m/s) multiplicatively to a grid (N_XPIX x N_FIBERS x 3)."
function refit_apply_sky_terms!(G::AbstractArray, theta::AbstractMatrix)
    u = ((1:size(G, 1)) .- 1024.5) ./ 1024
    for c in 1:3
        all(isfinite, view(theta, :, c)) || continue
        for f in 1:size(G, 2)
            ph = (f - 150.5) / 149.5
            G[:, f, c] .*= 1 .+ (theta[1, c] .+ theta[2, c] .* u .+ theta[3, c] * ph) ./ REFIT_C_MPS
        end
    end
    return G
end

# ------------------------------------------------------------------------------------------
# plate era: comb shape template (file written by scripts/refit/build_comb_shape_template.jl)
# ------------------------------------------------------------------------------------------
"Legendre polynomials P0..P6 at t."
@inline function refit_legendre7(t)
    P = zeros(7)
    P[1] = 1.0
    P[2] = t
    for k in 2:6
        P[k + 1] = ((2k - 1) * t * P[k] - (k - 1) * P[k - 1]) / k
    end
    return P
end

"""Comb shape template of one telescope and block: per (Legendre coefficient P0..P6, chip,
fiber) the median pixel-displacement field of the comb night solutions against a fixed reference
solution (`ref_lin`, `ref_nl`); t = (p - 1025)/960. P2..P6 carry the shape."""
struct RefitTemplate
    tele::String
    block::String
    T::Array{Float64, 3}
    ref_lin::Matrix{Float64}
    ref_nl::Matrix{Float64}
    fingerprint::String
end

function load_refit_template(path::AbstractString, tele::AbstractString, block::AbstractString)
    (isempty(path) || !isfile(path)) && return nothing
    h5open(path, "r") do f
        (haskey(f, tele) && haskey(f[tele], "T_$(block)")) || return nothing
        RefitTemplate(String(tele), String(block), read(f["$(tele)/T_$(block)"]),
            read(f["$(tele)/ref_lin"]), read(f["$(tele)/ref_nl"]), arclamp_file_fingerprint(path))
    end
end

"Reference solution wavelength and dλ/dpix of fiber f, chip c at (possibly fractional) pixel p."
function refit_ref_wave_disp(Tp::RefitTemplate, f, c, p)
    x = (p - 1024) / 2048
    cp = size(Tp.ref_nl, 2) ÷ 2 - 1
    if c == 2
        xt, dc = x, 1.0
    else
        o = c == 1 ? 0 : cp + 1
        xt, dc = cp == 0 ? (Tp.ref_nl[f, o + 1] + x, 1.0) :
                 (Tp.ref_nl[f, o + 1] + Tp.ref_nl[f, o + 2] * x, Tp.ref_nl[f, o + 2])
    end
    lin = view(Tp.ref_lin, f, :)
    lam = evalpoly(xt, lin)
    dl = 0.0
    for k in 2:length(lin)
        dl += (k - 1) * lin[k] * xt^(k - 2)
    end
    return lam, dl * dc / 2048
end

"Pixel where the reference solution of fiber f, chip c equals `lambda` (Newton from q0)."
function refit_ref_invert(Tp::RefitTemplate, f, c, lambda, q0)
    q = q0
    for _ in 1:30
        l, s = refit_ref_wave_disp(Tp, f, c, q)
        st = (l - lambda) / s
        q -= st
        abs(st) < 1e-9 && break
    end
    return q
end

"Template model wavelength of fiber f, chip c at reference-frame pixel pr; P = (zp_R, zp_G, zp_B, stretch)."
function refit_template_lambda(Tp::RefitTemplate, P, f, c, pr)
    t = (pr - 1025) / 960
    Lg = refit_legendre7(t)
    q = pr + P[c] + (Tp.T[2, c, f] + P[4]) * t
    for k in 3:7
        q += Tp.T[k, c, f] * Lg[k]
    end
    return refit_ref_wave_disp(Tp, f, c, q)[1]
end

"Pixel shift (px) of an arc state on chip c of fiber f at pixel p (reference-frame pixel = p + shift)."
function refit_state_shift_px(A::RefitArc, f, c, p, st)
    a, b = refit_arc_shift_AB(A, f, c, st)
    x = (p - 1024) / 2048
    bc = c == 1 ? A.TH[f, 2] : c == 3 ? A.TH[f, 4] : 1.0
    return N_XPIX * (a + b * x) / bc
end

"""Fit the template to the night's kept arc lines (lab wavelengths of list `ll`): per fiber 3 chip
zero points + 1 common stretch, the global position terms from the arc fit. Returns (P 4 x N, ok)."""
function refit_fit_template(Tp::RefitTemplate, A::RefitArc, ll::ArcLineList)
    P = fill(NaN, 4, N_FIBERS)
    ok = falses(N_FIBERS)
    L = A.lines
    jmap = Dict(s => j for (j, s) in enumerate(A.used))
    rowsf = [Int[] for _ in 1:N_FIBERS]
    for i in eachindex(L.id)
        L.keep[i] || continue
        haskey(jmap, L.seq[i]) || continue
        haskey(ll.labwave, L.id[i]) || continue
        push!(rowsf[L.fiber[i]], i)
    end
    for f in 1:N_FIBERS
        A.okf[f] || continue
        ii = rowsf[f]
        length(ii) < 30 && continue
        nr = length(ii)
        X = zeros(nr, 4)
        y = zeros(nr)
        w = zeros(nr)
        for (k, i) in enumerate(ii)
            c = L.chip[i]
            j = jmap[L.seq[i]]
            d, t = arcfit_seq_state(A.g, A.lay, j)
            st = (A.grp_u[j], d, t, A.tau_u[j])
            pr = L.pix[i] + refit_state_shift_px(A, f, c, L.pix[i], st)
            q = refit_ref_invert(Tp, f, c, ll.labwave[L.id[i]], pr)
            tt = (pr - 1025) / 960
            Lg = refit_legendre7(tt)
            yy = q - pr - Tp.T[2, c, f] * tt
            for kk in 3:7
                yy -= Tp.T[kk, c, f] * Lg[kk]
            end
            X[k, c] = 1.0
            X[k, 4] = tt
            y[k] = yy
            w[k] = 1.0 / (max(L.pixvar[i], 0.0) + 0.02^2)
        end
        all(any(X[:, c] .!= 0) for c in 1:3) || continue
        good = isfinite.(y)
        X, y, w = X[good, :], y[good], w[good]
        length(y) < 30 && continue
        keep = abs.(y .- median(y)) .< 3.0
        b = zeros(4)
        for _ in 1:4
            sw = sqrt.(w[keep])
            b = (X[keep, :] .* sw) \ (y[keep] .* sw)
            r = y .- X * b
            z = r .* sqrt.(w)
            s = 1.4826 * median(abs.(z[keep]))
            keep = (abs.(r) .< 1.3) .& (abs.(z) .< max(5s, 1e-6))
            count(keep) < 10 && break
        end
        P[:, f] .= b
        ok[f] = true
    end
    return P, ok
end

# ------------------------------------------------------------------------------------------
# stored solution: per fiber, per chip Legendre coefficients in t = (pix - 1024.5)/1023.5
# ------------------------------------------------------------------------------------------
"Legendre basis (N_XPIX x (deg+1)) on t = (pix - 1024.5)/1023.5."
function refit_coeff_basis(deg::Integer)
    t = ((1:N_XPIX) .- 1024.5) ./ 1023.5
    B = zeros(N_XPIX, deg + 1)
    B[:, 1] .= 1.0
    deg >= 1 && (B[:, 2] .= t)
    for k in 2:deg
        B[:, k + 1] .= ((2k - 1) .* t .* B[:, k] .- (k - 1) .* B[:, k - 1]) ./ k
    end
    return B
end

"""Fit the per-fiber, per-chip coefficients to a wavelength grid (N_XPIX x N_FIBERS x 3).
Returns (C (deg+1) x N_FIBERS x 3, max |residual| Å). Non-finite fiber-chips give NaN."""
function refit_grid_to_coeffs(G::AbstractArray{<:Real, 3}, B::AbstractMatrix)
    deg1 = size(B, 2)
    C = fill(NaN, deg1, N_FIBERS, N_CHIPS)
    F = qr(B)
    mx = 0.0
    for c in 1:N_CHIPS, f in 1:N_FIBERS
        y = view(G, :, f, c)
        all(isfinite, y) || continue
        cf = F \ Vector{Float64}(y)
        C[:, f, c] .= cf
        mx = max(mx, maximum(abs.(B * cf .- y)))
    end
    return C, mx
end

"Evaluate stored coefficients to a wavelength grid (N_XPIX x N_FIBERS x 3), comparable with chipWaveSoln."
function refit_coeffs_to_grid(C::AbstractArray{<:Real, 3})
    B = refit_coeff_basis(size(C, 1) - 1)
    G = fill(NaN, N_XPIX, size(C, 2), size(C, 3))
    for c in axes(C, 3), f in axes(C, 2)
        all(isfinite, view(C, :, f, c)) || continue
        G[:, f, c] .= B * C[:, f, c]
    end
    return G
end

"""
    refit_exposure_wave(fname, expnum) -> Array (N_XPIX, N_FIBERS, 3)

The refit wavelength solution of exposure `expnum` from a night product
`waveCalRefit_<tele>_<mjd>.h5` (directly comparable with pass 1's `chipWaveSoln`).
"""
function refit_exposure_wave(fname::AbstractString, expnum::Integer)
    h5open(fname, "r") do f
        en = read(f["exp_expnum"])
        k = findfirst(==(expnum), en)
        k === nothing && throw(ArgumentError("exposure $(expnum) not in $(fname)"))
        refit_coeffs_to_grid(f["wave_coeffs"][:, :, :, k])
    end
end

# ------------------------------------------------------------------------------------------
# pass-1 comparison
# ------------------------------------------------------------------------------------------
"Pass-1 adopted solution file and dataset for an exposure (read-only), or `nothing`."
function refit_pass1_source(apred_dir, wavecal_dir, tele, mjd, e, kind, wt)
    mdir = joinpath(apred_dir, string(mjd))
    es = _arc_expstr(e)
    if wt == "waveCalNightfpiDither" || wt == "waveCalNightskyDither"
        fn = joinpath(mdir, "$(wt)_$(tele)_$(mjd)_$(es)_$(kind).h5")
        return isfile(fn) ? (fn, "chipWaveSoln") : nothing
    elseif startswith(wt, "error_night_ave")
        fn = joinpath(wavecal_dir, "wavecalNightAve_$(tele)_$(mjd).h5")
        grp = endswith(wt, "fpi") ? "fpi" : "sky"
        return isfile(fn) ? (fn, "$(grp)/nightAve_wave_soln") : nothing
    end
    return nothing
end

"""Comparison of a new grid on the flat pixel grid against pass 1: dv[c, f] = flat mean over
REFIT_PXG of c (λ_new - λ_p1)/λ_p1 (m/s; the first-order RV prediction of HOLDOUT_RV/AV2), the
per-chip fiber-median zero point, the fiber-median shape rms (dv about its flat mean) and the
across-fiber robust sigma of dv."""
function refit_vs_pass1(Gnew_pxg::AbstractArray, Gp1_pxg::AbstractArray)
    dv = fill(NaN32, N_CHIPS, N_FIBERS)
    shape = fill(NaN, N_CHIPS, N_FIBERS)
    for f in 1:N_FIBERS, c in 1:N_CHIPS
        a = view(Gnew_pxg, :, f, c)
        b = view(Gp1_pxg, :, f, c)
        (all(isfinite, a) && all(isfinite, b)) || continue
        v = (a .- b) ./ b .* REFIT_C_MPS
        m = mean(v)
        dv[c, f] = Float32(m)
        shape[c, f] = sqrt(mean((v .- m) .^ 2))
    end
    zp = [_arc_nanmedian(Float64.(dv[c, :])) for c in 1:N_CHIPS]
    sh = [_arc_nanmedian(shape[c, :]) for c in 1:N_CHIPS]
    rs = [_arc_robstd(Float64.(dv[c, :])) for c in 1:N_CHIPS]
    return dv, zp, sh, rs
end
