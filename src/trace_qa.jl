# Trace quality: per-flat QA records, a gate, a deterministic choice of the night's trace,
# whole-slit (template-free) fiber identification, and the per-fiber trace bitmask.
#
# Background (2026-11-17 trace study, TRACE_ROBUST.md / TRACE_ROBUST_V2.md): the 0.4.0 trace
# problems were failures of the CHOICE of trace (last flat by file name, a fixed survey-wide
# fallback file) and of fiber IDENTIFICATION on frames that are not normal flats (faint,
# saturated, partially lit), plus a two-neighbour extrapolation of dark fibers that walks off
# across a dark 30-fiber bundle. Production's per-column Gauss-Hermite refinement of lit fibers
# on a good flat was not the problem and is kept.
#
# This file follows the calibration-selection pattern the newer stages converged on (FPI gate,
# fluxing-flat choice): measured per-exposure statistics with thresholds placed in a gap of the
# archive distribution, a deterministic choice by a stated rule, a named fallback, and the
# decision plus its statistics written next to the product and carried as bits into the 1D
# products.
#
# Behaviour switches. The defaults are the robust behaviour (PI decision 2026-11-18); the
# pre-existing behaviour is one switch away and is reproduced bit for bit:
#   make_traces_from_flats.jl --trace_identify robust|production   (fiber identification)
#   pipeline_2d_1d.jl         --trace_select   gated|last          (which flat of the night)
#   pipeline_2d_1d.jl         --trace_fallback none|fixed_file     (no own-night flat)
# With --trace_fallback none a night without an eligible flat of its own is not extracted
# (listed loudly) until the epoch templates (PR B) give it a template level. Nothing here needs
# a template or epoch file: the identification is template-free and uses only the static
# reference positions in data/traceCenterPositions_*.
# What is always on: the per-flat QA record, the traceMain QA/provenance metadata and the
# per-fiber `bitmsk_trace` in the 1D products.

using HDF5
using Statistics: median, mean, quantile
using Polynomials: fit

_medornan(v) = isempty(v) ? NaN : median(v)

# ---------------------------------------------------------------------------------------------
# Per-fiber trace bitmask (`bitmsk_trace`), stored as `TraceBits` = UInt16: the smallest type
# that holds the 10 bits defined here (2^0-2^9, five of them reserved for PR B) plus the
# per-exposure registration bit planned for PR C, with 5 bits to spare. Unsigned like the
# other compact flag words (`exposure_flags` UInt8, `wavecal_flags` UInt32). Written per chip
# into traceMain and ar1D* (length N_FIBERS) and stacked (N_CHIPS, N_FIBERS) into ar1Duni*;
# `trace_flags` (metadata) is the same type. Readers convert with `to_trace_bits`, so files
# written before this change (Int64 datasets) read the same.
# Night-level bits are replicated into every fiber, so a consumer needs only this array.
# Severity follows the wavelength-flag scheme (ARCLAMP_PRODUCTION_DESIGN.md §2.4,
# `wavecal_flags`): "bad" bits are proposed for the default mask (arMADGICS, PR C).
# Bits 1, 2, 3, 5, 6 are set only with epoch templates (PR B: src/slit_model.jl, src/trace_levels.jl).
# ---------------------------------------------------------------------------------------------
"Storage type of `bitmsk_trace` and `trace_flags`."
const TraceBits = UInt16
"Convert stored trace bits (any integer type, e.g. Int64 from files written before UInt16) to `TraceBits`."
to_trace_bits(x::Integer) = TraceBits(x)
to_trace_bits(x::AbstractArray{<:Integer}) = TraceBits.(x)
"info: trace from a gated flat of the night that is not the preferred kind (dome, not quartz)."
const TRACE_OTHER_FLAT = TraceBits(2^0)
"info: level T1, the epoch template registered on another calibration frame of the night (src/trace_levels.jl)."
const TRACE_REGISTERED_CAL = TraceBits(2^1)
"info: level T2, the epoch template registered on the stacked science frames of the night."
const TRACE_REGISTERED_SCI = TraceBits(2^2)
"bad: level T3, the epoch template at the epoch-median registration (nothing of the night registered)."
const TRACE_TEMPLATE_ONLY = TraceBits(2^3)
"bad: this fiber had no light on the trace source (median height < 5 % of the median fiber);
its position is not measured but placed from its lit neighbours."
const TRACE_FIBER_UNLIT = TraceBits(2^4)
"screen (per fiber): departure of a template-seeded flat fit from the registered template > 0.3 px."
const TRACE_FIBER_DELTA = TraceBits(2^5)
"info: the night is after the newest flat of an open-ended (unconfirmed) last epoch."
const TRACE_EPOCH_PROVISIONAL = TraceBits(2^6)
"bad: the night used the fixed survey-wide fallback file (`data/fallback_traces`), i.e. a
trace from another night and possibly another instrument epoch."
const TRACE_FIXED_FALLBACK = TraceBits(2^7)
"warn: the trace source did not pass every gate criterion (`trace_gate_reasons` says which)."
const TRACE_UNGATED = TraceBits(2^8)
"warn: no QA record exists for the trace source (products made before QA records existed)."
const TRACE_QA_MISSING = TraceBits(2^9)

"Per-fiber bits (not replicated night-level bits); excluded from the night-level `trace_flags`."
const TRACE_PER_FIBER_BITS = TRACE_FIBER_UNLIT | TRACE_FIBER_DELTA

"Bits proposed for the default mask downstream (arMADGICS PR C)."
const TRACE_BAD_BITS = TRACE_TEMPLATE_ONLY | TRACE_FIBER_UNLIT | TRACE_FIXED_FALLBACK

const TRACE_BIT_NAMES = Dict(
    TRACE_OTHER_FLAT => "TRACE_OTHER_FLAT", TRACE_REGISTERED_CAL => "TRACE_REGISTERED_CAL",
    TRACE_REGISTERED_SCI => "TRACE_REGISTERED_SCI", TRACE_TEMPLATE_ONLY => "TRACE_TEMPLATE_ONLY",
    TRACE_FIBER_UNLIT => "TRACE_FIBER_UNLIT", TRACE_FIBER_DELTA => "TRACE_FIBER_DELTA",
    TRACE_EPOCH_PROVISIONAL => "TRACE_EPOCH_PROVISIONAL",
    TRACE_FIXED_FALLBACK => "TRACE_FIXED_FALLBACK", TRACE_UNGATED => "TRACE_UNGATED",
    TRACE_QA_MISSING => "TRACE_QA_MISSING")

"Names of the bits set in `b`, for logs."
trace_bit_names(b::Integer) = [v for (k, v) in sort(collect(TRACE_BIT_NAMES)) if (b & k) != 0]

# Trace level of the night (traceMain metadata `trace_level`).
const TRACE_LEVEL_OWN_FLAT = 0        # T0: a flat of the night
const TRACE_LEVEL_REGISTERED_CAL = 1  # T1 (PR B)
const TRACE_LEVEL_REGISTERED_SCI = 2  # T2 (PR B)
const TRACE_LEVEL_TEMPLATE = 3        # T3 (PR B)
const TRACE_LEVEL_FIXED_FALLBACK = 4  # legacy: fixed survey-wide file
const TRACE_LEVEL_NONE = -1           # no trace (only with --trace_fallback none)

"Night-level bit implied by a trace level."
function trace_level_bits(level::Integer)
    level == TRACE_LEVEL_REGISTERED_CAL && return TRACE_REGISTERED_CAL
    level == TRACE_LEVEL_REGISTERED_SCI && return TRACE_REGISTERED_SCI
    level == TRACE_LEVEL_TEMPLATE && return TRACE_TEMPLATE_ONLY
    level == TRACE_LEVEL_FIXED_FALLBACK && return TRACE_FIXED_FALLBACK
    return zero(TraceBits)
end

# ---------------------------------------------------------------------------------------------
# Gate thresholds. Calibrated on pass-1 (dr21) flats, 2026-11-17; the distributions and the
# gaps the thresholds sit in are documented in the PR (validation folder 2026_11_17/trace_pr).
# ---------------------------------------------------------------------------------------------
"A fiber is lit if its median height exceeds this fraction of the median fiber (as production)."
const TRACE_LIT_FRAC = 0.05
"Minimum number of lit fibers for a flat to serve as a full trace source."
const TRACE_GATE_NLIT_MIN = 240
"""
Maximum robust RMS (px) of the per-fiber band centers about a per-fiber cubic in x.
Measures whether each fiber's trace is a smooth curve (defocused, sparse, faint or
mis-walked flats are not). Pass-1 traceMain flats of 680 no-regression nights (2,040
chip-flats): max 0.0059 (APO) / 0.0098 (LCO) px; refused frames of the wrong nights
0.04-0.24 px (the LCO 57817 focus sweep 0.04-0.08).
"""
const TRACE_GATE_FIT_RMS_MAX = 0.02
"""
Maximum fraction of lit fibers, in the worst band across the chip, whose trace center is
more than 1 px from the nearest fiber peak of the light (`trace_light_qa`, 9 bands). A block of fibers
on the wrong alias, or an outward column walk that twists half the slit, fails here even when
every fiber is individually smooth; so does a defocused or partially lit flat. (The spacing
pattern statistics `spacing_rms`/`spacing_nbad` are recorded but not gated: on faint early LCO
flats they fire on traces the light check accepts.) Pass-1 traceMain flats of 680 no-regression
nights: p99 0.013 / 0.003, max 0.017 (APO) / 0.003 (LCO); production traces of the wrong nights
0.05-1.0, their good flats <= 0.019.
"""
const TRACE_GATE_LIGHT_MISMATCH_MAX = 0.03
"""
Minimum fraction of lit fibers whose central position agrees (< 0.35 spacing) with the
whole-slit identification of the light (`trace_identify_slit`). A labelling error moves
whole blocks of fibers by a spacing.
"""
const TRACE_GATE_ID_MATCH_MIN = 0.98
"""
Soft criterion: fraction of chip pixels flagged fully saturated. LCO R carries ~0.018 on every
flat (chronic defects); the saturated APO 57969 exposure 9 has 0.23-0.24 on G/B.
"""
const TRACE_GATE_SAT_FRAC_MAX = 0.05

"Minimum alias margin (cost of the best alias over the best solution) to accept a labelling."
const TRACE_ID_MARGIN_MIN = 1.5
"Census: peaks must have at least this S/N (height over the combined-profile error)."
const TRACE_CENSUS_SNR_MIN = 5.0

# lamp-off thresholds (moved verbatim from make_traces_from_flats.jl): ~20 % of the typical
# value of the faintest chip and flat kind on lamp-on nights. flux_med is recorded only
# (see the note in make_traces_from_flats.jl).
lamp_flux68_thresh(tele) = tele == "apo" ? 40.0 : 10.0
lamp_fluxmed_thresh(tele) = tele == "apo" ? 16.0 : 10.0

"""
    trace_lamp_metrics(image_data, pix_bitmask, tele)

The lamp-off test of `make_traces_from_flats.jl`, unchanged: over pixels without
`bad_pix_bits`, the half 16-84 % range of the (non-NaN, non-zero) values must reach
`lamp_flux68_thresh(tele)`. Returns `(flux_med, flux_68p, thresh, lamp_ok, frac_good)`.
The single implementation used both for the decision and for the QA record.
"""
function trace_lamp_metrics(image_data, pix_bitmask, tele)
    good_pixels = ((pix_bitmask .& bad_pix_bits) .== 0)
    s = nanzeropercentile(image_data[good_pixels], percent_vec = [16, 50, 84])
    flux_med = s[2]
    flux_68p = 0.5 * (s[3] - s[1])
    thresh = lamp_flux68_thresh(tele)
    # `!(x < t)` reproduces the original `if flux_68p < thresh -> skip` exactly, NaN included
    return (flux_med = flux_med, flux_68p = flux_68p, thresh = thresh,
        lamp_ok = !(flux_68p < thresh), frac_good = sum(good_pixels) / length(good_pixels))
end

# ---------------------------------------------------------------------------------------------
# Center-band profile (production's algorithm, factored out of `trace_extract`)
# ---------------------------------------------------------------------------------------------
"""
    trace_center_band(image_data, noise_image, good_pixels; mid = 1025, n_center_cols = 100)

Collapse columns `mid ± n_center_cols` into one spatial profile: per-row outlier rejection
about the median (5 sigma) and an inverse-variance weighted mean. Returns
`(comb_fluxes, comb_errs, comb_ivars)` as length-2048 vectors (comb_errs = Inf where no
pixel survived). This is exactly the computation `trace_extract` uses to find fibers.
`nsigma` (default 5, production's value) is the per-row clipping threshold; `Inf` keeps every
good pixel (used by `slit_bands` for arc frames, whose light is in a few columns per row).
"""
function trace_center_band(image_data, noise_image, good_pixels; mid = 1025, n_center_cols = 100,
        nsigma = 5)
    r = (mid - n_center_cols):(mid + n_center_cols)
    cutout_fluxes = image_data[r, :]'
    cutout_errs = noise_image[r, :]'
    cutout_masks = good_pixels[r, :]'

    cutout_fluxes[.!cutout_masks] .= 0
    cutout_ivars = cutout_errs .^ -2
    cutout_ivars[.!cutout_masks] .= 0
    cutout_errs[.!cutout_masks] .= Inf

    med_flux = nanzeromedian(cutout_fluxes, 2)
    sigma_dists = (cutout_fluxes .- med_flux) .* (cutout_ivars .^ (0.5))
    good_sigmas = (abs.(sigma_dists) .< nsigma) .& (cutout_ivars .!= 0)

    cutout_fluxes[.!good_sigmas] .= 0
    cutout_ivars[.!good_sigmas] .= 0
    cutout_errs[.!good_sigmas] .= Inf

    comb_ivars = nansum(cutout_ivars, 2)
    good_ivars = comb_ivars .> 0
    comb_errs = zeros(size(comb_ivars)) .+ Inf
    comb_fluxes = zeros(size(comb_ivars))
    comb_fluxes[good_ivars] = nansum((cutout_fluxes .* cutout_ivars), 2)[good_ivars] ./
                              comb_ivars[good_ivars]
    comb_errs[good_ivars] = comb_ivars[good_ivars] .^ (-0.5)
    return comb_fluxes, comb_errs, comb_ivars
end

# ---------------------------------------------------------------------------------------------
# Template-free peak census of the center band
# ---------------------------------------------------------------------------------------------
"""
    trace_peak_census(comb_fluxes, comb_errs; lit_frac, snr_min, edge = 11)

Find fiber peaks in a center-band profile without any reference: local maxima of a lightly
smoothed profile, height above the local minimum within ±4 px, kept if the height is at least
`lit_frac` of the median bright-peak height and `snr_min` times the profile error; peaks
closer than 3 px are merged (the higher survives); positions refined by a parabola.

Also numbers the peaks into slots from their own median spacing (gaps of n spacings are
n-1 missing slots). For a flat on which every fiber is lit, `nslots == 300` fixes the labels
without any reference; this is the independent cross-check of `trace_identify_slit`.
"""
function trace_peak_census(comb_fluxes::AbstractArray, comb_errs::AbstractArray;
        lit_frac = TRACE_LIT_FRAC, snr_min = TRACE_CENSUS_SNR_MIN, edge = 11)
    comb_fluxes = vec(comb_fluxes)
    comb_errs = vec(comb_errs)
    n = length(comb_fluxes)
    p = [isfinite(comb_errs[i]) && isfinite(comb_fluxes[i]) ? Float64(comb_fluxes[i]) : 0.0
         for i in 1:n]
    w = exp.(-0.5 .* ((-3:3) ./ 0.8) .^ 2)
    w ./= sum(w)
    ps = [sum(w[k] * p[clamp(y + k - 4, 1, n)] for k in 1:7) for y in 1:n]
    cand = [y for y in max(edge, 2):min(n - edge, n - 1)
            if ps[y] > ps[y - 1] && ps[y] >= ps[y + 1]]
    h = [ps[y] - minimum(ps[max(1, y - 4):min(n, y + 4)]) for y in cand]
    snr = [h[i] / (isfinite(comb_errs[cand[i]]) ? comb_errs[cand[i]] : Inf) for i in eachindex(cand)]
    empty = (y = Float64[], height = Float64[], snr = Float64[], hmed = NaN, spacing = NaN,
        slots = Int[], nslots = 0, nmissing = 0, median_snr = NaN)
    isempty(h) && return empty
    q90 = quantile(h, 0.9)
    bright = h[h .> 0.3 * q90]
    hmed = isempty(bright) ? NaN : median(bright)
    keep = (h .>= lit_frac * hmed) .& (snr .>= snr_min)
    pk = cand[keep]
    hk = h[keep]
    sk = snr[keep]
    isempty(pk) && return empty
    # merge peaks closer than 3 px
    keep = trues(length(pk))
    for i in 2:length(pk)
        j = findlast(keep[1:(i - 1)])
        isnothing(j) && continue
        if pk[i] - pk[j] < 3
            hk[i] > hk[j] ? (keep[j] = false) : (keep[i] = false)
        end
    end
    pk = pk[keep]
    hk = hk[keep]
    sk = sk[keep]
    # sub-pixel position: parabola through the smoothed maximum
    y = Float64[]
    for q in pk
        a, b, c = ps[q - 1], ps[q], ps[q + 1]
        den = a - 2b + c
        push!(y, q + (den < 0 ? clamp(0.5 * (a - c) / den, -0.5, 0.5) : 0.0))
    end
    sp = length(y) > 1 ? median(diff(y)) : NaN
    slots = Int[1]
    for i in 2:length(y)
        push!(slots, slots[end] + max(1, round(Int, (y[i] - y[i - 1]) / sp)))
    end
    return (y = y, height = hk, snr = sk, hmed = hmed, spacing = sp, slots = slots,
        nslots = slots[end], nmissing = slots[end] - length(y), median_snr = median(sk))
end

# ---------------------------------------------------------------------------------------------
# Whole-slit identification against the reference fiber positions
# ---------------------------------------------------------------------------------------------
function _nearest_sorted(R::AbstractVector, v)
    i = searchsortedfirst(R, v)
    i == 1 && return 1
    i > length(R) && return length(R)
    return (v - R[i - 1]) <= (R[i] - v) ? i - 1 : i
end

function _slit_cost(y, R, a, b, cap2)
    c = 0.0
    for v in y
        # nearest reference fiber to v under v = a + b * R
        u = (v - a) / b
        j = _nearest_sorted(R, u)
        d = v - (a + b * R[j])
        c += min(d * d, cap2)
    end
    return c
end

"""
    trace_identify_slit(y, refpos; nfib = 300, dmax = 40.0, census_nslots = 0)

Label every detected peak (positions `y`, sorted) with a fiber index by registering the
reference fiber positions `refpos` (fiber 1..300 near the chip center; the same static file
production identifies against) to ALL peaks at once: one shift `a` and one scale `b`,
`y = a + b * refpos[f]`. The cost of a solution is the sum over peaks of the squared distance
to the nearest shifted reference fiber, capped at 0.4 fiber spacings, so a peak where the
solution expects no fiber (light beyond fiber 1 or 300, or a mis-registered bundle gap)
costs the cap. Never per-peak rounding: a solution that is half a spacing off would put whole
blocks of fibers on the wrong alias, which is exactly what this replaces.

Returns a NamedTuple with the fiber `label` of each peak (0 = unmatched), `a`, `b`, the
alias `margin` (lowest cost of any solution at least 0.6 spacings away, over the best cost),
`n_unmatched`, `n_dup` (two peaks on one fiber), `resid_rms`, `positions` (all `nfib`
fibers at the center band: matched peaks where lit, the registered reference corrected by
the residuals of the nearest matched neighbours where not), `matched` (Bool per fiber),
`census_agrees` (if the census numbered exactly `nfib` slots, whether its first slot is the
registration's fiber 1), `ok` and `reason`.
"""
function trace_identify_slit(y::AbstractVector, refpos::AbstractVector; nfib = 300,
        dmax = 40.0, dstep = 0.1, census_nslots = 0, census_slots = Int[],
        margin_min = TRACE_ID_MARGIN_MIN, max_unmatched_frac = 0.02)
    R = Float64.(refpos)
    sp = median(diff(R))
    cap2 = (0.4 * sp)^2
    floor2 = (0.05)^2 # per-peak cost floor so a perfect match has a finite margin
    fail(reason) = (label = zeros(Int, length(y)), a = NaN, b = NaN, margin = NaN,
        n_unmatched = length(y), n_outside = 0, n_dup = 0, resid_rms = NaN, positions = fill(NaN, nfib),
        matched = falses(nfib), census_agrees = false, ok = false, reason = reason,
        spacing = sp, cost = NaN)
    length(y) < 10 && return fail("id_too_few_peaks")
    ds = collect((-dmax):dstep:dmax)
    costs = [_slit_cost(y, R, d, 1.0, cap2) for d in ds]
    a = ds[argmin(costs)]
    b = 1.0
    label = zeros(Int, length(y))
    for _ in 1:4
        for (i, v) in enumerate(y)
            j = _nearest_sorted(R, (v - a) / b)
            label[i] = abs(v - (a + b * R[j]))^2 < cap2 ? j : 0
        end
        m = label .> 0
        sum(m) < 10 && return fail("id_too_few_matched")
        X = hcat(ones(sum(m)), R[label[m]])
        a, b = X \ y[m]
    end
    best = _slit_cost(y, R, a, b, cap2)
    # alias: the best solution at least 0.6 spacings away (same scale)
    alias = Inf
    for d in ds
        abs(d - a) >= 0.6 * sp || continue
        alias = min(alias, _slit_cost(y, R, d, b, cap2))
    end
    margin = (alias + floor2 * length(y)) / (best + floor2 * length(y))
    m = label .> 0
    n_unmatched = count(.!m)
    # light beyond the outermost fibers: the decisive evidence against a relabelled slit
    lo_edge = a + b * R[1] - 0.6 * b * sp
    hi_edge = a + b * R[end] + 0.6 * b * sp
    n_outside = count(v -> v < lo_edge || v > hi_edge, y)
    labs = label[m]
    n_dup = length(labs) - length(unique(labs))
    resid = y[m] .- (a .+ b .* R[labs])
    resid_rms = 1.4826 * median(abs.(resid .- median(resid)))
    # positions for every fiber: matched peak, else registered reference + interpolated residual
    positions = a .+ b .* R
    matched = falses(nfib)
    rf = zeros(nfib)
    for (i, j) in enumerate(label)
        j == 0 && continue
        positions[j] = y[i]
        matched[j] = true
        rf[j] = y[i] - (a + b * R[j])
    end
    mi = findall(matched)
    for f in findall(.!matched)
        lo = findlast(<(f), mi)
        hi = findfirst(>(f), mi)
        r = if !isnothing(lo) && !isnothing(hi)
            fl, fh = mi[lo], mi[hi]
            rf[fl] + (rf[fh] - rf[fl]) * (f - fl) / (fh - fl)
        elseif !isnothing(lo)
            rf[mi[lo]]
        elseif !isnothing(hi)
            rf[mi[hi]]
        else
            0.0
        end
        positions[f] += r
    end
    census_agrees = if census_nslots == nfib && length(census_slots) == length(y)
        # census slot s <-> fiber s: every matched peak must carry the same number
        all(label[i] == 0 || label[i] == census_slots[i] for i in eachindex(y))
    else
        true # census not decisive (dark fibers at an edge or unlit gaps): no evidence against
    end
    reason = if n_dup > 0
        "id_duplicate"
    elseif n_outside > 0
        "id_light_outside_slit"
    elseif n_unmatched > max(2, max_unmatched_frac * length(y))
        "id_unmatched_light"
    elseif margin < margin_min
        "id_ambiguous"
    elseif !census_agrees
        "id_census_disagrees"
    else
        "ok"
    end
    return (label = label, a = a, b = b, margin = margin, n_unmatched = n_unmatched,
        n_outside = n_outside, n_dup = n_dup, resid_rms = resid_rms, positions = positions, matched = matched,
        census_agrees = census_agrees, ok = reason == "ok", reason = reason, spacing = sp,
        cost = best)
end

# ---------------------------------------------------------------------------------------------
# Unlit fibers: placed from lit neighbours on BOTH sides along the reference spacing pattern
# ---------------------------------------------------------------------------------------------
"""
    fill_unlit_fibers_robust!(params, unlit, refpos; nnb = 8)

`params` is (n_x, n_fibers, 3) = (height, center, width) per column. For every unlit fiber and
column, fit `center = alpha + beta * refpos` to up to `nnb` nearest lit fibers below AND above
it (one side only at a slit edge) and evaluate at the fiber's reference position; the width is
interpolated linearly in reference position between the nearest lit fibers; the height is 0.
Replaces production's linear extrapolation from the two nearest lit fibers in index, which
walks off by several pixels across a dark 30-fiber bundle (fibers 151-180, APO 57969-58015).
"""
function fill_unlit_fibers_robust!(params::AbstractArray{<:Real, 3}, unlit::AbstractVector{<:Integer},
        refpos::AbstractVector; nnb = 8)
    nfib = size(params, 2)
    isempty(unlit) && return params
    lit = setdiff(1:nfib, unlit)
    isempty(lit) && return params
    R = Float64.(refpos)
    for f in unlit
        below = filter(<(f), lit)
        above = filter(>(f), lit)
        nb = vcat(below[max(1, end - nnb + 1):end], above[1:min(nnb, end)])
        flo = isempty(below) ? nothing : below[end]
        fhi = isempty(above) ? nothing : above[1]
        X = hcat(ones(length(nb)), R[nb])
        for x in axes(params, 1)
            c = params[x, nb, 2]
            ok = isfinite.(c)
            if count(ok) >= 2
                coef = X[ok, :] \ c[ok]
                params[x, f, 2] = coef[1] + coef[2] * R[f]
            end
            params[x, f, 3] = if !isnothing(flo) && !isnothing(fhi)
                params[x, flo, 3] +
                (params[x, fhi, 3] - params[x, flo, 3]) * (R[f] - R[flo]) / (R[fhi] - R[flo])
            elseif !isnothing(flo)
                params[x, flo, 3]
            else
                params[x, fhi, 3]
            end
            params[x, f, 1] = 0
        end
    end
    return params
end

"Reference fiber positions near the chip center (the file production identifies against)."
function trace_reference_positions(tele, chip; data_path = "./data/")
    mjd = tele == "apo" ? "60782" : "60768"
    h5open(joinpath(data_path, "traceCenterPositions_$(tele)_$(mjd)_$(chip).h5"), "r") do f
        fi = Int.(read(f["fiber_index"]))
        pos = Float64.(read(f["fiber_trace_pos_near_chip_center"]))
        o = sortperm(fi)
        return pos[o]
    end
end

# ---------------------------------------------------------------------------------------------
# QA metrics of a fitted trace (any identification mode)
# ---------------------------------------------------------------------------------------------
"""
    trace_params_qa(trace_params; nband = 32, lit_frac)

Metrics of a per-column trace (n_x, 300, 3):
- `nlit`, `lit` (median height > `lit_frac` of the median fiber);
- `fit_rms`: robust RMS of the per-fiber band-median centers (`nband` bands) about a per-fiber
  cubic in x, over lit fibers;
- `spacing_rms`, `spacing_nbad`: the fiber-to-fiber spacing pattern of each band minus its
  median over bands, after removing a running median over ±5 fiber pairs (smooth slit scale and
  curvature changes), robust RMS and count of deviations > 1 px, over adjacent lit pairs;
- `center_mid` (300): band-median centers at the central band; `width_med`, `height_med`.
"""
function trace_params_qa(trace_params::AbstractArray{<:Real, 3}; nband = 32,
        lit_frac = TRACE_LIT_FRAC)
    nx, nfib, _ = size(trace_params)
    hmed_f = [_medornan(filter(isfinite, view(trace_params, :, f, 1))) for f in 1:nfib]
    hpos = filter(v -> isfinite(v) && v > 0, hmed_f)
    ref = isempty(hpos) ? NaN : median(hpos)
    lit = [isfinite(hmed_f[f]) && hmed_f[f] > lit_frac * ref for f in 1:nfib]
    bw = nx ÷ nband
    xc = [(b - 0.5) * bw for b in 1:nband]
    u = (xc .- nx / 2) ./ (nx / 2)
    C = fill(NaN, nband, nfib)
    for b in 1:nband, f in 1:nfib
        v = filter(isfinite, view(trace_params, ((b - 1) * bw + 1):(b * bw), f, 2))
        isempty(v) || (C[b, f] = median(v))
    end
    res = Float64[]
    for f in findall(lit)
        g = isfinite.(C[:, f])
        count(g) >= 8 || continue
        pf = fit(u[g], C[g, f], 3)
        append!(res, C[g, f] .- pf.(u[g]))
    end
    fit_rms = isempty(res) ? NaN : 1.4826 * median(abs.(res .- median(res)))
    # spacing pattern stability along x
    bmid = nband ÷ 2
    pairs = [f for f in 1:(nfib - 1) if lit[f] && lit[f + 1]]
    devs = Float64[]
    function detrended_spacing(b)
        s = [C[b, f + 1] - C[b, f] for f in pairs]
        return [s[i] - median(s[max(1, i - 5):min(length(s), i + 5)]) for i in eachindex(s)]
    end
    if length(pairs) > 20
        S = reduce(hcat, [detrended_spacing(b) for b in 1:nband])   # pairs x bands
        sref = [_medornan(filter(isfinite, S[i, :])) for i in axes(S, 1)]
        for b in 1:nband
            append!(devs, filter(isfinite, S[:, b] .- sref))
        end
    end
    spacing_rms = isempty(devs) ? NaN : 1.4826 * median(abs.(devs .- median(devs)))
    spacing_nbad = count(abs.(devs) .> 1.0)
    wv = filter(isfinite, vec(trace_params[:, lit, 3]))
    hv = filter(isfinite, vec(trace_params[:, lit, 1]))
    return (nlit = count(lit), lit = lit, fit_rms = fit_rms, spacing_rms = spacing_rms,
        spacing_nbad = spacing_nbad, center_mid = C[bmid, :],
        width_med = isempty(wv) ? NaN : median(wv), height_med = isempty(hv) ? NaN : median(hv))
end

"""
    trace_light_qa(image_data, noise_image, good_pixels, trace_params, lit;
        xs = (128, 384, ..., 1025, ..., 1920), half = 32, tol = 1.0)

Does the trace sit on the light at every part of the chip? For bands of `2 half + 1` columns
centered at `xs` (every 256 columns plus the chip center, where production's center fit pins the
trace; LCO 57803 G in 0.4.0 was off by up to 2 px there only), find the fiber peaks of the light (`trace_peak_census`) and, for every lit
fiber, the distance from its trace center (band median) to the nearest peak. Returns
`light_offset_med` (median |distance|, px), `light_offset_p99`, `light_nbad` (fiber-bands
with |distance| > `tol` px) and `light_mismatch_frac` (the worst band's fraction of lit fibers
with |distance| > `tol`). A trace that is right at the chip center but walks onto the wrong
alias or twists towards the chip edges (LCO 60702 G/R in 0.4.0) fails here even though every
fiber's curve is individually smooth; so does a defocused or partially lit flat, whose light
does not resolve into fiber peaks.
"""
function trace_light_qa(image_data, noise_image, good_pixels, trace_params::AbstractArray{<:Real, 3},
        lit::AbstractVector{Bool}; xs = (128, 384, 640, 896, 1025, 1152, 1408, 1664, 1920), half = 32,
        tol = 1.0)
    nx = size(trace_params, 1)
    dall = Float64[]
    worst = 0.0
    nbad = 0
    nlit = count(lit)
    for x in xs
        cf, ce, _ = trace_center_band(image_data, noise_image, good_pixels; mid = x, n_center_cols = half)
        cen = trace_peak_census(cf, ce)
        cols = max(1, x - half):min(nx, x + half)
        nb = 0
        for f in findall(lit)
            c = _medornan(filter(isfinite, view(trace_params, cols, f, 2)))
            isfinite(c) || continue
            d = isempty(cen.y) ? Inf : minimum(abs.(cen.y .- c))
            push!(dall, d)
            d > tol && (nb += 1)
        end
        nbad += nb
        nlit > 0 && (worst = max(worst, nb / nlit))
    end
    fin = filter(isfinite, dall)
    return (light_offset_med = _medornan(fin),
        light_offset_p99 = isempty(fin) ? NaN : quantile(fin, 0.99),
        light_nbad = nbad, light_mismatch_frac = nlit == 0 ? NaN : worst)
end

"""
    trace_id_match(center_mid, lit, ident; tol = 0.35)

Fraction of lit fibers whose central-band center lies within `tol` fiber spacings of the
position the whole-slit identification of the light gives that fiber. A trace whose labels
are off for a block of fibers scores low.
"""
function trace_id_match(center_mid::AbstractVector, lit::AbstractVector{Bool}, ident; tol = 0.35)
    ident.ok || return NaN
    n = 0
    k = 0
    for f in eachindex(center_mid)
        (lit[f] && ident.matched[f] && isfinite(center_mid[f])) || continue
        n += 1
        abs(center_mid[f] - ident.positions[f]) < tol * ident.spacing && (k += 1)
    end
    return n == 0 ? NaN : k / n
end

# ---------------------------------------------------------------------------------------------
# The per-flat QA record
# ---------------------------------------------------------------------------------------------
"""
    flat_trace_qa(image_data, ivar_image, pix_bitmask, tele, chip; trace_params = nothing,
        outcome = "ok", refpos, extra = Dict())

The QA record of one flat and chip, as a `Dict{String,Any}` of scalars and strings (written
to `<flat_type>TraceQA_<tele>_<mjd>_<exp>_<chip>.h5`). Computed for EVERY flat, including the
ones the trace stage refuses: the refusal is the record's `outcome`, never a silent drop.
Includes the lamp metrics, the saturation fraction, the template-free census of the center
band, the whole-slit identification of the light, and, when a trace was fitted, its
`trace_params_qa` metrics and the agreement of its labels with the identification. The gate
verdict (`trace_gate`) is added by the caller.
"""
function flat_trace_qa(image_data, ivar_image, pix_bitmask, tele, chip; trace_params = nothing,
        outcome = "ok", refpos, mid = 1025, n_center_cols = 100,
        extra = Dict{String, Any}())
    qa = Dict{String, Any}()
    lm = trace_lamp_metrics(image_data, pix_bitmask, tele)
    qa["flux_med"] = lm.flux_med
    qa["flux_68p"] = lm.flux_68p
    qa["flux_68p_thresh"] = lm.thresh
    qa["lamp_ok"] = lm.lamp_ok
    qa["frac_good_pix"] = lm.frac_good
    good = ((pix_bitmask .& bad_pix_bits) .== 0) .& (ivar_image .> 0) .& isfinite.(image_data)
    qa["sat_frac"] = count((pix_bitmask .& bad_fully_saturated) .!= 0) / length(pix_bitmask)
    noise = 1 ./ sqrt.(ivar_image)
    cf, ce, _ = trace_center_band(image_data, noise, good; mid = mid, n_center_cols = n_center_cols)
    cen = trace_peak_census(cf, ce)
    qa["census_npeaks"] = length(cen.y)
    qa["census_nslots"] = cen.nslots
    qa["census_nmissing"] = cen.nmissing
    qa["census_spacing"] = cen.spacing
    qa["census_hmed"] = cen.hmed
    qa["census_median_snr"] = cen.median_snr
    id = trace_identify_slit(cen.y, refpos; census_nslots = cen.nslots, census_slots = cen.slots)
    qa["id_reason"] = id.reason
    qa["id_ok"] = id.ok
    qa["id_shift_px"] = id.a
    qa["id_scale"] = id.b
    qa["id_margin"] = id.margin
    qa["id_n_unmatched"] = id.n_unmatched
    qa["id_n_outside"] = id.n_outside
    qa["id_n_dup"] = id.n_dup
    qa["id_resid_rms"] = id.resid_rms
    qa["id_nlit"] = count(id.matched)
    if isnothing(trace_params) || isempty(trace_params)
        qa["has_trace"] = false
        for k in ("nlit", "spacing_nbad")
            qa[k] = -1
        end
        for k in ("fit_rms", "spacing_rms", "width_med", "height_med", "id_match_frac", "trace_shift_px",
            "light_offset_med", "light_offset_p99", "light_mismatch_frac")
            qa[k] = NaN
        end
        qa["light_nbad"] = -1
    else
        tq = trace_params_qa(trace_params)
        qa["has_trace"] = true
        qa["nlit"] = tq.nlit
        qa["fit_rms"] = tq.fit_rms
        qa["spacing_rms"] = tq.spacing_rms
        qa["spacing_nbad"] = tq.spacing_nbad
        qa["width_med"] = tq.width_med
        qa["height_med"] = tq.height_med
        qa["id_match_frac"] = trace_id_match(tq.center_mid, tq.lit, id)
        # rigid position of the trace at the chip center relative to the reference positions
        qa["trace_shift_px"] = _medornan(filter(isfinite, (tq.center_mid .- refpos)[tq.lit]))
        lq = trace_light_qa(image_data, noise, good, trace_params, tq.lit)
        for (k, v) in pairs(lq)
            qa[String(k)] = v
        end
    end
    qa["outcome"] = outcome
    qa["tele"] = String(tele)
    qa["chip"] = String(chip)
    for (k, v) in extra
        qa[k] = v
    end
    lvl, reasons = trace_gate(qa)
    qa["gate_level"] = lvl
    qa["gate_reasons"] = join(reasons, ";")
    return qa
end

"""
    trace_gate(qa) -> (gate_level, reasons)

`gate_level` 2 = passes every criterion; 1 = fails only soft criteria (saturation), usable
when no flat of the night passes fully, flagged `TRACE_UNGATED`; 0 = refused. `reasons` names
every failed criterion (empty when the level is 2).
"""
function trace_gate(qa::AbstractDict)
    hard = String[]
    soft = String[]
    g(k, d = NaN) = get(qa, k, d)
    get(qa, "outcome", "ok") == "ok" || push!(hard, "outcome:" * string(g("outcome")))
    g("lamp_ok", false) || push!(hard, "lamp_off")
    if get(qa, "has_trace", false)
        g("nlit", 0) >= TRACE_GATE_NLIT_MIN || push!(hard, "nlit")
        (isfinite(g("fit_rms")) && g("fit_rms") <= TRACE_GATE_FIT_RMS_MAX) || push!(hard, "fit_rms")
        (isfinite(g("light_mismatch_frac")) && g("light_mismatch_frac") <= TRACE_GATE_LIGHT_MISMATCH_MAX) ||
            push!(hard, "light_mismatch")
        (isfinite(g("id_match_frac")) && g("id_match_frac") >= TRACE_GATE_ID_MATCH_MIN) ||
            push!(hard, "id_match")
    else
        push!(hard, "no_trace")
    end
    g("id_ok", false) || push!(hard, "id:" * string(g("id_reason", "")))
    (isfinite(g("sat_frac")) && g("sat_frac") <= TRACE_GATE_SAT_FRAC_MAX) || push!(soft, "saturated")
    lvl = !isempty(hard) ? 0 : (!isempty(soft) ? 1 : 2)
    return lvl, vcat(hard, soft)
end

# ---------------------------------------------------------------------------------------------
# QA record I/O (one small file per flat and chip, written atomically)
# ---------------------------------------------------------------------------------------------
trace_qa_path(trace_dir, flat_type, tele, mjd, expnum, chip) = joinpath(trace_dir,
    "$(flat_type)_flats", "$(mjd)", "$(flat_type)TraceQA_$(tele)_$(mjd)_$(expnum)_$(chip).h5")

"Write a QA record to `path` (temporary file + rename, so readers never see a partial file)."
function write_trace_qa(path::AbstractString, qa::AbstractDict)
    mkpath(dirname(path))
    tmp = path * ".tmp$(getpid())"
    md = Dict{String, Any}(String(k) => v for (k, v) in qa)
    safe_jldsave(tmp, md; trace_qa_version = 1)
    mv(tmp, path; force = true)
    return path
end

"Read a QA record (the metadata group), or `nothing` if the file does not exist."
function read_trace_qa(path::AbstractString)
    isfile(path) || return nothing
    return read_metadata(path)
end

# ---------------------------------------------------------------------------------------------
# Choice of the night's trace
# ---------------------------------------------------------------------------------------------
"""
    TraceCandidate

One per-flat trace of the night: `flat_type` ("quartz"/"dome"), exposure number, the trace
file, and its QA record (`nothing` if none exists).
"""
struct TraceCandidate
    flat_type::String
    expnum::Int
    trace_file::String
    qa::Union{Nothing, Dict{String, Any}}
end

"""
    select_night_trace(cands; mode = "gated", type_order = ["quartz","dome"], sci_mid = nothing)

Choose the night's trace among `cands` (all flats of the night and chip that have a trace).

- `mode = "last"`: the pre-existing rule (`--trace_select last`), bit for bit: the first `flat_type` in `type_order`
  that has any trace, and of those the lexically last file name.
- `mode = "gated"` (default): among candidates with gate level >= 1, rank by (a) gate level (2 before
  1), (b) `type_order`, (c) |exposure - median science exposure of the night| (`sci_mid`;
  ignored if `nothing`), (d) `fit_rms`, (e) larger exposure number. Deterministic: no file
  order, no dependence on what else has run. Candidates without a QA record are not eligible.

Returns `(chosen::Union{Nothing,TraceCandidate}, info::Dict)`; `info` holds the rule, the
ranked candidate table and the chosen candidate's QA summary, for the traceMain metadata.
"""
function select_night_trace(cands::AbstractVector{TraceCandidate}; mode = "gated",
        type_order = ["quartz", "dome"], sci_mid = nothing)
    info = Dict{String, Any}("trace_select_rule" => mode, "trace_n_candidates" => length(cands))
    tord(c) = (i = findfirst(==(c.flat_type), type_order); isnothing(i) ? length(type_order) + 1 : i)
    qa_lvl(c) = isnothing(c.qa) ? -1 : Int(get(c.qa, "gate_level", 0))
    qa_rms(c) = isnothing(c.qa) ? Inf : (v = Float64(get(c.qa, "fit_rms", Inf)); isfinite(v) ? v : Inf)
    ranked = TraceCandidate[]
    chosen = nothing
    if mode == "last"
        for t in type_order
            tl = sort(filter(c -> c.flat_type == t, cands), by = c -> basename(c.trace_file))
            if !isempty(tl)
                chosen = tl[end]
                break
            end
        end
        ranked = sort(cands, by = c -> (tord(c), basename(c.trace_file)), rev = false)
    elseif mode == "gated"
        dist(c) = isnothing(sci_mid) ? 0.0 : abs(c.expnum - sci_mid)
        key(c) = (-max(qa_lvl(c), 0), tord(c), dist(c), qa_rms(c), -c.expnum)
        ranked = sort(cands, by = key)
        el = filter(c -> qa_lvl(c) >= 1, ranked)
        chosen = isempty(el) ? nothing : el[1]
    else
        error("unknown trace selection mode $(mode) (expected \"last\" or \"gated\")")
    end
    info["trace_cand_files"] = String[basename(c.trace_file) for c in ranked]
    info["trace_cand_gate_level"] = Int[qa_lvl(c) for c in ranked]
    info["trace_cand_gate_reasons"] = String[isnothing(c.qa) ? "qa_missing" :
                                       string(get(c.qa, "gate_reasons", "")) for c in ranked]
    info["trace_cand_fit_rms"] = Float64[qa_rms(c) for c in ranked]
    info["trace_sci_mid_exp"] = isnothing(sci_mid) ? -1.0 : Float64(sci_mid)
    return chosen, info
end

"Median exposure number of the night's object exposures (nothing if none)."
function night_science_mid_exposure(dfalmanac)
    hasproperty(dfalmanac, :image_type) || return nothing
    idx = findall(==("object"), String.(dfalmanac.image_type))
    isempty(idx) && return nothing
    return median(idx)
end

"""
    trace_fiber_bits(trace_params, level; other_flat = false, ungated = false,
        qa_missing = false, lit_frac)

The per-fiber `bitmsk_trace` (length n_fibers) of a night's trace: the level bit and the
night-level warn bits in every fiber, plus `TRACE_FIBER_UNLIT` where the trace source's
median height is below `lit_frac` of the median fiber. For the fixed fallback file the
"unlit" fibers are those of the file's own night, so `TRACE_FIBER_UNLIT` is not set there
(the whole night is already `TRACE_FIXED_FALLBACK`).

With templates (PR B): `fiber_delta` (per fiber, px; from a template-seeded flat fit) sets
`TRACE_FIBER_DELTA` above `delta_px` (0.3 px); `unlit_mask` (levels T1-T3: the fibers dark on
every flat of the epoch) sets `TRACE_FIBER_UNLIT`; `provisional` sets
`TRACE_EPOCH_PROVISIONAL` in every fiber.
"""
function trace_fiber_bits(trace_params::AbstractArray{<:Real, 3}, level::Integer;
        other_flat = false, ungated = false, qa_missing = false, lit_frac = TRACE_LIT_FRAC,
        fiber_delta = nothing, unlit_mask = nothing, provisional = false,
        delta_px = TRACE_FIBER_DELTA_PX)
    nfib = size(trace_params, 2)
    night = trace_level_bits(level)
    other_flat && (night |= TRACE_OTHER_FLAT)
    ungated && (night |= TRACE_UNGATED)
    qa_missing && (night |= TRACE_QA_MISSING)
    provisional && (night |= TRACE_EPOCH_PROVISIONAL)
    bits = fill(night, nfib)
    if !isnothing(fiber_delta)
        for f in 1:nfib
            (isfinite(fiber_delta[f]) && fiber_delta[f] > delta_px) && (bits[f] |= TRACE_FIBER_DELTA)
        end
    end
    if !isnothing(unlit_mask)
        for f in 1:nfib
            unlit_mask[f] && (bits[f] |= TRACE_FIBER_UNLIT)
        end
    end
    if level == TRACE_LEVEL_OWN_FLAT
        hmed_f = [_medornan(filter(isfinite, view(trace_params, :, f, 1))) for f in 1:nfib]
        hpos = filter(v -> isfinite(v) && v > 0, hmed_f)
        ref = isempty(hpos) ? NaN : median(hpos)
        for f in 1:nfib
            (isfinite(hmed_f[f]) && hmed_f[f] > lit_frac * ref) || (bits[f] |= TRACE_FIBER_UNLIT)
        end
    end
    return bits
end

"""
    trace_bits_from_metadata(trace_params, trace_metadata)

`bitmsk_trace` for a traceMain that predates the QA records (no stored `bitmsk_trace`):
the level from `trace_type` ("*_fallback" = fixed file), unlit fibers from the heights,
`TRACE_QA_MISSING` set.
"""
function trace_bits_from_metadata(trace_params, trace_metadata::AbstractDict)
    tt = string(get(trace_metadata, "trace_type", ""))
    level = endswith(tt, "_fallback") ? TRACE_LEVEL_FIXED_FALLBACK : TRACE_LEVEL_OWN_FLAT
    return trace_fiber_bits(trace_params, level; other_flat = (tt == "dome"), qa_missing = true)
end

# ---------------------------------------------------------------------------------------------
# traceMain: the night's trace for one chip
# ---------------------------------------------------------------------------------------------
"""
    night_trace_candidates(outdir, tele, mjd, chip; type_order = ["quartz","dome"])

Every per-flat trace file of the night and chip (`<type>_flats/<mjd>/<type>Trace_*.h5`) with
its QA record (or `nothing`).
"""
function night_trace_candidates(outdir, tele, mjd, chip; type_order = ["quartz", "dome"])
    cands = TraceCandidate[]
    for t in type_order
        d = joinpath(outdir, "$(t)_flats", "$(mjd)")
        isdir(d) || continue
        rx = Regex("^$(t)Trace_$(tele)_$(mjd)_(\\d+)_$(chip)\\.h5\$")
        for fn in sort(readdir(d))
            m = match(rx, fn)
            isnothing(m) && continue
            e = m.captures[1]
            qa = read_trace_qa(trace_qa_path(outdir, t, tele, mjd, e, chip))
            push!(cands, TraceCandidate(t, parse(Int, e), joinpath(d, fn), qa))
        end
    end
    return cands
end

"""
    make_traceMain(outdir, tele, mjd, chip; select = "gated", fallback = "none",
        fallback_file = "", type_order = ["quartz","dome"], almanac_file = "",
        checkpoint_mode = "commit_same", template_cal = nothing)

Write `apred/<mjd>/traceMain_<tele>_<mjd>_<chip>.h5`, the trace every exposure of the night
is extracted with: a copy of the chosen per-flat trace (or of `fallback_file`) plus
- metadata: `trace_type`, `trace_found_match`, `trace_orig_param_fname` (as before) and
  `trace_level`, `trace_flags` (night-level bits), `trace_select_rule`,
  `trace_fallback_rule`, the chosen flat's gate level, reasons and QA metrics
  (`trace_qa_*`), and the ranked candidate table (`trace_cand_*`);
- dataset `bitmsk_trace` (per fiber, `trace_fiber_bits`).

`select = "gated"` (default) uses `select_night_trace` and re-decides when the night's
candidates changed (e.g. dome traces made after the quartz pass of a bulk run);
`select = "last"` reproduces the previous choice exactly (and the previous checkpoint
behaviour). `fallback = "none"` (default) writes no traceMain when the night has no eligible
trace (and removes a stale one): the night's exposures are then not extracted, and the caller
logs them; `fallback_file` is not needed and not read. `fallback = "fixed_file"` (previous
behaviour) copies `fallback_file`. Returns the trace level.

`template_cal` (a `TraceTemplateCal` from `--trace_template_cal`, or `nothing` = off): its
identity for this (telescope, chip, night) (`trace_template_metadata`, src/trace_template_cal.jl)
is written into the metadata, and a checkpointed traceMain whose recorded `trace_template_id`
differs from the configured one (a traceMain without the key counts as "none") is recomputed,
under either selection rule.
"""
function make_traceMain(outdir, tele, mjd, chip; select = "gated", fallback = "none",
        fallback_file = "", type_order = ["quartz", "dome"], almanac_file = "",
        checkpoint_mode = "commit_same", template_cal = nothing, data_path = "./data/",
        fallback_kw = (;))
    select in ("last", "gated") || error("--trace_select must be last or gated, got $(select)")
    fb_levels, fb_fixed = parse_trace_fallback(fallback)
    if !isempty(fb_levels) && isnothing(template_cal)
        error("--trace_fallback $(fallback) needs the epoch templates (--trace_template_cal), which are off")
    end
    traceMain_fname = joinpath(outdir, "apred", "$(mjd)", "traceMain_$(tele)_$(mjd)_$(chip).h5")
    exists = check_file(traceMain_fname, mode = checkpoint_mode)
    # a traceMain made with another template (or without one) is recomputed
    if exists && template_needs_recompute(read_metadata(traceMain_fname), template_cal, tele, chip, mjd)
        exists = false
    end
    if select == "last" && exists
        return get(read_metadata(traceMain_fname), "trace_level", TRACE_LEVEL_OWN_FLAT)
    end
    cands = night_trace_candidates(outdir, tele, mjd, chip; type_order = type_order)
    sci_mid = nothing
    if select == "gated" && almanac_file != "" && isfile(almanac_file)
        sci_mid = try
            night_science_mid_exposure(read_almanac_exp_df(almanac_file, tele, mjd))
        catch
            nothing
        end
    end
    chosen, info = select_night_trace(cands; mode = select, type_order = type_order,
        sci_mid = sci_mid)
    # own-night levels T1 > T2 > T3 (src/trace_levels.jl) when no flat of the night is eligible
    fb = nothing
    if isnothing(chosen) && !isempty(fb_levels)
        fb = trace_fallback_chain(outdir, tele, mjd, chip; levels = fb_levels,
            template_cal = template_cal, almanac_file = almanac_file, data_path = data_path,
            checkpoint_mode = checkpoint_mode, fallback_kw...)
    end
    want = !isnothing(chosen) ? chosen.trace_file : !isnothing(fb) ? fb.source :
           (fb_fixed ? fallback_file : "")
    if exists && select == "gated"
        md = read_metadata(traceMain_fname)
        if get(md, "trace_select_rule", "") == "gated" &&
           get(md, "trace_orig_param_fname", "") == want &&
           get(md, "trace_fallback_rule", "") == fallback
            return get(md, "trace_level", TRACE_LEVEL_OWN_FLAT)
        end
    end
    if isnothing(chosen) && isnothing(fb) && !fb_fixed
        @warn "No trace for $(tele) $(mjd) chip $(chip): no flat of the night is eligible ($(length(cands)) candidate(s)) and --trace_fallback $(fallback) gave no trace; the night's exposures will not be extracted."
        isfile(traceMain_fname) && rm(traceMain_fname)
        return TRACE_LEVEL_NONE
    end
    md = Dict{String, Any}()
    if !isnothing(chosen)
        level = TRACE_LEVEL_OWN_FLAT
        src = chosen.trace_file
        md["trace_type"] = chosen.flat_type
        md["trace_found_match"] = true
        qa = chosen.qa
    elseif !isnothing(fb)
        level = fb.level
        src = fb.source
        md["trace_type"] = level == TRACE_LEVEL_REGISTERED_CAL ? "registered_cal" :
                           level == TRACE_LEVEL_REGISTERED_SCI ? "registered_sci" : "template"
        md["trace_found_match"] = false
        qa = fb.qa
        merge!(md, fb.info)
        md["trace_reg"] = reg_vector(fb.reg)
        md["trace_reg_source"] = level == TRACE_LEVEL_TEMPLATE ? "epoch_median" : string(qa["kind"])
    else
        isfile(fallback_file) || error("Fallback trace file $(fallback_file) does not exist")
        level = TRACE_LEVEL_FIXED_FALLBACK
        src = fallback_file
        md["trace_type"] = "quartz_fallback"
        md["trace_found_match"] = false
        qa = nothing
    end
    md["trace_orig_param_fname"] = src
    md["trace_level"] = level
    md["trace_fallback_rule"] = fallback
    gl = (isnothing(qa) || !isnothing(fb)) ? -1 : Int(get(qa, "gate_level", 0))
    md["trace_gate_level"] = gl
    md["trace_gate_reasons"] = isnothing(qa) ? (level == TRACE_LEVEL_OWN_FLAT ? "qa_missing" : "") :
                               !isnothing(fb) ? string(get(qa, "reg_gate_reasons", "")) :
                               string(get(qa, "gate_reasons", ""))
    for k in ("fit_rms", "light_mismatch_frac", "spacing_rms", "id_match_frac", "sat_frac", "id_margin",
        "id_shift_px", "trace_shift_px")
        md["trace_qa_" * k] = isnothing(qa) ? NaN : Float64(get(qa, k, NaN))
    end
    for k in ("nlit", "spacing_nbad")
        md["trace_qa_" * k] = isnothing(qa) ? -1 : Int(get(qa, k, -1))
    end
    md["trace_identify"] = isnothing(qa) ? "" : string(get(qa, "identify", ""))
    merge!(md, info)
    merge!(md, trace_template_metadata(template_cal, tele, chip, mjd))
    own = level == TRACE_LEVEL_OWN_FLAT
    mkpath(dirname(traceMain_fname))
    tmp = traceMain_fname * ".tmp$(getpid())"
    fiber_delta = nothing
    unlit_mask = nothing
    if isnothing(fb)
        trace_params, fiber_delta = h5open(src, "r") do f
            read(f["trace_params"]),
            haskey(f, "template_fiber_delta_px") ? read(f["template_fiber_delta_px"]) : nothing
        end
        cp(src, tmp, force = true)
    else
        # T1-T3: the registered template (no per-fiber delta); epoch-dark fibers flagged unlit
        trace_params = registered_trace_params(fb.tmpl, fb.reg)
        unlit_mask = .!fb.tmpl.lit
        isfile(tmp) && rm(tmp)
        safe_jldsave(tmp; trace_params = trace_params,
            regularized_trace_params = regularize_trace(trace_params), no_metadata = true)
    end
    provisional = (own && !isnothing(fiber_delta) || !isnothing(fb)) &&
                  trace_template_provisional(template_cal, tele, chip, parse(Int, string(mjd)))
    bits = trace_fiber_bits(trace_params, level;
        other_flat = own && md["trace_type"] != first(type_order),
        ungated = own && gl in (0, 1), qa_missing = own && gl == -1,
        fiber_delta = fiber_delta, unlit_mask = unlit_mask, provisional = provisional)
    md["trace_flags"] = reduce(|, bits; init = zero(TraceBits)) & ~TRACE_PER_FIBER_BITS
    md["trace_n_unlit"] = count(b -> (b & TRACE_FIBER_UNLIT) != 0, bits)
    md["trace_n_fiber_delta"] = count(b -> (b & TRACE_FIBER_DELTA) != 0, bits)
    chmod(tmp, 0o644)
    modify_saved_dict(tmp, "metadata", md)
    h5open(tmp, "r+") do f
        haskey(f, "bitmsk_trace") && delete_object(f, "bitmsk_trace")
        f["bitmsk_trace"] = bits
    end
    mv(tmp, traceMain_fname, force = true)
    if level != TRACE_LEVEL_OWN_FLAT || gl in (0, 1)
        @warn "traceMain $(tele) $(mjd) $(chip): level $(level), flags $(join(trace_bit_names(md["trace_flags"]), ",")) ($(md["trace_gate_reasons"]))"
    end
    return level
end

# ---------------------------------------------------------------------------------------------
# 1D products
# ---------------------------------------------------------------------------------------------
"""
    trace_bits_for_1d(traceFname, regularized_trace_params, trace_metadata) -> (bits, md)

The per-fiber `bitmsk_trace` an ar1D product of this traceMain carries: the stored dataset
(written by `make_traceMain`), or, for a traceMain made before that, `trace_bits_from_metadata`
plus the `trace_flags` / `trace_level` metadata it implies (`md`, empty otherwise: a new
traceMain's metadata already holds them).
"""
function trace_bits_for_1d(traceFname, regularized_trace_params, trace_metadata::AbstractDict)
    bits = h5open(traceFname, "r") do f
        haskey(f, "bitmsk_trace") ? to_trace_bits(read(f["bitmsk_trace"])) : nothing
    end
    md = Dict{String, Any}()
    if isnothing(bits)
        bits = trace_bits_from_metadata(regularized_trace_params, trace_metadata)
        md["trace_flags"] = reduce(|, bits; init = zero(TraceBits)) & ~TRACE_PER_FIBER_BITS
        md["trace_level"] = endswith(string(get(trace_metadata, "trace_type", "")), "_fallback") ?
                            TRACE_LEVEL_FIXED_FALLBACK : TRACE_LEVEL_OWN_FLAT
    end
    return bits, md
end

"""
    trace_exposure_summary(bitmsk_trace_stack, trace_level_chips) -> Dict

Exposure-level trace provenance of a chip-combined (ar1Duni*) product: `trace_flags` (OR of the
night-level bits over chips and fibers), `trace_level_chips` and `trace_level_max`.
"""
function trace_exposure_summary(bitmsk_trace_stack::AbstractMatrix{<:Integer},
        trace_level_chips::AbstractVector{<:Integer})
    return Dict{String, Any}(
        "trace_flags" => reduce(|, to_trace_bits(bitmsk_trace_stack); init = zero(TraceBits)) &
                         ~TRACE_PER_FIBER_BITS,
        "trace_level_chips" => collect(Int, trace_level_chips),
        "trace_level_max" => maximum(trace_level_chips))
end
