# Arc-lamp (ThAr + UNe) wavelength calibration, stage 1: frozen line list, local-noise peak
# finder, association, joint night fit with measured per-sequence positions, and arc-only
# outlier rejection of sequences. Writes one per-night product (`arclampWaveCalNight_<tele>_<mjd>.h5`)
# that NOTHING in the pipeline reads yet. It runs only when `pipeline_2d_1d.jl --arclamp_wavecal
# true` is given (default false), so production results are unchanged.
#
# Provenance of the methods (offline studies, 2026-11-16/17):
#   * finder: the "N=8" local-noise matched-filter finder (`get_arclamp_peaks_local`). It works in
#     detector pixels only and uses the production `fit_gauss_and_bias` for the centroids.
#   * association: catalog scan + quadratic self-calibration (the `associate` of the studies).
#   * joint night fit (Gauss-Newton with a Schur complement) and position clustering: the best
#     arc fit of the ARC_VERDICT2 study, validated there on 1,563 nights.
#   * sequence rejection: arc data only (arcs come first in the wavelength chain; nothing from
#     the sky or FPI stages is read). Study 2026_11_18/arc_gate replaced the earlier gate against
#     the night's object-exposure solutions.
#
# Positions are MEASURED from each sequence's own lines. The recorded dither (DITHPIX) is used for one thing only: to pair the
# ThAr and UNe exposures of an almanac arc sequence that sit at the same commanded position.
# A wrong label there is caught by the measured ThAr-vs-UNe position check (`ARCSEQ_LAMP_MISMATCH`).
#
# Model, per fiber f (green-chip-relative coordinate x = (pix - 1024)/2048):
#   x'_R = a_R + b_R x,  x'_G = x,  x'_B = a_B + b_B x          (chip transform)
#   lambda = p0 + p1 z + p2 z^2 + p3 z^3 + p4 z^4,  z = x'_c + shift(s, f, c)
#   shift(s, f, c) = d_s + t_s phi(f)                         (sequence s; s = 1 is the reference)
#                  + [position 2] (alpha_c + beta_c x)          (second measured position cluster)
#                  + tau_s (u_c + v_c x)                        (linear in time, per chip)
# alpha_G = u_G = 0 (degenerate with d_s). The chip scales b_R, b_B are either fitted per fiber
# (no calibration supplied) or held fixed at the comb-measured values (`chip_scales`, the
# N_FIBERS x 2 matrix returned by the chip-transform calibration lookup of PR #415).

using HDF5
using CSV
using DataFrames
using LinearAlgebra
using Statistics: median, mean, quantile

# ------------------------------------------------------------------------------------------
# small helpers (prefixed to avoid clashing with other files of the package)
# ------------------------------------------------------------------------------------------
_arc_finite(x) = filter(isfinite, x)
_arc_nanmedian(x) = (y = _arc_finite(x); isempty(y) ? NaN : median(y))
_arc_nanmean(x) = (y = _arc_finite(x); isempty(y) ? NaN : mean(y))
"Half the 16-84 percentile range of the finite entries (NaN if fewer than 3)."
function _arc_robstd(x)
    y = _arc_finite(x)
    length(y) < 3 && return NaN
    q = quantile(y, [0.16, 0.84])
    return (q[2] - q[1]) / 2
end
_arc_ximport(pix) = (pix .- (N_XPIX ÷ 2)) ./ N_XPIX
_arc_expstr(e::Integer) = lpad(string(e % 10000), 4, '0')

"FNV-1a 64-bit hash of a file's bytes (16 hex digits); identifies the exact line list used."
function arclamp_file_fingerprint(path::AbstractString)
    h = 0xcbf29ce484222325
    for b in read(path)
        h = (h ⊻ UInt64(b)) * 0x00000100000001b3
    end
    return string(h, base = 16, pad = 16)
end

# ------------------------------------------------------------------------------------------
# local-noise arc peak finder (pixels only; no wavelength enters)
# ------------------------------------------------------------------------------------------
const ARCPK_BIT_PARTIAL_SAT = bad_partial_saturated
const ARCPK_SAT_BITS = bad_fully_saturated

"""Default parameters of the local-noise arc peak finder (all overridable by keyword).
They are generic pixel-space constants: detection threshold, matched-filter width (about the
APOGEE LSF sigma), background/noise windows, plausible LSF sigma range, and sanity cuts."""
const ARCPK_DEFAULTS = (
    nsig = 8.0,            # detection: matched-filter S/N >= nsig against the local noise
    kern_sigma = 1.1,      # px, matched-filter gaussian
    kern_hw = 4,
    bg_hw = 30, noise_hw = 60,
    edge = 16,             # px from the array ends excluded (production n_offset + 1)
    sig_lo = 0.6, sig_hi = 2.2,
    max_shift = 1.5,       # px, |fitted centre - detection pixel|
    max_cen_err = 0.5,     # px, 1-sigma centroid error
    max_masked_core = 0,   # masked pixels allowed within +-2 px of the centre
    min_fit_snr = 0.0,
    dup_sep = 0.6,         # px
    sat_pad = 4            # px, exclusion around a fully saturated run
)

const ARCPK_F_SIGMA = 1
const ARCPK_F_CLAMP = 2
const ARCPK_F_SHIFT = 4
const ARCPK_F_COV = 8
const ARCPK_F_MASKED = 16
const ARCPK_F_LOWSNR = 32
const ARCPK_F_DUP = 64
const ARCPK_F_SATNB = 128

"""Running statistic `f` of `y` over pixels where `use` is true, window half-width `hw`,
evaluated every `step` pixels and linearly interpolated; grid points with fewer than `min_n`
samples are filled from the nearest valid grid point."""
function arcpk_running_stat(y::AbstractVector{Float64}, use::AbstractVector{Bool}, hw::Int, f;
        step::Int = 4, min_n::Int = 8)
    n = length(y)
    grid = collect(1:step:n)
    grid[end] != n && push!(grid, n)
    vals = fill(NaN, length(grid))
    buf = Float64[]
    for (k, g) in enumerate(grid)
        empty!(buf)
        @inbounds for i in max(1, g - hw):min(n, g + hw)
            use[i] && push!(buf, y[i])
        end
        length(buf) >= min_n && (vals[k] = f(buf))
    end
    all(isnan, vals) && return fill(NaN, n)
    good = findall(!isnan, vals)
    for k in eachindex(vals)
        isnan(vals[k]) || continue
        j = good[argmin(abs.(good .- k))]
        vals[k] = vals[j]
    end
    out = Vector{Float64}(undef, n)
    k = 1
    @inbounds for i in 1:n
        while k < length(grid) - 1 && grid[k + 1] < i
            k += 1
        end
        g0, g1 = grid[k], grid[k + 1]
        t = g1 == g0 ? 0.0 : (i - g0) / (g1 - g0)
        out[i] = (1 - t) * vals[k] + t * vals[k + 1]
    end
    return out
end

"One-sided (lower-half) robust sigma: median minus the 15.87th percentile."
arcpk_lower_sigma(v) = median(v) - quantile(v, 0.1587)

"""Local background (iterated running median with lines masked) and local noise (ivar-based
sigma times a local robust excess factor that is never below 1) of one 1D arc spectrum."""
function arcpk_background_noise(flux, ivar, good; bg_hw = 30, noise_hw = 60)
    n = length(flux)
    sig_ivar = [good[i] && ivar[i] > 0 ? 1 / sqrt(ivar[i]) : Inf for i in 1:n]
    bg1 = arcpk_running_stat(flux, good, bg_hw, median)
    z1 = (flux .- bg1) ./ sig_ivar
    line = falses(n)
    for i in 1:n
        if good[i] && z1[i] > 3
            line[max(1, i - 3):min(n, i + 3)] .= true
        end
    end
    use = good .& .!line
    bg = arcpk_running_stat(flux, use, bg_hw, median)
    any(isnan, bg) && (bg = bg1)
    any(isnan, bg) && (bg = zeros(n))
    z = (flux .- bg) ./ sig_ivar
    s = arcpk_running_stat(z, use .& isfinite.(z), noise_hw, arcpk_lower_sigma; min_n = 20)
    s[isnan.(s)] .= 1.0
    excess = max.(s, 1.0)
    return bg, sig_ivar .* excess, excess, line
end

"Post-fit sanity cuts on the candidate fits; returns a flag vector (0 = accepted)."
function arcpk_candidate_flags(params, cov, det_pix, n_masked_core, near_sat; p = ARCPK_DEFAULTS)
    n = size(params, 1)
    fl = zeros(Int, n)
    for k in 1:n
        s = params[k, 3]
        (s < p.sig_lo || s > p.sig_hi) && (fl[k] |= ARCPK_F_SIGMA)
        (s <= 0.1 + 1e-9 || s >= 3.0 - 1e-9) && (fl[k] |= ARCPK_F_CLAMP)
        abs(params[k, 2] - det_pix[k]) > p.max_shift && (fl[k] |= ARCPK_F_SHIFT)
        v = cov[k, 2, 2]
        if !(all(isfinite, @view cov[k, :, :])) || !(v > 0) || sqrt(v) > p.max_cen_err ||
           !(cov[k, 1, 1] > 0)
            fl[k] |= ARCPK_F_COV
        end
        n_masked_core[k] > p.max_masked_core && (fl[k] |= ARCPK_F_MASKED)
        if p.min_fit_snr > 0 && cov[k, 1, 1] > 0
            params[k, 1] / sqrt(cov[k, 1, 1]) < p.min_fit_snr && (fl[k] |= ARCPK_F_LOWSNR)
        end
        near_sat[k] && (fl[k] |= ARCPK_F_SATNB)
    end
    ok = findall(==(0), fl)
    o = ok[sortperm(params[ok, 2])]
    for j in 2:length(o)
        a, b = o[j - 1], o[j]
        (fl[a] != 0 || fl[b] != 0) && continue
        if abs(params[b, 2] - params[a, 2]) < p.dup_sep
            if params[a, 1] >= params[b, 1]
                fl[b] |= ARCPK_F_DUP
            else
                fl[a] |= ARCPK_F_DUP
            end
        end
    end
    return fl
end

"Contiguous runs of fully saturated pixels: vector of (lo, hi, centre)."
function arcpk_saturated_runs(mask)
    n = length(mask)
    runs = Tuple{Int, Int, Float64}[]
    i = 1
    while i <= n
        if (mask[i] & ARCPK_SAT_BITS) != 0
            j = i
            while j < n && (mask[j + 1] & ARCPK_SAT_BITS) != 0
                j += 1
            end
            push!(runs, (i, j, (i + j) / 2))
            i = j + 1
        else
            i += 1
        end
    end
    return runs
end

"""
    get_arclamp_peaks_local(flux, ivar, mask; kw...) -> (inds, params, cov, info)

Local-noise arc peak finder for one fiber's 1D spectrum. Detects local maxima of a matched-filter
S/N (against a local background and a local noise) above `nsig`, fits them with the production
gaussian + bias model (`fit_gauss_and_bias`, neighbour subtraction), applies sanity cuts, and
returns the ACCEPTED peaks sorted by centre. `params[k, :]` = (integrated height, centre in
1-based pixels, sigma in pixels, bias), as the production finder. Fully saturated lines are
flagged (`info.sat`) and never fitted. Nothing here knows a wavelength.
"""
function get_arclamp_peaks_local(flux::AbstractVector, ivar::AbstractVector,
        mask::AbstractVector; kw...)
    p = merge(ARCPK_DEFAULTS, (; kw...))
    n = length(flux)
    flux = Float64.(flux)
    ivar = Float64.(ivar)
    good = ((mask .& bad_pix_bits) .== 0) .& (ivar .> 0) .& isfinite.(flux) .& isfinite.(ivar)
    sat = (mask .& ARCPK_SAT_BITS) .!= 0
    good .&= .!sat
    psat = (mask .& ARCPK_BIT_PARTIAL_SAT) .!= 0
    empty = ([], zeros(0, 4), zeros(0, 4, 4),
        (; snr = Float64[], det_pix = Int[], cand_params = zeros(0, 4),
            cand_cov = zeros(0, 4, 4), cand_snr = Float64[], cand_det = Int[],
            cand_nmask = Int[], cand_nearsat = Bool[], cand_flags = Int[],
            cand_npsat = Int[], sat = Tuple{Int, Int, Float64}[], noise_excess = Float64[]))
    count(good) < 100 && return empty

    bg, sig, excess, _ = arcpk_background_noise(flux, ivar, good;
        bg_hw = p.bg_hw, noise_hw = p.noise_hw)
    w = [good[i] && isfinite(sig[i]) ? 1 / sig[i]^2 : 0.0 for i in 1:n]
    d = flux .- bg

    kr = (-p.kern_hw):(p.kern_hw)
    g = exp.(-0.5 .* (kr ./ p.kern_sigma) .^ 2)
    snr = zeros(n)
    @inbounds for i in 1:n
        num = 0.0
        den = 0.0
        for (jj, j) in enumerate(kr)
            k = i + j
            (1 <= k <= n) || continue
            num += w[k] * g[jj] * d[k]
            den += w[k] * g[jj]^2
        end
        snr[i] = den > 0 ? num / sqrt(den) : 0.0
    end

    runs = arcpk_saturated_runs(mask)
    satzone = falses(n)
    for (lo, hi, _) in runs
        satzone[max(1, lo - p.sat_pad):min(n, hi + p.sat_pad)] .= true
    end

    det = Int[]
    @inbounds for i in (p.edge):(n - p.edge + 1)
        snr[i] >= p.nsig || continue
        (snr[i] >= snr[i - 1] && snr[i] > snr[i + 1]) || continue
        push!(det, i)
    end
    isempty(det) &&
        return (empty[1], empty[2], empty[3], merge(empty[4], (; sat = runs, noise_excess = excess)))
    near_sat = [satzone[i] for i in det]

    n_offset = 15
    n_pad = 20
    det_fit = clamp.(det, n_offset + 1, n - n_offset - 1)
    offset_inds = -4:4
    best_model_offset_inds = (-n_offset):n_offset
    fit_inds = det_fit .+ offset_inds'
    best_model_fit_inds = det_fit .+ best_model_offset_inds'
    fluxz = copy(flux)
    ivarz = copy(ivar)
    fluxz[.!good] .= 0.0
    ivarz[.!good] .= 0.0
    fit_fluxes = fluxz[fit_inds]
    fg = zeros(Float64, length(det), 4)
    fg[:, 1] .= max.(0.01, vec(maximum(fit_fluxes, dims = 2)))
    fg[:, 2] .= det_fit
    fg[:, 3] .= 1.3
    fg[:, 4] .= 0
    fg[:, 1] .*= (2π)^0.5 .* fg[:, 3]
    flux_pad = zeros(n + 2n_pad)
    ivar_pad = zeros(n + 2n_pad)
    flux_pad[(1 + n_pad):(end - n_pad)] .= fluxz
    ivar_pad[(1 + n_pad):(end - n_pad)] .= ivarz
    fg[:, 2] .+= n_pad
    params, cov = fit_gauss_and_bias(flux_pad, ivar_pad, fg, fit_inds, best_model_fit_inds,
        offset_inds; n_iter = 10, return_cov = true, use_first_guess_heights = true,
        max_center_move = 3, min_widths = 0.1, max_widths = 3.0, n_pad = n_pad,
        n_offset = n_offset)
    params[:, 2] .-= n_pad

    nmask = [count(!, good[max(1, i - 2):min(n, i + 2)]) for i in det]
    fl = arcpk_candidate_flags(params, cov, det, nmask, near_sat; p = p)
    acc = findall(==(0), fl)
    acc = acc[sortperm(params[acc, 2])]
    info = (; snr = snr[det[acc]], det_pix = det[acc], cand_params = params, cand_cov = cov,
        cand_snr = snr[det], cand_det = det, cand_nmask = nmask, cand_nearsat = near_sat,
        cand_flags = fl, cand_npsat = [count(psat[max(1, i - 2):min(n, i + 2)]) for i in det],
        sat = runs, noise_excess = excess)
    return collect(1:length(acc)), params[acc, :], cov[acc, :, :], info
end

"""
    arclamp_exposure_peaks(fnames_by_chip; nsig = 8.0) -> NamedTuple of vectors

Run `get_arclamp_peaks_local` on every fiber of the three chips of one arc exposure
(`fnames_by_chip` = the three ar1Dcal paths, R, G, B; opened read-only). Returns the accepted
peaks as flat vectors: pix, pixvar (centroid variance), height, sigma, fiber, chip (1..3).
"""
function arclamp_exposure_peaks(fnames_by_chip::AbstractVector{<:AbstractString}; nsig = 8.0)
    pix = Float64[]
    pv = Float64[]
    hgt = Float64[]
    sg = Float64[]
    fib = Int[]
    chp = Int[]
    for (ci, fn) in enumerate(fnames_by_chip)
        isfile(fn) || continue
        flux, ivar, mask = h5open(fn, "r") do f
            read(f["flux_1d"]), read(f["ivar_1d"]), read(f["mask_1d"])
        end
        for fb in 1:size(flux, 2)
            _, prm, cov, _ = get_arclamp_peaks_local(view(flux, :, fb), view(ivar, :, fb),
                view(mask, :, fb); nsig = nsig)
            for k in 1:size(prm, 1)
                x = prm[k, 2]
                v = cov[k, 2, 2]
                (isfinite(x) && 1 <= x <= N_XPIX && v > 0) || continue
                push!(pix, x)
                push!(pv, v)
                push!(hgt, prm[k, 1])
                push!(sg, prm[k, 3])
                push!(fib, fb)
                push!(chp, ci)
            end
        end
    end
    return (; pix, pixvar = pv, height = hgt, sigma = sg, fiber = fib, chip = chp)
end

# ------------------------------------------------------------------------------------------
# frozen line list
# ------------------------------------------------------------------------------------------
const ARCLAMP_LAMPS = ("ThAr", "UNe")

struct ArcLineSlice
    wave::Vector{Float64}       # association wavelength, sorted
    store::Vector{Bool}         # true: lab-matched feature, kept in the line table
    lab_id::Vector{Int}         # line id (0 = guard feature without a lab match)
end

"""
    ArcLineList

The frozen arc line list (`data/arclamp_linelist_v1.csv`, built by
`scripts/cal/arclamp_linelist/build_arclamp_linelist.jl`). Every row is one recurring arc feature
of one lamp on one chip. Features without a laboratory match stay in as association guards.
`fit_ids` is the set of line ids fitted (the `selection` column), at their laboratory wavelengths.
"""
struct ArcLineList
    path::String
    version::String
    fingerprint::String
    selection::String
    df::DataFrame
    slices::Dict{Tuple{Int, Int}, ArcLineSlice}
    seed::Dict{Int, Vector{Float64}}     # per chip: quartic in ximport(pix) -> wave (seed only)
    labwave::Dict{Int, Float64}
    lamp_of_id::Dict{Int, Int}
    fit_ids::Set{Int}
end

const ARCLAMP_LINELIST_DEFAULT = joinpath(dirname(@__DIR__), "data", "arclamp_linelist_v1.csv")

"""
    load_arclamp_linelist(path = ARCLAMP_LINELIST_DEFAULT; selection = "use_default")

Read the frozen line list. `selection` names the Bool column that decides which lab lines are
fitted:
  * "use_default"  : shape-clean lab Th/U/Ar lines (93 catalog lines), without the DRP re-derived (blended) lines
                     and without neon (the stage-1 default, pending PI confirmation);
  * "prelim_classB": the September lab class-B Th/U/Ar lines present in the catalog (70);
  * "use_with_blended": use_default plus the 3 DRP re-derived (blended) catalog features;
  * "use_with_neon": use_default plus every lab-matched neon feature (at lab wavelength).
Any other Bool column of the file can be named.

`neon_offset_A`: optional Dict lamp ("ThAr", "UNe") => offset in Å added to the laboratory
wavelength of the NEON lines only (a measured per-lamp offset for the telescope being reduced,
from a calibration file; design option L4 (c)). It has an effect only if the selection fits
neon. The selection name recorded in products includes the offsets.
"""
function load_arclamp_linelist(path::AbstractString = ARCLAMP_LINELIST_DEFAULT;
        selection::AbstractString = "use_default", neon_offset_A = nothing)
    isfile(path) || throw(ArgumentError("arc line list not found: $(path)"))
    df = CSV.read(path, DataFrame; comment = "#")
    hasproperty(df, Symbol(selection)) ||
        throw(ArgumentError("$(path): no selection column '$(selection)'"))
    sel = Bool.(df[!, Symbol(selection)])
    store = Bool.(df.store)
    lab_id = Int.(df.lab_id)
    any(sel .& .!store) &&
        throw(ArgumentError("$(path): selection '$(selection)' contains rows with store = false"))
    slices = Dict{Tuple{Int, Int}, ArcLineSlice}()
    for (li, ln) in enumerate(ARCLAMP_LAMPS), c in 1:N_CHIPS
        m = findall((df.lamp .== ln) .& (df.chip_idx .== c))
        o = m[sortperm(df.assoc_wave[m])]
        slices[(li, c)] = ArcLineSlice(Float64.(df.assoc_wave[o]), store[o], lab_id[o])
    end
    seed = Dict{Int, Vector{Float64}}()
    for c in 1:N_CHIPS
        m = findall((df.chip_idx .== c) .& isfinite.(df.med_pix))
        A = positional_poly_mat(_arc_ximport(Float64.(df.med_pix[m])), porder = 4)
        seed[c] = A \ Float64.(df.assoc_wave[m])
    end
    labwave = Dict{Int, Float64}()
    lamp_of_id = Dict{Int, Int}()
    isneon = Bool.(df.neon)
    for i in 1:nrow(df)
        lab_id[i] > 0 || continue
        w = Float64(df.lab_wave[i])
        if isneon[i] && neon_offset_A !== nothing
            w += Float64(get(neon_offset_A, String(df.lamp[i]), 0.0))
        end
        labwave[lab_id[i]] = w
        lamp_of_id[lab_id[i]] = df.lamp[i] == "ThAr" ? 1 : 2
    end
    fit_ids = Set(lab_id[sel .& (lab_id .> 0)])
    version = splitext(basename(path))[1]
    selname = neon_offset_A === nothing ? String(selection) :
              String(selection) * "+neon_offset(" *
              join(["$(k)=$(neon_offset_A[k])" for k in sort(collect(keys(neon_offset_A)))], ",") * ")"
    return ArcLineList(String(path), version, arclamp_file_fingerprint(path), selname,
        df, slices, seed, labwave, lamp_of_id, fit_ids)
end

"Seed wavelength of pixels on chip c (quartic fitted to the list's median pixel positions)."
arclamp_rough_lambda(ll::ArcLineList, c, pix) =
    positional_poly_mat(_arc_ximport(pix), porder = 4) * ll.seed[c]

"Per-chip dispersion |dλ/dpix| (Å/px) at the chip centre, from the list's seed polynomial."
function arclamp_seed_dispersion(ll::ArcLineList)
    d = zeros(N_CHIPS)
    for c in 1:N_CHIPS
        l = arclamp_rough_lambda(ll, c, [1023.5, 1024.5])
        d[c] = abs(l[2] - l[1])
    end
    return Tuple(d)
end

# ------------------------------------------------------------------------------------------
# association
# ------------------------------------------------------------------------------------------
"""
    arclamp_associate(lam0, xt, wave; tol = 1.0, ...) -> Vector{Int}

Associate peaks with seed wavelengths `lam0` (and positions `xt`) to the sorted catalog
wavelengths `wave`: a global offset scan (±`scan_halfwidth` Å), then iterated quadratic
self-calibration in `xt` with sigma clipping, one-to-one assignment within `tol` Å, and a final
outlier cut. Returns the catalog index per peak (0 = not associated).
"""
function arclamp_associate(lam0, xt, wave::AbstractVector{Float64}; tol = 1.0,
        scan_halfwidth = 15.0, scan_step = 0.05, max_iter = 5, nsigma = 7.0,
        scatter_floor = 0.05)
    n = length(lam0)
    n == 0 && return zeros(Int, 0)
    isempty(wave) && return zeros(Int, n)
    bestoff, bestcount, bestcost = 0.0, -1, Inf
    for off in (-scan_halfwidth):scan_step:scan_halfwidth
        c = 0
        cost = 0.0
        @inbounds for i in 1:n
            l = lam0[i] + off
            j = searchsortedfirst(wave, l)
            d = Inf
            j > 1 && (d = min(d, abs(wave[j - 1] - l)))
            j <= length(wave) && (d = min(d, abs(wave[j] - l)))
            if d < tol
                c += 1
                cost += d^2
            end
        end
        if c > bestcount || (c == bestcount && cost < bestcost)
            bestcount = c
            bestcost = cost
            bestoff = off
        end
    end
    function assign(lams)
        a = zeros(Int, n)
        dist = fill(Inf, n)
        for i in 1:n
            l = lams[i]
            j = searchsortedfirst(wave, l)
            best, bd = 0, tol
            for jj in max(1, j - 1):min(length(wave), j)
                dd = abs(wave[jj] - l)
                dd < bd && (bd = dd; best = jj)
            end
            a[i] = best
            dist[i] = bd
        end
        owner = Dict{Int, Int}()
        for i in 1:n
            a[i] == 0 && continue
            o = get(owner, a[i], 0)
            if o == 0 || dist[i] < dist[o]
                o != 0 && (a[o] = 0)
                owner[a[i]] = i
            else
                a[i] = 0
            end
        end
        return a
    end
    lamcur = lam0 .+ bestoff
    a = assign(lamcur)
    prev = copy(a)
    prev2 = fill(-1, n)
    for _ in 1:max_iter
        m = findall(!=(0), a)
        length(m) < 5 && break
        y = [wave[a[i]] - lam0[i] for i in m]
        X = positional_poly_mat(xt[m], porder = 2)
        keep = trues(length(m))
        coef = X \ y
        for _ in 1:4
            coef = X[keep, :] \ y[keep]
            r = y .- X * coef
            s = max(_arc_robstd(r[keep]), scatter_floor)
            nk = abs.(r .- _arc_nanmedian(r[keep])) .< nsigma * s
            (count(nk) < 4 || all(nk .== keep)) && break
            keep = nk
        end
        lamcur = lam0 .+ positional_poly_mat(xt, porder = 2) * coef
        a = assign(lamcur)
        if all(a .== prev)
            break
        elseif all(a .== prev2)
            for i in 1:n
                (a[i] != prev[i]) && (a[i] = 0)
            end
            break
        end
        prev2 = prev
        prev = copy(a)
    end
    m = findall(!=(0), a)
    if length(m) >= 8
        r = [wave[a[i]] - lamcur[i] for i in m]
        s = max(_arc_robstd(r), scatter_floor)
        med = _arc_nanmedian(r)
        for (k, i) in enumerate(m)
            abs(r[k] - med) > nsigma * s && (a[i] = 0)
        end
    end
    return a
end

"""
    arclamp_sequence_lines(peaks_thar, peaks_une, ll::ArcLineList; tol = 1.0)

Associate the peaks of one arc sequence (one ThAr and one UNe exposure, from
`arclamp_exposure_peaks`) per fiber and chip. Returns a line table (NamedTuple of vectors:
fiber, chip, lamp (1 ThAr, 2 UNe), id, pix, pixvar, height) of the lab-matched features.
"""
function arclamp_sequence_lines(peaks_thar, peaks_une, ll::ArcLineList; tol = 1.0)
    rows = (fiber = Int16[], chip = Int8[], lamp = Int8[], id = Int16[], pix = Float32[],
        pixvar = Float32[], height = Float32[])
    for (li, P) in ((1, peaks_thar), (2, peaks_une))
        byf = [Int[] for _ in 1:N_FIBERS]
        for i in eachindex(P.pix)
            push!(byf[P.fiber[i]], i)
        end
        for f in 1:N_FIBERS, c in 1:N_CHIPS
            m = [i for i in byf[f] if P.chip[i] == c]
            isempty(m) && continue
            pix = Float64.(P.pix[m])
            cs = ll.slices[(li, c)]
            a = arclamp_associate(arclamp_rough_lambda(ll, c, pix), _arc_ximport(pix), cs.wave;
                tol = tol)
            for (k, i) in enumerate(m)
                a[k] == 0 && continue
                cs.store[a[k]] || continue
                id = cs.lab_id[a[k]]
                id > 0 || continue
                push!(rows.fiber, f)
                push!(rows.chip, c)
                push!(rows.lamp, li)
                push!(rows.id, id)
                push!(rows.pix, P.pix[i])
                push!(rows.pixvar, P.pixvar[i])
                push!(rows.height, P.height[i])
            end
        end
    end
    return rows
end

# ------------------------------------------------------------------------------------------
# joint night fit
# ------------------------------------------------------------------------------------------
const ARCFIT_NLOC = 9          # a_R b_R a_B b_B p0..p4
"Fiber coordinate of the per-sequence slit tilt."
arcfit_phi(f) = (f - 150.5) / 150.0

"""Index layout of the night-global parameters: per-sequence shift and tilt (sequence 1 is the
reference), the second-position per-chip zero point (R, B) and stretch (R, G, B), and the
linear-in-time per-chip zero point (R, B) and stretch (R, G, B). 0 = absent."""
struct ArcGlobalLayout
    S::Int
    has2::Bool
    timeok::Bool
    idx_d::Vector{Int}
    idx_t::Vector{Int}
    idx_alpha::NTuple{3, Int}
    idx_beta::NTuple{3, Int}
    idx_u::NTuple{3, Int}
    idx_v::NTuple{3, Int}
    n::Int
end
function ArcGlobalLayout(S, grp, tau; force_notime = false, force_nopos = false)
    k = 0
    idd = zeros(Int, S)
    idt = zeros(Int, S)
    for s in 2:S
        idd[s] = (k += 1)
        idt[s] = (k += 1)
    end
    has2 = !force_nopos && any(==(1), grp) && any(==(2), grp)
    ia = (0, 0, 0)
    ib = (0, 0, 0)
    if has2
        ia = (k + 1, 0, k + 2)
        k += 2
        ib = (k + 1, k + 2, k + 3)
        k += 3
    end
    # the time terms are identifiable only if one position group has >= 2 sequences > 1 h apart
    timeok = false
    if !force_notime
        for g in (1, 2)
            ts = tau[grp .== g]
            length(ts) >= 2 && (maximum(ts) - minimum(ts)) * 12 > 1.0 && (timeok = true)
        end
    end
    iu = (0, 0, 0)
    iv = (0, 0, 0)
    if timeok
        iu = (k + 1, 0, k + 2)
        k += 2
        iv = (k + 1, k + 2, k + 3)
        k += 3
    end
    return ArcGlobalLayout(S, has2, timeok, idd, idt, ia, ib, iu, iv, k)
end

struct ArcJointData
    seq::Vector{Int}
    fiber::Vector{Int}
    chip::Vector{Int}
    x::Vector{Float64}
    y::Vector{Float64}
    w::Vector{Float64}
    rows_of_fiber::Vector{Vector{Int}}
end
"""Fit rows: `use[i]` and finite `y[i]`. Weight 1/((σ_pix² + floor²)·disp_c²)."""
function ArcJointData(L, y::Vector{Float64}, use::AbstractVector{Bool};
        floor_px = 0.02, disp_chip = (0.2359, 0.2828, 0.3265))
    sel = findall(i -> use[i] && isfinite(y[i]), eachindex(y))
    seq = Int.(L.seq[sel])
    fib = Int.(L.fiber[sel])
    ch = Int.(L.chip[sel])
    x = _arc_ximport(Float64.(L.pix[sel]))
    w = [1.0 / ((max(Float64(L.pixvar[i]), 0.0) + floor_px^2) * disp_chip[L.chip[i]]^2)
         for i in sel]
    rof = [Int[] for _ in 1:N_FIBERS]
    for (k, f) in enumerate(fib)
        push!(rof[f], k)
    end
    return ArcJointData(seq, fib, ch, x, y[sel], w, rof), sel
end

@inline function _arcfit_shift!(gd, lay::ArcGlobalLayout, g, s, grp, tau, f, c, x)
    fill!(gd, 0.0)
    sh = 0.0
    if s > 1
        sh += g[lay.idx_d[s]] + g[lay.idx_t[s]] * arcfit_phi(f)
        gd[lay.idx_d[s]] = 1.0
        gd[lay.idx_t[s]] = arcfit_phi(f)
    end
    if lay.has2 && grp[s] == 2
        ia = lay.idx_alpha[c]
        ib = lay.idx_beta[c]
        ia > 0 && (sh += g[ia]; gd[ia] = 1.0)
        ib > 0 && (sh += g[ib] * x; gd[ib] = x)
    end
    if lay.timeok
        iu = lay.idx_u[c]
        iv = lay.idx_v[c]
        t = tau[s]
        iu > 0 && (sh += g[iu] * t; gd[iu] = t)
        iv > 0 && (sh += g[iv] * x * t; gd[iv] = x * t)
    end
    return sh
end
@inline function _arcfit_xt(th, c, x)
    c == 2 && return x
    c == 1 && return th[1] + th[2] * x
    return th[3] + th[4] * x
end
@inline _arcfit_poly4(p0, p1, p2, p3, p4, z) = p0 + z * (p1 + z * (p2 + z * (p3 + z * p4)))
@inline _arcfit_dpoly4(p1, p2, p3, p4, z) = p1 + z * (2p2 + z * (3p3 + z * 4p4))

function _arcfit_fiber_resid(D::ArcJointData, rows, th, g, lay, grp, tau, f)
    r = zeros(length(rows))
    gd = zeros(max(lay.n, 1))
    for (k, i) in enumerate(rows)
        c = D.chip[i]
        x = D.x[i]
        xs = _arcfit_xt(th, c, x) + _arcfit_shift!(gd, lay, g, D.seq[i], grp, tau, f, c, x)
        r[k] = D.y[i] - _arcfit_poly4(th[5], th[6], th[7], th[8], th[9], xs)
    end
    return r
end

"Initial chip offsets (a_R, b_R, a_B, b_B) for the Gauss-Newton seed (seed only)."
arcfit_offsets0(tele) = tele == "apo" ? [-1.070, 1.0, 1.0755, 1.0] : [-1.07295, 1.0, 1.076, 1.0]

"""Gauss-Newton joint fit on the rows with `keep`. With `chip_scales` (N_FIBERS x 2: b_R, b_B)
the chip scales are held fixed. Returns (TH :: N_FIBERS x 9, g, ok-fibers)."""
function arcfit_gauss_newton!(D::ArcJointData, keep::AbstractVector{Bool}, lay::ArcGlobalLayout,
        grp, tau, tele; theta0 = nothing, g0 = nothing, maxit = 25, minrows = 40, tol = 1e-9,
        chip_scales = nothing)
    n = lay.n
    TH = theta0 === nothing ? zeros(N_FIBERS, ARCFIT_NLOC) : copy(theta0)
    g = g0 === nothing ? zeros(n) : copy(g0)
    okf = falses(N_FIBERS)
    o0 = arcfit_offsets0(tele)
    gd = zeros(max(n, 1))
    fixed = chip_scales !== nothing
    if theta0 === nothing
        for f in 1:N_FIBERS
            rows = [i for i in D.rows_of_fiber[f] if keep[i]]
            r1 = [i for i in rows if D.seq[i] == 1]
            length(r1) < 15 && (r1 = rows)
            length(r1) < 15 && continue
            th = vcat(o0, zeros(5))
            if fixed
                th[2] = chip_scales[f, 1]
                th[4] = chip_scales[f, 2]
            end
            xs = [_arcfit_xt(th, D.chip[i], D.x[i]) for i in r1]
            A = positional_poly_mat(xs, porder = 4)
            sw = sqrt.(D.w[r1])
            th[5:9] .= (A .* sw) \ (D.y[r1] .* sw)
            TH[f, :] .= th
            okf[f] = true
        end
        if n > 0
            for s in 2:(lay.S)
                v = Float64[]
                for f in 1:N_FIBERS
                    okf[f] || continue
                    th = TH[f, :]
                    for i in D.rows_of_fiber[f]
                        (keep[i] && D.seq[i] == s) || continue
                        xs = _arcfit_xt(th, D.chip[i], D.x[i])
                        push!(v, (D.y[i] - _arcfit_poly4(th[5:9]..., xs)) /
                                 _arcfit_dpoly4(th[6:9]..., xs))
                    end
                end
                isempty(v) || (g[lay.idx_d[s]] = median(v))
            end
        end
    else
        for f in 1:N_FIBERS
            okf[f] = all(isfinite, TH[f, :])
        end
    end
    function total_loss(TH, g)
        Ls = 0.0
        for f in 1:N_FIBERS
            okf[f] || continue
            rows = [i for i in D.rows_of_fiber[f] if keep[i]]
            r = _arcfit_fiber_resid(D, rows, TH[f, :], g, lay, grp, tau, f)
            Ls += sum(D.w[rows] .* r .^ 2)
        end
        return Ls
    end
    Lold = total_loss(TH, g)
    for it in 1:maxit
        Ssch = zeros(n, n)
        rhs = zeros(n)
        Haa_inv = Vector{Matrix{Float64}}(undef, N_FIBERS)
        Habs = Vector{Matrix{Float64}}(undef, N_FIBERS)
        gas = Vector{Vector{Float64}}(undef, N_FIBERS)
        for f in 1:N_FIBERS
            okf[f] || continue
            rows = [i for i in D.rows_of_fiber[f] if keep[i]]
            if length(rows) < minrows
                okf[f] = false
                continue
            end
            th = TH[f, :]
            m = length(rows)
            A = zeros(m, ARCFIT_NLOC)
            B = zeros(m, n)
            r = zeros(m)
            w = D.w[rows]
            for (k, i) in enumerate(rows)
                c = D.chip[i]
                x = D.x[i]
                xs = _arcfit_xt(th, c, x) + _arcfit_shift!(gd, lay, g, D.seq[i], grp, tau, f, c, x)
                lp = _arcfit_dpoly4(th[6], th[7], th[8], th[9], xs)
                r[k] = D.y[i] - _arcfit_poly4(th[5], th[6], th[7], th[8], th[9], xs)
                if c == 1
                    A[k, 1] = lp
                    A[k, 2] = lp * x
                elseif c == 3
                    A[k, 3] = lp
                    A[k, 4] = lp * x
                end
                z = 1.0
                for j in 0:4
                    A[k, 5 + j] = z
                    z *= xs
                end
                n > 0 && for j in 1:n
                    B[k, j] = lp * gd[j]
                end
            end
            # a chip with no rows: freeze its transform; fixed chip scales: freeze b_R, b_B
            for (cc, cols) in ((1, (1, 2)), (3, (3, 4)))
                if !any(D.chip[i] == cc for i in rows)
                    A[:, cols[1]] .= 0
                    A[:, cols[2]] .= 0
                end
            end
            if fixed
                A[:, 2] .= 0
                A[:, 4] .= 0
            end
            Aw = A .* w
            Haa = A' * Aw
            for j in 1:ARCFIT_NLOC
                Haa[j, j] == 0 && (Haa[j, j] = 1.0)
            end
            Hi = try
                inv(Symmetric(Haa))
            catch
                okf[f] = false
                continue
            end
            ga = Aw' * r
            Haa_inv[f] = Hi
            gas[f] = ga
            if n > 0
                Hab = Aw' * B
                Hbb = B' * (B .* w)
                gb = (B .* w)' * r
                Habs[f] = Hab
                Ssch .+= Hbb .- Hab' * Hi * Hab
                rhs .+= gb .- Hab' * (Hi * ga)
            end
        end
        dg = zeros(n)
        if n > 0
            for j in 1:n
                Ssch[j, j] += 1e-12 * max(Ssch[j, j], 1e-30)
            end
            dg = try
                Symmetric(Ssch) \ rhs
            catch
                pinv(Ssch) * rhs
            end
        end
        dTH = zeros(N_FIBERS, ARCFIT_NLOC)
        for f in 1:N_FIBERS
            okf[f] || continue
            dTH[f, :] .= Haa_inv[f] * (n > 0 ? gas[f] .- Habs[f] * dg : gas[f])
        end
        step = 1.0
        accepted = false
        Lnew = Lold
        for _ in 1:6
            THn = TH .+ step .* dTH
            gn = g .+ step .* dg
            Lnew = total_loss(THn, gn)
            if Lnew <= Lold * (1 + 1e-12)
                TH = THn
                g = gn
                accepted = true
                break
            end
            step /= 2
        end
        accepted || break
        rel = (Lold - Lnew) / max(Lold, 1e-30)
        Lold = Lnew
        rel < tol && it > 2 && break
    end
    TH[.!okf, :] .= NaN
    return TH, g, okf
end

"""Joint fit with rejection: |r| < `resid_cut` Å, then rounds of `clipk`-sigma clipping of
r·√w per fiber. Returns (TH, g, okf, keep, resid)."""
function arcfit_joint(D::ArcJointData, lay, grp, tau, tele; theta0 = nothing, g0 = nothing,
        chip_scales = nothing, resid_cut = 0.35, clipk = 5.0)
    keep = fill(true, length(D.y))
    TH, g, okf = arcfit_gauss_newton!(D, keep, lay, grp, tau, tele; theta0 = theta0, g0 = g0,
        chip_scales = chip_scales)
    for round in 1:3
        changed = false
        for f in 1:N_FIBERS
            rows = D.rows_of_fiber[f]
            isempty(rows) && continue
            okf[f] || continue
            r = _arcfit_fiber_resid(D, rows, TH[f, :], g, lay, grp, tau, f)
            if round == 1
                nk = abs.(r) .< resid_cut
            else
                z = r .* sqrt.(D.w[rows])
                kk = [keep[i] for i in rows]
                s = _arc_robstd(z[kk])
                (isnan(s) || s <= 0) && continue
                nk = kk .& (abs.(z .- _arc_nanmedian(z[kk])) .< clipk * s)
            end
            for (k, i) in enumerate(rows)
                if keep[i] != nk[k]
                    keep[i] = nk[k]
                    changed = true
                end
            end
        end
        changed || (round > 1 && break)
        TH0 = copy(TH)
        for f in 1:N_FIBERS
            okf[f] || (TH0[f, :] .= NaN)
        end
        TH, g, okf = arcfit_gauss_newton!(D, keep, lay, grp, tau, tele; theta0 = TH0, g0 = g,
            chip_scales = chip_scales)
    end
    resid = fill(NaN, length(D.y))
    for f in 1:N_FIBERS
        okf[f] || continue
        rows = D.rows_of_fiber[f]
        resid[rows] .= _arcfit_fiber_resid(D, rows, TH[f, :], g, lay, grp, tau, f)
    end
    return (; TH, g, okf, keep, resid)
end

"Per-sequence (shift, tilt) in x units (0, 0 for the reference sequence)."
arcfit_seq_state(g, lay::ArcGlobalLayout, s) =
    s == 1 ? (0.0, 0.0) : (g[lay.idx_d[s]], g[lay.idx_t[s]])

"""
    arcfit_eval(th, g, lay, f, c, px; grpv = 1, d = 0.0, t = 0.0, tau = 0.0)

Wavelength of fiber `f`, chip `c` at pixels `px` for a given state: position group `grpv`,
sequence shift `d` and tilt `t` (x units) and time `tau` (units of 12 h from the night's mean
sequence time).
"""
function arcfit_eval(th, g, lay::ArcGlobalLayout, f, c, px; grpv = 1, d = 0.0, t = 0.0,
        tau = 0.0)
    out = similar(px, Float64)
    for (k, p) in enumerate(px)
        x = (p - 1024) / 2048
        xs = _arcfit_xt(th, c, x) + d + t * arcfit_phi(f)
        if lay.has2 && grpv == 2
            ia = lay.idx_alpha[c]
            ib = lay.idx_beta[c]
            ia > 0 && (xs += g[ia])
            ib > 0 && (xs += g[ib] * x)
        end
        if lay.timeok
            iu = lay.idx_u[c]
            iv = lay.idx_v[c]
            iu > 0 && (xs += g[iu] * tau)
            iv > 0 && (xs += g[iv] * x * tau)
        end
        out[k] = _arcfit_poly4(th[5], th[6], th[7], th[8], th[9], xs)
    end
    return out
end

# ------------------------------------------------------------------------------------------
# measured positions
# ------------------------------------------------------------------------------------------
"""
    arclamp_position_clusters(shift_px; outlier_px = 3.0, group_px = 0.25)

Measured per-sequence shifts (px, relative to the reference sequence) -> (posout, cluster id,
number of clusters). A sequence more than `outlier_px` from the night median is a position
outlier (cluster 0). The rest are single-linkage clustered with gap `group_px`.
"""
function arclamp_position_clusters(shift_px::AbstractVector; outlier_px = 3.0, group_px = 0.25)
    S = length(shift_px)
    posout = abs.(shift_px .- median(shift_px)) .> outlier_px
    keepi = findall(.!posout)
    cl = zeros(Int, S)
    isempty(keepi) && return posout, cl, 0
    o = keepi[sortperm(shift_px[keepi])]
    ncl = 1
    cl[o[1]] = 1
    for k in 2:length(o)
        shift_px[o[k]] - shift_px[o[k - 1]] > group_px && (ncl += 1)
        cl[o[k]] = ncl
    end
    return posout, cl, ncl
end

# ------------------------------------------------------------------------------------------
# arc-only sequence rejection
# ------------------------------------------------------------------------------------------
# Arcs come first in the wavelength chain, so a sequence may be judged only against the night's
# other arc sequences. Every sequence keeps its own measured shift and tilt in the joint model,
# so a pure displacement is absorbed and is NOT a reason to reject (study 2026_11_18/arc_gate:
# a 0.1 px displacement of one sequence moves the night solution at the other sequences' states
# by < 0.03 m/s). What is rejected is a sequence whose lines do not fit the shared night model
# as the others do:
#   * quality: larger line-residual scatter, or fewer lines surviving the clipping;
#   * shape: its binned residual pattern differs from the night's median pattern beyond a
#     per-chip zero point and stretch (pixel direction) or tilt (fiber direction);
#   * chip offset: its per-chip residual zero points differ from the night's by a
#     chip-differential amount. The model has no per-sequence per-chip zero point, so this part
#     leaks into the other sequences' solutions (study: 600-1500 m/s per px, chip-dependent).
const ARCSTAT_NXB = 8                  # pixel bins per chip
const ARCSTAT_NFB = 10                 # fiber bins per chip

"""
    arclamp_sequence_stats(D, F, nseq, disp_chip) -> Vector of NamedTuples

Per-sequence statistics about a joint fit (rows of `D`, fit `F`, sequences 1..nseq of the fit):
`n_rows`, `n_keep`, `keep_frac`, `rsig_px` (robust scatter of the kept residuals, px), and, from
the binned median residuals (8 pixel bins x 10 fiber bins per chip, bins with >= 10 kept rows)
minus the night's median binned pattern over the sequences:
  * `chip_zp_px`, `chip_stretch_px` (per chip, px): straight line in x through the difference;
  * `chipoff_px`: the largest |zero point of R or B minus that of G| (chip-differential offset);
  * `shape_px`: RMS of the difference after the per-chip line (pixel bins) and after a per-chip
    line in fiber coordinate (fiber bins).
For a single sequence the differences are zero by construction.
"""
function arclamp_sequence_stats(D::ArcJointData, F, nseq::Integer, disp_chip)
    xe = range(-0.5, 0.5, length = ARCSTAT_NXB + 1)
    fe = range(-1.0, 1.0, length = ARCSTAT_NFB + 1)
    xc = collect((xe[1:(end - 1)] .+ xe[2:end]) ./ 2)
    fc = collect((fe[1:(end - 1)] .+ fe[2:end]) ./ 2)
    binof(v, e, nb) = clamp(searchsortedlast(e, v), 1, nb)
    Bx = fill(NaN, nseq, N_CHIPS, ARCSTAT_NXB)
    Bf = fill(NaN, nseq, N_CHIPS, ARCSTAT_NFB)
    base = map(1:nseq) do j
        rows = [i for i in eachindex(D.seq)
                if D.seq[i] == j && F.okf[D.fiber[i]] && isfinite(F.resid[i])]
        kept = [i for i in rows if F.keep[i]]
        rpx = [F.resid[i] / disp_chip[D.chip[i]] for i in kept]
        for c in 1:N_CHIPS
            m = [D.chip[i] == c for i in kept]
            x = D.x[kept[m]]
            ph = arcfit_phi.(D.fiber[kept[m]])
            r = rpx[m]
            bx = [binof(v, xe, ARCSTAT_NXB) for v in x]
            bf = [binof(v, fe, ARCSTAT_NFB) for v in ph]
            for b in 1:ARCSTAT_NXB
                sel = bx .== b
                count(sel) >= 10 && (Bx[j, c, b] = median(r[sel]))
            end
            for b in 1:ARCSTAT_NFB
                sel = bf .== b
                count(sel) >= 10 && (Bf[j, c, b] = median(r[sel]))
            end
        end
        (n_rows = length(rows), n_keep = length(kept),
         keep_frac = isempty(rows) ? NaN : length(kept) / length(rows), rsig_px = _arc_robstd(rpx))
    end
    refx = [_arc_nanmedian(Bx[:, c, b]) for c in 1:N_CHIPS, b in 1:ARCSTAT_NXB]
    reff = [_arc_nanmedian(Bf[:, c, b]) for c in 1:N_CHIPS, b in 1:ARCSTAT_NFB]
    out = map(1:nseq) do j
        zp = fill(NaN, N_CHIPS)
        st = fill(NaN, N_CHIPS)
        res = Float64[]
        for c in 1:N_CHIPS
            for (B, ref, cc, setzp) in ((Bx, refx, xc, true), (Bf, reff, fc, false))
                d = B[j, c, :] .- ref[c, :]
                ok = isfinite.(d)
                count(ok) >= 3 || continue
                A = hcat(ones(count(ok)), cc[ok])
                cf = A \ d[ok]
                append!(res, d[ok] .- A * cf)
                if setzp
                    zp[c] = cf[1]
                    st[c] = cf[2]
                end
            end
        end
        offs = [abs(zp[c] - zp[2]) for c in (1, 3) if isfinite(zp[c]) && isfinite(zp[2])]
        merge(base[j], (shape_px = isempty(res) ? NaN : sqrt(mean(res .^ 2)),
            chipoff_px = isempty(offs) ? NaN : maximum(offs), chip_zp_px = zp, chip_stretch_px = st))
    end
    return out
end

"""Default thresholds of the arc-only sequence rejection (all are options of
`fit_arclamp_night`, recorded in the product metadata; values chosen on the 2026_11_18/arc_gate
sample, pending PI confirmation).

Relative test (nights with at least `min_seq_compare` sequences): for a statistic `q` of sequence
s, with the night median `m` and robust spread `σ` = max(1.4826 MAD, `q_scale_floor`) over the
night's sequences, s is an outlier when its excess (q − m; m − q for the kept fraction `keep`)
exceeds both `reject_k` σ and the absolute excess floor `q_floor`. Statistics: `rsig`, `keep`,
`shape` and `chipoff` of `arclamp_sequence_stats`.

Absolute test (every night): `rsig` > `rsig_max_px` or `keep` < `keep_min` (or no kept rows)
marks a sequence as failing the absolute ceilings."""
const ARCREJ_DEFAULTS = (reject_k = 5.0, min_seq_compare = 3,
    rsig_scale_floor_px = 0.001, rsig_floor_px = 0.01,
    keep_scale_floor = 0.0005, keep_floor = 0.02,
    shape_scale_floor_px = 0.0005, shape_floor_px = 0.005,
    chipoff_scale_floor_px = 0.001, chipoff_floor_px = 0.01,
    rsig_max_px = 0.12, keep_min = 0.90)

"""
    arclamp_sequence_reject(stats; opts...) -> (quality, shape, chipoff, absfail, score, compared)

Judge a night's sequences from their `arclamp_sequence_stats` (rule in `ARCREJ_DEFAULTS`).
`quality[s]`: relative outlier in residual scatter or kept fraction, or no kept rows; `shape[s]`
/ `chipoff[s]`: relative outlier in residual shape / chip-differential offset (all three only
when `compared`, i.e. at least `min_seq_compare` sequences). `absfail[s]`: fails an absolute
ceiling. `score[s]`: the largest ratio of a statistic's excess to its threshold (> 1 = fails
something; the largest score marks the worst sequence). `fit_arclamp_night` removes failing
sequences one at a time, worst first, refitting after each.
"""
function arclamp_sequence_reject(stats::AbstractVector; reject_k = ARCREJ_DEFAULTS.reject_k,
        min_seq_compare = ARCREJ_DEFAULTS.min_seq_compare,
        rsig_scale_floor_px = ARCREJ_DEFAULTS.rsig_scale_floor_px,
        rsig_floor_px = ARCREJ_DEFAULTS.rsig_floor_px,
        keep_scale_floor = ARCREJ_DEFAULTS.keep_scale_floor, keep_floor = ARCREJ_DEFAULTS.keep_floor,
        shape_scale_floor_px = ARCREJ_DEFAULTS.shape_scale_floor_px,
        shape_floor_px = ARCREJ_DEFAULTS.shape_floor_px,
        chipoff_scale_floor_px = ARCREJ_DEFAULTS.chipoff_scale_floor_px,
        chipoff_floor_px = ARCREJ_DEFAULTS.chipoff_floor_px,
        rsig_max_px = ARCREJ_DEFAULTS.rsig_max_px, keep_min = ARCREJ_DEFAULTS.keep_min)
    S = length(stats)
    rs = [Float64(x.rsig_px) for x in stats]
    kf = [Float64(x.keep_frac) for x in stats]
    sh = [Float64(x.shape_px) for x in stats]
    co = [Float64(x.chipoff_px) for x in stats]
    empty = [x.n_keep == 0 || !isfinite(x.rsig_px) || !isfinite(x.keep_frac) for x in stats]
    compared = S >= min_seq_compare
    function rel(q, sfloor, afloor, sgn)
        compared || return zeros(S)
        m = _arc_nanmedian(q)
        y = _arc_finite(q)
        sc = max(length(y) >= 2 ? 1.4826 * median(abs.(y .- m)) : 0.0, sfloor)
        thr = max(reject_k * sc, afloor)
        return [isfinite(q[s]) ? sgn * (q[s] - m) / thr : 0.0 for s in 1:S]
    end
    r_rs = rel(rs, rsig_scale_floor_px, rsig_floor_px, 1)
    r_kf = rel(kf, keep_scale_floor, keep_floor, -1)
    r_sh = rel(sh, shape_scale_floor_px, shape_floor_px, 1)
    r_co = rel(co, chipoff_scale_floor_px, chipoff_floor_px, 1)
    quality = (r_rs .> 1) .| (r_kf .> 1) .| (compared .& empty)
    shape = (r_sh .> 1) .& .!quality
    chipoff = (r_co .> 1) .& .!quality .& .!shape
    a_rs = [isfinite(x) ? x / rsig_max_px : Inf for x in rs]
    a_kf = [isfinite(x) ? (1 - x) / (1 - keep_min) : Inf for x in kf]
    absfail = (a_rs .> 1) .| (a_kf .> 1) .| empty
    score = [maximum((r_rs[s], r_kf[s], r_sh[s], r_co[s], a_rs[s], a_kf[s])) for s in 1:S]
    return quality, shape, chipoff, absfail, score, compared
end

# ------------------------------------------------------------------------------------------
# the night fit
# ------------------------------------------------------------------------------------------
# per-sequence flag bits (UInt16; the meanings are also written to the product metadata).
# Bit positions are never reused: 2 and 4 (the retired gate B against the object exposures'
# solutions) are reserved.
const ARCSEQ_POS_OUTLIER = 0x0001     # > outlier_px from the night median: excluded
#     0x0002 reserved (was ARCSEQ_GATE_FAIL, retired)
#     0x0004 reserved (was ARCSEQ_UNGATED, retired)
const ARCSEQ_LAMP_MISMATCH = 0x0008   # ThAr and UNe rows disagree in position by > 0.1 px
const ARCSEQ_FEW_LINES = 0x0010       # < 1000 kept rows in the sequence
const ARCSEQ_QUALITY = 0x0020         # rejected: residual scatter / kept fraction / ceiling
const ARCSEQ_SHAPE = 0x0040           # rejected: residual-shape outlier
const ARCSEQ_NOT_COMPARED = 0x0080    # night has < min_seq_compare sequences: not judged
const ARCSEQ_ABS_FAIL = 0x0100     # fails an absolute ceiling (rsig_max_px, keep_min)
const ARCSEQ_CHIPOFF = 0x0200      # rejected: chip-differential offset outlier
const ARCSEQ_REJECTED = ARCSEQ_QUALITY | ARCSEQ_SHAPE | ARCSEQ_CHIPOFF   # mask: rejected
# night flag bits (UInt16). 8 and 16 (retired gate B) are reserved.
const ARCNIGHT_MULTIPOS = 0x0001    # > 2 measured position clusters: no shared position terms
const ARCNIGHT_NO_TIME_TERMS = 0x0002  # time terms not identifiable
const ARCNIGHT_FEW_FIBERS = 0x0004  # < 250 good fibers
#     0x0008 reserved (was ARCNIGHT_ALL_GATED, retired)
#     0x0010 reserved (was ARCNIGHT_NO_GATE, retired)
const ARCNIGHT_CHIPSCALE_FREE = 0x0020  # chip scales fitted (no chip-transform calibration used)
const ARCNIGHT_POS_OUTLIER = 0x0040 # at least one sequence was a position outlier
const ARCNIGHT_SEQ_REJECTED = 0x0080  # at least one sequence rejected by the arc-only rule
const ARCNIGHT_FEW_SEQ = 0x0100     # < min_seq_compare sequences: rejection not possible
const ARCNIGHT_ABS_FAIL = 0x0200    # sequences in the final fit fail the absolute ceilings
"Bit meanings of `seq_flags` (written to the product metadata)."
const ARCLAMP_SEQFLAG_BITS = (
    "1: position outlier (> position_outlier_px from the night median), excluded",
    "2: reserved (retired gate B: gate failure)",
    "4: reserved (retired gate B: ungated)",
    "8: ThAr and UNe positions disagree by > lamp_mismatch_px",
    "16: fewer than 1000 kept rows",
    "32: rejected (line-residual scatter, kept fraction, no rows, or absolute ceiling)",
    "64: rejected, residual-shape outlier (pattern differs from the night's median pattern)",
    "128: not judged (night has fewer than min_seq_compare sequences)",
    "256: fails an absolute ceiling (reject_rsig_max_px, reject_keep_min)",
    "512: rejected, chip-differential offset outlier among the night's sequences")
"Bit meanings of `night_flags` (written to the product metadata)."
const ARCLAMP_NIGHTFLAG_BITS = (
    "1: more than 2 measured position clusters, no shared position terms",
    "2: time terms not identifiable",
    "4: fewer than 250 good fibers",
    "8: reserved (retired gate B: all gated sequences failed)",
    "16: reserved (retired gate B: no object reference)",
    "32: chip scales fitted (no chip-transform calibration used)",
    "64: at least one position outlier",
    "128: at least one sequence rejected by the arc-only rule",
    "256: fewer than min_seq_compare sequences, rejection not possible",
    "512: sequences in the final fit fail the absolute ceilings (night solution suspect)")

"""
    fit_arclamp_night(lines, tseq, tele, ll; chip_scales = nothing, reject = (;), ...)

Fit one night's arc sequences jointly, from arc data only.

  * `lines`: per sequence s, a line table from `arclamp_sequence_lines`.
  * `tseq`: per sequence, the mean mid-exposure MJD of its ThAr and UNe exposures.
  * `chip_scales`: N_FIBERS x 2 (b_R, b_B) to hold the chip scales fixed, or `nothing`.
  * `reject`: overrides of the rejection thresholds (keys of `ARCREJ_DEFAULTS`).

Stage 1 (shift + tilt per sequence only) measures each sequence's position; the positions are
clustered and position outliers excluded. The joint fit of the remaining sequences gives each
sequence's residual statistics (`arclamp_sequence_stats`), judged by `arclamp_sequence_reject`:
 1. relative outliers among the night's sequences are removed one at a time, worst first, with a
    refit after each removal (kept in the tables, flagged);
 2. if sequences fail the absolute ceilings, the worst is removed and the night refitted, as long
    as at least 2 sequences remain; if that does not leave a set that passes, these absolute
    removals are undone and the night is flagged (`ARCNIGHT_ABS_FAIL`) instead.
Nights with fewer than `min_seq_compare` sequences are flagged, never trimmed. Nothing from later pipeline
stages (sky, FPI) is used. Returns a NamedTuple with everything the night product stores.
"""
function fit_arclamp_night(lines::AbstractVector, tseq::AbstractVector{<:Real}, tele::AbstractString,
        ll::ArcLineList; chip_scales = nothing, group_px = 0.25, outlier_px = 3.0,
        floor_px = 0.02, disp_chip = arclamp_seed_dispersion(ll), lamp_mismatch_px = 0.1,
        reject::NamedTuple = NamedTuple())
    unknown = setdiff(keys(reject), keys(ARCREJ_DEFAULTS))
    isempty(unknown) || throw(ArgumentError("unknown rejection options: $(collect(unknown))"))
    ropt = merge(ARCREJ_DEFAULTS, reject)
    S = length(lines)
    S >= 1 || throw(ArgumentError("no arc sequences"))
    L = (seq = vcat([fill(Int8(k), length(lines[k].fiber)) for k in 1:S]...),
        fiber = vcat([x.fiber for x in lines]...), chip = vcat([x.chip for x in lines]...),
        lamp = vcat([x.lamp for x in lines]...), id = vcat([x.id for x in lines]...),
        pix = vcat([x.pix for x in lines]...), pixvar = vcat([x.pixvar for x in lines]...),
        height = vcat([x.height for x in lines]...))
    tseq = Float64.(tseq)
    any(isnan, tseq) && (tseq = [isnan(x) ? _arc_nanmean(tseq) : x for x in tseq])
    all(isnan, tseq) && (tseq = collect(1.0:S) ./ 24)
    ylab = [get(ll.labwave, Int(i), NaN) for i in L.id]
    mainuse = [Int(i) in ll.fit_ids for i in L.id]
    tmean = mean(tseq)
    tau = (tseq .- tmean) .* 2.0

    # stage 1: one shift + tilt per sequence; measured positions
    D1, _ = ArcJointData(L, ylab, mainuse; floor_px = floor_px, disp_chip = disp_chip)
    lay1 = ArcGlobalLayout(S, ones(Int, S), tau; force_notime = true, force_nopos = true)
    F1 = arcfit_joint(D1, lay1, ones(Int, S), tau, tele; chip_scales = chip_scales)
    dpx1 = [s == 1 ? 0.0 : F1.g[lay1.idx_d[s]] * N_XPIX for s in 1:S]
    posout, igrp, ncl = arclamp_position_clusters(dpx1; outlier_px = outlier_px,
        group_px = group_px)
    grp = ncl == 2 ? [igrp[s] == 2 ? 2 : 1 for s in 1:S] : ones(Int, S)

    function fit_subset(useseq)
        mp = zeros(Int, S)
        k = 0
        for s in 1:S
            useseq[s] && (mp[s] = (k += 1))
        end
        rowuse = [mainuse[i] && useseq[L.seq[i]] for i in eachindex(L.id)]
        Lm = merge(L, (seq = Int8[mp[s] > 0 ? mp[s] : 0 for s in L.seq],))
        Dm, selm = ArcJointData(Lm, ylab, rowuse; floor_px = floor_px, disp_chip = disp_chip)
        g2 = grp[useseq]
        ta = tau[useseq]
        lay = ArcGlobalLayout(k, g2, ta)
        Fm = arcfit_joint(Dm, lay, g2, ta, tele; chip_scales = chip_scales)
        return (; F = Fm, D = Dm, lay, grp = g2, tau = ta, sidx = findall(useseq), sel = selm)
    end

    full_all = fit_subset(.!posout)

    # arc-only rejection: each sequence against the night's other sequences (arc data only)
    seq_rsig_px = fill(NaN, S)
    seq_keep_frac = fill(NaN, S)
    seq_shape_px = fill(NaN, S)
    seq_chipoff_px = fill(NaN, S)
    seq_score = fill(NaN, S)
    seq_chip_zp_px = fill(NaN, S, N_CHIPS)
    seq_chip_stretch_px = fill(NaN, S, N_CHIPS)
    seq_absfail = falses(S)
    rej_quality = falses(S)
    rej_shape = falses(S)
    rej_chipoff = falses(S)
    rej_round = zeros(Int8, S)          # removal round (0 = not removed)
    function judge(FF)
        st = arclamp_sequence_stats(FF.D, FF.F, length(FF.sidx), disp_chip)
        q, sh, co, ab, sc, cmp = arclamp_sequence_reject(st; ropt...)
        for (j, s) in enumerate(FF.sidx)
            seq_rsig_px[s] = st[j].rsig_px
            seq_keep_frac[s] = st[j].keep_frac
            seq_shape_px[s] = st[j].shape_px
            seq_chipoff_px[s] = st[j].chipoff_px
            seq_score[s] = sc[j]
            seq_chip_zp_px[s, :] .= st[j].chip_zp_px
            seq_chip_stretch_px[s, :] .= st[j].chip_stretch_px
            seq_absfail[s] = ab[j]
        end
        return st, q, sh, co, ab, sc, cmp
    end
    cur = .!posout
    compared = count(cur) >= ropt.min_seq_compare
    full = full_all
    snap = nothing                      # state before the first absolute-only removal
    rnd = 0
    abs_fail_night = false
    while true
        st, q, sh, co, ab, sc, cmp = judge(full)
        rel = q .| sh .| co
        if cmp && any(rel)
            # worst relative outlier first
            jj = findall(rel)
            j = jj[argmax(sc[jj])]
        elseif any(ab) && count(cur) >= ropt.min_seq_compare
            snap === nothing && (snap = (cur = copy(cur), full = full, rq = copy(rej_quality),
                rs = copy(rej_shape), rc = copy(rej_chipoff), rr = copy(rej_round)))
            jj = findall(ab)
            j = jj[argmax(sc[jj])]
        else
            if any(ab) && snap !== nothing
                # absolute removals did not leave a passing set: undo them, flag the night
                cur, full = snap.cur, snap.full
                rej_quality .= snap.rq
                rej_shape .= snap.rs
                rej_chipoff .= snap.rc
                rej_round .= snap.rr
                judge(full)
                abs_fail_night = true
            elseif any(ab)
                abs_fail_night = true
            end
            break
        end
        s = full.sidx[j]
        rnd += 1
        rej_round[s] = rnd
        if q[j] || ab[j]
            rej_quality[s] = true
        elseif sh[j]
            rej_shape[s] = true
        else
            rej_chipoff[s] = true
        end
        cur[s] = false
        full = fit_subset(cur)
    end
    rejected = rej_quality .| rej_shape .| rej_chipoff
    F = full.F
    lay = full.lay

    # per-row residuals (all stored rows, also lines outside the fit set)
    resid_all = fill(NaN32, length(L.id))
    for i in eachindex(L.id)
        s = Int(L.seq[i])
        j = findfirst(==(s), full.sidx)
        j === nothing && continue
        f = Int(L.fiber[i])
        F.okf[f] || continue
        d, t = arcfit_seq_state(F.g, lay, j)
        m = arcfit_eval(view(F.TH, f, :), F.g, lay, f, Int(L.chip[i]), [Float64(L.pix[i])];
            grpv = full.grp[j], d = d, t = t, tau = full.tau[j])[1]
        resid_all[i] = Float32(ylab[i] - m)
    end
    keep_main = falses(length(L.id))
    keep_main[full.sel[F.keep]] .= true

    # per-sequence QA
    seqflags = zeros(UInt16, S)
    lamp_off_px = fill(NaN, S)
    nrows_seq = zeros(Int, S)
    for s in 1:S
        posout[s] && (seqflags[s] |= ARCSEQ_POS_OUTLIER)
        rej_quality[s] && (seqflags[s] |= ARCSEQ_QUALITY)
        rej_shape[s] && (seqflags[s] |= ARCSEQ_SHAPE)
        rej_chipoff[s] && (seqflags[s] |= ARCSEQ_CHIPOFF)
        (compared || posout[s]) || (seqflags[s] |= ARCSEQ_NOT_COMPARED)
        seq_absfail[s] && (seqflags[s] |= ARCSEQ_ABS_FAIL)
        r = [i for i in eachindex(L.id) if L.seq[i] == s && keep_main[i]]
        nrows_seq[s] = length(r)
        nrows_seq[s] < 1000 && (seqflags[s] |= ARCSEQ_FEW_LINES)
        rt = [resid_all[i] / disp_chip[L.chip[i]] for i in r if L.lamp[i] == 1]
        ru = [resid_all[i] / disp_chip[L.chip[i]] for i in r if L.lamp[i] == 2]
        lamp_off_px[s] = _arc_nanmedian(rt) - _arc_nanmedian(ru)
        abs(lamp_off_px[s]) > lamp_mismatch_px && (seqflags[s] |= ARCSEQ_LAMP_MISMATCH)
    end
    # final-fit per-sequence state, in px (NaN for sequences not in the final fit)
    seq_shift_px = fill(NaN, S)
    seq_tilt_px = fill(NaN, S)
    for (j, s) in enumerate(full.sidx)
        d, t = arcfit_seq_state(F.g, lay, j)
        seq_shift_px[s] = d * N_XPIX
        seq_tilt_px[s] = t * N_XPIX
    end
    gpx(ix) = ntuple(c -> ix[c] > 0 ? F.g[ix[c]] * N_XPIX : NaN, 3)
    nflags = 0x0000
    ncl > 2 && (nflags |= ARCNIGHT_MULTIPOS)
    lay.timeok || (nflags |= ARCNIGHT_NO_TIME_TERMS)
    count(F.okf) < 250 && (nflags |= ARCNIGHT_FEW_FIBERS)
    any(rejected) && (nflags |= ARCNIGHT_SEQ_REJECTED)
    compared || (nflags |= ARCNIGHT_FEW_SEQ)
    abs_fail_night && (nflags |= ARCNIGHT_ABS_FAIL)
    chip_scales === nothing && (nflags |= ARCNIGHT_CHIPSCALE_FREE)
    any(posout) && (nflags |= ARCNIGHT_POS_OUTLIER)

    # per-chip residual summary of the fitted rows (mÅ)
    chip_med_mA = fill(NaN, N_CHIPS)
    chip_rsig_mA = fill(NaN, N_CHIPS)
    for c in 1:N_CHIPS
        r = [resid_all[i] for i in eachindex(L.id) if keep_main[i] && L.chip[i] == c]
        chip_med_mA[c] = 1000 * _arc_nanmedian(r)
        chip_rsig_mA[c] = 1000 * _arc_robstd(r)
    end
    wrms = sqrt(_arc_nanmean((F.resid[F.keep] .^ 2) .* full.D.w[F.keep]))

    return (; S, L, tseq, tau, tmean, dpx1, posout, igrp, n_clusters = ncl, grp,
        rej_quality, rej_shape, rej_chipoff, rej_round, compared, seq_rsig_px, seq_keep_frac,
        seq_shape_px, seq_chipoff_px,
        seq_score, seq_absfail, seq_chip_zp_px, seq_chip_stretch_px, reject_opts = ropt,
        used_seq = full.sidx, TH = F.TH, okf = F.okf, g = F.g, lay, grp_used = full.grp,
        tau_used = full.tau, resid = resid_all, keep_main, seqflags, lamp_off_px, nrows_seq,
        seq_shift_px, seq_tilt_px, pos2_zp_px = gpx(lay.idx_alpha), pos2_stretch_px = gpx(lay.idx_beta),
        time_zp_px = gpx(lay.idx_u), time_stretch_px = gpx(lay.idx_v), night_flags = nflags,
        chip_med_mA, chip_rsig_mA, wrms, n_rows = length(full.D.y), n_keep = count(F.keep),
        disp_chip = collect(disp_chip), TH_stage1 = F1.TH, g_stage1 = F1.g)
end

"""Wavelength grid (N_XPIX x N_FIBERS x 3) of a night fit at the state of final-fit sequence
`j` (default: the reference sequence)."""
function arclamp_night_grid(R; j = 1)
    G = fill(NaN, N_XPIX, N_FIBERS, N_CHIPS)
    d, t = arcfit_seq_state(R.g, R.lay, j)
    px = collect(1.0:N_XPIX)
    for f in 1:N_FIBERS
        R.okf[f] || continue
        for c in 1:N_CHIPS
            G[:, f, c] .= arcfit_eval(view(R.TH, f, :), R.g, R.lay, f, c, px;
                grpv = R.grp_used[j], d = d, t = t, tau = R.tau_used[j])
        end
    end
    return G
end

# ------------------------------------------------------------------------------------------
# I/O: arc sequences from the almanac, object references, the night product
# ------------------------------------------------------------------------------------------
"""
    arclamp_night_sequences(almanac_file, tele, mjd, apred_mjd_dir)

Arc sequences of a night: the almanac's `sequences/arclamps` blocks (or contiguous runs of arc
exposures when the almanac has none), split into (ThAr, UNe) pairs that share a commanded dither
label (the label pairs the two lamps; positions are measured afterwards). A pair is kept only if
the ar1Dcal files of both exposures exist for all three chips. Returns a vector of NamedTuples
(thar, une, dithpix_label).
"""
function arclamp_night_sequences(almanac_file::AbstractString, tele::AbstractString,
        mjd::Integer, apred_mjd_dir::AbstractString)
    df, blocks = h5open(almanac_file, "r") do f
        g = f["raw/$(tele)/$(mjd)"]
        df = DataFrame(read(g["exposures"]))
        blocks = (haskey(g, "sequences") && haskey(g["sequences"], "arclamps")) ?
                 Matrix{Int}(read(g["sequences/arclamps"])) : zeros(Int, 2, 0)
        df, blocks
    end
    exps = Int.(df.exposure)
    isarc = String.(df.image_type) .== "arclamp"
    if size(blocks, 2) == 0
        # contiguous runs of arc exposures
        runs = Tuple{Int, Int}[]
        i = 1
        while i <= length(exps)
            if isarc[i]
                j = i
                while j < length(exps) && isarc[j + 1]
                    j += 1
                end
                push!(runs, (exps[i], exps[j]))
                i = j + 1
            else
                i += 1
            end
        end
        blocks = isempty(runs) ? zeros(Int, 2, 0) : hcat([[a, b] for (a, b) in runs]...)
    end
    idx = Dict(exps[i] => i for i in eachindex(exps))
    has1d(e) = all(isfile(joinpath(apred_mjd_dir,
                       "ar1Dcal_$(tele)_$(mjd)_$(_arc_expstr(e))_$(chip)_arclamp.h5"))
                   for chip in CHIP_LIST)
    out = NamedTuple{(:thar, :une, :dithpix_label), Tuple{Int, Int, Float64}}[]
    for k in 1:size(blocks, 2)
        thar = Int[]
        une = Int[]
        for e in blocks[1, k]:blocks[2, k]
            haskey(idx, e) || continue
            i = idx[e]
            isarc[i] || continue
            if df.lamp_thar[i] == 1
                has1d(e) && push!(thar, e)
            elseif df.lamp_une[i] == 1
                has1d(e) && push!(une, e)
            end
        end
        (isempty(thar) || isempty(une)) && continue
        lab(e) = round(Float64(df.dithered_pixels[idx[e]]), digits = 3)
        for dl in unique(lab.(vcat(thar, une)))
            th = [e for e in thar if lab(e) == dl]
            un = [e for e in une if lab(e) == dl]
            (isempty(th) || isempty(un)) && continue
            push!(out, (thar = th[1], une = un[1], dithpix_label = dl))
        end
    end
    return out
end

"Mid-exposure MJD stored in an ar1Dcal file's metadata (NaN if absent)."
function arclamp_mjd_mid(fn::AbstractString)
    isfile(fn) || return NaN
    h5open(fn, "r") do f
        (haskey(f, "metadata") && haskey(f["metadata"], "mjd_mid_exposure")) ?
        Float64(read(f["metadata/mjd_mid_exposure"])) : NaN
    end
end

"""
    get_and_save_arclamp_night_wavecal(tele, mjd; apred_dir, almanac_file, outname = nothing,
        linelist = load_arclamp_linelist(), chip_transform_cal = nothing, chip_scales = nothing,
        chip_scale_source = nothing, nsig = 8.0, checkpoint_mode = "commit_same",
        disp_chip = nothing, reject = (;))

Stage-1 arc night product. Finds the night's arc sequences, runs the local-noise finder on the
first ThAr and first UNe exposure of each, associates with the frozen line list, fits the night
jointly with measured positions and arc-only sequence rejection, and
writes `arclampWaveCalNight_<tele>_<mjd>.h5` (default: in `apred_dir/<mjd>/`). Inputs are opened
read-only. Returns the output path, or `nothing` if the night has no usable arc sequence.
Chip scales: `chip_scales` (N_FIBERS x 2) if given, else `chip_scales(chip_transform_cal, tele,
mjd)` (the pipeline passes its loaded `--chip_transform_cal`), else fitted freely; the source is
recorded in the product (`chip_scale_source`). `reject` overrides thresholds of the arc-only sequence rejection (keys of `ARCREJ_DEFAULTS`);
the thresholds used are recorded in the product metadata. Only arc data and static calibration
inputs are read: this stage runs before the sky and FPI wavelength stages.
Nothing in the pipeline reads this product yet.
"""
function get_and_save_arclamp_night_wavecal(tele::AbstractString, mjd::Integer;
        apred_dir::AbstractString, almanac_file::AbstractString, outname = nothing,
        linelist::ArcLineList = load_arclamp_linelist(), chip_transform_cal = nothing,
        chip_scales = nothing, chip_scale_source = nothing,
        nsig = 8.0, checkpoint_mode = "commit_same", disp_chip = nothing,
        reject::NamedTuple = NamedTuple())
    # chip scales: an explicit matrix wins; else the chip-transform calibration's block for this
    # (tele, mjd); else fitted freely (flagged ARCNIGHT_CHIPSCALE_FREE)
    if chip_scales === nothing && chip_transform_cal !== nothing
        chip_scales = ApogeeReduction.chip_scales(chip_transform_cal, tele, mjd)
        chip_scale_source = something(chip_scale_source, chip_scales === nothing ?
            "free (no block for this night in $(chip_transform_id(chip_transform_cal)))" :
            "chip_transform_cal $(chip_transform_id(chip_transform_cal))")
    end
    chip_scale_source = something(chip_scale_source, chip_scales === nothing ? "free" : "fixed (given)")
    mdir = joinpath(apred_dir, string(mjd))
    outname = something(outname, joinpath(mdir, "arclampWaveCalNight_$(tele)_$(mjd).h5"))
    check_file(outname, mode = checkpoint_mode) && return outname
    seqs = arclamp_night_sequences(almanac_file, tele, mjd, mdir)
    isempty(seqs) && return nothing
    ar1d(e) = [joinpath(mdir, "ar1Dcal_$(tele)_$(mjd)_$(_arc_expstr(e))_$(chip)_arclamp.h5")
               for chip in CHIP_LIST]
    peaks = Dict{Int, Any}()
    for s in seqs, e in (s.thar, s.une)
        haskey(peaks, e) || (peaks[e] = arclamp_exposure_peaks(ar1d(e); nsig = nsig))
    end
    lines = [arclamp_sequence_lines(peaks[s.thar], peaks[s.une], linelist) for s in seqs]
    tth = [arclamp_mjd_mid(ar1d(s.thar)[1]) for s in seqs]
    tun = [arclamp_mjd_mid(ar1d(s.une)[1]) for s in seqs]
    tseq = [_arc_nanmean([a, b]) for (a, b) in zip(tth, tun)]
    R = fit_arclamp_night(lines, tseq, tele, linelist; chip_scales = chip_scales, reject = reject,
        disp_chip = something(disp_chip, arclamp_seed_dispersion(linelist)))
    write_arclamp_night_wavecal(outname, R, seqs, linelist; tele = tele, mjd = mjd,
        chip_scales = chip_scales, chip_scale_source = chip_scale_source, nsig = nsig)
    return outname
end

"""Write the night product (see `get_and_save_arclamp_night_wavecal`). Datasets:
  fiber_params (N_FIBERS x 9: a_R b_R a_B b_B p0..p4), fiber_ok, the global terms
  (global_params with its decoded per-sequence and per-chip parts), per-sequence tables (seq_*),
  per-sequence rejection statistics, flags, the per-chip residual summary, and the line table (lines_*).
  Wavelength grids are not stored (about 15 MB per night); `arclamp_night_wavelength` evaluates
  the solution from the stored parameters at any sequence state."""
function write_arclamp_night_wavecal(outname, R, seqs, ll::ArcLineList; tele, mjd,
        chip_scales, chip_scale_source, nsig)
    S = R.S
    used = falses(S)
    used[R.used_seq] .= true
    metadata = Dict{String, Any}(
        "product" => "arclampWaveCalNight", "stage" => 1, "tele" => String(tele),
        "mjd" => Int(mjd), "linelist_version" => ll.version,
        "linelist_fingerprint" => ll.fingerprint, "linelist_selection" => ll.selection,
        "n_fit_lines" => length(ll.fit_ids), "chip_scale_source" => String(chip_scale_source),
        "finder" => "get_arclamp_peaks_local", "finder_nsig" => Float64(nsig),
        "weight_floor_px" => 0.02, "group_px" => 0.25, "position_outlier_px" => 3.0,
        "sequence_rejection" => "arc-only: outliers among the night's sequences (arclamp_sequence_reject)",
        "inputs" => "arc ar1Dcal files, almanac, line list, chip scales; no sky/FPI products",
        "seq_flag_bits" => collect(String, ARCLAMP_SEQFLAG_BITS),
        "night_flag_bits" => collect(String, ARCLAMP_NIGHTFLAG_BITS),
        "dither_label_use" => "pairs ThAr with UNe within an almanac arc block only",
        "read_by" => "nothing (stage 1)")
    for (k, v) in pairs(R.reject_opts)
        metadata["reject_" * replace(String(k), "reject_" => "")] = v
    end
    safe_jldsave(outname, metadata;
        fiber_params = R.TH, fiber_ok = Vector{Bool}(R.okf),
        chip_scales = chip_scales === nothing ? fill(NaN, N_FIBERS, 2) : Matrix{Float64}(chip_scales),
        seq_thar_exp = [s.thar for s in seqs], seq_une_exp = [s.une for s in seqs],
        seq_dithpix_label = [s.dithpix_label for s in seqs], seq_mjd_mid = R.tseq,
        seq_tau = R.tau, seq_used = Vector{Bool}(used), seq_flags = R.seqflags,
        seq_stage1_shift_px = R.dpx1, seq_pos_cluster = R.igrp, seq_pos_group = R.grp,
        seq_shift_px = R.seq_shift_px, seq_tilt_px = R.seq_tilt_px,
        seq_rsig_px = R.seq_rsig_px, seq_keep_frac = R.seq_keep_frac, seq_shape_px = R.seq_shape_px,
        seq_chipoff_px = R.seq_chipoff_px,
        seq_chip_zp_px = R.seq_chip_zp_px, seq_chip_stretch_px = R.seq_chip_stretch_px,
        seq_reject_score = R.seq_score, seq_reject_round = R.rej_round,
        night_compared = R.compared, seq_lamp_offset_px = R.lamp_off_px,
        seq_n_rows = R.nrows_seq,
        n_pos_clusters = R.n_clusters, has_pos2 = R.lay.has2, has_time_terms = R.lay.timeok,
        pos2_zp_px = collect(R.pos2_zp_px), pos2_stretch_px = collect(R.pos2_stretch_px),
        time_zp_px = collect(R.time_zp_px), time_stretch_px = collect(R.time_stretch_px),
        tau_mjd_mean = R.tmean, global_params = R.g, night_flags = R.night_flags,
        n_good_fibers = count(R.okf), n_rows = R.n_rows, n_keep = R.n_keep,
        weighted_rms = R.wrms, chip_resid_median_mA = R.chip_med_mA,
        chip_resid_rsig_mA = R.chip_rsig_mA, disp_chip = R.disp_chip,
        lines_seq = Int8.(R.L.seq), lines_fiber = Int16.(R.L.fiber), lines_chip = Int8.(R.L.chip),
        lines_lamp = Int8.(R.L.lamp), lines_id = Int16.(R.L.id), lines_pix = Float32.(R.L.pix),
        lines_pixvar = Float32.(R.L.pixvar), lines_height = Float32.(R.L.height),
        lines_resid = Float32.(R.resid), lines_keep = Vector{Bool}(R.keep_main))
    return outname
end

"""
    arclamp_night_wavelength(fname; seq = nothing, px = 1.0:N_XPIX) -> Array (px, N_FIBERS, 3)

Evaluate a stored night product at the state of almanac-order sequence `seq` (default: the
reference sequence, i.e. the first sequence used in the final fit). Sequences not used in the
final fit (position outliers, rejected sequences) have no fitted state and raise an error.
"""
function arclamp_night_wavelength(fname::AbstractString; seq = nothing, px = 1.0:N_XPIX)
    TH, okf, g, used, grp, tau, has2 = h5open(fname, "r") do f
        read(f["fiber_params"]), Bool.(read(f["fiber_ok"])), read(f["global_params"]),
        findall(Bool.(read(f["seq_used"]))), read(f["seq_pos_group"]), read(f["seq_tau"]),
        Bool(read(f["has_pos2"]))
    end
    s = seq === nothing ? used[1] : seq
    j = findfirst(==(s), used)
    j === nothing && throw(ArgumentError("sequence $(s) is not in the final fit of $(fname)"))
    grp_u = has2 ? Int.(grp[used]) : ones(Int, length(used))
    lay = ArcGlobalLayout(length(used), grp_u, tau[used])
    lay.n == length(g) ||
        throw(ArgumentError("$(fname): stored global terms do not match the rebuilt layout"))
    d, t = arcfit_seq_state(g, lay, j)
    pxv = collect(Float64, px)
    G = fill(NaN, length(pxv), N_FIBERS, N_CHIPS)
    for f in 1:N_FIBERS
        okf[f] || continue
        for c in 1:N_CHIPS
            G[:, f, c] .= arcfit_eval(view(TH, f, :), g, lay, f, c, pxv; grpv = grp_u[j], d = d,
                t = t, tau = tau[used][j])
        end
    end
    return G
end
