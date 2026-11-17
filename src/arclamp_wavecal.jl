# Arc-lamp (ThAr + UNe) wavelength calibration, stage 1: frozen line list, local-noise peak
# finder, association, joint night fit with measured per-sequence positions, and the 0.02 px
# sequence-position gate. Writes one per-night product (`arclampWaveCalNight_<tele>_<mjd>.h5`)
# that NOTHING in the pipeline reads yet. It runs only when `pipeline_2d_1d.jl --arclamp_wavecal
# true` is given (default false), so production results are unchanged.
#
# Provenance of the methods (offline studies, 2026-11-16/17):
#   * finder: the "N=8" local-noise matched-filter finder (`get_arclamp_peaks_local`). It works in
#     detector pixels only and uses the production `fit_gauss_and_bias` for the centroids.
#   * association: catalog scan + quadratic self-calibration (the `associate` of the studies).
#   * joint night fit (Gauss-Newton with a Schur complement), position clustering and gate B:
#     the best arc fit of the ARC_VERDICT2 study, validated there on 1,563 nights.
#
# Positions are MEASURED from each sequence's own lines (and, for the gate, from the object
# exposures' own solutions). The recorded dither (DITHPIX) is used for one thing only: to pair the
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
# positions and the gate
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

"""
    arclamp_gate(offset_px; gate_px = 0.02) -> (dev_px, gated, flagged)

Gate B: `offset_px[s]` is sequence s's measured offset (px) from the nearest-in-time object
exposure at the same measured position (NaN if none). The night median is removed; a sequence
is flagged when |deviation| > `gate_px`. Ungated sequences (NaN) are never flagged.
"""
function arclamp_gate(offset_px::AbstractVector; gate_px = 0.02)
    dev = offset_px .- _arc_nanmedian(offset_px)
    gated = isfinite.(dev)
    flagged = gated .& (abs.(dev) .> gate_px)
    return dev, gated, flagged
end

"""Object-exposure reference for the gate: a (pixel window x fibers 100:200) slice of the
exposure's chip-G wavelength solution, with its mid-exposure MJD."""
struct ArcObjRef
    expnum::Int
    mjd_mid::Float64
    kind::String                     # "fpi" or "sky" (which per-exposure solution file)
    px0::Int                         # first pixel of the slice
    fibers::UnitRange{Int}
    wave::Matrix{Float64}            # (pixels, fibers) chip G
end

const ARCGATE_FIBERS = 100:200
const ARCGATE_PX = collect(900.0:25.0:1150.0)

"""Measured offset (px) of arc state (`th`, `g`, ...) from an object exposure's chip-G
solution: median over fibers 100-200 and pixels 900-1150 of (λ_arc − λ_obj)/(dλ_obj/dpix)."""
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

# ------------------------------------------------------------------------------------------
# the night fit
# ------------------------------------------------------------------------------------------
# per-sequence flag bits
const ARCSEQ_POS_OUTLIER = 1        # > outlier_px from the night median: excluded
const ARCSEQ_GATE_FAIL = 2          # |gate deviation| > gate_px: excluded from the final fit
const ARCSEQ_UNGATED = 4            # no same-position object exposure to gate against
const ARCSEQ_LAMP_MISMATCH = 8      # ThAr and UNe rows disagree in position by > 0.1 px
const ARCSEQ_FEW_LINES = 16         # < 1000 fitted rows in the sequence
# night flag bits
const ARCNIGHT_MULTIPOS = 1         # > 2 measured position clusters: no shared position terms
const ARCNIGHT_NO_TIME_TERMS = 2    # time terms not identifiable
const ARCNIGHT_FEW_FIBERS = 4       # < 250 good fibers
const ARCNIGHT_ALL_GATED = 8        # every gated sequence failed the gate
const ARCNIGHT_NO_GATE = 16         # no sequence could be gated (no object references)
const ARCNIGHT_CHIPSCALE_FREE = 32  # chip scales fitted (no chip-transform calibration given)
const ARCNIGHT_POS_OUTLIER = 64     # at least one sequence was a position outlier

"""
    fit_arclamp_night(lines, tseq, tele, ll; objrefs = ArcObjRef[], chip_scales = nothing, ...)

Fit one night's arc sequences jointly.

  * `lines`: per sequence s, a line table from `arclamp_sequence_lines`.
  * `tseq`: per sequence, the mean mid-exposure MJD of its ThAr and UNe exposures.
  * `objrefs`: object-exposure references for the gate (may be empty).
  * `chip_scales`: N_FIBERS x 2 (b_R, b_B) to hold the chip scales fixed, or `nothing`.

Stage 1 (shift + tilt per sequence only) measures each sequence's position; the positions are
clustered; gate B is evaluated against the object references; the final fit uses the
non-outlier, non-flagged sequences with the second-position and time terms when identifiable.
Returns a NamedTuple with everything the night product stores.
"""
function fit_arclamp_night(lines::AbstractVector, tseq::AbstractVector{<:Real}, tele::AbstractString,
        ll::ArcLineList; objrefs::AbstractVector{ArcObjRef} = ArcObjRef[],
        chip_scales = nothing, gate_px = 0.02, group_px = 0.25, outlier_px = 3.0,
        floor_px = 0.02, disp_chip = arclamp_seed_dispersion(ll), lamp_mismatch_px = 0.1)
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

    # gate B against the object references (measured positions only)
    function offsets(FF)
        off = fill(NaN, S)
        obj = zeros(Int, S)
        for (j, s) in enumerate(FF.sidx)
            d, t = arcfit_seq_state(FF.F.g, FF.lay, j)
            best = 0
            bt = Inf
            bo = NaN
            for (ie, R) in enumerate(objrefs)
                o = arclamp_offset_to_ref(FF.F.TH, FF.F.okf, FF.F.g, FF.lay, FF.grp[j], d, t,
                    FF.tau[j], R)
                (isfinite(o) && abs(o) < group_px) || continue
                dt = abs(R.mjd_mid - tseq[s])
                if dt < bt
                    bt = dt
                    best = ie
                    bo = o
                end
            end
            best == 0 && continue
            off[s] = bo
            obj[s] = objrefs[best].expnum
        end
        return off, obj
    end
    off_all, gobj = offsets(full_all)
    gdev, gated, flag = arclamp_gate(off_all; gate_px = gate_px)
    # flagged sequences are refitted out; if nothing would remain, the all-sequence fit is kept
    useseq = .!flag .& .!posout
    all_gated = any(gated) && !any(gated .& .!flag)
    full = (count(flag) > 0 && count(useseq) >= 1) ? fit_subset(useseq) : full_all
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
    seqflags = zeros(Int, S)
    lamp_off_px = fill(NaN, S)
    nrows_seq = zeros(Int, S)
    for s in 1:S
        posout[s] && (seqflags[s] |= ARCSEQ_POS_OUTLIER)
        flag[s] && (seqflags[s] |= ARCSEQ_GATE_FAIL)
        gated[s] || (seqflags[s] |= ARCSEQ_UNGATED)
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
    nflags = 0
    ncl > 2 && (nflags |= ARCNIGHT_MULTIPOS)
    lay.timeok || (nflags |= ARCNIGHT_NO_TIME_TERMS)
    count(F.okf) < 250 && (nflags |= ARCNIGHT_FEW_FIBERS)
    all_gated && (nflags |= ARCNIGHT_ALL_GATED)
    any(gated) || (nflags |= ARCNIGHT_NO_GATE)
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
        gate_offset_px = off_all, gate_dev_px = gdev, gate_obj = gobj, gated, gate_flag = flag,
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
    arclamp_object_refs(apred_mjd_dir, tele, mjd)

Gate references: for every object exposure with a per-exposure wavelength solution
(`waveCalNightfpiDither_*_object.h5`, else `waveCalNightskyDither_*_object.h5`), a chip-G slice
of its solution (pixels 880-1170, fibers 100-200) and its mid-exposure MJD. Opened read-only.
"""
function arclamp_object_refs(apred_mjd_dir::AbstractString, tele::AbstractString, mjd::Integer)
    isdir(apred_mjd_dir) || return ArcObjRef[]
    fs = readdir(apred_mjd_dir)
    found = Dict{Int, Tuple{String, String}}()
    for (kind, pre) in (("sky", "waveCalNightskyDither"), ("fpi", "waveCalNightfpiDither"))
        rx = Regex("^$(pre)_$(tele)_$(mjd)_(\\d{4})_object\\.h5\$")
        for fn in fs
            m = match(rx, fn)
            m === nothing && continue
            found[parse(Int, m.captures[1])] = (kind, joinpath(apred_mjd_dir, fn))   # fpi wins
        end
    end
    px0, px1 = 880, 1170
    refs = ArcObjRef[]
    for e in sort(collect(keys(found)))
        kind, fn = found[e]
        W = try
            h5open(fn, "r") do f
                Float64.(f["chipWaveSoln"][px0:px1, ARCGATE_FIBERS, 2])
            end
        catch
            continue
        end
        t = arclamp_mjd_mid(joinpath(apred_mjd_dir,
            "ar1Dcal_$(tele)_$(mjd)_$(_arc_expstr(e))_$(FIRST_CHIP)_object.h5"))
        push!(refs, ArcObjRef(e, t, kind, px0, ARCGATE_FIBERS, W))
    end
    return refs
end

"""
    get_and_save_arclamp_night_wavecal(tele, mjd; apred_dir, almanac_file, outname = nothing,
        linelist = load_arclamp_linelist(), chip_scales = nothing, chip_scale_source = "free",
        nsig = 8.0, checkpoint_mode = "commit_same", disp_chip = nothing)

Stage-1 arc night product. Finds the night's arc sequences, runs the local-noise finder on the
first ThAr and first UNe exposure of each, associates with the frozen line list, fits the night
jointly with measured positions and gate B (against the night's object-exposure solutions), and
writes `arclampWaveCalNight_<tele>_<mjd>.h5` (default: in `apred_dir/<mjd>/`). Inputs are opened
read-only. Returns the output path, or `nothing` if the night has no usable arc sequence.
Nothing in the pipeline reads this product yet.
"""
function get_and_save_arclamp_night_wavecal(tele::AbstractString, mjd::Integer;
        apred_dir::AbstractString, almanac_file::AbstractString, outname = nothing,
        linelist::ArcLineList = load_arclamp_linelist(), chip_scales = nothing,
        chip_scale_source::AbstractString = chip_scales === nothing ? "free" : "fixed",
        nsig = 8.0, checkpoint_mode = "commit_same", disp_chip = nothing)
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
    refs = arclamp_object_refs(mdir, tele, mjd)
    R = fit_arclamp_night(lines, tseq, tele, linelist; objrefs = refs, chip_scales = chip_scales,
        disp_chip = something(disp_chip, arclamp_seed_dispersion(linelist)))
    write_arclamp_night_wavecal(outname, R, seqs, linelist; tele = tele, mjd = mjd,
        chip_scales = chip_scales, chip_scale_source = chip_scale_source, nsig = nsig,
        n_objrefs = length(refs))
    return outname
end

"""Write the night product (see `get_and_save_arclamp_night_wavecal`). Datasets:
  fiber_params (N_FIBERS x 9: a_R b_R a_B b_B p0..p4), fiber_ok, the global terms
  (global_params with its decoded per-sequence and per-chip parts), per-sequence tables (seq_*),
  gate quantities, flags, the per-chip residual summary, and the line table (lines_*).
  Wavelength grids are not stored (about 15 MB per night); `arclamp_night_wavelength` evaluates
  the solution from the stored parameters at any sequence state."""
function write_arclamp_night_wavecal(outname, R, seqs, ll::ArcLineList; tele, mjd,
        chip_scales, chip_scale_source, nsig, n_objrefs)
    S = R.S
    used = falses(S)
    used[R.used_seq] .= true
    metadata = Dict{String, Any}(
        "product" => "arclampWaveCalNight", "stage" => 1, "tele" => String(tele),
        "mjd" => Int(mjd), "linelist_version" => ll.version,
        "linelist_fingerprint" => ll.fingerprint, "linelist_selection" => ll.selection,
        "n_fit_lines" => length(ll.fit_ids), "chip_scale_source" => String(chip_scale_source),
        "finder" => "get_arclamp_peaks_local", "finder_nsig" => Float64(nsig),
        "weight_floor_px" => 0.02, "gate" => "B", "gate_px" => 0.02, "group_px" => 0.25,
        "position_outlier_px" => 3.0, "n_object_refs" => n_objrefs,
        "dither_label_use" => "pairs ThAr with UNe within an almanac arc block only",
        "read_by" => "nothing (stage 1)")
    safe_jldsave(outname, metadata;
        fiber_params = R.TH, fiber_ok = Vector{Bool}(R.okf),
        chip_scales = chip_scales === nothing ? fill(NaN, N_FIBERS, 2) : Matrix{Float64}(chip_scales),
        seq_thar_exp = [s.thar for s in seqs], seq_une_exp = [s.une for s in seqs],
        seq_dithpix_label = [s.dithpix_label for s in seqs], seq_mjd_mid = R.tseq,
        seq_tau = R.tau, seq_used = Vector{Bool}(used), seq_flags = R.seqflags,
        seq_stage1_shift_px = R.dpx1, seq_pos_cluster = R.igrp, seq_pos_group = R.grp,
        seq_shift_px = R.seq_shift_px, seq_tilt_px = R.seq_tilt_px,
        seq_gate_offset_px = R.gate_offset_px, seq_gate_dev_px = R.gate_dev_px,
        seq_gate_obj = R.gate_obj, seq_lamp_offset_px = R.lamp_off_px,
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
final fit (position outliers, gate failures) have no fitted state and raise an error.
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
