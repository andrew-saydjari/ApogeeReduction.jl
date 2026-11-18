# ---------------------------------------------------------------------------------------------
# Per-exposure registration of the night trace to the exposure's own light (branch trace-reg).
#
# The night's trace comes from one flat. The slit moves during the night (APO: ~0.05-0.1 px between
# the evening science frames and the end-of-night flat), so before extraction the night trace can be
# moved rigidly onto the light of the exposure being extracted:
#
#     y_f(x) -> y_f(x) + a + b u + c v,   u = (x - 1024) / 1024,  v = (y_f(x) - 1024) / 1024
#
# a = shift (px), b = rotation (px at the chip edge in x), c = scale (px at the slit edge). The
# parameters are measured from the science light itself: 32 bands of 65 columns (5-sigma clipped
# column means, as the T2 level), every fiber's light center by the production Gauss-Hermite
# refinement (all 300 fibers modelled jointly, started at the night trace), then a weighted robust
# fit over the fiber-bands with S/N >= `snr_min` (stars and sky alike; faint and empty fibers drop
# out by S/N). Nothing recorded by the instrument (dither, focus, any mechanism readback) is used.
# Half splits (even/odd fibers) give an empirical noise estimate that the gate checks.
#
# OFF by default (`pipeline_2d_1d.jl --trace_register off`). When on, a refused registration falls
# back to the night trace and sets TRACE_EXPOSURE_REG_FAILED; an applied one sets
# TRACE_EXPOSURE_REGISTERED. Both are exposure-level bits replicated into every fiber of
# `bitmsk_trace`, and the parameters are written into the 1D metadata (`trace_expreg_*`).
# ---------------------------------------------------------------------------------------------
using LinearAlgebra: det, inv, diag

"Minimum S/N (band height / its error) of a fiber-band used by the exposure registration."
const TRACE_EXPREG_SNR_MIN = 20.0
"Minimum number of fiber-bands (of 32 x 300) with S/N >= TRACE_EXPREG_SNR_MIN."
const TRACE_EXPREG_NMIN = 200
"Minimum number of distinct fibers contributing."
const TRACE_EXPREG_NFIB_MIN = 20
"""
Largest accepted |shift| (px). Measured on 4,496 science chip-frames of 52 nights (trace-reg
study, 2026-11-18, pass-1 traceMain): p99 |a| = 0.25 px; the 7 frames above 0.5 px are all LCO 57802
(commissioning). A larger value means a wrong night trace or a frame that is not a normal science
exposure, which registration must not paper over.
"""
const TRACE_EXPREG_AMAX = 0.5
"Largest accepted |rotation| and |scale| (px at the chip / slit edge)."
const TRACE_EXPREG_BMAX = 0.3
"Largest accepted half-split disagreement |a_even - a_odd| / 2 (px): an empirical noise check.
Same study: median 0.0006-0.0026 px, p99 0.010 px; 9 of 4,496 chip-frames above 0.02 px."
const TRACE_EXPREG_HALF_MAX = 0.02
"Weight floor (px) for a fiber-band: fibers do not move together at the 0.01 px level, so no single bright fiber dominates."
const TRACE_EXPREG_SIG_FLOOR = 0.01

"""
    exposure_trace_offsets(image, ivar, pix_bitmask, trace_params, P; nb = 32, half = 32,
        nsigma = 5, maxmove = 1.0) -> NamedTuple of vectors (b, f, u, v, r, sig, snr)

Per band (`nb` bands of `2 half + 1` columns) and fiber: the light center (Gauss-Hermite
refinement of all fibers jointly, started at the trace averaged over the band's columns) minus
that trace position (`r`), its formal error (`sig`) and the band S/N (`snr`).
"""
function exposure_trace_offsets(image_data, ivar_image, pix_bitmask, trace_params, P::SlitProfile;
        nb = 32, half = 32, nsigma = 5, maxmove = 1.0)
    nx = size(trace_params, 1)
    good = ((pix_bitmask .& bad_pix_bits) .== 0) .& (ivar_image .> 0) .& isfinite.(image_data)
    xc, prof, err = slit_bands(image_data, 1 ./ sqrt.(ivar_image), good; nb = nb, half = half,
        nsigma = nsigma)
    nf = size(trace_params, 2)
    B = Int[]; F = Int[]; U = Float64[]; V = Float64[]; R = Float64[]; S = Float64[]; SN = Float64[]
    for b in 1:nb
        c = clamp(round(Int, xc[b]), half + 1, nx - half)
        cols = (c - half):(c + half)
        mu0 = [mean(view(trace_params, cols, f, 2)) for f in 1:nf]
        w0 = [mean(view(trace_params, cols, f, 3)) for f in 1:nf]
        all(isfinite, mu0) && all(isfinite, w0) || continue
        _, h, _, _, he = slit_heights_chi2(P, prof[b, :], _band_weights(err[b, :]), mu0, w0)
        np, cv = slit_gh_refine(P, prof[b, :], err[b, :], mu0, w0, h; max_center_move = maxmove,
            minw = 0.6, maxw = 1.6, n_iter = 10, return_cov = true)
        for f in 1:nf
            v2 = cv[f, 2, 2]
            (isfinite(v2) && v2 > 0 && isfinite(np[f, 2]) && isfinite(h[f]) && isfinite(he[f]) && he[f] > 0) || continue
            push!(B, b); push!(F, f); push!(U, (c - 1024) / 1024); push!(V, (mu0[f] - 1024) / 1024)
            push!(R, np[f, 2] - mu0[f]); push!(S, sqrt(v2)); push!(SN, h[f] / he[f])
        end
    end
    return (b = B, f = F, u = U, v = V, r = R, sig = S, snr = SN)
end

# weighted least squares with 4-sigma clipping on normalized residuals
function _expreg_wfit(A, y, w)
    k = trues(length(y))
    c = zeros(size(A, 2))
    for _ in 1:4
        sw = sqrt.(w[k])
        c = (A[k, :] .* sw) \ (y[k] .* sw)
        z = (y .- A * c) .* sqrt.(w)
        s = 1.4826 * median(abs.(z[k]))
        k .= abs.(z) .< 4 * max(s, 1e-6)
        count(k) > size(A, 2) || (k .= true; break)
    end
    sw = sqrt.(w[k])
    Ak = A[k, :] .* sw
    z = (y[k] .- A[k, :] * c) .* sw
    chi2nu = sum(abs2, z) / max(1, count(k) - size(A, 2))
    F = Ak' * Ak
    e = det(F) > 0 ? sqrt.(abs.(diag(inv(F)) .* max(chi2nu, 1.0))) : fill(Inf, size(A, 2))
    return c, e, k
end

"Design matrix columns of the registration model (:shift, :shift_rot or :shift_rot_scale)."
function _expreg_design(model::Symbol, u, v)
    model == :shift && return reshape(ones(length(u)), :, 1)
    model == :shift_rot && return hcat(ones(length(u)), u)
    model == :shift_rot_scale && return hcat(ones(length(u)), u, v)
    error("exposure registration model must be :shift, :shift_rot or :shift_rot_scale, got $(model)")
end

"""
    fit_exposure_registration(m; model = :shift_rot, snr_min, sig_floor, maxabs = 0.6,
        nmin, nfib_min, amax, bmax, half_max) -> NamedTuple

Fit the rigid model to the offsets `m` of `exposure_trace_offsets`. Returns `outcome` ("ok" or
the refusal: "too_few", "too_large", "inconsistent_halves", "nonfinite"), `a`, `b`, `c` (absent
terms are 0), their formal errors, `n` fiber-bands and `nfib` fibers used, `rms` (robust rms of
the residuals, px), `half_a` (|a_even - a_odd| / 2 over even/odd fibers) and `model`.
"""
function fit_exposure_registration(m; model::Symbol = :shift_rot, snr_min = TRACE_EXPREG_SNR_MIN,
        sig_floor = TRACE_EXPREG_SIG_FLOOR, maxabs = 0.6, nmin = TRACE_EXPREG_NMIN,
        nfib_min = TRACE_EXPREG_NFIB_MIN, amax = TRACE_EXPREG_AMAX, bmax = TRACE_EXPREG_BMAX,
        half_max = TRACE_EXPREG_HALF_MAX)
    sel = findall((m.snr .>= snr_min) .& (abs.(m.r) .< maxabs))
    nfib = length(unique(m.f[sel]))
    res(outcome; a = 0.0, b = 0.0, c = 0.0, ea = NaN, eb = NaN, ec = NaN, rms = NaN, half_a = NaN) =
        (outcome = outcome, model = model, a = a, b = b, c = c, a_err = ea, b_err = eb, c_err = ec,
         n = length(sel), nfib = nfib, rms = rms, half_a = half_a)
    (length(sel) >= nmin && nfib >= nfib_min) || return res("too_few")
    u = m.u[sel]; v = m.v[sel]; y = m.r[sel]
    w = 1 ./ (m.sig[sel] .^ 2 .+ sig_floor^2)
    A = _expreg_design(model, u, v)
    co, er, k = _expreg_wfit(A, y, w)
    np = length(co)
    a = co[1]; b = np >= 2 ? co[2] : 0.0; c = np >= 3 ? co[3] : 0.0
    ea = er[1]; eb = np >= 2 ? er[2] : 0.0; ec = np >= 3 ? er[3] : 0.0
    rms = 1.4826 * median(abs.((y .- A * co)[k]))
    fe = findall(iseven.(m.f[sel])); fo = findall(isodd.(m.f[sel]))
    half_a = if length(fe) >= nmin ÷ 2 && length(fo) >= nmin ÷ 2
        ce, _, _ = _expreg_wfit(A[fe, :], y[fe], w[fe])
        cO, _, _ = _expreg_wfit(A[fo, :], y[fo], w[fo])
        abs(ce[1] - cO[1]) / 2
    else
        NaN
    end
    kw = (a = a, b = b, c = c, ea = ea, eb = eb, ec = ec, rms = rms, half_a = half_a)
    all(isfinite, (a, b, c)) || return res("nonfinite"; kw...)
    (abs(a) <= amax && abs(b) <= bmax && abs(c) <= bmax) || return res("too_large"; kw...)
    (isfinite(half_a) && half_a <= half_max) || return res("inconsistent_halves"; kw...)
    return res("ok"; kw...)
end

"""
    register_exposure_trace(image, ivar, pix_bitmask, trace_params, P; model = :shift_rot, kw...)

Measure (`exposure_trace_offsets`) and fit (`fit_exposure_registration`) the registration of
`trace_params` (the night trace, n_x x n_fibers x 3) to one exposure. Keywords go to the fit.
"""
function register_exposure_trace(image_data, ivar_image, pix_bitmask, trace_params, P::SlitProfile;
        model::Symbol = :shift_rot, nb = 32, half = 32, kw...)
    m = exposure_trace_offsets(image_data, ivar_image, pix_bitmask, trace_params, P; nb = nb, half = half)
    return fit_exposure_registration(m; model = model, kw...)
end

"""
    apply_exposure_registration(trace_params, reg) -> trace_params (copy)

The night trace moved by the registration: centers `y + a + b u + c v` (v from the unmoved
center); heights and widths unchanged. A refused registration (`reg.outcome != "ok"`) returns an
unchanged copy.
"""
function apply_exposure_registration(trace_params::AbstractArray{<:Real, 3}, reg)
    T = copy(trace_params)
    reg.outcome == "ok" || return T
    for x in axes(T, 1), f in axes(T, 2)
        y = trace_params[x, f, 2]
        T[x, f, 2] = y + reg.a + reg.b * (x - 1024) / 1024 + reg.c * (y - 1024) / 1024
    end
    return T
end

"Metadata block of an exposure registration (written into the 1D products)."
function exposure_registration_metadata(reg)
    return Dict{String, Any}("trace_expreg_outcome" => String(reg.outcome),
        "trace_expreg_model" => String(reg.model),
        "trace_expreg_a" => Float64(reg.a), "trace_expreg_b" => Float64(reg.b),
        "trace_expreg_c" => Float64(reg.c), "trace_expreg_a_err" => Float64(reg.a_err),
        "trace_expreg_b_err" => Float64(reg.b_err), "trace_expreg_c_err" => Float64(reg.c_err),
        "trace_expreg_n" => Int(reg.n), "trace_expreg_nfib" => Int(reg.nfib),
        "trace_expreg_rms" => Float64(reg.rms), "trace_expreg_half_a" => Float64(reg.half_a))
end

"Exposure-level trace bit of a registration result (`nothing` = registration off)."
exposure_registration_bits(reg) = isnothing(reg) ? zero(TraceBits) :
                                   (reg.outcome == "ok" ? TRACE_EXPOSURE_REGISTERED : TRACE_EXPOSURE_REG_FAILED)
