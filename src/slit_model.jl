# Slit model: per-epoch trace templates registered to a frame (PR B of the robust-trace design,
# TRACE_ROBUST_V2.md §4-§7; production port of the study's prototype `RobustTrace.jl`).
#
# Model of the fiber centers on a frame:
#     y_f(x) = T_f(x) + d(u) + s(u) * (T_f(x) - 1024) + delta_f(x),    u = (x - 1024) / 1024
#   T_f(x)    epoch template (a calibration product, one per telescope, chip and epoch)
#   d(u)      shift, quadratic in u (offset, rotation, bow);  s(u) scale, linear in u
#   delta_f   per-fiber departure, measured only on lit fibers of a flat (T0)
#
# Three uses:
#   * T0, a flat of the night: `trace_extract_template` registers the template (whole-slit chi^2
#     with free per-fiber heights), checks the fiber labels (throughput fingerprint and light at
#     the slit edges; inside an epoch the label shift is fixed at 0 and a k != 0 preference is a
#     refusal, never a relabel), then runs production's per-column Gauss-Hermite refinement
#     (`fit_gaussians`) started from the registered template IN EVERY COLUMN. There is no outward
#     column walk, so one bad column cannot propagate (the LCO 60702 G/R twist).
#   * T1/T2, another frame of the night (dome flat, arcs, FPI, stacked science frames):
#     `register_frame` in registration-only mode (k fixed at 0, |d - d_epoch| <= 3 px); the
#     trace is the registered template with no per-fiber delta.
#   * T3, no usable frame: the template at the epoch-median registration.
#
# Profile shape and refinement are production's (`load_trace_hyperparams` tables,
# `pdf_func_mult`, `fit_gaussians` in traceExtract_GH.jl); band profiles are PR A's
# `trace_center_band`. Nothing here reads or writes files.

using LinearAlgebra: SymTridiagonal, I
using Statistics: median, mean, quantile, std, cor
using Polynomials: fit, Polynomial

# ---------------------------------------------------------------------------------------------
# Thresholds. Values fixed by the design reports are cited; the others are switches (keyword
# arguments) and listed in the PR as open choices.
# ---------------------------------------------------------------------------------------------
"Per-fiber departure from the registered template above which `TRACE_FIBER_DELTA` is set (px; V2 §1.4)."
const TRACE_FIBER_DELTA_PX = 0.3
"Registration-only frames: |d - d_epoch| search half-range (px; V2 §5.2 T1)."
const TRACE_REG_DMAX = 3.0
"Registration-only gate: robust RMS of band-refined centers about a per-fiber cubic (px; study s16/s19 gate)."
const TRACE_REG_FIT_RMS_MAX = 0.08
"Registration-only gate: minimum number of fibers with light in the bands (study s16 gate)."
const TRACE_REG_NLIT_MIN = 100
"Full registration (T0) search half-range around the epoch-median shift (px; V2 §6.3, prototype dmax)."
const TRACE_TEMPLATE_DMAX = 40.0
"Label shifts tested by the identification (prototype kmax)."
const TRACE_TEMPLATE_KMAX = 3
"Throughput-fingerprint margin above which it decides the label (prototype; else chi^2 / edge light)."
const TRACE_TEMPLATE_TP_MARGIN = 0.05

# ---------------------------------------------------------------------------------------------
# Templates: `TraceTemplate` (src/trace_template_cal.jl) holds, per fiber, centers and widths as
# polynomials in u = (x - 1024) / 1024 (`center_poly`, `width_poly`, evaluated by
# `template_centers` / `template_widths`), the relative throughput `height_rel` (the
# identification fingerprint), `lit` (false = dark on every flat of the epoch) and the
# epoch-median registration (`reg_d_median`, `reg_s_median`).
# ---------------------------------------------------------------------------------------------
nfibers(t) = size(t.center_poly, 2)

_polyval(c::AbstractVector, u) = (v = 0.0; for k in length(c):-1:1; v = v * u + c[k]; end; v)

"""
    TemplateRegistration

`d(u) = d[1] + d[2] u + d[3] u^2`, `s(u) = s[1] + s[2] u`. Applied to the template:
`y_f(x) = T_f(x) + d(u) + s(u) (T_f(x) - 1024)`.
"""
struct TemplateRegistration
    d::Vector{Float64}
    s::Vector{Float64}
end
TemplateRegistration(v::AbstractVector) = TemplateRegistration(Float64[v[1], v[2], v[3]], Float64[v[4], v[5]])
reg_vector(r::TemplateRegistration) = vcat(r.d, r.s)
_finite0(v) = isfinite(v) ? Float64(v) : 0.0
"The epoch-median registration of a template (constant shift and scale; level T3)."
epoch_median_registration(t) = TemplateRegistration([_finite0(t.reg_d_median), 0.0, 0.0],
    [_finite0(t.reg_s_median), 0.0])

"Registered template centers at columns `xs`: (length(xs), n_fibers)."
function registered_centers(t, r::TemplateRegistration, xs::AbstractVector)
    T = template_centers(t, collect(Float64, xs))
    for (i, x) in enumerate(xs)
        u = (x - 1024) / 1024
        T[i, :] .= T[i, :] .+ _polyval(r.d, u) .+ _polyval(r.s, u) .* (T[i, :] .- 1024)
    end
    return T
end

"The template identity string recorded in products (as `trace_template_id`)."
template_id(t) = "$(t.tele)_$(t.chip)_$(t.epoch):$(t.fingerprint)"

# ---------------------------------------------------------------------------------------------
# Production profiles
# ---------------------------------------------------------------------------------------------
"The production Gauss-Hermite profile tables of one telescope and chip (`load_trace_hyperparams`)."
struct SlitProfile
    x_prof_min::Float64
    x_prof_max_ind::Int
    n_sub::Int
    min_prof_fib::Int
    max_prof_fib::Int
    all_y_prof::Matrix{Float64}
    all_y_prof_deriv::Matrix{Float64}
end
SlitProfile(hp::Tuple) = SlitProfile(hp[2], hp[3], hp[4], hp[5], hp[6], hp[7], hp[8])
slit_profile(tele, chip; data_path = "./data/") = SlitProfile(get_default_trace_hyperparams(tele, chip;
    profile_path = data_path, plot_path = tempdir()))

_hp(P::SlitProfile) = (P.x_prof_min, P.x_prof_max_ind, P.n_sub, P.min_prof_fib, P.all_y_prof,
    P.all_y_prof_deriv)
_clampfib(P::SlitProfile, fibs) = clamp.(fibs, P.min_prof_fib, P.max_prof_fib)

"""
Pixel-integrated unit profiles (production `pdf_func_mult`) of fibers `fibs` at centers `mu`,
widths `w`, on pixels `y0[f]:y0[f]+2half`. Returns `(A, y0)`, `A` is (n, 2half+1).
"""
function slit_unit_profiles(P::SlitProfile, mu, w, fibs; half = 6)
    y0 = round.(Int, mu) .- half
    X = y0 .+ (0:(2half))'
    A = pdf_func_mult(X, mu, w, _clampfib(P, fibs), _hp(P)...)
    return A, y0
end

# ---------------------------------------------------------------------------------------------
# Band profiles
# ---------------------------------------------------------------------------------------------
"""
    slit_bands(image_data, noise_image, good_pixels; nb = 32, half = 8, nsigma = 5)

`nb` spatial profiles at band centers `xc = 32, 96, ..., 2016` (`(b - 1) * 2048/nb + 1024/nb`),
each the collapse of columns `xc ± half` by `trace_center_band` (PR A; outlier-clipped,
inverse-variance weighted). Returns `(xc, prof, err)`, `prof` and `err` are (nb, 2048); `err` is
`Inf` where no pixel survived and `prof` is 0 there.
"""
function slit_bands(image_data, noise_image, good_pixels; nb = 32, half = 8, nsigma = 5)
    nx = size(image_data, 1)
    w = nx ÷ nb
    xc = [(b - 1) * w + w ÷ 2 for b in 1:nb]
    ny = size(image_data, 2)
    prof = zeros(nb, ny)
    err = fill(Inf, nb, ny)
    for b in 1:nb
        c = clamp(xc[b], half + 1, nx - half)
        cf, ce, _ = trace_center_band(image_data, noise_image, good_pixels; mid = c,
            n_center_cols = half, nsigma = nsigma)
        prof[b, :] .= vec(cf)
        err[b, :] .= vec(ce)
    end
    return Float64.(xc), prof, err
end

"""
Combine the band profiles of several frames (e.g. the night's science frames) by the per-pixel
median; the error is 1.2533 x the rms error over frames / sqrt(n).
"""
function stack_bands(profs::AbstractVector, errs::AbstractVector)
    n = length(profs)
    P = similar(profs[1])
    E = similar(errs[1])
    for i in eachindex(P)
        v = [p[i] for (p, e) in zip(profs, errs) if isfinite(e[i])]
        ev = [e[i] for e in errs if isfinite(e[i])]
        if isempty(v)
            P[i] = 0.0
            E[i] = Inf
        else
            P[i] = median(v)
            E[i] = 1.2533 * sqrt(mean(ev .^ 2)) / sqrt(length(v))
        end
    end
    return P, E
end

# ---------------------------------------------------------------------------------------------
# Linear heights and chi^2 of a set of fiber positions against one band profile
# ---------------------------------------------------------------------------------------------
function _band_weights(e)
    wt = zeros(length(e))
    for y in eachindex(e)
        (isfinite(e[y]) && e[y] > 0) && (wt[y] = 1 / e[y]^2)
    end
    return wt
end

"""
    slit_heights_chi2(P, p, wt, mu, w; fibs, half = 6) -> (chi2, h, use, model, herr)

Least-squares heights of fibers at `mu` (widths `w`, profile of fiber `fibs[f]`) against the
band profile `p` with weights `wt` (adjacent-fiber overlaps solved jointly, tridiagonal), and
the chi^2 over the WHOLE profile (light where the model has no fiber counts against it).
Fibers whose window leaves the chip are not used (height 0). `herr` is the height error
from the diagonal of the normal matrix (neighbour covariance ignored).
"""
function slit_heights_chi2(P::SlitProfile, p, wt, mu, w; fibs = 1:length(mu), half = 6)
    n = length(mu)
    ny = length(p)
    use = [isfinite(mu[f]) && isfinite(w[f]) && w[f] > 0 &&
           round(Int, mu[f]) - half >= 1 && round(Int, mu[f]) + half <= ny for f in 1:n]
    mus = ifelse.(use, mu, 1024.0)
    ws = ifelse.(use, w, 1.0)
    A, y0 = slit_unit_profiles(P, mus, ws, collect(fibs); half = half)
    A[.!use, :] .= 0
    A[.!isfinite.(A)] .= 0
    m = 2half + 1
    dg = zeros(n)
    off = zeros(max(n - 1, 0))
    rhs = zeros(n)
    @inbounds for f in 1:n
        if !use[f]
            dg[f] = 1.0
            continue
        end
        for i in 1:m
            y = y0[f] + i - 1
            dg[f] += wt[y] * A[f, i]^2
            rhs[f] += wt[y] * A[f, i] * p[y]
        end
        if f < n && use[f + 1]
            for i in 1:m
                y = y0[f] + i - 1
                j = y - y0[f + 1] + 1
                (1 <= j <= m) || continue
                off[f] += wt[y] * A[f, i] * A[f + 1, j]
            end
        end
    end
    dg .+= 1e-12
    h = n == 0 ? Float64[] : SymTridiagonal(dg, off) \ rhs
    h[.!use] .= 0
    model = zeros(ny)
    @inbounds for f in 1:n
        use[f] || continue
        for i in 1:m
            model[y0[f] + i - 1] += h[f] * A[f, i]
        end
    end
    chi2 = 0.0
    @inbounds for y in 1:ny
        chi2 += wt[y] * (p[y] - model[y])^2
    end
    herr = [use[f] ? 1 / sqrt(dg[f]) : Inf for f in 1:n]
    return chi2, h, use, model, herr
end

# ---------------------------------------------------------------------------------------------
# Registration of template positions to one band profile
# ---------------------------------------------------------------------------------------------
"""
    slit_register_band(P, p, e, T, W; dmax, dstep = 0.25, dcenter = 0.0, scenter = 0.0,
        srange = 0.004, alias_sep = 4.0)

Register template centers `T` (one band, all fibers) to profile `p`: coarse grid in the shift
`d` over `dcenter ± dmax` at the scales `scenter .+ scoarse`, then three local grids in (d, s). Returns
`(d, s, chi2, margin)`, where `margin` is the lowest chi^2 of the coarse grid at least
`alias_sep` px from the solution, over the best chi^2.
"""
function slit_register_band(P::SlitProfile, p, e, T, W; dmax = TRACE_TEMPLATE_DMAX, dstep = 0.25,
        dcenter = 0.0, scenter = 0.0, srange = 0.004, alias_sep = 4.0, fibs = 1:length(T),
        scoarse = (0.0,))
    wt = _band_weights(e)
    cost(d, s) = first(slit_heights_chi2(P, p, wt, T .+ d .+ s .* (T .- 1024), W; fibs = fibs))
    ds = collect((dcenter - dmax):dstep:(dcenter + dmax))
    # coarse grid over d at each scale of `scoarse` (offsets from scenter): a slit whose scale
    # differs from the template's by a few 1e-3 is off by px at the slit ends, which can make a
    # wrong alias win a grid run at a single scale (LCO 57802, commissioning)
    sgrid = [scenter + o for o in scoarse if abs(o) <= srange + 1e-12]
    chis = fill(Inf, length(ds))
    sbest = fill(scenter, length(ds))
    for sv in sgrid, (i, d) in enumerate(ds)
        c = cost(d, sv)
        if c < chis[i]
            chis[i] = c
            sbest[i] = sv
        end
    end
    k = argmin(chis)
    c0, d0, s0 = chis[k], ds[k], sbest[k]
    for (dd, ss, ns) in ((0.1, 0.001, 4), (0.02, 0.0002, 5), (0.005, 0.00005, 4))
        for s in (s0 - ns * ss):ss:(s0 + ns * ss), d in (d0 - 5dd):dd:(d0 + 5dd)
            abs(s - scenter) <= srange + 1e-12 || continue
            abs(d - dcenter) <= dmax + 1e-9 || continue
            c = cost(d, s)
            if c < c0
                c0, d0, s0 = c, d, s
            end
        end
    end
    far = abs.(ds .- d0) .>= alias_sep
    alias = any(far) ? minimum(chis[far]) : Inf
    return (d = d0, s = s0, chi2 = c0, margin = alias / max(c0, eps()))
end

"extend a template row by `ne` slots at each end with the mean end spacing"
function _extend_row(T, ne)
    s1 = mean(diff(T[1:6]))
    s2 = mean(diff(T[(end - 5):end]))
    return vcat([T[1] - s1 * k for k in ne:-1:1], T, [T[end] + s2 * k for k in 1:ne])
end

_safecorr(a, b) = (length(a) > 10 && std(a) > 0 && std(b) > 0) ? cor(a, b) : NaN

"""
    slit_label_test(P, p, e, T, W, href; d, s, kmax = 3, tp_margin = 0.05)

Which fiber is fiber 1? With the registered positions extended by `kmax` slots at each end,
compare every label shift k in -kmax:kmax by (i) the throughput fingerprint (correlation of
log heights with the template's `href`) and (ii) chi^2 with only the 300 slots of that
labelling (light left outside the slit counts against it). Returns `(kbest, ktp, kchi,
tpmargin, edge_light)`; the fingerprint decides when its margin is at least `tp_margin`.
"""
function slit_label_test(P::SlitProfile, p, e, T, W, href; d, s, kmax = TRACE_TEMPLATE_KMAX,
        tp_margin = TRACE_TEMPLATE_TP_MARGIN)
    nf = length(T)
    kmax == 0 && return (kbest = 0, ktp = 0, kchi = 0, tpmargin = NaN, edge_light = NaN)
    wt = _band_weights(e)
    Te = _extend_row(T, kmax)
    We = vcat(fill(W[1], kmax), W, fill(W[end], kmax))
    fibe = vcat(fill(1, kmax), 1:nf, fill(nf, kmax))
    mue = Te .+ d .+ s .* (Te .- 1024)
    _, he, _, _ = slit_heights_chi2(P, p, wt, mue, We; fibs = fibe)
    hpos = max.(he, 0.0)
    hmed = median(hpos[(kmax + 1):(kmax + nf)])
    ks = (-kmax):kmax
    tp = fill(NaN, length(ks))
    chik = fill(Inf, length(ks))
    edge = fill(NaN, length(ks))
    g = isfinite.(href) .& (href .> 0.02)
    for (i, k) in enumerate(ks)
        sl = (kmax + 1 + k):(kmax + nf + k)
        chik[i] = first(slit_heights_chi2(P, p, wt, mue[sl], We[sl]; fibs = 1:nf))
        out = setdiff(1:length(mue), sl)
        edge[i] = sum(hpos[out]) / max(hmed, eps())
        a = log.(max.(hpos[sl][g], 1e-3 * max(hmed, eps())))
        tp[i] = _safecorr(a, log.(href[g]))
    end
    kchi = ks[argmin(chik)]
    if any(isfinite, tp)
        o = sortperm(replace(tp, NaN => -Inf), rev = true)
        ktp = ks[o[1]]
        tpm = tp[o[1]] - tp[o[2]]
    else
        ktp, tpm = kchi, NaN
    end
    kbest = (isfinite(tpm) && tpm >= tp_margin) ? ktp : kchi
    return (kbest = kbest, ktp = ktp, kchi = kchi, tpmargin = tpm, edge_light = edge[kmax + 1])
end

# ---------------------------------------------------------------------------------------------
# Band refinement with production's Gauss-Hermite fitter
# ---------------------------------------------------------------------------------------------
"""
    slit_gh_refine(P, p, e, mu0, w0, h0; fibs, max_center_move, minw, maxw, n_iter)

Production's `fit_gaussians` on one spatial profile (a band, or one column), all fibers at once
with neighbours subtracted, started at centers `mu0`, widths `w0`, heights `h0`. Centers move at
most `max_center_move` px, widths stay in `[minw, maxw] .* w0`. Returns (n, 3) parameters
(height, center, width), and the (n, 3, 3) covariances when `return_cov`.
"""
function slit_gh_refine(P::SlitProfile, p, e, mu0, w0, h0; fibs = 1:length(mu0),
        max_center_move = 2, minw = 0.5, maxw = 2.0, n_iter = 10, return_cov = false)
    ny = length(p)
    guess = hcat(max.(ifelse.(isfinite.(h0), h0, 0.01), 0.01), clamp.(mu0, 11, ny - 11), w0)
    offset_inds = -4:4
    yc = round.(Int, guess[:, 2])
    fit_inds = yc .+ offset_inds'
    bm_inds = yc .+ (-10:10)'
    return fit_gaussians(p, e, copy(guess), fit_inds, bm_inds, offset_inds,
        _clampfib(P, collect(fibs)), P.x_prof_min, P.x_prof_max_ind, P.n_sub, P.min_prof_fib,
        P.all_y_prof, P.all_y_prof_deriv; n_iter = n_iter, return_cov = return_cov,
        use_first_guess_heights = true, max_center_move = max_center_move,
        min_widths = minw .* guess[:, 3], max_widths = maxw .* guess[:, 3])
end

# ---------------------------------------------------------------------------------------------
# Registration of a template to a frame (all bands)
# ---------------------------------------------------------------------------------------------
"""
    register_frame(P, tmpl, xc, prof, err; mode = :full, dmax, kmax, lit_frac, refine = true)

Register the epoch template to the band profiles of a frame.

1. whole-slit registration at the central band(s), searching `d_epoch ± dmax` (`d_epoch` the
   epoch-median shift; `mode = :full` default dmax 40 px, `:register_only` 3 px);
2. (`:full` only) label test: `kbest != 0` is a refusal (`outcome = "id_template_k"`), never a
   relabel;
3. per-band registration within ±2 px of the central solution, `d(u)` quadratic and `s(u)`
   linear fitted with 4-sigma clipping;
4. (`refine`) every fiber of every band refined by `fit_gaussians` from the registered
   position; QA: `fit_rms` (robust RMS of the refined centers minus the registration about a
   per-fiber cubic in u, lit fibers; bands where the fit hit a bound and was reset are excluded),
   `nlit` (median band height > `lit_frac` of the median fiber and median band height S/N >=
   `snr_min`, PR A's census threshold), band heights.

Returns a NamedTuple: `reg::TemplateRegistration`, `outcome` ("ok" or a named refusal),
central `d`, `s`, `margin`, label test fields, `fit_rms`, `nlit`, `lit` (Bool per fiber),
`H` (nb x 300 band heights), `MU` (nb x 300 refined centers), `MU0` (registered centers).
"""
function register_frame(P::SlitProfile, tmpl, xc, prof, err; mode = :full,
        dmax = (mode == :full ? TRACE_TEMPLATE_DMAX : TRACE_REG_DMAX),
        kmax = (mode == :full ? TRACE_TEMPLATE_KMAX : 0), lit_frac = TRACE_LIT_FRAC,
        refine = true, tp_margin = TRACE_TEMPLATE_TP_MARGIN, snr_min = TRACE_CENSUS_SNR_MIN)
    mode in (:full, :register_only) || error("register_frame: mode must be :full or :register_only")
    nb = length(xc)
    nf = nfibers(tmpl)
    r0 = epoch_median_registration(tmpl)
    Tb = template_centers(tmpl, collect(Float64, xc))   # nb x nf
    Wb = template_widths(tmpl, collect(Float64, xc))
    u = (xc .- 1024) ./ 1024
    # 1) central registration (two central bands averaged)
    mid = clamp(searchsortedfirst(xc, 1024.0), 2, nb)
    cb = (mid - 1):mid
    pc = vec(mean(prof[cb, :], dims = 1))
    ec = vec(sqrt.(mean(err[cb, :] .^ 2, dims = 1)) ./ sqrt(length(cb)))
    Tc = vec(mean(Tb[cb, :], dims = 1))
    Wc = vec(mean(Wb[cb, :], dims = 1))
    uc = mean(u[cb])
    dcen = _polyval(r0.d, uc)
    scen = _polyval(r0.s, uc)
    rc = slit_register_band(P, pc, ec, Tc, Wc; dmax = dmax, dcenter = dcen, scenter = scen,
        scoarse = mode == :full ? (-0.004:0.001:0.004) : (0.0,))
    # 2) labels
    lt = slit_label_test(P, pc, ec, Tc, Wc, tmpl.height_rel; d = rc.d, s = rc.s, kmax = kmax,
        tp_margin = tp_margin)
    base = (d = rc.d, s = rc.s, chi2 = rc.chi2, margin = rc.margin, kbest = lt.kbest,
        ktp = lt.ktp, kchi = lt.kchi, tpmargin = lt.tpmargin, edge_light = lt.edge_light,
        mode = mode)
    empty = (reg = TemplateRegistration([rc.d, 0.0, 0.0], [rc.s, 0.0]), fit_rms = NaN, nlit = 0,
        lit = falses(nf), H = zeros(nb, nf), MU = fill(NaN, nb, nf), MU0 = fill(NaN, nb, nf),
        xc = collect(Float64, xc))
    if lt.kbest != 0
        return merge(base, empty, (outcome = "id_template_k",))
    end
    # 3) per-band registration around the central solution
    db = zeros(nb)
    sb = zeros(nb)
    for b in 1:nb
        r = slit_register_band(P, prof[b, :], err[b, :], Tb[b, :], Wb[b, :]; dmax = 2.0,
            dcenter = rc.d, scenter = rc.s, srange = 0.002)
        db[b] = r.d
        sb[b] = r.s
    end
    keep = trues(nb)
    pd = fit(u, db, 2)
    ps = fit(u, sb, 1)
    for _ in 1:2
        rd = db .- pd.(u)
        σ = 1.4826 * median(abs.(rd[keep])) + 0.02
        keep .= abs.(rd) .< 4σ
        count(keep) >= 4 || (keep .= true; break)
        pd = fit(u[keep], db[keep], 2)
        ps = fit(u[keep], sb[keep], 1)
    end
    cpad(c, n) = (v = zeros(n); cc = collect(Float64, c); v[1:min(n, length(cc))] .= cc[1:min(n, length(cc))]; v)
    reg = TemplateRegistration(cpad(pd.coeffs, 3), cpad(ps.coeffs, 2))
    MU0 = similar(Tb)
    for b in 1:nb
        MU0[b, :] .= Tb[b, :] .+ _polyval(reg.d, u[b]) .+ _polyval(reg.s, u[b]) .* (Tb[b, :] .- 1024)
    end
    H = zeros(nb, nf)
    SN = zeros(nb, nf)
    for b in 1:nb
        _, h, _, _, he = slit_heights_chi2(P, prof[b, :], _band_weights(err[b, :]), MU0[b, :], Wb[b, :])
        H[b, :] .= h
        SN[b, :] .= h ./ he
    end
    # lit: median band height above `lit_frac` of the median fiber AND median band S/N >= snr_min
    hmed = vec(median(H, dims = 1))
    snmed = vec(median(SN, dims = 1))
    hpos = filter(v -> isfinite(v) && v > 0, hmed)
    ref = isempty(hpos) ? NaN : median(hpos)
    lit = (hmed .> lit_frac * ref) .& isfinite.(hmed) .& (snmed .>= snr_min)
    fit_rms = NaN
    MU = fill(NaN, nb, nf)
    if refine
        for b in 1:nb
            np, cv = slit_gh_refine(P, prof[b, :], err[b, :], MU0[b, :], Wb[b, :], H[b, :];
                max_center_move = 1.5, minw = 0.6, maxw = 1.6, n_iter = 6, return_cov = true)
            # `fit_gaussians` resets a fiber that hit a bound to its first guess and marks its
            # covariance NaN: such a band is a failed measurement, not a perfect one
            failed = [!isfinite(cv[f, 2, 2]) for f in 1:nf]
            MU[b, :] .= ifelse.(failed, NaN, np[:, 2])
            H[b, :] .= np[:, 1]
        end
        res = Float64[]
        for f in findall(lit)
            r = MU[:, f] .- MU0[:, f]
            ok = isfinite.(r) .& (abs.(r) .< 1.49) .& (H[:, f] .> lit_frac * ref)
            count(ok) >= 6 || continue
            k = copy(ok)
            pf = fit(u[k], r[k], 3)
            for _ in 1:2
                rr = r .- pf.(u)
                σ = 1.4826 * median(abs.(rr[k])) + 0.005
                k .= ok .& (abs.(rr) .< 4σ)
                count(k) >= 6 || break
                pf = fit(u[k], r[k], 3)
            end
            append!(res, (r .- pf.(u))[ok])
        end
        fit_rms = isempty(res) ? NaN : 1.4826 * median(abs.(res))
    end
    return merge(base, (reg = reg, fit_rms = fit_rms, nlit = count(lit), lit = collect(lit), H = H,
        MU = MU, MU0 = MU0, xc = collect(Float64, xc), outcome = "ok"))
end

"""
    registration_gate(rf; fit_rms_max = TRACE_REG_FIT_RMS_MAX, nlit_min = TRACE_REG_NLIT_MIN)
        -> (pass::Bool, reasons::Vector{String})

Gate of a registration-only frame (levels T1/T2): label shift 0, `fit_rms` below
`fit_rms_max` and at least `nlit_min` fibers with light (the study's alternative-frame gate,
V2 §5.1-5.2).
"""
function registration_gate(rf; fit_rms_max = TRACE_REG_FIT_RMS_MAX, nlit_min = TRACE_REG_NLIT_MIN)
    reasons = String[]
    rf.outcome == "ok" || push!(reasons, "outcome:" * rf.outcome)
    (isfinite(rf.fit_rms) && rf.fit_rms <= fit_rms_max) || push!(reasons, "fit_rms")
    rf.nlit >= nlit_min || push!(reasons, "nlit")
    return isempty(reasons), reasons
end

# ---------------------------------------------------------------------------------------------
# Trace parameters from a template
# ---------------------------------------------------------------------------------------------
"""
    registered_trace_params(tmpl, reg; nx = 2048, heights = tmpl.height_rel)

(nx, n_fibers, 3) trace parameters of the registered template with no per-fiber delta: heights
(default the template's relative throughput; extraction does not use them), registered centers,
template widths. The trace of levels T1, T2 (night registration) and T3 (epoch median).
"""
function registered_trace_params(tmpl, reg::TemplateRegistration; nx = N_XPIX,
        heights = tmpl.height_rel)
    nf = nfibers(tmpl)
    tp = zeros(nx, nf, 3)
    tp[:, :, 1] .= reshape(heights, 1, nf)
    tp[:, :, 2] .= registered_centers(tmpl, reg, 1:nx)
    tp[:, :, 3] .= template_widths(tmpl, collect(Float64, 1:nx))
    return tp
end

"""
    template_fiber_delta(trace_params, tmpl, reg, lit; nband = 32)

Per-fiber departure of a fitted trace from the registered template: maximum over `nband`
x-bands of |band median of (center - registered template center)|, lit fibers only (NaN for
unlit fibers).
"""
function template_fiber_delta(trace_params::AbstractArray{<:Real, 3}, tmpl,
        reg::TemplateRegistration, lit::AbstractVector{Bool}; nband = 32)
    nx, nf, _ = size(trace_params)
    R = registered_centers(tmpl, reg, 1:nx)
    bw = nx ÷ nband
    delta = fill(NaN, nf)
    for f in findall(lit)
        m = 0.0
        for b in 1:nband
            cols = ((b - 1) * bw + 1):(b * bw)
            v = filter(isfinite, view(trace_params, cols, f, 2) .- view(R, cols, f))
            isempty(v) || (m = max(m, abs(median(v))))
        end
        delta[f] = m
    end
    return delta
end

"""
    trace_extract_template(image_data, ivar_image, tmpl, P; good_pixels, lit_frac,
        nband = 32, band_half = 8, dmax, kmax, max_center_move = 2, verbose = false)

Level T0 with a template: the trace of a flat.

1. `register_frame` (`:full`): whole-slit registration around the epoch-median shift, label
   test (k != 0 refuses with `outcome = "id_template_k"`), per-band registration, band heights;
2. production's per-column Gauss-Hermite refinement (`fit_gaussians`, the per-column settings
   of `trace_extract`: 10 iterations, centers within 2 px, widths 0.5-2 x the first guess) of
   the lit fibers, started in EVERY column from the registered template (centers), the
   template widths and the band heights interpolated in x: no outward walk, so a bad column
   cannot propagate into its neighbours;
3. unlit fibers (median band height < `lit_frac` of the median fiber) sit at the registered
   template position with height 0 (`TRACE_FIBER_UNLIT` downstream);
4. `delta` per fiber (`template_fiber_delta`; `TRACE_FIBER_DELTA` above 0.3 px).

`seed = :template_delta` (default) adds to every column's start each lit fiber's smooth
departure from the registered template, measured by a first pass of the same fit on every
`seed_step`-th column (band medians, robust cubic in u); `seed = :template` starts every column
at the registered template itself. Why: `fit_gaussians` resets a fiber whose center or width hits
its bound to the first guess, so with a template start a fiber whose pattern departs from the
template by a large fraction of a pixel (LCO 60702 R: median 0.6 px) gets columns reset to the
template, i.e. jumps; with the smooth departure in the seed those resets are harmless. Neither
seed uses a neighbouring column's fit (no walk).

Returns `(trace_params, trace_param_covs, info)`; on a refusal the arrays are empty and
`info.outcome` names the reason. `info` also holds the registration and label statistics.
"""
function trace_extract_template(image_data, ivar_image, tmpl, P::SlitProfile;
        good_pixels = trues(size(image_data)), lit_frac = TRACE_LIT_FRAC, nband = 32,
        band_half = 8, dmax = TRACE_TEMPLATE_DMAX, kmax = TRACE_TEMPLATE_KMAX,
        max_center_move = 2, tp_margin = TRACE_TEMPLATE_TP_MARGIN, seed = :template_delta,
        seed_step = 8, verbose = false)
    good = good_pixels .& (ivar_image .> 0) .& isfinite.(image_data)
    noise = 1 ./ sqrt.(ivar_image)
    xc, prof, err = slit_bands(image_data, noise, good; nb = nband, half = band_half)
    rf = register_frame(P, tmpl, xc, prof, err; mode = :full, dmax = dmax, kmax = kmax,
        lit_frac = lit_frac, refine = false, tp_margin = tp_margin)
    nf = nfibers(tmpl)
    if rf.outcome != "ok"
        return [], [], merge(rf, (delta = fill(NaN, nf),))
    end
    nx = size(image_data, 1)
    lit = rf.lit
    fl = findall(lit)
    # per-column inputs: flux with bad pixels zeroed and Inf errors (as production)
    flux = copy(image_data)
    errs = copy(noise)
    flux[.!good] .= 0
    errs[.!good] .= Inf
    xcs = rf.xc
    MU0 = registered_centers(tmpl, rf.reg, 1:nx)
    W0 = template_widths(tmpl, collect(Float64, 1:nx))
    function hguess(x)   # band heights interpolated linearly in x (flat beyond the outer bands)
        j = clamp(searchsortedlast(xcs, x), 1, length(xcs) - 1)
        t = clamp((x - xcs[j]) / (xcs[j + 1] - xcs[j]), 0.0, 1.0)
        return (1 - t) .* rf.H[j, :] .+ t .* rf.H[j + 1, :]
    end
    function fit_columns!(params, covs, cols, start)
        for x in cols
            if !isempty(fl)
                np, cv = slit_gh_refine(P, view(flux, x, :), view(errs, x, :), start[x, fl],
                    W0[x, fl], hguess(x)[fl]; fibs = fl, max_center_move = max_center_move,
                    minw = 0.5, maxw = 2.0, n_iter = 10, return_cov = true)
                params[x, fl, :] .= np
                isnothing(covs) || (covs[x, fl, :, :] .= cv)
            end
        end
    end
    # Seed of every column: the registered template, plus (seed = :template_delta, default) each lit
    # fiber's smooth departure from it measured by a first pass on every `seed_step`-th column
    # (band medians, robust cubic in u). Either way no column is started from another column's fit.
    START = copy(MU0)
    seed_delta = zeros(nf)
    if seed == :template_delta && !isempty(fl)
        cols1 = collect((seed_step ÷ 2):seed_step:nx)
        p1 = fill(NaN, nx, nf, 3)
        fit_columns!(p1, nothing, cols1, MU0)
        uu = (cols1 .- 1024) ./ 1024
        ub = Float64[]
        nbs = max(8, length(cols1) ÷ 8)
        edges = round.(Int, range(0, length(cols1), length = nbs + 1))
        for b in 1:nbs
            push!(ub, mean(uu[(edges[b] + 1):edges[b + 1]]))
        end
        for f in fl
            r = p1[cols1, f, 2] .- MU0[cols1, f]
            rb = [_medornan(filter(isfinite, r[(edges[b] + 1):edges[b + 1]])) for b in 1:nbs]
            ok = isfinite.(rb)
            count(ok) >= 6 || continue
            pf = fit(ub[ok], rb[ok], 3)
            k = copy(ok)
            for _ in 1:2
                rr = rb .- pf.(ub)
                σ = 1.4826 * median(abs.(rr[k])) + 0.005
                k .= ok .& (abs.(rr) .< 4σ)
                count(k) >= 6 || break
                pf = fit(ub[k], rb[k], 3)
            end
            seed_delta[f] = maximum(abs.(pf.(ub)))
            for x in 1:nx
                START[x, f] += pf((x - 1024) / 1024)
            end
        end
    elseif seed != :template
        error("trace_extract_template: seed must be :template_delta or :template")
    end
    params = zeros(nx, nf, 3)
    covs = zeros(nx, nf, 3, 3)
    fit_columns!(params, covs, 1:nx, START)
    ul = findall(.!lit)
    params[:, ul, 1] .= 0
    params[:, ul, 2] .= MU0[:, ul]
    params[:, ul, 3] .= W0[:, ul]
    delta = template_fiber_delta(params, tmpl, rf.reg, lit)
    verbose && @info "trace_extract_template: d=$(round(rf.d, digits = 3)) s=$(round(rf.s, digits = 5)) nlit=$(rf.nlit) max delta $(round(maximum(filter(isfinite, delta); init = 0.0), digits = 3))"
    return params, covs, merge(rf, (delta = delta,))
end

"""
    template_qa_fields(info) -> Dict

The template-registration fields of a QA record (prefixed `tmpl_`), from the `info` of
`trace_extract_template` or `register_frame`.
"""
function template_qa_fields(info; tmpl = nothing)
    q = Dict{String, Any}()
    q["tmpl_outcome"] = String(info.outcome)
    q["tmpl_mode"] = String(Symbol(info.mode))
    q["tmpl_d_center"] = Float64(info.d)
    q["tmpl_s_center"] = Float64(info.s)
    q["tmpl_margin"] = Float64(info.margin)
    q["tmpl_kbest"] = Int(info.kbest)
    q["tmpl_ktp"] = Int(info.ktp)
    q["tmpl_kchi"] = Int(info.kchi)
    q["tmpl_tpmargin"] = Float64(info.tpmargin)
    q["tmpl_edge_light"] = Float64(info.edge_light)
    q["tmpl_reg"] = reg_vector(info.reg)
    q["tmpl_nlit"] = Int(info.nlit)
    q["tmpl_fit_rms"] = Float64(info.fit_rms)
    if hasproperty(info, :delta)
        dv = filter(isfinite, info.delta)
        q["tmpl_delta_max"] = isempty(dv) ? NaN : maximum(dv)
        q["tmpl_delta_med"] = isempty(dv) ? NaN : median(dv)
        q["tmpl_n_fiber_delta"] = count(v -> isfinite(v) && v > TRACE_FIBER_DELTA_PX, info.delta)
    end
    if !isnothing(tmpl)
        q["tmpl_id"] = template_id(tmpl)
    end
    return q
end

"""
    template_free_crosscheck(image_data, ivar_image, trace_params, P, refpos; good_pixels,
        lit, mid, n_center_cols, tol_px = 0.01)

Cross-check of a template-seeded flat fit against the template-free path (TRACE_ROBUST_V2 §7.1):
PR A's whole-slit identification of the light in the center band (`trace_peak_census` +
`trace_identify_slit`) and a Gauss-Hermite fit of that band started at the identified positions,
compared with the template-seeded trace averaged over the same columns. Returns a Dict of QA
fields: `tmpl_xcheck_id_reason`, `tmpl_xcheck_label_agree` (fraction of lit, identified fibers
within 0.35 spacings: a label disagreement moves whole fibers by a spacing),
`tmpl_xcheck_med_abs_px` / `tmpl_xcheck_p90_abs_px` (|template-free - template| over those
fibers) and `tmpl_xcheck_flag` = labels disagree, or identification refused, or median > `tol_px`
(0.01 px, the design's tolerance). A QA flag only; it does not gate.
"""
function template_free_crosscheck(image_data, ivar_image, trace_params, P::SlitProfile, refpos;
        good_pixels = trues(size(image_data)), lit = trues(size(trace_params, 2)),
        mid = min(1025, size(image_data, 1) ÷ 2 + 1), n_center_cols = min(100, size(image_data, 1) ÷ 2 - 1),
        tol_px = 0.01)
    q = Dict{String, Any}()
    good = good_pixels .& (ivar_image .> 0) .& isfinite.(image_data)
    cf, ce, _ = trace_center_band(image_data, 1 ./ sqrt.(ivar_image), good; mid = mid, n_center_cols = n_center_cols)
    cf = vec(cf); ce = vec(ce)
    cen = trace_peak_census(cf, ce)
    nf = size(trace_params, 2)
    id = trace_identify_slit(cen.y, refpos; nfib = nf, census_nslots = cen.nslots, census_slots = cen.slots)
    q["tmpl_xcheck_id_reason"] = id.reason
    cols = max(1, mid - n_center_cols):min(size(trace_params, 1), mid + n_center_cols)
    tc = [_medornan(filter(isfinite, view(trace_params, cols, f, 2))) for f in 1:nf]
    tw = [_medornan(filter(isfinite, view(trace_params, cols, f, 3))) for f in 1:nf]
    if !id.ok
        q["tmpl_xcheck_label_agree"] = NaN
        q["tmpl_xcheck_med_abs_px"] = NaN
        q["tmpl_xcheck_p90_abs_px"] = NaN
        q["tmpl_xcheck_flag"] = true
        return q
    end
    start = clamp.(id.positions, 11, length(cf) - 11)
    _, h, _, _, _ = slit_heights_chi2(P, cf, _band_weights(ce), start, ifelse.(isfinite.(tw), tw, 1.0))
    np, cv = slit_gh_refine(P, cf, ce, start, ifelse.(isfinite.(tw), tw, 1.0), h; max_center_move = 2,
        minw = 0.6, maxw = 1.6, n_iter = 10, return_cov = true)
    use = [lit[f] && id.matched[f] && isfinite(cv[f, 2, 2]) && isfinite(tc[f]) for f in 1:nf]
    d = np[use, 2] .- tc[use]
    agree = isempty(d) ? NaN : mean(abs.(d) .< 0.35 * id.spacing)
    ok = abs.(d) .< 0.35 * id.spacing
    da = abs.(d[ok])
    q["tmpl_xcheck_label_agree"] = agree
    q["tmpl_xcheck_med_abs_px"] = isempty(da) ? NaN : median(da)
    q["tmpl_xcheck_p90_abs_px"] = isempty(da) ? NaN : quantile(da, 0.9)
    q["tmpl_xcheck_flag"] = !(agree == 1.0) || !(q["tmpl_xcheck_med_abs_px"] <= tol_px)
    return q
end
