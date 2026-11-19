using LinearAlgebra: Diagonal, diag, qr, eigen, eigvals, Symmetric

# Gauss-Newton / Levenberg-Marquardt solves with analytic Jacobians for the polynomial wavelength fits.
#
# The per-fiber sky wavelength fit, the per-fiber sky/FPI dither fit and the per-exposure dither of the
# FPI comb fit are least-squares problems whose model is a polynomial in a transformed pixel
# coordinate: lambda = sum_k c_k xt^k, with xt = a_0 + a_1 x + ... (chip transform) or a dither
# polynomial of xt. Their Jacobians are exact and cheap, so Gauss-Newton with Levenberg-Marquardt
# damping converges to the minimum to numerical precision in a handful of iterations. The result
# does not depend on the Optim/LineSearches version or on finite-difference step choices.
# If Gauss-Newton does not converge within WAVECAL_GN_MAXIT iterations or ends with a non-finite or
# higher objective than its start, the fit falls back to LBFGS (`wavecal_lbfgs`) and the fallback
# is recorded in a per-fit flag.

"""Iteration cap of the Gauss-Newton / Levenberg-Marquardt solves; at the cap the fit falls back to LBFGS."""
const WAVECAL_GN_MAXIT = 100

"""Solver flags (UInt8 bits) stored per fiber in the wavelength products."""
const WAVE_SOLVER_LBFGS_FALLBACK = 0x01
"""Chip R transform held at its calibration / default value (too few lines on R, or too few lines overall)."""
const WAVE_SOLVER_HOLD_R = 0x02
"""Chip B transform held at its calibration / default value."""
const WAVE_SOLVER_HOLD_B = 0x04

"""Fewest lines on a chip for fitting that chip's transform in the per-fiber sky fit."""
const WAVE_MIN_LINES_CHIP_TRANSFORM = 2

"""
Rule of the chip-transform hold in `fit_fiber_wavecal_gn` (the review decision on it is a change of
this one default):
- `:wide`   (as first proposed): hold a chip whose line count is below `WAVE_MIN_LINES_CHIP_TRANSFORM`,
  and hold both chips when the fiber has no more lines than free parameters;
- `:narrow`: only the per-chip rule;
- `:none`:   no hold (every chip transform free; measurement only).
"""
const WAVE_CHIP_HOLD_RULE = Ref(:wide)

"""
    lm_solve(rj, theta0, active; maxit = WAVECAL_GN_MAXIT) -> (theta, ssr, ssr0, iterations, converged)

Levenberg-Marquardt minimisation of `sum(abs2, r)` where `r, J = rj(theta)` returns the residual
`y - f(theta)` and the model Jacobian `df/dtheta`. Only the parameters flagged in `active` move.
`converged` is false only if the iteration cap is hit.
"""
function lm_solve(rj, theta0::AbstractVector{Float64}, active::AbstractVector{Bool};
        maxit = WAVECAL_GN_MAXIT)
    th = copy(theta0)
    r, J = rj(th)
    ssr0 = sum(abs2, r)
    ssr = ssr0
    lam = 0.0
    it = 0
    any(active) || return th, ssr, ssr0, 0, true
    step = zeros(length(th))
    for outer it in 1:maxit
        Ja = J[:, active]
        A = Ja' * Ja
        g = Ja' * r
        D = Diagonal(max.(diag(A), eps()))
        accepted = false
        rn, Jn, ssrn = r, J, ssr
        for _ in 1:40
            sa = try
                (A + lam * D) \ g
            catch
                nothing
            end
            if !isnothing(sa) && all(isfinite, sa)
                fill!(step, 0.0)
                step[active] .= sa
                rn, Jn = rj(th .+ step)
                ssrn = sum(abs2, rn)
                if isfinite(ssrn) && ssrn <= ssr
                    accepted = true
                    break
                end
            end
            lam = lam == 0 ? 1e-8 : 10 * lam
        end
        # no downhill step even with heavy damping: at the minimum to numerical precision
        accepted || return th, ssr, ssr0, it, true
        th .+= step
        small = maximum(abs.(step[active]) ./ (abs.(th[active]) .+ 1e-10)) < 1e-13
        dssr = ssr - ssrn
        r, J, ssr = rn, Jn, ssrn
        lam /= 10
        (small || dssr <= 1e-15 * ssr) && return th, ssr, ssr0, it, true
    end
    return th, ssr, ssr0, it, false
end

# chip / dither transform with the semantics of `transform_x_chips` (porder 0 = offset + x)
_tr(x, a) = length(a) == 1 ? a[1] + x : evalpoly(x, a)
_dtr_da(x, a, j) = length(a) == 1 ? 1.0 : x^(j - 1)        # d xt / d a_j (j = 1-based)
_dpoly(c, x) = sum((k - 1) * c[k] * x^(k - 2) for k in 2:length(c); init = 0.0)

"""
    sky_fiber_rj(theta, x, y, chip, wporder, cporder) -> (r, J)

Residual and model Jacobian of the per-fiber sky model: theta = [c_0..c_wporder, a_R (cporder+1),
a_B (cporder+1)]; xt = transform of x with a_R (chip 1), identity (chip 2), a_B (chip 3);
lambda = sum_k c_k xt^k.
"""
function sky_fiber_rj(theta, x, y, chip, wporder, cporder)
    nl = wporder + 1
    na = cporder + 1
    c = theta[1:nl]
    aR = theta[(nl + 1):(nl + na)]
    aB = theta[(nl + na + 1):(nl + 2na)]
    aG = cporder == 0 ? [0.0] : vcat(0.0, 1.0, zeros(cporder - 1))
    n = length(x)
    xt = similar(x, Float64)
    for i in 1:n
        xt[i] = _tr(x[i], chip[i] == 1 ? aR : chip[i] == 3 ? aB : aG)
    end
    r = y .- [evalpoly(t, c) for t in xt]
    J = zeros(n, nl + 2na)
    for i in 1:n
        for k in 1:nl
            J[i, k] = xt[i]^(k - 1)
        end
        if chip[i] == 1 || chip[i] == 3
            d = _dpoly(c, xt[i])
            a = chip[i] == 1 ? aR : aB
            off = chip[i] == 1 ? nl : nl + na
            for j in 1:na
                J[i, off + j] = d * _dtr_da(x[i], a, j)
            end
        end
    end
    return r, J
end

"""
    fit_fiber_wavecal_gn(x, y, chip, chipPolyParams; wporder, cporder, hold_params = nothing,
        min_lines_chip = WAVE_MIN_LINES_CHIP_TRANSFORM, maxit = WAVECAL_GN_MAXIT)
        -> (nlParamsOpt, linParamsOpt, linResid, flag)

Per-fiber chip-transform + polynomial wavelength fit by Gauss-Newton. `chipPolyParams` (N_CHIPS x
cporder+1) is the start. A chip's transform is fitted only if that chip has at least
`min_lines_chip` lines AND (rule `:wide`, see `WAVE_CHIP_HOLD_RULE`) the fiber has more lines than
free parameters; otherwise it is held at
`hold_params` (same layout as the nonlinear parameter vector: [a_R; a_B]; default = the start).
Falls back to LBFGS (over the free chip parameters) if Gauss-Newton does not converge.
"""
function fit_fiber_wavecal_gn(x, y, chip, chipPolyParams; wporder = 2, cporder = 1,
        hold_params = nothing, min_lines_chip = WAVE_MIN_LINES_CHIP_TRANSFORM, maxit = WAVECAL_GN_MAXIT)
    na = cporder + 1
    nl = wporder + 1
    p0 = ChipPolyParams2Params(chipPolyParams)
    hold = isnothing(hold_params) ? copy(p0) : collect(Float64, hold_params)
    flag = 0x00
    nR = count(==(1), chip)
    nB = count(==(3), chip)
    rule = WAVE_CHIP_HOLD_RULE[]
    rule in (:wide, :narrow, :none) || throw(ArgumentError("WAVE_CHIP_HOLD_RULE = $(rule)"))
    freeR = nR >= min_lines_chip || rule == :none
    freeB = nB >= min_lines_chip || rule == :none
    # wide rule: more lines than free parameters, else hold both chip transforms
    if rule == :wide && length(x) <= nl + na * (freeR + freeB)
        freeR = false
        freeB = false
    end
    freeR || (flag |= WAVE_SOLVER_HOLD_R)
    freeB || (flag |= WAVE_SOLVER_HOLD_B)
    pstart = copy(p0)
    freeR || (pstart[1:na] .= hold[1:na])
    freeB || (pstart[(na + 1):(2na)] .= hold[(na + 1):(2na)])
    cpp = copy(chipPolyParams)
    if length(x) == 0
        linResid, linParamsOpt = nonlinear_loss_fit!(cpp, pstart, x, y, chip;
            wporder = wporder, cporder = cporder, returnL2only = false)
        return pstart, linParamsOpt, linResid, flag
    end
    # polynomial at the start transform (linear), then the joint solve
    _, c0 = nonlinear_loss_fit!(cpp, pstart, x, y, chip; wporder = wporder, cporder = cporder,
        returnL2only = false)
    theta0 = vcat(c0, pstart)
    active = vcat(trues(nl), fill(freeR, na), fill(freeB, na))
    th, ssr, ssr0, _, conv = lm_solve(t -> sky_fiber_rj(t, x, y, chip, wporder, cporder), theta0, active;
        maxit = maxit)
    pOpt = th[(nl + 1):end]
    if !conv || !isfinite(ssr) || ssr > ssr0
        flag |= WAVE_SOLVER_LBFGS_FALLBACK
        free = vcat(fill(freeR, na), fill(freeB, na))
        function loss(pf)
            p = copy(pstart)
            p[free] .= pf
            nonlinear_loss_fit!(cpp, p, x, y, chip; wporder = wporder, cporder = cporder, returnL2only = true)
        end
        pOpt = copy(pstart)
        if any(free)
            pOpt[free] .= Optim.minimizer(wavecal_lbfgs(loss, pstart[free]))
        end
    end
    linResid, linParamsOpt = nonlinear_loss_fit!(cpp, pOpt, x, y, chip;
        wporder = wporder, cporder = cporder, returnL2only = false)
    return pOpt, linParamsOpt, linResid, flag
end

"""
    dither_rj(p, xt, y, linparam) -> (r, J)

Residual and Jacobian of the dither model lambda = poly(transform(xt; p); linparam).
"""
function dither_rj(p, xt, y, linparam)
    n = length(xt)
    xp = [_tr(t, p) for t in xt]
    r = y .- [evalpoly(t, linparam) for t in xp]
    J = zeros(n, length(p))
    for i in 1:n
        d = _dpoly(linparam, xp[i])
        for j in eachindex(p)
            J[i, j] = d * _dtr_da(xt[i], p, j)
        end
    end
    return r, J
end

"""
    fit_dither_gn(xt, y, linparam, p0; wporder, maxit = WAVECAL_GN_MAXIT) -> (params, flag)

Dither polynomial fit (`dither_nonlinear_loss_fit` model) by Gauss-Newton from `p0`, with the LBFGS
fallback (flag `WAVE_SOLVER_LBFGS_FALLBACK`).
"""
function fit_dither_gn(xt, y, linparam, p0; wporder = 2, maxit = WAVECAL_GN_MAXIT)
    length(xt) == 0 && return collect(Float64, p0), 0x00
    lin = collect(Float64, linparam[1:(wporder + 1)])
    th, ssr, ssr0, _, conv = lm_solve(p -> dither_rj(p, xt, y, lin), collect(Float64, p0), trues(length(p0));
        maxit = maxit)
    if !conv || !isfinite(ssr) || ssr > ssr0
        loss(p) = dither_nonlinear_loss_fit(p, xt, y, linparam; wporder = wporder, returnL2only = true)
        return Optim.minimizer(wavecal_lbfgs(loss, collect(Float64, p0))), WAVE_SOLVER_LBFGS_FALLBACK
    end
    return th, 0x00
end

# ------------------------------------------------------------------------------------------
# FPI comb fit: chip-transform fit and cavity (L, eps) fit
# ------------------------------------------------------------------------------------------

"""Solver of the two remaining comb fits (chip transform per fiber, cavity per night): `:gn`
(Gauss-Newton / Levenberg-Marquardt with analytic Jacobians) or `:lbfgs` (previous fits);
`:shadow` runs both on identical inputs, uses LBFGS and records both in `COMB_SOLVER_LOG`
(validation only)."""
const WAVECAL_COMB_SOLVER = Ref(:gn)
"""Validation records of `:shadow` runs (one NamedTuple per fit)."""
const COMB_SOLVER_LOG = NamedTuple[]

"""
    comb_fiber_rj(theta, x, y, chip, expi, dither, wporder) -> (r, J)

Residual and model Jacobian of one fiber in the comb chip-transform fit (`comb_exp_nonlinear_loss_fit!`
with `chip_only`, cporder 1): theta = [c_0..c_wporder, a_R, b_R, a_B, b_B]; xt = a_c + b_c x for
chips R (1) and B (3), x for G; xt' = dither[e, :] polynomial of xt for exposure e (held);
lambda = sum_k c_k xt'^k.
"""
function comb_fiber_rj(theta, x, y, chip, expi, dither, wporder)
    nl = wporder + 1
    c = theta[1:nl]
    n = length(x)
    r = similar(x, Float64)
    J = zeros(n, nl + 4)
    for i in 1:n
        ch = chip[i]
        a, b = ch == 1 ? (theta[nl + 1], theta[nl + 2]) : ch == 3 ? (theta[nl + 3], theta[nl + 4]) : (0.0, 1.0)
        xt = a + b * x[i]
        d = view(dither, expi[i], :)
        xp = evalpoly(xt, d)
        dd = _dpoly(d, xt)
        r[i] = y[i] - evalpoly(xp, c)
        for k in 1:nl
            J[i, k] = xp^(k - 1)
        end
        if ch == 1 || ch == 3
            g = _dpoly(c, xp) * dd
            off = ch == 1 ? nl : nl + 2
            J[i, off + 1] = g
            J[i, off + 2] = g * x[i]
        end
    end
    return r, J
end

"""
    fit_comb_chip_gn(x, y, chip, expi, dither, p0; wporder, maxit = WAVECAL_GN_MAXIT)
        -> (pOpt, linParamsOpt, converged)

Gauss-Newton fit of one fiber's chip transforms (cporder 1: a_R, b_R, a_B, b_B) jointly with its
polynomial, per-exposure dither held: the same minimum as the LBFGS variable-projection fit in
`comb_exp_get_and_save_fpi_wavecal`. A chip without lines keeps its start values (no data, no step).
"""
function fit_comb_chip_gn(x, y, chip, expi, dither, p0; wporder = 4, maxit = WAVECAL_GN_MAXIT)
    nl = wporder + 1
    xt = similar(x, Float64)
    for i in eachindex(x)
        ch = chip[i]
        a, b = ch == 1 ? (p0[1], p0[2]) : ch == 3 ? (p0[3], p0[4]) : (0.0, 1.0)
        xt[i] = evalpoly(a + b * x[i], view(dither, expi[i], :))
    end
    c0 = positional_poly_mat(xt, porder = wporder) \ y
    theta0 = vcat(c0, collect(Float64, p0))
    active = vcat(trues(nl), fill(any(==(1), chip), 2), fill(any(==(3), chip), 2))
    th, ssr, ssr0, _, conv = lm_solve(t -> comb_fiber_rj(t, x, y, chip, expi, dither, wporder),
        theta0, active; maxit = maxit)
    ok = conv && isfinite(ssr) && ssr <= ssr0
    return th[(nl + 1):end], th[1:nl], ok
end

"""
    cavity_schur(Jf, K) -> (S, KtK)

Contribution of one fiber to the cavity information with that fiber's own solution free:
S = K~' K~ with K~ the cavity Jacobian columns `K` (rows = the fiber's peaks) projected off the
column space of the fiber's model Jacobian `Jf` (polynomial, chip transforms, dithers), and the
unprojected K' K for scale.
"""
function cavity_schur(Jf, K)
    KtK = K' * K
    size(Jf, 1) == 0 && return zeros(2, 2), KtK
    Q = Matrix(qr(Jf).Q)[:, 1:min(size(Jf)...)]
    Kp = K .- Q * (Q' * K)
    return Kp' * Kp, KtK
end

"""Relative eigenvalue below which a direction of the cavity information (after projecting off
the free per-fiber solutions) counts as not identifiable and is not moved."""
const CAVITY_SCHUR_RTOL = 1e-6

"""
    fit_cavity_gn(peak_waves, peak_ints, good, cavity, m_offset, S, KtK; maxit = 20)
        -> (cavity, m_offset, nheld)

Gauss-Newton fit of the cavity size L and the fractional order offset eps to the measured FPI peak
wavelengths (lambda = 2 L / (m + eps)), moving only along the directions that the comb can
identify with the per-fiber solutions free. `S` and `KtK` are the summed `cavity_schur` matrices
(parameter order [L, eps]); directions of `S` with eigenvalue below `CAVITY_SCHUR_RTOL` times the
matching scale of `KtK` are held at their start. With every per-fiber solution free both
directions are degenerate (L is a common velocity zero point, eps zero point + stretch), so the
cavity stays at its seed exactly; `nheld` is the number of held directions.
"""
function fit_cavity_gn(peak_waves, peak_ints, good, cavity, m_offset, S, KtK; maxit = 20)
    pw = peak_waves[good]; m = peak_ints[good]
    th = [Float64(cavity), Float64(m_offset)]
    # identifiable subspace from the Schur information, scaled by the unprojected one
    sc = sqrt.(max.(diag(KtK), eps()))
    Sn = S ./ (sc * sc')
    E = eigen(Symmetric(Sn))
    keep = E.values .> CAVITY_SCHUR_RTOL
    nheld = count(.!keep)
    nheld == 2 && return th[1], th[2], 2
    P = E.vectors[:, keep]                 # identifiable directions (scaled parameters)
    model(t) = 2 .* t[1] ./ (m .+ t[2])
    ssr = sum(abs2, pw .- model(th))
    for _ in 1:maxit
        r = pw .- model(th)
        Jc = hcat(2 ./ (m .+ th[2]), -2 .* th[1] ./ (m .+ th[2]) .^ 2) ./ sc'
        Jp = Jc * P
        step = (P * (Jp \ r)) ./ sc        # scaled parameters u = theta .* sc
        all(isfinite, step) || break
        tn = th .+ step
        ssrn = sum(abs2, pw .- model(tn))
        ssrn <= ssr || break
        th, dssr, ssr = tn, ssr - ssrn, ssrn
        dssr <= 1e-15 * ssr && break
    end
    return th[1], th[2], nheld
end

"""
    comb_fiber_cavity_info(x, chip, expi, nl4, lin, dither, m, cavity, m_offset, wporder) -> (S, KtK)

`cavity_schur` for one fiber of the comb fit at its fitted state: the fiber's free parameters are
its polynomial, chip transforms (cporder 1) and the dithers of exposures 2..n (dporder from
`dither`); the cavity columns are d lambda / d L and d lambda / d eps of lambda = 2 L / (m + eps).
"""
function comb_fiber_cavity_info(x, chip, expi, nl4, lin, dither, m, cavity, m_offset, wporder)
    theta = vcat(lin, nl4)
    _, Jc = comb_fiber_rj(theta, x, zeros(length(x)), chip, expi, dither, wporder)
    keepc = vec(any(!=(0.0), Jc, dims = 1))
    nd = size(dither, 2)
    ne = size(dither, 1)
    Jd = zeros(length(x), nd * max(ne - 1, 0))
    for i in eachindex(x)
        e = expi[i]
        e == 1 && continue
        ch = chip[i]
        a, b = ch == 1 ? (nl4[1], nl4[2]) : ch == 3 ? (nl4[3], nl4[4]) : (0.0, 1.0)
        xt = a + b * x[i]
        g = _dpoly(lin, evalpoly(xt, view(dither, e, :)))
        for j in 1:nd
            Jd[i, (e - 2) * nd + j] = g * xt^(j - 1)
        end
    end
    keepd = vec(any(!=(0.0), Jd, dims = 1))
    Jf = hcat(Jc[:, keepc], Jd[:, keepd])
    K = hcat(2 ./ (m .+ m_offset), -2 .* cavity ./ (m .+ m_offset) .^ 2)
    return cavity_schur(Jf, K)
end

# wavelength (exposure-1 frame) of a comb fiber solution at a few pixels per chip, for comparisons
function _comb_probe_waves(nl4, lin)
    xs = ((1.0, 512.0, 1024.0, 1536.0, 2048.0) .- (N_XPIX ÷ 2)) ./ N_XPIX
    w = Float64[]
    for ch in 1:3, x in xs
        a, b = ch == 1 ? (nl4[1], nl4[2]) : ch == 3 ? (nl4[3], nl4[4]) : (0.0, 1.0)
        push!(w, evalpoly(a + b * x, lin))
    end
    return w
end

"""
    comb_chip_fit(loss, p0, x, y, chip, expi, dither, wporder; fiber, iter) -> pOpt

The comb chip-transform fit of one fiber with the solver of `WAVECAL_COMB_SOLVER[]`; Gauss-Newton
falls back to LBFGS if it does not converge.
"""
function comb_chip_fit(loss, p0, x, y, chip, expi, dither, wporder; fiber = 0, iter = 0)
    mode = WAVECAL_COMB_SOLVER[]
    mode == :lbfgs && return Optim.minimizer(wavecal_lbfgs(loss, p0))
    t_gn = @elapsed (pg, cg, ok) = fit_comb_chip_gn(x, y, chip, expi, dither, p0; wporder = wporder)
    if mode == :gn
        return ok ? pg : Optim.minimizer(wavecal_lbfgs(loss, p0))
    end
    # :shadow: both, plus Gauss-Newton from two perturbed starts; LBFGS result used
    t_lb = @elapsed pl = Optim.minimizer(wavecal_lbfgs(loss, p0))
    linL = positional_poly_mat(_comb_xt(x, chip, expi, dither, pl), porder = wporder) \ y
    C = 299792458.0
    wl = _comb_probe_waves(pl, linL); wg = _comb_probe_waves(pg, cg)
    dvs = Float64[]
    for s in (1, -1)
        pp = collect(Float64, p0) .+ s .* [2e-3, 1e-4, -2e-3, -1e-4]
        pp2, cp2, _ = fit_comb_chip_gn(x, y, chip, expi, dither, pp; wporder = wporder)
        push!(dvs, maximum(abs.(C .* (_comb_probe_waves(pp2, cp2) .- wg) ./ wg)))
    end
    push!(COMB_SOLVER_LOG, (kind = :chip, fiber = fiber, iter = iter, n = length(x),
        ssr_lbfgs = loss(pl), ssr_gn = loss(pg), ssr_start = loss(p0), gn_ok = ok,
        dv_max = maximum(abs.(C .* (wg .- wl) ./ wl)), dv_chip = [maximum(abs.(C .* (wg[k] .- wl[k]) ./ wl[k])) for k in (1:5, 6:10, 11:15)],
        dv_perturb = maximum(dvs), t_gn = t_gn, t_lbfgs = t_lb))
    return pl
end

function _comb_xt(x, chip, expi, dither, p)
    xt = similar(x, Float64)
    for i in eachindex(x)
        ch = chip[i]
        a, b = ch == 1 ? (p[1], p[2]) : ch == 3 ? (p[3], p[4]) : (0.0, 1.0)
        xt[i] = evalpoly(a + b * x[i], view(dither, expi[i], :))
    end
    return xt
end

"""
    comb_cavity_fit(loss, e0, mults, peak_waves, peak_ints, good, cavity, m_offset, S, KtK; iter)
        -> (cavity, m_offset)

The comb's cavity / order-offset fit with the solver of `WAVECAL_COMB_SOLVER[]` (`fit_cavity_gn`
for `:gn`).
"""
function comb_cavity_fit(loss, e0, mults, peak_waves, peak_ints, good, cavity, m_offset, S, KtK; iter = 0)
    mode = WAVECAL_COMB_SOLVER[]
    lb() = (o = Optim.minimizer(wavecal_lbfgs(loss, e0)) .* mults; (cavity + o[1], o[2]))
    mode == :lbfgs && return lb()
    t_gn = @elapsed (Lg, eg, nheld) = fit_cavity_gn(peak_waves, peak_ints, good, cavity, m_offset, S, KtK)
    mode == :gn && return Lg, eg
    t_lb = @elapsed (Ll, el) = lb()
    sc = sqrt.(max.(diag(KtK), eps()))
    ev = eigvals(Symmetric(S ./ (sc * sc')))
    # undamped, unguarded 2-parameter least squares from the same start (what a plain GN step would do)
    pw = peak_waves[good]; m = peak_ints[good]
    Ju = hcat(2 ./ (m .+ m_offset), -2 .* cavity ./ (m .+ m_offset) .^ 2)
    du = (Ju' * Ju) \ (Ju' * (pw .- 2 .* cavity ./ (m .+ m_offset)))
    mmed = median(m); C = 299792458.0
    push!(COMB_SOLVER_LOG, (kind = :cavity, iter = iter, n = length(pw), cavity0 = cavity, eps0 = m_offset,
        L_lbfgs = Ll, eps_lbfgs = el, L_gn = Lg, eps_gn = eg, nheld = nheld, schur_eig = ev,
        dv_L_lbfgs = C * (Ll - cavity) / cavity, dv_eps_lbfgs = -C * (el - m_offset) / mmed,
        dv_L_plain = C * du[1] / cavity, dv_eps_plain = -C * du[2] / mmed,
        ssr_lbfgs = loss([(Ll - cavity), el] ./ mults), ssr_seed = loss(e0), t_gn = t_gn, t_lbfgs = t_lb))
    return Ll, el
end
