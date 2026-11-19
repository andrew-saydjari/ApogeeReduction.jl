using LinearAlgebra: Diagonal, diag

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

"""Switch for the chip-transform hold in `fit_fiber_wavecal_gn` (on in production; off only to
measure the solver change separately from the hold)."""
const WAVE_CHIP_TRANSFORM_HOLD = Ref(true)

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
`min_lines_chip` lines AND the fiber has more lines than free parameters; otherwise it is held at
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
    freeR = nR >= min_lines_chip || !WAVE_CHIP_TRANSFORM_HOLD[]
    freeB = nB >= min_lines_chip || !WAVE_CHIP_TRANSFORM_HOLD[]
    # more lines than free parameters, else hold both chip transforms
    if WAVE_CHIP_TRANSFORM_HOLD[] && length(x) <= nl + na * (freeR + freeB)
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
