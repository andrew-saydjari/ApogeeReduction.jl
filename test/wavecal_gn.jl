using ApogeeReduction: sky_fiber_rj, dither_rj, fit_fiber_wavecal_gn, fit_dither_gn, lm_solve,
                       WAVE_SOLVER_LBFGS_FALLBACK, WAVE_SOLVER_HOLD_R, WAVE_SOLVER_HOLD_B,
                       nonlinear_loss_fit!, dither_nonlinear_loss_fit, positional_poly_mat, transform_x_chips
using Random

@testset "Gauss-Newton wavelength fits" begin
    rng = MersenneTwister(3)
    lin = [16000.0, -580.0, 3.0, 0.5, 0.1]
    lam(x) = (positional_poly_mat([x], porder = 4) * lin)[1]
    # lines on three chips; x is green-chip-relative, R/B offsets ~ -/+1.07
    xR = collect(range(-0.45, 0.45, length = 6)); xG = collect(range(-0.4, 0.4, length = 6)); xB = collect(range(-0.45, 0.4, length = 6))
    x = vcat(xR, xG, xB); chip = vcat(fill(1, 6), fill(2, 6), fill(3, 6))
    aR, aB = -1.0712, 1.0753

    # analytic Jacobians against central finite differences
    for cporder in (0, 1)
        na = cporder + 1
        th = vcat(lin, cporder == 0 ? [aR] : [aR, 1.0002], cporder == 0 ? [aB] : [aB, 0.9993])
        y = zeros(length(x))
        r0, J = sky_fiber_rj(th, x, y, chip, 4, cporder)
        for j in eachindex(th)
            h = 1e-6 * max(1.0, abs(th[j]))
            tp = copy(th); tp[j] += h; tm = copy(th); tm[j] -= h
            fd = -(sky_fiber_rj(tp, x, y, chip, 4, cporder)[1] .- sky_fiber_rj(tm, x, y, chip, 4, cporder)[1]) ./ (2h)
            @test maximum(abs.(fd .- J[:, j])) <= 1e-5 * max(1.0, maximum(abs.(J[:, j])))
        end
    end
    p = [2.0e-4, 1.0 + 3.0e-5]
    xt = collect(range(-0.45, 0.45, length = 12))
    _, Jd = dither_rj(p, xt, zeros(12), lin)
    for j in 1:2
        h = 1e-8; pp = copy(p); pp[j] += h; pm = copy(p); pm[j] -= h
        fd = -(dither_rj(pp, xt, zeros(12), lin)[1] .- dither_rj(pm, xt, zeros(12), lin)[1]) ./ (2h)
        @test maximum(abs.(fd .- Jd[:, j])) <= 1e-4 * maximum(abs.(Jd[:, j]))
    end

    # noiseless per-fiber sky fit (cporder 0): recovers the injected chip offsets from a start 1 px off
    xtrue = [c == 1 ? aR + xi : c == 3 ? aB + xi : xi for (xi, c) in zip(x, chip)]
    y = lam.(xtrue)
    cpp0 = zeros(3, 1); cpp0[1, 1] = aR + 5e-4; cpp0[3, 1] = aB - 5e-4
    nlp, lp, res, flag = fit_fiber_wavecal_gn(x, y, chip, copy(cpp0); wporder = 4, cporder = 0)
    @test flag == 0x00
    @test abs(nlp[1] - aR) < 1e-9 && abs(nlp[2] - aB) < 1e-9
    @test maximum(abs.(res)) < 1e-6
    # the LBFGS fallback path (forced with maxit = 0) gives the same answer, flagged
    nlf, lpf, resf, flagf = fit_fiber_wavecal_gn(x, y, chip, copy(cpp0); wporder = 4, cporder = 0, maxit = 0)
    @test flagf & WAVE_SOLVER_LBFGS_FALLBACK != 0
    @test abs(nlf[1] - aR) < 1e-6 && abs(nlf[2] - aB) < 1e-6

    # chip-transform hold: one line on R -> R held at the hold value, flagged; B still fitted
    keep = vcat([1], 7:18)
    hold = [-1.0700, 1.0760]
    nlh, _, _, flagh = fit_fiber_wavecal_gn(x[keep], y[keep], chip[keep], copy(cpp0); wporder = 4, cporder = 0,
        hold_params = hold)
    @test flagh & WAVE_SOLVER_HOLD_R != 0 && flagh & WAVE_SOLVER_HOLD_B == 0
    @test nlh[1] == hold[1]
    @test nlh[2] != hold[2] && abs(nlh[2] - aB) < 1e-4   # fitted (the held R line pulls it slightly)
    # too few lines overall (7 lines, 5 polynomial + 2 offset parameters) -> both held
    k7 = [1, 2, 7, 8, 13, 14, 15]
    nl7, _, _, flag7 = fit_fiber_wavecal_gn(x[k7], y[k7], chip[k7], copy(cpp0); wporder = 4, cporder = 0, hold_params = hold)
    @test flag7 & WAVE_SOLVER_HOLD_R != 0 && flag7 & WAVE_SOLVER_HOLD_B != 0
    @test nl7 == hold

    # dither fit: GN recovers offset and scale; forced fallback agrees and is flagged
    yd = lam.(p[1] .+ p[2] .* xt)
    pg, fg = fit_dither_gn(xt, yd, lin, [0.0, 1.0]; wporder = 4)
    @test fg == 0x00 && abs(pg[1] - p[1]) < 1e-10 && abs(pg[2] - p[2]) < 1e-9
    pl, fl = fit_dither_gn(xt, yd, lin, [0.0, 1.0]; wporder = 4, maxit = 0)
    @test fl == WAVE_SOLVER_LBFGS_FALLBACK && abs(pl[1] - p[1]) < 1e-7
    # start independence (noisy data): perturbed starts converge to the same point
    yn = yd .+ 0.01 .* randn(rng, length(yd))
    pa, _ = fit_dither_gn(xt, yn, lin, [0.0, 1.0]; wporder = 4)
    for _ in 1:5
        pb, _ = fit_dither_gn(xt, yn, lin, [1e-3 * (2rand(rng) - 1), 1 + 1e-3 * (2rand(rng) - 1)]; wporder = 4)
        @test maximum(abs.(pa .- pb)) < 1e-10
    end
end
