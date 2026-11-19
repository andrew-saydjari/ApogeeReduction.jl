using ApogeeReduction: sky_fiber_rj, dither_rj, fit_fiber_wavecal_gn, fit_dither_gn, lm_solve,
                       WAVE_SOLVER_LBFGS_FALLBACK, WAVE_SOLVER_HOLD_R, WAVE_SOLVER_HOLD_B, WAVE_CHIP_HOLD_RULE,
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
    # the other hold rules: narrow = per-chip rule only (7 lines, 2 per chip -> both free);
    # none = never held
    @test WAVE_CHIP_HOLD_RULE[] == :wide
    try
        WAVE_CHIP_HOLD_RULE[] = :narrow
        _, _, _, f7n = fit_fiber_wavecal_gn(x[k7], y[k7], chip[k7], copy(cpp0); wporder = 4, cporder = 0, hold_params = hold)
        @test f7n & (WAVE_SOLVER_HOLD_R | WAVE_SOLVER_HOLD_B) == 0
        _, _, _, fhn = fit_fiber_wavecal_gn(x[keep], y[keep], chip[keep], copy(cpp0); wporder = 4, cporder = 0, hold_params = hold)
        @test fhn & WAVE_SOLVER_HOLD_R != 0 && fhn & WAVE_SOLVER_HOLD_B == 0
        WAVE_CHIP_HOLD_RULE[] = :none
        nln, _, _, fhx = fit_fiber_wavecal_gn(x[keep], y[keep], chip[keep], copy(cpp0); wporder = 4, cporder = 0, hold_params = hold)
        @test fhx & (WAVE_SOLVER_HOLD_R | WAVE_SOLVER_HOLD_B) == 0 && nln[1] != hold[1]
        WAVE_CHIP_HOLD_RULE[] = :bogus
        @test_throws ArgumentError fit_fiber_wavecal_gn(x, y, chip, copy(cpp0); wporder = 4, cporder = 0)
    finally
        WAVE_CHIP_HOLD_RULE[] = :wide
    end

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

using ApogeeReduction: comb_fiber_rj, fit_comb_chip_gn, fit_cavity_gn, cavity_schur,
                       comb_fiber_cavity_info, comb_exp_nonlinear_loss_fit!, WAVECAL_COMB_SOLVER,
                       comb_exp_get_and_save_fpi_wavecal, FPI_COMB_FIT_KWARGS, N_FIBERS
using HDF5
isdefined(@__MODULE__, :SyntheticWavecal) || include("synthetic_wavecal.jl")
using .SyntheticWavecal: SyntheticWavecal as SW

@testset "FPI comb: Gauss-Newton chip-transform and cavity fits" begin
    rng = MersenneTwister(5)
    lin = [16000.0, -580.0, 3.0, 0.5, 0.1]
    p4 = [-1.0712, 1.0011, 1.0761, 0.9878]
    dither = [0.0 1.0 0.0; 2e-4 (1 + 3e-5) 1e-6]
    n = 240
    x = rand(rng, n) .- 0.5
    chip = rand(rng, 1:3, n)
    expi = rand(rng, 1:2, n)
    model(c, p) = [evalpoly(evalpoly((chip[i] == 1 ? p[1] + p[2] * x[i] : chip[i] == 3 ? p[3] + p[4] * x[i] : x[i]),
        dither[expi[i], :]), c) for i in 1:n]
    y = model(lin, p4) .+ 1e-4 .* randn(rng, n)
    # analytic Jacobian vs central finite differences
    th = vcat(lin, p4) .+ 1e-4 .* randn(rng, 9)
    r0, J = comb_fiber_rj(th, x, y, chip, expi, dither, 4)
    Jfd = similar(J)
    for k in 1:9
        h = 1e-6 * max(1.0, abs(th[k]))
        tp = copy(th); tp[k] += h; tm = copy(th); tm[k] -= h
        Jfd[:, k] = -(comb_fiber_rj(tp, x, y, chip, expi, dither, 4)[1] .- comb_fiber_rj(tm, x, y, chip, expi, dither, 4)[1]) ./ (2h)
    end
    @test maximum(abs.(J .- Jfd) ./ (abs.(Jfd) .+ 1.0)) < 1e-5
    # recovery from a perturbed start, and the same minimum as the LBFGS variable-projection fit
    p0 = p4 .+ [1e-3, 1e-4, -1e-3, -1e-4]
    pg, cg, ok = fit_comb_chip_gn(x, y, chip, expi, dither, p0; wporder = 4)
    @test ok
    @test maximum(abs.(pg .- p4)) < 1e-5
    cpp = [0.0 1.0; 0.0 1.0; 0.0 1.0]
    loss(p) = comb_exp_nonlinear_loss_fit!(copy(cpp), dither, p, x, y, chip, expi;
        wporder = 4, cporder = 1, returnL2only = true, dporder = 2, chip_only = true)
    pl = ApogeeReduction.Optim.minimizer(ApogeeReduction.wavecal_lbfgs(loss, p0))
    @test loss(pg) <= loss(pl) * (1 + 1e-9)

    # cavity: lambda = 2 L / (m + eps) on 3 "fibers", each with its own free quartic in x
    L0, e0 = SW.CAV_APO, SW.MOFF_APO
    Ltrue, etrue = L0 * (1 + 3e-8), e0 + 0.02            # truth differs from the seed
    S = zeros(2, 2); KtK = zeros(2, 2)
    pw_all = Float64[]; m_all = Float64[]
    for f in 1:3
        m = collect(4400.0:4900.0)
        lam_true = 2Ltrue ./ (m .+ etrue)
        xf = (lam_true .- 16000) ./ -580 .+ 1e-4 * f       # the fiber's pixel coordinate
        # the per-fiber solution is fitted to the SEED expectation (free quartic), as in the comb loop
        yseed = 2L0 ./ (m .+ e0)
        A = hcat([xf .^ k for k in 0:4]...)
        c = A \ yseed
        pw = yseed .- (yseed .- A * c)                     # measured = expected - residual
        Jf = A
        K = hcat(2 ./ (m .+ e0), -2L0 ./ (m .+ e0) .^ 2)
        Sf, Kf = cavity_schur(Jf, K)
        S .+= Sf; KtK .+= Kf
        append!(pw_all, pw); append!(m_all, m)
    end
    good = trues(length(pw_all))
    Lg, eg, nheld = fit_cavity_gn(pw_all, m_all, good, L0, e0, S, KtK)
    @test nheld == 2                       # both directions degenerate with free per-fiber solutions
    @test Lg == L0 && eg == e0             # the cavity stays at its seed exactly
    # identifiable configuration (per-fiber solutions held: S = K'K): the fit recovers the truth
    pw_t = 2Ltrue ./ (m_all .+ etrue)
    Lt, et, nh = fit_cavity_gn(pw_t, m_all, good, L0, e0, KtK, KtK)
    @test nh == 0
    @test abs(Lt - Ltrue) / Ltrue < 1e-10 && abs(et - etrue) < 1e-6
end

@testset "FPI comb fit: Gauss-Newton vs LBFGS on synthetic peaks" begin
    lin, nlp = SW.night_params()
    out = Dict{Symbol, Any}()
    for mode in (:lbfgs, :gn)
        mktempdir() do dir
            f1 = SW.write_fpi_peaks(dir, "0001", fill(10_000, N_FIBERS))
            f2 = SW.write_fpi_peaks(dir, "0002", fill(10_000, N_FIBERS); d = [1.0e-4, 1.0])
            try
                WAVECAL_COMB_SOLVER[] = mode
                o = comb_exp_get_and_save_fpi_wavecal([f1, f2], lin, nlp; FPI_COMB_FIT_KWARGS...,
                    verbose = false, checkpoint_mode = "clobber")
                out[mode] = (o[4], h5open(h -> (read(h["fpi_cavity_size"]), read(h["fpi_m0_offset"])), o[1]))
            finally
                WAVECAL_COMB_SOLVER[] = :gn
            end
        end
    end
    wl, wg = out[:lbfgs][1], out[:gn][1]
    dv = 299792458.0 .* abs.(wg .- wl) ./ wl
    @test maximum(dv) < 1.0                                # m/s
    @test out[:gn][2] == (SW.CAV_APO, SW.MOFF_APO)         # cavity held at the seed
end
