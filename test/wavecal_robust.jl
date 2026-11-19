using ApogeeReduction: get_sky_wavecal, chip_offset_hold_values, fit_poly_robust,
                       fit_poly_percentile_clip, interpolate_wave_params, positional_poly_mat,
                       WAVE_CHIP_HOLD_R, WAVE_CHIP_HOLD_B, WAVE_MIN_LINES_CHIP_TRANSFORM,
                       WAVE_SMOOTH_MIN_LINES, N_FIBERS, N_CHIPS
using Random, Statistics

@testset "per-exposure sky fit: chip offset held below 2 lines (fix 1)" begin
    @test WAVE_MIN_LINES_CHIP_TRANSFORM == 2
    @test isnothing(chip_offset_hold_values(nothing, "skyLinePeaks_apo_60000_0001_R_object.h5", 0))
    lin = [16000.0, -580.0, 3.0, 0.5, 0.1]
    aR, aB = -1.0712, 1.0761                     # true offsets
    nslot = 12
    ux = fill(NaN, nslot, N_FIBERS); fw = fill(NaN, nslot, N_FIBERS); ci = zeros(Int, nslot, N_FIBERS)
    # fiber 1: 1 R line, 5 G, 5 B; fiber 2: 2 R lines, 5 G, 5 B; fiber 3: 5 R, 5 G, 1 B; fiber 4: 0 lines
    function put!(k, nR, nG, nB)
        j = 0
        for (c, n, a) in ((1, nR, aR), (2, nG, 0.0), (3, nB, aB)), t in 1:n
            j += 1
            x = -0.4 + 0.8 * (t - 0.5) / max(n, 1)
            ux[j, k] = x; ci[j, k] = c
            fw[j, k] = evalpoly(a + x, lin)
        end
    end
    put!(1, 1, 5, 5); put!(2, 2, 5, 5); put!(3, 5, 5, 1)
    # a corrupted single R line on fiber 1 (wrong by 1 Angstrom): a free offset would chase it
    fw[1, 1] += 1.0
    cpp0 = zeros(N_FIBERS, N_CHIPS, 1)
    cpp0[:, 1, 1] .= -1.070; cpp0[:, 3, 1] .= 1.0755        # start values
    hold = repeat([-1.0711 1.0760], N_FIBERS)              # "calibration" values
    linP, nlP, res, fl = get_sky_wavecal(ux, fw, ci, cpp0; cporder = 0, wporder = 4, hold_params = hold)
    @test fl[1] == WAVE_CHIP_HOLD_R && nlP[1, 1] == hold[1, 1]       # 1 R line: held at the calibration
    @test fl[2] == 0x00 && abs(nlP[2, 1] - aR) < 1e-6 && abs(nlP[2, 2] - aB) < 1e-6   # fitted
    @test fl[3] == WAVE_CHIP_HOLD_B && nlP[3, 2] == hold[3, 2] && abs(nlP[3, 1] - aR) < 1e-4   # the held B line pulls R slightly
    @test fl[4] == (WAVE_CHIP_HOLD_R | WAVE_CHIP_HOLD_B) && nlP[4, :] == hold[4, :]
    # without calibration values a held offset stays at the start
    _, nlP0, _, fl0 = get_sky_wavecal(ux, fw, ci, cpp0; cporder = 0, wporder = 4)
    @test fl0 == fl && nlP0[1, 1] == -1.070
end

@testset "bounded-influence smoothing fit (fix 3)" begin
    rng = MersenneTwister(7)
    x = collect(1.0:300.0)
    truth = 1e-3 .+ 2e-6 .* (x .- 150) .- 3e-9 .* (x .- 150) .^ 2
    sig = 1e-6
    y = truth .+ sig .* randn(rng, 300)
    f, k = fit_poly_robust(x, y, 2)
    @test all(k) || count(k) >= 295
    @test maximum(abs.(f.(x) .- truth)) < 5sig
    # one fiber moving out from the curve: the curve changes continuously and boundedly
    yref = f.(x)
    prev = yref
    maxjump = 0.0; maxdev = 0.0
    for d in vcat(0:0.25:20, 50, 1e2, 1e3, 1e6)
        yy = copy(y); yy[300] += d * sig
        g, _ = fit_poly_robust(x, yy, 2)
        cur = g.(x)
        maxjump = max(maxjump, maximum(abs.(cur .- prev)))
        maxdev = max(maxdev, maximum(abs.(cur .- yref)))
        prev = cur
        if d >= 50
            @test maximum(abs.(cur .- yref)) < 1e-3 * sig     # far outlier: no influence
        end
    end
    @test maxdev < 1.0 * sig      # never pulls the curve by more than the per-fiber noise
    @test maxjump < 0.2 * sig     # no switch between two solutions
    # non-finite values and keep = false carry no weight
    yy = copy(y); yy[10] = NaN; yy[20] = 1.0
    keep = trues(300); keep[20] = false
    g, kk = fit_poly_robust(x, yy, 2; keep = keep)
    @test !kk[10] && !kk[20]
    @test maximum(abs.(g.(x) .- truth)) < 5sig
    # percentile clip path still available and honours keep throughout its iterations
    h, kh = fit_poly_percentile_clip(x, yy, 2; keep = keep)
    @test !kh[20]
end

@testset "few-line fibers out of the slit smoothing (fix 2)" begin
    @test WAVE_SMOOTH_MIN_LINES == 10
    rng = MersenneTwister(11)
    fib = collect(1:N_FIBERS)
    nl = hcat(-1.07 .+ 1e-6 .* fib .+ 1e-6 .* randn(rng, N_FIBERS),
        1.0755 .+ 1e-6 .* randn(rng, N_FIBERS))
    ln = hcat(16000.0 .+ 0.01 .* fib, fill(-580.0, N_FIBERS), fill(3.0, N_FIBERS),
        fill(0.5, N_FIBERS), fill(0.1, N_FIBERS)) .+ 1e-6 .* randn(rng, N_FIBERS, 5)
    nlines = fill(30, N_FIBERS)
    # run-away chip offset on an edge fiber with few lines
    nlb = copy(nl); nlb[1, 1] += 3e-4; nlines[1] = 4
    ex = nlines .< WAVE_SMOOTH_MIN_LINES
    a_nl, a_lin = interpolate_wave_params(fib, nlb, ln; n_lines = nlines, smooth_exclude = ex)
    nlnan = copy(nl); nlnan[1, :] .= NaN
    lnnan = copy(ln); lnnan[1, :] .= NaN
    b_nl, b_lin = interpolate_wave_params(fib, nlnan, lnnan; n_lines = nlines)
    @test a_nl ≈ b_nl      # the excluded fiber has no influence at all
    @test a_lin[2:end, :] ≈ b_lin[2:end, :]
    @test isfinite(a_nl[1, 1]) && abs(a_nl[1, 1] - nl[1, 1]) < 1e-5   # and gets the smooth value
end
