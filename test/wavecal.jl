@testset "wavecal" begin
    @testset "get_fpi_exp_m0s empty in_exp (A5)" begin
        # 3 fibers, 2 exposures. Fiber 1 has peaks in both exposures, fiber 2
        # only in exposure 1, fiber 3 has NO ingested FPI peaks at all.
        # (Pre-fix, exp_m0s was an Int array and the peakless entries were
        # assigned NaN -> InexactError, killing the night's FPI wavecal.)
        n_fnames = 2
        # peak-major layout: rows = peaks, cols = fibers (0 = unused row)
        fpi_line_expInt = [1 1 0
                           1 0 0
                           2 0 0
                           2 0 0]
        peak_ints = Float64[1204 1210 0
                            1203 0 0
                            1205 0 0
                            1202 0 0]

        exp_m0s = ApogeeReduction.get_fpi_exp_m0s(fpi_line_expInt, peak_ints, n_fnames)

        @test size(exp_m0s) == (2, 3)
        @test eltype(exp_m0s) == Float64
        @test exp_m0s[1, 1] == 1203.0
        @test exp_m0s[2, 1] == 1202.0
        @test exp_m0s[1, 2] == 1210.0
        @test isnan(exp_m0s[2, 2])   # fiber with zero peaks in exposure 2
        @test all(isnan.(exp_m0s[:, 3]))  # fiber with zero peaks everywhere

        # downstream median-m0 computation must skip the NaNs and not throw
        med_m0 = round(ApogeeReduction.nanmedian(
            ApogeeReduction.nanmedian(exp_m0s, 2), 1)[1, 1])
        @test med_m0 == round(median([1203.0 + 1210.0, 2 * 1202.0]) / 2)

        # degenerate case: no fiber has any peaks -> all NaN, still no throw
        exp_m0s_empty = ApogeeReduction.get_fpi_exp_m0s(zeros(Int, 4, 3), peak_ints, n_fnames)
        @test all(isnan.(exp_m0s_empty))
        med_empty = ApogeeReduction.nanmedian(
            ApogeeReduction.nanmedian(exp_m0s_empty, 2), 1)[1, 1]
        @test isnan(med_empty)

        # NaN peak_ints (bad peaks) propagate to NaN m0 instead of throwing
        peak_ints_nan = copy(peak_ints)
        peak_ints_nan[1, 2] = NaN
        exp_m0s_nan = ApogeeReduction.get_fpi_exp_m0s(fpi_line_expInt, peak_ints_nan, n_fnames)
        @test isnan(exp_m0s_nan[1, 2])
        @test exp_m0s_nan[1, 1] == 1203.0
    end

    @testset "interpolate_wave_params: no slit-edge bias, outliers rejected" begin
        # Per-fiber wavelength parameters whose true fiber dependence is CURVED
        # (quadratic in fiber), plus noise, a few gross outlier fibers and one
        # NaN fiber. The pre-fix regulariser fitted linParams[:, 3:end] and all
        # nlParams LINEARLY in fiber, which leaves a systematic bias at the slit
        # edges (fibers 1-30 and 271-300) of ~ curvature * (u^2 - <u^2>).
        rng = MersenneTwister(20261117)
        nfib = ApogeeReduction.N_FIBERS
        fibInds = collect(1:nfib)
        u = (fibInds .- (nfib + 1) / 2) ./ ((nfib - 1) / 2)   # -1 .. 1
        edge = vcat(1:30, (nfib - 29):nfib)
        centre = 101:200
        # (offset, slope, curvature, noise) per column; curvature = 25 noise
        lin_spec = [(16000.0, 0.5, 0.25, 0.01), (130.0, 0.02, 0.05, 0.002),
            (-3.0, 0.002, 0.005, 2e-4), (0.5, 1e-3, 2.5e-3, 1e-4), (-0.2, 5e-4, 1e-3, 4e-5)]
        nl_spec = [(-1.070, 2e-4, 5e-4, 2e-5), (1.0755, -1e-4, 5e-4, 2e-5)]
        truth(s) = s[1] .+ s[2] .* u .+ s[3] .* u .^ 2
        lin_true = reduce(hcat, truth.(lin_spec))
        nl_true = reduce(hcat, truth.(nl_spec))
        lin_noise = [s[4] for s in lin_spec]'
        nl_noise = [s[4] for s in nl_spec]'
        linParams = lin_true .+ lin_noise .* randn(rng, nfib, length(lin_spec))
        nlParams = nl_true .+ nl_noise .* randn(rng, nfib, length(nl_spec))
        outliers = [3, 17, 150, 222, 296]
        for (k, f) in enumerate(outliers)
            sgn = isodd(k) ? 1 : -1
            linParams[f, :] .+= sgn .* 200 .* vec(lin_noise)
            nlParams[f, :] .-= sgn .* 200 .* vec(nl_noise)
        end
        linParams[60, :] .= NaN
        nlParams[60, :] .= NaN

        interp_nl, interp_lin = ApogeeReduction.interpolate_wave_params(
            fibInds, nlParams, linParams)

        @test size(interp_lin) == size(linParams)
        @test size(interp_nl) == size(nlParams)
        @test all(isfinite, interp_lin)
        @test all(isfinite, interp_nl)

        for (P, T, N) in ((interp_lin, lin_true, vec(lin_noise)),
            (interp_nl, nl_true, vec(nl_noise)))
            for j in axes(P, 2)
                r = (P[:, j] .- T[:, j]) ./ N[j]   # in units of the per-fiber noise
                # no slit-edge bias relative to the centre (the linear fit gave
                # ~ 0.55 * curvature = ~14 noise units here)
                @test abs(mean(r[edge]) - mean(r[centre])) < 1.0
                # outlier and NaN fibers are replaced by the smooth solution
                @test maximum(abs.(r[vcat(outliers, 60)])) < 3.0
                # the regularised solution is smooth: no fiber far off the truth
                @test maximum(abs.(r)) < 5.0
            end
        end
    end
end
