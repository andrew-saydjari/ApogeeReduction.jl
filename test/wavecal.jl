using HDF5

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

    @testset "interpolate_wave_params: per-MTP-block steps, joint across parameters" begin
        # Curved fiber dependence as above, plus a constant step per 30-fiber MTP
        # block that is shared (with fixed per-parameter signs, i.e. correlated)
        # by every non-constant parameter, plus noise, gross outliers, a NaN fiber
        # and one block with only 4 finite fibers. A smooth-only regulariser leaves
        # each block off by ~ its step; the block steps must recover them, and the
        # nearly empty block must fall back to the smooth solution.
        rng = MersenneTwister(20261118)
        nfib = ApogeeReduction.N_FIBERS
        fibInds = collect(1:nfib)
        u = (fibInds .- (nfib + 1) / 2) ./ ((nfib - 1) / 2)
        blk = (fibInds .- 1) .÷ 30 .+ 1
        lin_spec = [(16000.0, 0.5, 0.25, 0.01), (130.0, 0.02, 0.05, 0.002),
            (-3.0, 0.002, 0.005, 2e-4), (0.5, 1e-3, 2.5e-3, 1e-4), (-0.2, 5e-4, 1e-3, 4e-5)]
        nl_spec = [(-1.070, 2e-4, 5e-4, 2e-5), (1.0755, -1e-4, 5e-4, 2e-5)]
        truth(s) = s[1] .+ s[2] .* u .+ s[3] .* u .^ 2
        lin_noise = [s[4] for s in lin_spec]
        nl_noise = [s[4] for s in nl_spec]
        # block step amplitudes in units of the per-fiber noise; one pattern for all
        # parameters, with per-parameter signs (correlated, as in the real fits)
        step_units = [6.0, -5.0, 4.0, -6.0, 5.0, -4.0, 6.0, -5.0, 3.0, -6.0]
        lin_sign = [0.0, 1.0, -1.0, 1.0, -1.0]   # the per-MTP constant term has no step here
        nl_sign = [1.0, -1.0]
        lin_true = reduce(hcat, truth.(lin_spec)) .+
                   step_units[blk] .* (lin_sign .* lin_noise)'
        nl_true = reduce(hcat, truth.(nl_spec)) .+ step_units[blk] .* (nl_sign .* nl_noise)'
        linParams = lin_true .+ lin_noise' .* randn(rng, nfib, length(lin_spec))
        nlParams = nl_true .+ nl_noise' .* randn(rng, nfib, length(nl_spec))
        outliers = [3, 17, 150, 222, 296]
        for (k, f) in enumerate(outliers)
            sgn = isodd(k) ? 1 : -1
            linParams[f, :] .+= sgn .* 200 .* lin_noise
            nlParams[f, :] .-= sgn .* 200 .* nl_noise
        end
        linParams[60, :] .= NaN
        nlParams[60, :] .= NaN
        sparse_block = 7                     # fibers 181-210: only 4 finite fibers
        dead = setdiff(findall(blk .== sparse_block), [183, 191, 199, 207])
        linParams[dead, :] .= NaN
        nlParams[dead, :] .= NaN

        interp_nl, interp_lin = ApogeeReduction.interpolate_wave_params(
            fibInds, nlParams, linParams)
        smooth_nl, smooth_lin = ApogeeReduction.interpolate_wave_params(
            fibInds, nlParams, linParams; block_steps = false)
        @test all(isfinite, interp_lin) && all(isfinite, interp_nl)

        # block-mean error in noise units, per non-constant parameter and block
        function block_errors(P_nl, P_lin)
            E = Float64[]
            for (P, T, N, cols) in ((P_nl, nl_true, nl_noise, 1:2), (P_lin, lin_true, lin_noise, 2:5))
                for j in cols, b in 1:10
                    b == sparse_block && continue
                    push!(E, mean((P[blk .== b, j] .- T[blk .== b, j]) ./ N[j]))
                end
            end
            E
        end
        E_new = block_errors(interp_nl, interp_lin)
        E_smooth = block_errors(smooth_nl, smooth_lin)
        # the steps are recovered: every block mean within 1 noise unit of the truth
        @test maximum(abs.(E_new)) < 1.0
        # and the smooth-only (pre-block-step) regulariser fails that same check
        @test maximum(abs.(E_smooth)) > 3.0
        # the nearly empty block keeps the smooth solution (no step from 4 fibers)
        in_sparse = blk .== sparse_block
        @test interp_nl[in_sparse, :] == smooth_nl[in_sparse, :]
        @test interp_lin[in_sparse, 2:end] == smooth_lin[in_sparse, 2:end]
        # fibers with too few sky lines do not contribute to the steps: a block whose
        # fibers all have 6 lines keeps the smooth solution, the others are unchanged
        n_lines = fill(18, nfib)
        n_lines[blk .== 2] .= 6
        lines_nl, lines_lin = ApogeeReduction.interpolate_wave_params(
            fibInds, nlParams, linParams; n_lines = n_lines)
        @test lines_nl[blk .== 2, :] == smooth_nl[blk .== 2, :]
        @test lines_lin[blk .== 2, 2:end] == smooth_lin[blk .== 2, 2:end]
        @test lines_nl[blk .!= 2, :] == interp_nl[blk .!= 2, :]
        # outlier and NaN fibers get their block's solution, not their own values
        for (P, T, N, cols) in ((interp_nl, nl_true, nl_noise, 1:2),
            (interp_lin, lin_true, lin_noise, 2:5))
            for j in cols
                r = (P[:, j] .- T[:, j]) ./ N[j]
                @test maximum(abs.(r[vcat(outliers, 60)])) < 3.0
                @test maximum(abs.(r[.!in_sparse])) < 5.0
            end
        end
        # the block step itself: joint fiber selection, nothing below min_fibers
        R = randn(rng, 30, 3) .+ [1.0 -2.0 3.0]
        R[5, 2] = 50.0                       # outlier in ONE parameter drops the whole fiber
        R[9, :] .= NaN
        st = ApogeeReduction.mtp_block_steps(R)
        kept = setdiff(1:30, [5, 9])
        @test st ≈ vec(mean(R[kept, :], dims = 1)) atol = 1e-12
        @test isnothing(ApogeeReduction.mtp_block_steps(R[1:9, :]; min_fibers = 10))
    end
end

# ---------------------------------------------------------------------------
# FPI night-average grid (comb_exp_get_and_save_fpi_wavecal).
#
# The FPI night grid -- saved to waveCalFPI, returned, and stored by
# fpi_medwavecal_skyline_dither as `fpi/nightAve_wave_soln` -- must be the
# fitted night solution at the reference (first-exposure) dither. It used to be
# the LAST FPI exposure's dithered grid, because the per-exposure loop that
# writes the waveCalNightfpiDither files overwrote the same buffer.

# Independent evaluation of the wavecal model (not using the package helpers):
# x = (pix - 1024)/2048; chip G untransformed, R/B registered by nlParams
# [R0, R1, B0, B1] (cporder 1) or [R0, B0] (cporder 0); optional dither poly
# on the registered x; lambda = sum lin[k] xt^(k-1).
function _fpi_test_grid(lin, nl; dith = nothing)
    nf = size(lin, 1)
    cp = size(nl, 2) ÷ 2 - 1
    x = ((1:2048) .- 1024) ./ 2048
    W = zeros(2048, nf, 3)
    for f in 1:nf, c in 1:3
        xt = if c == 2
            collect(x)
        else
            p = c == 1 ? nl[f, 1:(cp + 1)] : nl[f, (cp + 2):end]
            cp == 0 ? p[1] .+ x : p[1] .+ p[2] .* x
        end
        if !isnothing(dith)
            d = dith[f, :]
            xt = sum(d[k] .* xt .^ (k - 1) for k in eachindex(d))
        end
        W[:, f, c] .= evalpoly.(xt, Ref(lin[f, :]))
    end
    return W
end

# Write fpiPeaks_{R,G,B} files for one exposure whose FPI lines (lambda_m =
# 2 cav / (m + m_off)) fall where `lin`/`nl` put them, shifted by `shift_pix`.
function _write_fpipeaks_model(fname_R, lin, nl; shift_pix = 0.0,
        cav = 3.7363125e7, m_off = 0.225)
    W = _fpi_test_grid(lin, nl)
    for (c, chip) in enumerate(("R", "G", "B"))
        cols = Vector{Vector{Float64}}(undef, size(lin, 1))
        for f in axes(lin, 1)
            lam = W[:, f, c]               # decreasing with pixel
            mlo = ceil(Int, 2cav / lam[1] - m_off)
            mhi = floor(Int, 2cav / lam[end] - m_off)
            xs = Float64[]
            for m in mlo:mhi
                lm = 2cav / (m + m_off)
                i = findlast(>=(lm), lam)
                (isnothing(i) || i >= 2048) && continue
                push!(xs, i + (lam[i] - lm) / (lam[i] - lam[i + 1]) + shift_pix)
            end
            cols[f] = filter(x -> 3 <= x <= 2046, xs)
        end
        n_peaks = maximum(length, cols)
        mat = fill(NaN, n_peaks, 4, size(lin, 1))
        cov = fill(NaN, n_peaks, 4, 4, size(lin, 1))
        for f in axes(lin, 1), (p, xp) in enumerate(cols[f])
            mat[p, :, f] .= (800.0, xp, 1.25, 9.0)
            cov[p, 2, 2, f] = 1e-4
        end
        h5open(replace(fname_R, "_R_" => "_$(chip)_"), "w") do fh
            write(fh, "fpi_line_mat", mat)
            write(fh, "fpi_line_cov_mat", cov)
        end
    end
    return fname_R
end

@testset "FPI night-average grid is the reference-dither solution" begin
    dir = mktempdir()
    nf = ApogeeReduction.N_FIBERS
    # seed (sky-like, cporder 0) and truth (FPI, cporder 1), typical APO values
    lin = repeat([16150.35 -579.42 -43.24 1.05 -0.015], nf)
    nl_sky = repeat([-1.0702 1.0754], nf)
    nl_true = repeat([-1.0702 1.0 1.0754 1.0], nf)
    # exposure 1 = reference; exposure 2 is dithered by 0.5 pix
    f1 = _write_fpipeaks_model(joinpath(dir, "fpiPeaks_apo_60000_0009_R_arclamp.h5"),
        lin, nl_true; shift_pix = 0.0)
    f2 = _write_fpipeaks_model(joinpath(dir, "fpiPeaks_apo_60000_0010_R_arclamp.h5"),
        lin, nl_true; shift_pix = 0.5)

    out = ApogeeReduction.comb_exp_get_and_save_fpi_wavecal(
        [f1, f2], lin, nl_sky; cporder = 1, wporder = 4, dporder = 2,
        n_sigma = 4, max_ang_sigma = 0.2, max_iter = 1, verbose = false,
        checkpoint_mode = "clobber")
    @test !isnothing(out)
    outname, night_lin, night_nl, night_grid = out

    saved_grid = h5open(f -> read(f["chipWaveSoln"]), outname)
    dith1 = joinpath(dir, "waveCalNightfpiDither_apo_60000_0009_arclamp.h5")
    dith2 = joinpath(dir, "waveCalNightfpiDither_apo_60000_0010_arclamp.h5")
    exp1_grid, exp2_grid, exp2_dith = h5open(f -> read(f["chipWaveSoln"]), dith1),
    h5open(f -> read(f["chipWaveSoln"]), dith2), h5open(f -> read(f["ditherParams"]), dith2)

    ref_grid = _fpi_test_grid(night_lin, night_nl)
    # what is returned (and so stored as nightAve_wave_soln) is what is saved ...
    @test saved_grid == night_grid
    # ... and it is the fitted night solution at the reference dither
    @test maximum(abs, night_grid .- ref_grid) < 1e-8
    @test maximum(abs, night_grid .- exp1_grid) < 1e-8
    # the dithered exposure's own grid is still its parameters + its dither ...
    @test maximum(abs, exp2_grid .- _fpi_test_grid(night_lin, night_nl; dith = exp2_dith)) < 1e-8
    # ... it really is dithered (~0.5 pix x ~0.28 A/pix) ...
    @test 0.05 < ApogeeReduction.nanmedian(abs.(exp2_grid .- ref_grid)) < 0.3
    # ... and the night grid is NOT the last exposure's dithered grid (pre-fix bug)
    @test ApogeeReduction.nanmedian(abs.(night_grid .- exp2_grid)) > 0.05
end
