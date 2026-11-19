using ApogeeReduction: dither_fit_mode, interpolate_bad_fibers!, relthrpt_unusable_for_peaks,
                       wavecal_fiber_quality, get_and_save_sky_dither_per_fiber,
                       comb_exp_get_and_save_fpi_wavecal, FPI_COMB_FIT_KWARGS,
                       DITHER_MIN_LINES, DITHER_FLAG_INTERP_NOLINES, DITHER_FLAG_INTERP_FEWLINES,
                       DITHER_FLAG_INTERP_LOWTHRPT, WAVECAL_BIT_NO_EXPOSURE_DITHER,
                       WAVECAL_BIT_FLAGS_MISSING, RELTHRPT_BROKEN_BIT, RELTHRPT_WARN_BIT, N_FIBERS
using HDF5
using ApogeeReduction: ApogeeReduction
isdefined(@__MODULE__, :SyntheticWavecal) || include("synthetic_wavecal.jl")
using .SyntheticWavecal: SyntheticWavecal as SW

@testset "per-fiber dither: minimum line count" begin
    @test DITHER_MIN_LINES == 6
    @test dither_fit_mode(0, false) == (:interp, DITHER_FLAG_INTERP_NOLINES)
    @test dither_fit_mode(1, false) == (:interp, DITHER_FLAG_INTERP_FEWLINES)
    @test dither_fit_mode(5, false) == (:interp, DITHER_FLAG_INTERP_FEWLINES)
    @test dither_fit_mode(6, false) == (:fit, 0x00)
    @test dither_fit_mode(20, false) == (:fit, 0x00)
    @test dither_fit_mode(20, true) == (:interp, DITHER_FLAG_INTERP_LOWTHRPT)   # unusable fiber, even with lines
    @test dither_fit_mode(0, true) == (:interp, DITHER_FLAG_INTERP_LOWTHRPT)

    # neighbour interpolation (straight line through the two nearest good fibers)
    p = Float64[1 10; 2 20; 3 30; 4 40; 5 50]
    good = [true, false, true, false, true]
    q = interpolate_bad_fibers!(p .* [1, 0, 1, 0, 1], good)
    @test q ≈ p
    q1 = interpolate_bad_fibers!(Float64[7 1; 0 0], [true, false])
    @test q1[2, :] == [7.0, 1.0]

    # relative-throughput mask read from the ar1Dcal file beside the skyLinePeaks file
    mktempdir() do d
        bm = zeros(Int, N_FIBERS)
        bm[7] = RELTHRPT_BROKEN_BIT | RELTHRPT_WARN_BIT
        bm[9] = RELTHRPT_WARN_BIT                      # warning only: still usable
        h5open(joinpath(d, "ar1Dcal_apo_60000_0001_R_object.h5"), "w") do f
            f["bitmsk_relthrpt"] = bm
        end
        u = relthrpt_unusable_for_peaks(joinpath(d, "skyLinePeaks_apo_60000_0001_R_object.h5"))
        @test findall(u) == [7]
        @test !any(relthrpt_unusable_for_peaks(joinpath(d, "skyLinePeaks_apo_60000_0002_R_object.h5")))
    end
end

# Every writer of a waveCalNight*Dither product flags: the per-fiber sky dither (dporder 1), the
# per-fiber dither against the FPI solution (dporder 2) and the FPI comb fit's per-FPI-exposure
# products. Synthetic data with a known dither, fibers with 0, 3, 5, 6 and 12 lines.
@testset "every dither writer records line count and flag" begin
    nl = fill(12, N_FIBERS)
    nl[50] = 0; nl[100] = 3; nl[101] = 5; nl[150] = 6
    dtrue = [2.0e-4, 1.0 + 3.0e-5]
    lin, nlp = SW.night_params()
    mktempdir() do dir
        fsky = SW.write_sky_peaks(dir, "0007", nl; d = dtrue)
        bm = zeros(Int, N_FIBERS); bm[200] = RELTHRPT_BROKEN_BIT
        h5open(joinpath(dir, "ar1Dcal_apo_60000_0007_R_object.h5"), "w") do f
            f["bitmsk_relthrpt"] = bm
        end
        for (wavetype, dporder) in (("sky", 1), ("fpi", 2))
            out = get_and_save_sky_dither_per_fiber(fsky, lin, nlp; dporder = dporder,
                wavetype = wavetype, checkpoint_mode = "clobber")
            @test basename(out) == "waveCalNight$(wavetype)Dither_apo_60000_0007_object.h5"
            fl, n, dp = h5open(f -> (read(f["dither_flag"]), read(f["dither_nlines"]),
                    read(f["ditherParams"])), out)
            @test length(fl) == N_FIBERS && length(n) == N_FIBERS
            @test fl[50] == DITHER_FLAG_INTERP_NOLINES && n[50] == 0
            @test fl[100] == DITHER_FLAG_INTERP_FEWLINES && n[100] == 3
            @test fl[101] == DITHER_FLAG_INTERP_FEWLINES && n[101] == 5
            @test fl[150] == 0x00 && n[150] == 6
            @test fl[200] == DITHER_FLAG_INTERP_LOWTHRPT
            @test count(!=(0x00), fl) == 4
            # fitted fibers recover the dither; interpolated ones get their neighbours' value
            # (every fiber has the same true dither here)
            @test abs(dp[150, 1] - dtrue[1]) < 1e-6
            @test abs(dp[100, 1] - dp[99, 1]) < 1e-7
            @test abs(dp[200, 1] - dtrue[1]) < 1e-6
            # the reinterpolation stage copies the flags (README: bitmsk_wavecal)
            bits, nq = wavecal_fiber_quality(out)
            @test bits == Int.(fl) && nq == Int.(n)
        end
        @test wavecal_fiber_quality(nothing) == (fill(WAVECAL_BIT_NO_EXPOSURE_DITHER, N_FIBERS), fill(-1, N_FIBERS))
        old = joinpath(dir, "waveCalNightskyDither_apo_60000_0099_object.h5")
        h5open(f -> (f["chipWaveSoln"] = zeros(2, 2)), old, "w")
        @test first(wavecal_fiber_quality(old)) == fill(WAVECAL_BIT_FLAGS_MISSING, N_FIBERS)

        # FPI comb fit: two FPI exposures, the second dithered; fiber 20 dark, fiber 40 keeps only
        # 3 peaks in the second exposure
        np1 = fill(10_000, N_FIBERS)
        np2 = copy(np1); np2[40] = 3
        lit = trues(N_FIBERS); lit[20] = false
        f1 = SW.write_fpi_peaks(dir, "0001", np1; lit = lit)
        f2 = SW.write_fpi_peaks(dir, "0002", np2; lit = lit, d = [1.0e-4, 1.0])
        ret = comb_exp_get_and_save_fpi_wavecal([f1, f2], lin, nlp; FPI_COMB_FIT_KWARGS...,
            verbose = false, checkpoint_mode = "clobber")
        @test !isnothing(ret)
        for (e, f) in enumerate((f1, f2))
            out = replace(replace(f, "fpiPeaks" => "waveCalNightfpiDither"), "_R_" => "_")
            fl, n, dp = h5open(h -> (read(h["dither_flag"]), read(h["dither_nlines"]),
                    read(h["ditherParams"])), out)
            @test fl[20] == DITHER_FLAG_INTERP_NOLINES && n[20] == 0
            @test minimum(n[setdiff(1:N_FIBERS, [20, 40])]) >= DITHER_MIN_LINES   # counted per exposure
            if e == 2
                # 1 peak per chip: the comb's per-chip outlier mask may drop them, so 0-3 peaks used
                @test fl[40] in (DITHER_FLAG_INTERP_FEWLINES, DITHER_FLAG_INTERP_NOLINES) && n[40] <= 3
                @test dp[40, :] ≈ (dp[39, :] .+ dp[41, :]) ./ 2
            end
        end
    end
end

# every writer of a waveCalNight*Dither product (per-fiber sky/FPI dither and the FPI comb fit's
# per-exposure products) passes the per-fiber line count and flag (source check: a new writer
# without them fails here)
@testset "dither products carry dither_nlines and dither_flag (source)" begin
    src = readlines(joinpath(pkgdir(ApogeeReduction), "src", "wavecal.jl"))
    writers = findall(l -> occursin("=> \"waveCalNight", l), src)
    @test length(writers) >= 2
    for i in writers
        j = findnext(l -> occursin("safe_jldsave(", l), src, i)
        @test !isnothing(j)
        k = findnext(l -> occursin("no_metadata", l), src, j)
        call = join(src[j:k], " ")
        @test occursin("dither_nlines", call)
        @test occursin("dither_flag", call)
    end
end
