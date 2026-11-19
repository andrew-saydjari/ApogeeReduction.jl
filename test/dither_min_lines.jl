using ApogeeReduction: dither_fit_mode, fit_fiber_dither, relthrpt_unusable_for_peaks,
                       DITHER_MIN_LINES, DITHER_MIN_LINES_FULL, DITHER_FLAG_INTERP_NOLINES,
                       DITHER_FLAG_INTERP_FEWLINES, DITHER_FLAG_OFFSET_ONLY, DITHER_FLAG_INTERP_LOWTHRPT,
                       RELTHRPT_BROKEN_BIT, RELTHRPT_WARN_BIT, positional_poly_mat, N_FIBERS
using HDF5
using ApogeeReduction: ApogeeReduction

@testset "per-fiber dither: minimum line counts" begin
    @test dither_fit_mode(0, false) == (:interp, DITHER_FLAG_INTERP_NOLINES)
    @test dither_fit_mode(1, false) == (:interp, DITHER_FLAG_INTERP_FEWLINES)
    @test dither_fit_mode(DITHER_MIN_LINES - 1, false) == (:interp, DITHER_FLAG_INTERP_FEWLINES)
    @test dither_fit_mode(DITHER_MIN_LINES, false) == (:offset, DITHER_FLAG_OFFSET_ONLY)
    @test dither_fit_mode(DITHER_MIN_LINES_FULL - 1, false) == (:offset, DITHER_FLAG_OFFSET_ONLY)
    @test dither_fit_mode(DITHER_MIN_LINES_FULL, false) == (:full, 0x00)
    @test dither_fit_mode(20, false) == (:full, 0x00)
    @test dither_fit_mode(20, true) == (:interp, DITHER_FLAG_INTERP_LOWTHRPT)   # unusable fiber, even with lines
    @test dither_fit_mode(4, false; dporder = 0) == (:full, 0x00)               # nothing but the offset anyway

    # synthetic fiber: quartic wavelength solution in the normalised pixel coordinate
    lin = [16000.0, -580.0, 3.0, 0.5, 0.1]
    lam(x) = (positional_poly_mat([x], porder = 4) * lin)[1]
    d_true = [2.0e-4, 1.0 + 3.0e-5]
    xt = collect(range(-0.45, 0.45, length = 18))
    y = lam.(d_true[1] .+ d_true[2] .* xt)
    p = fit_fiber_dither(xt, y, lin, [0.0, 1.0]; wporder = 4)
    @test abs(p[1] - d_true[1]) < 1e-8
    @test abs(p[2] - d_true[2]) < 1e-7
    # offset only: the scale stays at its start value, the offset absorbs the shift
    y0 = lam.(d_true[1] .+ xt[1:4])
    po = fit_fiber_dither(xt[1:4], y0, lin, [0.0, 1.0]; wporder = 4, offset_only = true)
    @test po[2] == 1.0
    @test abs(po[1] - d_true[1]) < 1e-8

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

# every writer of a waveCalNight*Dither product (per-fiber sky/FPI dither and the FPI comb fit's
# per-exposure products) must write the per-fiber line count and flag
@testset "dither products carry dither_nlines and dither_flag" begin
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
