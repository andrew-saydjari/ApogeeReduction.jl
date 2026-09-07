@testset "ap3D" begin
    # note: this whole test is implicitely done with gain = 1.0 for all pixels, because using a 
    # non-unity would complicate the test while adding little of substance.
    detector_dims = (2560, 2048)
    readVarMat = 25 * ones(Float64, detector_dims)

    # ------------------------------------------------------------
    # build fake image
    # ------------------------------------------------------------
    # not guarenteed to be reproducible across julia versions
    rng = MersenneTwister(101)

    # pretend the pixels saturate at 90% of digital saturation value
    saturationMat = 2^16 * ones(Float64, detector_dims) * 0.9
    # set 40,000 (~1%) of pixels to saturation at a low value
    low_saturation_pixels = rand(rng, eachindex(saturationMat), 40000)
    saturationMat[low_saturation_pixels] .= 2^16 * 0.1

    n_reads = 11
    n_diffs = n_reads - 1

    flux_per_reads = vcat(
        10 .^ (range(start = log10(10), stop = log10(4000), length = 2048)), #real pixels
        fill(0.0, 2560 - 2048) # reference array
    )
    # should this have read noise?
    dcounts = (randn(rng, (detector_dims..., n_diffs)) .* (flux_per_reads .^ 0.5)) .+ flux_per_reads

    # pepper with cosmic rays.
    n_crs = 10000
    cr_count = 1e4
    crs = zeros(size(dcounts))
    crs[rand(rng, eachindex(dcounts), n_crs)] .= cr_count
    crs[2049:end, :, :] .= 0 # no CRs in reference array
    dcounts .+= crs
    true_cr_mask = sum(crs .> 0, dims = 3) .> 0

    # without noise, cosmic rays, or saturation effects
    true_im = ones(Float32, detector_dims) .* flux_per_reads

    datacube = zeros(Float32, (detector_dims..., n_reads))
    datacube[:, :, (begin + 1):end] .+= cumsum(dcounts, dims = 3)

    # saturate pixels that are above saturationMat, and record the real last unsaturated read for 
    # each pixel
    true_last_unsaturated = Matrix{Int}(undef, size(saturationMat))
    for I in CartesianIndices(true_last_unsaturated)
        last_unsaturated = findlast(datacube[I, :] .< saturationMat[I])
        # subtract 1 because of diffs vs reads
        true_last_unsaturated[I] = last_unsaturated - 1
        datacube[I, (last_unsaturated + 1):n_reads] .= saturationMat[I]
    end

    datacube .+= randn(rng, (detector_dims..., n_reads)) .* (readVarMat .^ 0.5)

    # ------------------------------------------------------------
    # fit the synthetic datacube
    # ------------------------------------------------------------

    # compute the last unsaturated read for each pixel.
    last_unsaturated = ApogeeReduction.get_last_unsaturated_read(datacube, saturationMat)

    dimages = ApogeeReduction.diffify_datacube!(datacube, last_unsaturated)
    # try to identify any cosmic rays
    not_cosmic_ray = ApogeeReduction.outlier_mask(dimages, last_unsaturated)
    CRimage = sum(.!not_cosmic_ray, dims = 3)[:, :, 1]

    dimage, ivarimage, chisqimage = ApogeeReduction.sutr_wood(
        dimages, ones(size(saturationMat)), readVarMat, last_unsaturated, not_cosmic_ray)

    # ------------------------------------------------------------
    # chop off the reference array for the tests
    # ------------------------------------------------------------
    dimage = dimage[1:2048, :]
    ivarimage = ivarimage[1:2048, :]
    chisqimage = chisqimage[1:2048, :]
    CRimage = CRimage[1:2048, :]
    true_im = true_im[1:2048, :]
    true_cr_mask = true_cr_mask[1:2048, :, 1]

    flux_diffs = (dimage .- true_im)
    zscores = flux_diffs .* sqrt.(ivarimage)
    high_flux_mask = (true_im .> 1000)

    saturated = true_last_unsaturated .< n_diffs

    # were the last unsaturated reads correctly identified?
    @test mean(last_unsaturated .== true_last_unsaturated) > 0.997
    saturated = true_last_unsaturated .< n_diffs
    let l = true_last_unsaturated[saturated]
        # some error comes from read noise, some from the fudge-factor in the saturation flagging
        @test mean((l .- 1) .<= last_unsaturated[saturated] .<= (l .+ 1)) == 1
    end

    saturated = saturated[1:2048, 1:2048]

    # nearly everything flagged as CR is actually CR
    # as the number of reads increases, this quantity should approach one.
    @test mean(true_cr_mask[CRimage .> 0]) > 0.97

    # nearly all CR pixels are flagged
    @test mean(CRimage[true_cr_mask]) > 0.99

    # mean zscore should be 0, but it's biased at low fluxes
    @test isapprox(mean(zscores), 0.0, atol = 0.05)
    @test isapprox(mean(zscores[saturated]), 0.0, atol = 0.05)
    # at higher fluxes, the mean zscore should be unbiased
    @test isapprox(mean(zscores[high_flux_mask]), 0.0, atol = 3e-3)
    # std(zscore) should be 1
    @test isapprox(std(zscores), 1, atol = 0.007) #0.001
    # are the fluxes correct?
    @test isapprox(mean(dimage ./ true_im), 1, atol = 3e-3) #1e-3
    # less biased at high fluxes
    @test isapprox(mean((dimage ./ true_im)[high_flux_mask]), 1, atol = 3e-4)

    # Same tests, but only for CR pixels, and with looser tolerances because there are fewer
    @test isapprox(mean(zscores[true_cr_mask]), 0.0, atol = 0.05)
    @test isapprox(mean(zscores[true_cr_mask .& high_flux_mask]), 0.0, atol = 1.6e-2) # 1e-2
    @test isapprox(std(zscores[true_cr_mask]), 1, atol = 0.03) #TODO not passing?
    @test isapprox(mean((dimage ./ true_im)[true_cr_mask]), 1, atol = 0.1)
    @test isapprox(mean((dimage ./ true_im)[true_cr_mask .& high_flux_mask]), 1, atol = 1e-3)

    # these are specific to this random seed.  Some may be implementation-specific.
    # it is reasonable to delete any that seem more implementation-specific
    flux_mean_z = mean(zscores, dims = 2)
    @test isapprox(flux_mean_z[end], -0.04463, atol = 0.08)
    flux_std_z = std(zscores, dims = 2)
    @test isapprox(flux_std_z[end], 1.03251, atol = 0.06)
    flux_mean = mean(dimage, dims = 2)
end

@testset "dcs" begin
    # A2 regression: the DCS ivar photon term must be dimage ./ gainMat (DN^2
    # Poisson variance), not gainMat ./ dimage (which made the ivar
    # read-noise-only and inflated DCS SNRs by 10-1000x).
    rng = MersenneTwister(202)
    npix = (200, 200)
    gain = 1.9
    readVar = 25.0
    gainMat = fill(gain, npix)
    readVarMat = fill(readVar, npix)

    n_reads = 4
    ndiffs = n_reads - 1
    flux_e_per_read = 1000.0 # electrons per read: photon noise dominates read noise

    # accumulated electron counts with independent (Gaussian-approximated)
    # Poisson increments, converted to DN, plus per-read read noise
    ncube = zeros(npix..., n_reads)
    for r in 2:n_reads
        ncube[:, :, r] .= ncube[:, :, r - 1] .+ flux_e_per_read .+
                          sqrt(flux_e_per_read) .* randn(rng, npix)
    end
    dcube = ncube ./ gain .+ sqrt(readVar) .* randn(rng, (npix..., n_reads))

    dimage, ivarimage, chisqimage, CRimage = ApogeeReduction.dcs(dcube, gainMat, readVarMat)

    # exact analytic contract: per-read flux and Poisson+read ivar
    dtot = dcube[:, :, end] .- dcube[:, :, 1]
    @test dimage ≈ dtot ./ ndiffs
    @test ivarimage ≈ ndiffs^2 ./ (2 * readVar .+ max.(dtot, 0) ./ gain)
    @test all(ivarimage .> 0) && all(isfinite.(ivarimage))
    @test all(chisqimage .== 0)
    @test all(CRimage .== 0)

    # statistical calibration: z-scores of the measured per-read flux against
    # the truth must be ~unit normal. With the pre-fix (inverted) photon term
    # the reported ivar is read-noise-only here and std(z) comes out ~4.2.
    truth = flux_e_per_read / gain
    z = (dimage .- truth) .* sqrt.(ivarimage)
    @test isapprox(mean(z), 0.0, atol = 0.02)
    @test isapprox(std(z), 1.0, atol = 0.03)

    # dimage <= 0 guard: photon term dropped, ivar = read-noise-only, positive
    dcube_neg = zeros(2, 2, 2)
    dcube_neg[:, :, 2] .= [-50.0 -1e-3; 0.0 10.0]
    d2, iv2, _, _ = ApogeeReduction.dcs(dcube_neg, fill(gain, 2, 2), fill(readVar, 2, 2))
    @test d2 == dcube_neg[:, :, 2]
    @test iv2[1, 1] == iv2[1, 2] == iv2[2, 1] == 1 / (2 * readVar)
    @test iv2[2, 2] == 1 / (2 * readVar + 10.0 / gain)
    @test all(iv2 .> 0) && all(isfinite.(iv2))
end

@testset "detector calibration guards" begin
    # Regression test for the defect described in
    #   .../2026_09_06/lco_gain_investigation/LCO_ERROR_REPORT.md
    #   .../2026_09_06/lco_gain_recovery/LCO_CALIB_RECOVERY.md
    # From 2025-06-11 to 2026-09-06 every LCO reduction ran its 2D error model on
    # byte-identical copies of APO's gain and read-noise maps. The pre-existing
    # @warn-plus-fallback only fires on a MISSING file, so the placeholder copies
    # defeated it silently for 15 months. These tests assert that the copied-file
    # situation now raises.

    writefits(path, dat) = begin
        f = FITS(path, "w")
        write(f, dat)
        close(f)
    end

    rng = MersenneTwister(4242)

    # ------------------------------------------------------------------
    # unit-level: the guard itself, on small arrays (fast)
    # ------------------------------------------------------------------
    mktempdir() do dir
        n = 40
        apo = 1.80 .+ 0.05 .* randn(rng, n, n)
        writefits(joinpath(dir, "gain_apo_R.fits"), apo)

        # (1) an exact copy must raise, from either telescope's point of view
        writefits(joinpath(dir, "gain_lco_R.fits"), copy(apo))
        lcoPath = joinpath(dir, "gain_lco_R.fits")
        err = try
            ApogeeReduction.assert_calib_map_telescope_specific(
                lcoPath, "gain", "lco", "R", copy(apo))
            nothing
        catch e
            e
        end
        @test err isa ErrorException
        # the message must name BOTH files and be actionable
        @test occursin("gain_lco_R.fits", err.msg)
        @test occursin("gain_apo_R.fits", err.msg)
        @test occursin("md5sum", err.msg)
        # symmetric: loading apo against the lco copy raises too
        @test_throws ErrorException ApogeeReduction.assert_calib_map_telescope_specific(
            joinpath(dir, "gain_apo_R.fits"), "gain", "apo", "R", copy(apo))

        # (2) NOT defeatable by perturbing a few pixels. This is the real hazard:
        #     the deployed pass_clean APO files differ from their pass4 parents by
        #     as little as ONE pixel out of 4.16e6, so an exact-hash check alone
        #     would be defeated by any post-copy cleaning pass.
        nearcopy = copy(apo)
        nearcopy[1:8] .+= 1.0   # 8/1600 = 0.5% of pixels differ
        writefits(joinpath(dir, "gain_lco_R.fits"), nearcopy)
        @test_throws ErrorException ApogeeReduction.assert_calib_map_telescope_specific(
            lcoPath, "gain", "lco", "R", nearcopy)

        # (3) NaNs in identical positions still count as identical (isequal, not ==)
        nanapo = copy(apo)
        nanapo[3, 3] = NaN
        writefits(joinpath(dir, "gain_apo_R.fits"), nanapo)
        writefits(joinpath(dir, "gain_lco_R.fits"), copy(nanapo))
        @test_throws ErrorException ApogeeReduction.assert_calib_map_telescope_specific(
            lcoPath, "gain", "lco", "R", copy(nanapo))

        # (4) genuinely different maps must pass (LCO really does have a different gain)
        writefits(joinpath(dir, "gain_apo_R.fits"), apo)
        lco = 2.70 .+ 0.05 .* randn(rng, n, n)
        writefits(joinpath(dir, "gain_lco_R.fits"), lco)
        @test isnothing(ApogeeReduction.assert_calib_map_telescope_specific(
            lcoPath, "gain", "lco", "R", lco))

        # (5) a different chip's file is not compared against
        writefits(joinpath(dir, "gain_apo_G.fits"), lco)
        @test isnothing(ApogeeReduction.assert_calib_map_telescope_specific(
            lcoPath, "gain", "lco", "R", lco))
    end

    # ------------------------------------------------------------------
    # plausibility trip-wire
    # ------------------------------------------------------------------
    for med in (1.45, 1.81, 2.56, 2.70)  # real APO and LCO gains
        @test isnothing(ApogeeReduction.assert_calib_map_plausible(
            "p", "gain", "lco", "R", med,
            ApogeeReduction.GAIN_MEDIAN_BOUNDS_E_PER_DN, "e-/DN"))
    end
    for med in (1e-3, 0.0, 1e4, NaN)     # units error / reciprocal / corrupt
        @test_throws ErrorException ApogeeReduction.assert_calib_map_plausible(
            "p", "gain", "lco", "R", med,
            ApogeeReduction.GAIN_MEDIAN_BOUNDS_E_PER_DN, "e-/DN")
    end
    for med in (2.93, 5.24, 8.01, 11.04) # real APO and LCO read noise, DN
        @test isnothing(ApogeeReduction.assert_calib_map_plausible(
            "p", "rdnoise", "lco", "R", med,
            ApogeeReduction.READ_NOISE_MEDIAN_BOUNDS_DN, "DN"))
    end
    for med in (0.01, 1e6, NaN)
        @test_throws ErrorException ApogeeReduction.assert_calib_map_plausible(
            "p", "rdnoise", "lco", "R", med,
            ApogeeReduction.READ_NOISE_MEDIAN_BOUNDS_DN, "DN")
    end

    # ------------------------------------------------------------------
    # end-to-end: the guard is actually wired into the load path, at the real
    # 2040x2040 map size the loaders require. This is the test that would have
    # caught the 15-month defect.
    # ------------------------------------------------------------------
    mktempdir() do dir
        caldir = dir * "/"  # loaders build paths by string concatenation
        m = 2040
        gapo = 1.80 .+ 0.05 .* randn(rng, m, m)
        rapo = 11.0 .+ 0.4 .* randn(rng, m, m)
        writefits(joinpath(dir, "gain_apo_R.fits"), gapo)
        writefits(joinpath(dir, "rdnoise_apo_R.fits"), rapo)

        # the exact 2026-09-06 failure: LCO files are copies of APO's
        writefits(joinpath(dir, "gain_lco_R.fits"), copy(gapo))
        writefits(joinpath(dir, "rdnoise_lco_R.fits"), copy(rapo))
        @test_throws ErrorException ApogeeReduction.load_gain_maps(caldir, "lco", "R")
        @test_throws ErrorException ApogeeReduction.load_read_var_maps(caldir, "lco", "R")

        # with the real (recovered) LCO values in place, both load cleanly and
        # carry the LCO numbers, not APO's
        glco = 2.556 .+ 0.05 .* randn(rng, m, m)
        rlco = 3.877 .+ 0.3 .* randn(rng, m, m)
        writefits(joinpath(dir, "gain_lco_R.fits"), glco)
        writefits(joinpath(dir, "rdnoise_lco_R.fits"), rlco)
        gd = ApogeeReduction.load_gain_maps(caldir, "lco", "R")
        rd = ApogeeReduction.load_read_var_maps(caldir, "lco", "R")
        @test size(gd["R"]) == (2560, 2048)
        @test gd["R"][5:2044, 5:2044] == glco
        @test isapprox(median(gd["R"]), 2.556, atol = 0.02)
        @test isapprox(sqrt(median(rd["R"])), 3.877, atol = 0.05)
        # APO still loads, and is unchanged by all of this
        gda = ApogeeReduction.load_gain_maps(caldir, "apo", "R")
        @test gda["R"][5:2044, 5:2044] == gapo
    end
end
