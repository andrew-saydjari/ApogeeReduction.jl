@testset "ap3D" begin
    # not guarenteed to be reproducible across julia versions
    rng = MersenneTwister(101)

    detector_dims = (2560, 2048)
    gainMat = 1.9 * ones(Float64, detector_dims)
    readVarMat = 25 * ones(Float64, detector_dims)
    # pretend the pixels saturate at 90% of digital saturation value
    saturationMat = 2^16 * ones(Float64, detector_dims) * 0.9

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
    true_im = ones(Float32, detector_dims) .* flux_per_reads ./ gainMat

    datacube = zeros(Float32, (detector_dims..., n_reads))
    datacube[:, :, (begin + 1):end] .+= cumsum(dcounts, dims = 3)
    datacube .+= randn(rng, (detector_dims..., n_reads)) .* (readVarMat .^ 0.5) .* gainMat
    datacube ./= gainMat

    # compute the last unsaturated read for each pixel.
    last_unsaturated = ApogeeReduction.get_last_unsaturated_read(datacube, saturationMat)

    dimages = ApogeeReduction.diffify_datacube!(datacube, last_unsaturated)
    # try to identify any cosmic rays
    not_cosmic_ray = ApogeeReduction.outlier_mask(dimages, last_unsaturated)
    CRimage = sum(.!not_cosmic_ray, dims = 3)[:, :, 1]

    dimage, ivarimage, chisqimage = ApogeeReduction.sutr_wood(
        dimages, gainMat, readVarMat, last_unsaturated, not_cosmic_ray)

    # chop off the reference array for the tests
    dimage = dimage[1:2048, :]
    ivarimage = ivarimage[1:2048, :]
    chisqimage = chisqimage[1:2048, :]
    CRimage = CRimage[1:2048, :]
    last_unsaturated = last_unsaturated[1:2048, :]
    true_im = true_im[1:2048, :]
    true_cr_mask = true_cr_mask[1:2048, :, 1]

    flux_diffs = (dimage .- true_im)
    zscores = flux_diffs .* sqrt.(ivarimage)

    high_flux_mask = (true_im .> 1000)

    @test all(last_unsaturated .== n_diffs)

    # nearly everything flagged as CR is actually CR
    # as the number of reads increases, this quantity should approach one.
    @test mean(true_cr_mask[CRimage .> 0]) > 0.97

    # nearly all CR pixels are flagged
    @test mean(CRimage[true_cr_mask])≈1 rtol=2e-3

    # mean zscore should be 0, but it's biased at low fluxes
    @test isapprox(mean(zscores), 0.0, atol = 0.05)
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
