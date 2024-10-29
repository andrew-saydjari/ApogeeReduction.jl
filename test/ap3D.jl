@testset "ap3D" begin
    rng = MersenneTwister(101) #should we switch to stableRNG for my peace of mind?

    detector_dims = (2560, 2048) #TODO name better
    gainMat = 1.9 * ones(Float32, detector_dims)
    readVarMat = 25 * ones(Float32, detector_dims)

    n_reads = 20
    n_diffs = n_reads - 1

    flux_per_reads = 10 .^ (range(start = log10(0.01), stop = log10(10000), length = 2560))
    dcounts = (randn(rng, (detector_dims..., n_diffs)) .* (flux_per_reads .^ 0.5)) .+
              flux_per_reads

    # pepper with cosmic rays. These diffs should be excluded
    n_crs = 1000
    cr_count = 1e6
    crs = zeros(size(dcounts))
    crs[rand(rng, eachindex(dcounts), n_crs)] .= cr_count
    dcounts .+= crs
    cr_mask = sum(crs .> 0, dims = 3) .> 0

    true_im = ones(Float32, detector_dims) .* flux_per_reads

    outdat = zeros(Float32, (detector_dims..., n_reads))
    outdat[:, :, (begin + 1):end] .+= cumsum(dcounts, dims = 3)
    outdat .+= randn(rng, (detector_dims..., n_reads)) .* (readVarMat .^ 0.5)
    outdat ./= gainMat

    @time dimage, ivarimage, chisqimage = ApogeeReduction.sutr_wood!(
        outdat, gainMat, readVarMat)

    flux_diffs = (dimage .- true_im)
    zscores = flux_diffs .* sqrt.(ivarimage)

    # mean zscore should be 0, but it's biased at low fluxes
    @test isapprox(mean(zscores), 0.0, atol = 0.05)
    # at higher fluxes, the mean zscore should be unbiased
    @test isapprox(mean(zscores[(end - 500):end, :]), 0.0, atol = 0.003)
    # std(zscore) should be 1
    @test isapprox(std(zscores), 1, atol = 0.001)

    # were CRs handled well?
    @test isapprox(mean(zscores[cr_mask]), 0.0, atol = 0.05)
    @test isapprox(std(zscores[cr_mask]), 1, atol = 0.001)

    # these are specific to this implementation and random seed
    flux_mean_z = mean(zscores, dims = 2)
    @test isapprox(flux_mean_z[end], -0.04463, atol = 0.08)
    flux_std_z = std(zscores, dims = 2)
    @test isapprox(flux_std_z[end], 1.03251, atol = 0.06)
    flux_mean = mean(dimage, dims = 2)
    @test isapprox(flux_mean[end], 9998.98, atol = 3)
    flux_err_mean = mean(ivarimage .^ (-0.5), dims = 2)
    @test isapprox(flux_err_mean[begin], 0.20749, atol = 0.003)
    @test isapprox(flux_err_mean[end], 22.9450, atol = 0.003)
end
