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
    cr_count = 1e5
    dcounts[rand(eachindex(dcounts), 100)] .= cr_count

    true_im = ones(Float32, detector_dims) .* flux_per_reads

    outdat = zeros(Float32, (detector_dims..., n_reads))
    outdat[:, :, (begin + 1):end] .+= cumsum(dcounts, dims = 3)
    outdat .+= randn(rng, (detector_dims..., n_reads)) .* (readVarMat .^ 0.5)
    outdat ./= gainMat

    @time dimage, ivarimage, chisqimage = ApogeeReduction.sutr_tb!(
        outdat, gainMat, readVarMat)

    flux_diffs = (dimage .- true_im)
    zscores = flux_diffs .* sqrt.(ivarimage)

    full_mean_z = mean(zscores)
    full_std_z = std(zscores)
    flux_mean_z = mean(zscores, dims = 2)
    flux_std_z = std(zscores, dims = 2)
    flux_mean = mean(dimage, dims = 2)
    flux_diff_mean = mean(flux_diffs, dims = 2)
    flux_err_mean = mean(ivarimage .^ (-0.5), dims = 2)

    #expected outputs when n_reads=20, seed=101
    @test isapprox(full_mean_z, -0.05, atol = 0.001)
    @test isapprox(full_std_z, 1.0427, atol = 0.001)
    @test isapprox(flux_mean_z[begin], -0.1298, atol = 0.06)
    @test isapprox(flux_mean_z[end], -0.04463, atol = 0.08)
    @test isapprox(flux_std_z[begin], 1.0104, atol = 0.01)
    @test isapprox(flux_std_z[end], 1.03251, atol = 0.06)
    @test isapprox(flux_mean[begin], -0.002658, atol = 0.01)
    @test isapprox(flux_mean[end], 9998.98, atol = 3)
    @test isapprox(flux_err_mean[begin], 0.20749, atol = 0.003)
    @test isapprox(flux_err_mean[end], 23.57064, atol = 0.003)
end
