@testset "ap3D" begin
    rng = MersenneTwister(101) #should we switch to stableRNG for my peace of mind?

    gainMat = 1.9 * ones(Float32, 2560, 2048)
    readVarMat = 25 * ones(Float32, 2560, 2048)

    n_reads = 20
    n_diffs = n_reads - 1

    flux_per_reads = 10 .^ (LinRange(log10(0.01), log10(10000), 2560))
    dcounts = (randn(rng, (2560, 2048, n_diffs)) .* (flux_per_reads .^ 0.5)) .+
              flux_per_reads

    true_im = ones(Float32, 2560, 2048) .* flux_per_reads

    outdat = zeros(Float32, (2560, 2048, n_reads))
    outdat[begin:end, begin:end, (begin + 1):end] .+= cumsum(dcounts, dims = 3)
    outdat .+= randn(rng, (2560, 2048, n_reads)) .* (readVarMat .^ 0.5)
    outdat ./= gainMat

    dimage, ivarimage, chisqimage = sutr_tb(outdat, gainMat, readVarMat)
    #    dimage, ivarimage, chisqimage = dcs(outdat,gainMat,readVarMat)

    flux_diffs = (dimage .- true_im)
    zscores = flux_diffs ./ (ivarimage .^ (-0.5))

    full_mean_z = mean(zscores)
    full_std_z = std(zscores)
    flux_mean_z = mean(zscores, dims = 2)
    flux_std_z = std(zscores, dims = 2)
    flux_mean = mean(dimage, dims = 2)
    flux_diff_mean = mean(flux_diffs, dims = 2)
    flux_err_mean = mean(ivarimage .^ (-0.5), dims = 2)

    #expected outputs when n_reads=20, seed=101, extract_method=sutr_tb
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
