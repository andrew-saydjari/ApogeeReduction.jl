using ApogeeReduction
using Test
using Base
using Random

src_dir = "../"
include(src_dir * "src/ap3D.jl")

using LinearAlgebra
BLAS.set_num_threads(1)
using FITSIO, HDF5, FileIO, JLD2
using DataFrames, EllipsisNotation, StatsBase

@testset "ApogeeReduction.jl" begin
    # Write your tests here.

    rng = MersenneTwister(101)

    gainMat = 1.9*ones(Float32,2560,2048)
    readVarMat = 25*ones(Float32,2560,2048)

    n_reads = 20
    n_diffs = n_reads-1

    flux_per_reads = 10 .^ (LinRange(log10(0.01),log10(10000),2560))
    dcounts = (randn(rng,(2560,2048,n_diffs)) .* (flux_per_reads .^ 0.5)) .+ flux_per_reads

    true_im = ones(Float32,2560,2048) .* flux_per_reads

    outdat = zeros(Float32,(2560,2048,n_reads))
    outdat[begin:end,begin:end,begin+1:end] .+= cumsum(dcounts,dims=3)
    outdat .+= randn(rng,(2560,2048,n_reads)) .* (readVarMat .^ 0.5)
    outdat ./= gainMat

    dimage, ivarimage, chisqimage = sutr_tb(outdat,gainMat,readVarMat)
#    dimage, ivarimage, chisqimage = dcs(outdat,gainMat,readVarMat)

    flux_diffs = (dimage .- true_im)
    zscores = flux_diffs ./ (ivarimage .^ (-0.5))

    full_mean_z = mean(zscores)
    full_std_z = std(zscores)
    flux_mean_z = mean(zscores,dims=2)
    flux_std_z = std(zscores,dims=2)
    flux_mean = mean(dimage,dims=2)
    flux_diff_mean = mean(flux_diffs,dims=2)
    flux_err_mean = mean(ivarimage .^ (-0.5),dims=2)

    #expected outputs when n_reads=20, seed=101, extract_method=sutr_tb
    #and the flux_per_reads go from 0.01 to 10000 in 2560 log10 steps
    @test isapprox(full_mean_z,-0.05,atol=0.001)
    @test isapprox(full_std_z,1.0427,atol=0.001)
    @test isapprox(flux_mean_z[begin],-0.1298,atol=0.001)
    @test isapprox(flux_mean_z[end],-0.04463,atol=0.001)
    @test isapprox(flux_std_z[begin],1.0104,atol=0.001)
    @test isapprox(flux_std_z[end],1.03251,atol=0.001)
    @test isapprox(flux_mean[begin],-0.002658,atol=0.00001)
    @test isapprox(flux_mean[end],9998.98,atol=0.1)
    @test isapprox(flux_err_mean[begin],0.20749,atol=0.0001)
    @test isapprox(flux_err_mean[end],23.57064,atol=0.001)

#    println("Flux Bottom Mean ",flux_mean[begin:begin+3])
#    println("Flux Bottom Diff Mean ",flux_diff_mean[begin:begin+3])
#    println("Flux Bottom Mean Err ",flux_err_mean[begin:begin+3])
#    println("Flux Top Mean ",flux_mean[end-3:end])
#    println("Flux Top Diff Mean ",flux_diff_mean[end-3:end])
#    println("Flux Top Mean Err ",flux_err_mean[end-3:end],"\n")
#    println("Z Score Mean ",full_mean_z)
#    println("Z Score Std ",full_std_z)
#    println("Z Score Bottom Mean ",flux_mean_z[begin:begin+3])
#    println("Z Score Bottom Std ",flux_std_z[begin:begin+3])
#    println("Z Score Top Mean ",flux_mean_z[end-3:end])
#    println("Z Score Top Std ",flux_std_z[end-3:end])

end


