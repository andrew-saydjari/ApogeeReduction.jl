# Handling the 3D data cube
using LinearAlgebra: SymTridiagonal, Diagonal
using Statistics: mean

function apz2cube(fname)
    f = FITS(fname)
    hdr = read_header(f[2])
    avg_dcounts = read(f[2])
    cubedat = zeros(Float32, size(avg_dcounts)..., length(f) - 2) #use float bc of NaN
    cubedat[:, :, 1] .= read(f[3])
    for i in 2:(length(f) - 2)
        cubedat[:, :, i] .= read(f[i + 2]) .+ avg_dcounts .+ cubedat[:, :, i - 1]
    end
    close(f)
    return cubedat, hdr
end

# feeling against this sort of subtraction, but we do see why one might want to do it in the darks
function refcorr(dcubedat)
    # subtracts reference array with proper read ordering
    dcubedat_out = copy(dcubedat[1:2048, 1:2048, ..])
    dcubedat_out[1:512, :, ..] .-= dcubedat[2049:end, ..]
    dcubedat_out[513:1024, :, ..] .-= dcubedat[end:-1:2049, ..]
    dcubedat_out[1025:1536, :, ..] .-= dcubedat[2049:end, ..]
    dcubedat_out[1537:2048, :, ..] .-= dcubedat[end:-1:2049, ..]
    return dcubedat_out
end

function refarray_zpt!(dcubedat)
    mean_ref_vec = mean(dcubedat[2049:end, :, ..], dims = (1, 2))
    dcubedat .-= mean_ref_vec
    mean_ref_val = mean(dcubedat[2049:end, ..])
    dcubedat .-= mean_ref_val
    return
end

function vert_ref_edge_corr_amp!(dcubedat_out)
    dcubedat_out[1:512, ..] .-= mean([mean(dcubedat_out[1:512, 1:4, ..]),
        mean(dcubedat_out[1:512, (end - 3):(end - 1), ..])])
    dcubedat_out[513:1024, ..] .-= mean([mean(dcubedat_out[513:1024, 1:4, ..]),
        mean(dcubedat_out[513:1024, (end - 3):(end - 1), ..])])
    dcubedat_out[1025:1536, ..] .-= mean([mean(dcubedat_out[1025:1536, 1:4, ..]),
        mean(dcubedat_out[1025:1536, (end - 3):(end - 1), ..])])
    dcubedat_out[1537:2048, ..] .-= mean([mean(dcubedat_out[1537:2048, 1:4, ..]),
        mean(dcubedat_out[1537:2048, (end - 3):(end - 1), ..])])
    dcubedat_out[2049:end, ..] .-= mean([mean(dcubedat_out[2049:end, 1:4, ..]),
        mean(dcubedat_out[2049:end, (end - 3):(end - 1), ..])])
    return
end

function dcs(dcubedat, gainMat, readVarMat; firstind = 1)
    ndiffs = size(dcubedat, 3) - firstind
    dimage = gainMat .* (dcubedat[:, :, end] .- dcubedat[:, :, firstind])
    # bad to use measured flux as the photon noise
    ivarimage = 1 ./ (2 .* readVarMat .+ abs.(dimage))
    return dimage ./ ndiffs, (ndiffs .^ 2) .* ivarimage, nothing # no chisq, just mirroring sutr_tb
end

"""
    TODO

# Arguments
- `datacube` has shape (npix_x,npix_y,n_reads)
- `gainMat`: The gain for each pixel (npix_x,npix_y)
- `read_var_mat`: the read noise (as a variance) for each pixel (npix_x,npix_y)

# Keyword Arguments
- firstind TODO
- good_diffs, if provided, should contain booleans and have shape (npix_x, npix_y, n_reads-1). If
  not provided, all diffs will be used.

!!! note
    This assumes that all images are sequential (ie separated by the same read time).

This is loosely based on the algorithm described in 
[Brandt 2024](https://doi.org/10.1088/1538-3873/ad38d9), but with significant simplifications arises 
from the fact that all reads are sequential and that we don't need to group reads to save bandwidth,
 unlike JWST.
"""
function sutr(
        dcubedat, gainMat, readVarMat; firstind = 1, good_diffs = nothing, n_iter = 2)
    dimages = gainMat .* diff(dcubedat[:, :, firstind:end], dims = 3)

    # TODO use good_diffs
    if isnothing(good_diffs)
        good_diffs = ones(Bool, size(dimages))
    end
    ndiffs = size(dimages, 3)

    # this the starting guess, it will get refined to approach the maximum likelyhood estimate
    rate = mean(dimages; dims = 3)
    rate[rate .< 0] .= 0

    # working variables allocate here and reuse inside loop
    ones_vec = ones(ndiffs)
    α = 0.0
    β = 0.0
    c_prime = ones(ndiffs - 1)
    d_prime = ones(ndiffs)
    x = ones(ndiffs)
    residuals = ones(ndiffs)

    # returned arrays
    ivars = Matrix{Float64}(undef, size(dimages)[1:2])
    chi2s = Matrix{Float64}(undef, size(dimages)[1:2])
    for pixel_ind in CartesianIndices(dimages[:, :, 1])
        β = -readVarMat[pixel_ind]

        for _ in 1:n_iter
            α = 2 * readVarMat[pixel_ind] + rate[pixel_ind]
            solve_tridiagonal_toeplitz!(x, c_prime, d_prime, α, β, ones_vec)
            ivars[pixel_ind] = sum(x)

            solve_tridiagonal_toeplitz!(x, c_prime, d_prime, α, β, dimages[pixel_ind, :])
            rate[pixel_ind] = sum(x) / ivars[pixel_ind]
        end

        residuals .= dimages[pixel_ind, :] .- rate[pixel_ind]
        solve_tridiagonal_toeplitz!(x, c_prime, d_prime, α, β, residuals)
        chi2s[pixel_ind] = residuals' * x
    end
    rate, ivars, chi2s
end

"""
    solve_tridiagonal_toeplitz!(x, c_prime, d_prime, α, β, b)

Computes C^{-1}b for a tridiagonal Toeplitz matrix C, with diagonal value α and off-diagonal value 
β.  When forward modelling photon counts given reads at a steady interal (sequential), assuming the 
Poisson errors can be approximated as Gaussian, the read noise is constant, and flux is unchanging, 
the covariance matrix, C,  corresponding to the data  is a tridiagonal Toeplitz matrix. 

In computing the maximum likelihood estimate of the rate, the χ², and the estimated error of the 
MLE, we need to solve the equation Cx = b. This function solves this equation given the above 
assumptions about C using the Thomas algorithm.

Arguments:
- `x` is the output vector
- `c_prime` and `d_prime` are working vectors (notation borrowed from wikipedia)
- `α` is the diagonal value of the tridiagonal Toeplitz matrix
- `β` is the off-diagonal value of the tridiagonal Toeplitz matrix
- `b` is the right-hand side of the equation
"""
function solve_tridiagonal_toeplitz!(x, c_prime, d_prime, α, β, b)
    n = length(b)

    # Forward elimination
    c_prime[1] = β / α
    d_prime[1] = b[1] / α
    for i in 2:(n - 1)
        denominator = α - β * c_prime[i - 1]
        c_prime[i] = β / denominator
        d_prime[i] = (b[i] - β * d_prime[i - 1]) / denominator
    end
    d_prime[n] = (b[n] - β * d_prime[n - 1]) / (α - β * c_prime[n - 1])

    # Back substitution
    x[n] = d_prime[n]
    for i in (n - 1):-1:1
        x[i] = d_prime[i] - c_prime[i] * x[i + 1]
    end
end