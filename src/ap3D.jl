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
Created by Kevin McKinnon on Aug 20, 2024, slightly refactored by Adam Wheeler.

Based on Tim Brandt's [SUTR python code](https://github.com/t-brandt/fitramp) described in 
[Brandt 2024](https://doi.org/10.1088/1538-3873/ad38d9).

This assumes that all images are sequential (ie separated by one read time).

# Arguments
- datacube has shape (npix_x,npix_y,n_reads)
- gainMat TODO
- read_var_mat has shape (npix_x,npix_y)

# Keyword Arguments
- firstind TODO
- good_diffs, if provided, should contain booleans and have shape (npix_x, npix_y, n_reads-1). If
  not provided, all diffs will be used.
"""
function sutr_reference(dcubedat, gainMat, readVarMat; firstind = 1, good_diffs = nothing)
    dimages = gainMat .* diff(dcubedat[:, :, firstind:end], dims = 3)

    if isnothing(good_diffs)
        good_diffs = ones(Bool, size(dimages))
    end

    ndiffs = size(dimages, 3)

    alpha_readnoise = 2
    beta_readnoise = -1

    alpha_phnoise = ones(ndiffs)
    beta_phnoise = zeros(ndiffs - 1)

    final_countrates = zeros(Float64, size(dcubedat, 1), size(dcubedat, 2))
    final_vars = zeros(Float64, size(dcubedat, 1), size(dcubedat, 2))
    final_chisqs = zeros(Float64, size(dcubedat, 1), size(dcubedat, 2))

    #slice the datacube to analyze each row sequentially to keep runtime down
    for s_ind in 1:size(dimages, 1)
        diffs = transpose(dimages[s_ind, :, :]) # shape = (ndiffs, npix)
        read_var = readVarMat[s_ind, :] # shape = (npix)
        read_var = reshape(read_var, (1, size(read_var, 1)))
        diffs2use = transpose(good_diffs[s_ind, :, :]) # shape = (ndiffs, npix)
        d = diffs .* diffs2use # diffs with bad values set to 0 

        # Repeat the process after getting a good first guess to remove a bias (as 
        # discussed in the paper).
        for r_ind in 1:2 # DO NOT CHANGE THIS!
            if r_ind == 1
                #initialize first guess of countrates with mean
                countrateguess = sum(d, dims = 1) ./ sum(diffs2use, dims = 1)
                countrateguess .*= (countrateguess .> 0)
            else
                countrateguess = (final_countrates[s_ind, :] .*
                                  (final_countrates[s_ind, :] .> 0))'
            end

            # Elements of the covariance matrix
            alpha = countrateguess .* alpha_phnoise .+ read_var .* alpha_readnoise
            beta = countrateguess .* beta_phnoise .+ read_var .* beta_readnoise

            # rescale the covariance matrix to a determinant of order 1 to
            # avoid possible overflow/underflow.  The uncertainty and chi
            # squared value will need to be scaled back later.
            scale = exp.(mean(log.(alpha), dims = 1))
            alpha ./= scale
            beta ./= scale

            ndiffs, npix = size(alpha)
            # Mask resultant differences that should be ignored.  This is half
            # of what we need to do to mask these resultant differences; the
            # rest comes later.
            beta .= beta .* diffs2use[2:end, :] .* diffs2use[1:(end - 1), :]

            # All definitions and formulas here are in the paper.
            theta = ones(Float64, ndiffs + 1, npix)
            theta[2, :] .= alpha[1, :]
            for i in 3:(ndiffs + 1)
                theta[i, :] .= alpha[i - 1, :] .* theta[i - 1, :] .-
                               beta[i - 2, :] .^ 2 .* theta[i - 2, :]
            end

            phi = ones(Float64, ndiffs + 1, npix)
            phi[ndiffs, :] .= alpha[ndiffs, :]
            for i in (ndiffs - 1):-1:1
                phi[i, :] .= alpha[i, :] .* phi[i + 1, :] .-
                             beta[i, :] .^ 2 .* phi[i + 2, :]
            end

            sgn = ones(Float64, ndiffs, npix)
            sgn[1:2:end, :] .= -1

            Phi = zeros(Float64, ndiffs, npix)
            for i in (ndiffs - 1):-1:1
                Phi[i, :] .= Phi[i + 1, :] .* beta[i, :] .+
                             sgn[i + 1, :] .* beta[i, :] .* phi[i + 2, :]
            end

            # This one is defined later in the paper and is used for jump
            # detection and pedestal fitting.

            PhiD = zeros(Float64, ndiffs, npix)
            for i in (ndiffs - 1):-1:1
                PhiD[i, :] .= (PhiD[i + 1, :] .+
                               sgn[i + 1, :] .* d[i + 1, :] .* phi[i + 2, :]) .* beta[i, :]
            end

            Theta = zeros(Float64, ndiffs, npix)
            Theta[1, :] .= -theta[1, :]
            for i in 2:(ndiffs - 1)
                Theta[i, :] .= Theta[i - 1, :] .* beta[i - 1, :] .+ sgn[i, :] .* theta[i, :]
            end

            ThetaD = zeros(Float64, ndiffs + 1, npix)
            ThetaD[2, :] .= -d[1, :] .* theta[1, :]
            for i in 2:(ndiffs - 1)
                ThetaD[i + 1, :] .= beta[i - 1, :] .* ThetaD[i, :] .+
                                    sgn[i, :] .* d[i, :] .* theta[i, :]
            end

            beta_extended = ones(Float64, ndiffs, npix)
            beta_extended[2:end, :] .= beta

            # C' and B' in the paper

            theta_ndiffs = theta[ndiffs + 1, :]
            theta_ndiffs = reshape(theta_ndiffs, (1, size(theta_ndiffs, 1)))
            dC = sgn ./ theta_ndiffs .*
                 (phi[2:end, :] .* Theta .+ theta[1:(end - 1), :] .* Phi)
            dC .*= diffs2use

            dB = sgn ./ theta_ndiffs .*
                 (phi[2:end, :] .* ThetaD[2:end, :] .+ theta[1:(end - 1), :] .* PhiD)

            # {\cal A}, {\cal B}, {\cal C} in the paper
            A = 2 * sum(
                d .* sgn ./ theta_ndiffs .* beta_extended .* phi[2:end, :] .*
                ThetaD[1:(end - 1), :],
                dims = 1)
            A .+= sum(
                d .^ 2 .* theta[1:(end - 1), :] .* phi[2:end, :] ./ theta_ndiffs, dims = 1)

            B = sum(d .* dC, dims = 1)
            C = sum(dC, dims = 1)

            countrate = B ./ C
            chisq = (A .- B .^ 2 ./ C) ./ scale
            var = scale ./ C
            weights = dC ./ C

            #use first countrate measurement to improve alpha,beta definitions
            #and extract better, unbiased countrate
            final_countrates[s_ind, :] .= countrate[begin, :]
            final_vars[s_ind, :] .= var[begin, :]
            final_chisqs[s_ind, :] = chisq[begin, :]
        end
    end

    return final_countrates, final_vars .^ (-1), final_chisqs
end

"""
assume all timesteps are the same.
"""
function sutr(
        dcubedat, gainMat, readVarMat; firstind = 1, good_diffs = nothing, n_iter = 2)
    dimages = gainMat .* diff(dcubedat[:, :, firstind:end], dims = 3)
    if isnothing(good_diffs)
        good_diffs = ones(Bool, size(dimages))
    end
    ndiffs = size(dimages, 3)

    #ignore good_diffs for now TODO fix
    # this the starting guess, it will get refined to approach the maximum likelyhood estimate
    rate = mean(dimages; dims = 3)
    rate[rate .< 0] .= 0

    # working arrays allocate here and reuse inside loop
    α = ones(ndiffs)
    β = ones(ndiffs - 1)
    C = SymTridiagonal(α, β)
    # C \ b uses the LDLt factorization
    ones_vec = ones(ndiffs)

    # returned arrays
    ivars = Matrix{Float64}(undef, size(dimages)[1:2])
    chi2s = Matrix{Float64}(undef, size(dimages)[1:2])
    for pixel_ind in CartesianIndices(dimages[:, :, 1])
        β .= -readVarMat[pixel_ind]

        for _ in 1:n_iter
            α .= 2 * readVarMat[pixel_ind] + rate[pixel_ind]
            ivars[pixel_ind] = ones_vec' * (C \ ones_vec)
            @views rate[pixel_ind] = (ones_vec' * (C \ dimages[pixel_ind, :])) /
                                     ivars[pixel_ind]
        end

        residuals = dimages[pixel_ind, :] .- rate[pixel_ind]
        chi2s[pixel_ind] = residuals' * (C \ residuals)
    end
    rate, ivars, chi2s
end

"""
assume all timesteps are the same.
"""
function sutr_eigen(
        dcubedat, gainMat, readVarMat; firstind = 1, good_diffs = nothing, n_iter = 2)
    dimages = gainMat .* diff(dcubedat[:, :, firstind:end], dims = 3)
    if isnothing(good_diffs)
        good_diffs = ones(Bool, size(dimages))
    end
    ndiffs = size(dimages, 3)

    #ignore good_diffs for now TODO fix
    # this the starting guess, it will get refined to approach the maximum likelyhood estimate
    rate = mean(dimages; dims = 3)
    rate[rate .< 0] .= 0

    # working arrays
    ones_vec = ones(ndiffs)
    invΛ = Diagonal(ones(ndiffs))
    Q = @. sin((1:ndiffs) * (1:ndiffs)' * π / (ndiffs + 1))
    Q ./= sqrt.(sum(Q .^ 2, dims = 2)) #TODO do this analytically?
    invC = Q * invΛ * Q

    # returned arrays
    ivars = Matrix{Float64}(undef, size(dimages)[1:2])
    chi2s = Matrix{Float64}(undef, size(dimages)[1:2])

    for pixel_ind in CartesianIndices(dimages[:, :, 1])
        β = -readVarMat[pixel_ind]
        for _ in 1:n_iter
            α = 2 * readVarMat[pixel_ind] + rate[pixel_ind]
            @. invΛ.diag = (α + 2β * cos((1:ndiffs) * π / (ndiffs + 1)))^(-1)
            invC .= Q * invΛ * Q
            @views rate[pixel_ind] = (ones_vec' * (invC * dimages[pixel_ind, :])) /
                                     sum(invC)
        end

        ivars[pixel_ind] = sum(invC)

        residuals = dimages[pixel_ind, :] .- rate[pixel_ind]
        chi2s[pixel_ind] = residuals' * invC * residuals
    end
    rate, ivars, chi2s
end

function sutr_thomas(
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