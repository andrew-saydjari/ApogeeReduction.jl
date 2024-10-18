# Handling the 3D data cube
using LinearAlgebra: SymTridiagonal, Diagonal
using Statistics: mean
using TimerOutputs

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
    return dimage ./ ndiffs, (ndiffs .^ 2) .* ivarimage, nothing # no chisq, just mirroring sutr_tb!
end

"""
    sutr_tb!(datacube, gainMat, readVarMat; firstind = 1, good_diffs = nothing, n_repeat = 2)

# Arguments
- `datacube` has shape (npix_x,npix_y,n_reads)
- `gainMat`: The gain for each pixel (npix_x,npix_y)
- `read_var_mat`: the read noise (as a variance) for each pixel (npix_x,npix_y)

# Keyword Arguments
- firstind 
- good_diffs, if provided, should contain booleans and have shape (npix_x, npix_y, n_reads-1). If
  not provided, all diffs will be used.

!!! warning
    Andrew Saydjari reports a bug in this code's error bars and chi2 values 2024_10_17.
    This mutates datacube. The difference images are written to datacute[:, :, firstindex+1:end]

!!! note
    This assumes that all images are sequential (ie separated by the same read time).

Created by Kevin McKinnon and sped up by Adam Wheeler. 
Based on [Tim Brandt's SUTR python code](https://github.com/t-brandt/fitramp).
"""
function sutr_tb!(
        datacube, gainMat, readVarMat; firstind = 1, good_diffs = nothing, n_repeat = 2)
    # First version by Kevin McKinnon on Aug 20, 2024 
    # based on Tim Brandt SUTR python code (https://github.com/t-brandt/fitramp)
    # datacube has shape (npix_x,npix_y,n_reads)
    # read_var_mat has shape (npix_x,npix_y)
    # good_diffs is boolean and has shape (npix_x,npix_y,n_reads-1)

    # assumes all images are sequential (ie separated by one read time)

    # construct the differences images in place, overwriting datacube
    for i in size(datacube, 3):-1:(firstind + 1)
        @views datacube[:, :, i] .= gainMat .* (datacube[:, :, i] .- datacube[:, :, i - 1])
    end
    # this view is to minimize indexing headaches
    dimages = view(datacube, :, :, (firstind + 1):size(datacube, 3))

    if isnothing(good_diffs)
        good_diffs = ones(Bool, size(dimages))
    end

    rates = zeros(Float64, size(datacube, 1), size(datacube, 2))
    final_vars = zeros(Float64, size(datacube, 1), size(datacube, 2))
    final_chisqs = zeros(Float64, size(datacube, 1), size(datacube, 2))

    npix, ndiffs = size(dimages)[[1, 3]]

    diffs = zeros(Float64, npix, ndiffs)
    read_var = zeros(Float64, npix)
    diffs2use = zeros(Bool, npix, ndiffs)
    d = zeros(Float64, npix, ndiffs)

    theta = ones(Float64, npix, ndiffs + 1)
    phi = ones(Float64, npix, ndiffs + 1)
    Phi = zeros(Float64, npix, ndiffs)
    #PhiD = zeros(Float64, npix, ndiffs)
    Theta = zeros(Float64, npix, ndiffs)
    ThetaD = zeros(Float64, npix, ndiffs + 1)
    scale = Vector{Float64}(undef, npix)
    beta = Matrix{Float64}(undef, npix, ndiffs - 1)
    beta_extended = ones(npix, ndiffs)
    dC = zeros(Float64, npix, ndiffs)
    A = Vector{Float64}(undef, npix)
    B = Vector{Float64}(undef, npix)
    C = Vector{Float64}(undef, npix)

    # define this one here
    sgn = ones(Float64, npix, ndiffs)
    sgn[:, 1:2:end] .= -1

    #slice the datacube to analyze each row sequentially to keep runtime down
    # THIS LOOP IS ENTIRELY NONALLOCATING SO EDIT WITH CARE
    for c_ind in axes(dimages, 2)
        # these involve more copying than is really necessary
        @views diffs .= dimages[:, c_ind, :]
        @views read_var .= readVarMat[:, c_ind]
        @views diffs2use .= good_diffs[:, c_ind, :] # shape = (ndiffs, npix)
        d .= (diffs .* diffs2use)

        # this is done in a loop to be nonallocating
        for i in axes(rates, 1)
            @views rates[i, c_ind] = sum(d[i, :]) / sum(diffs2use[i, :])
        end

        for _ in 1:n_repeat
            # scale is the same at alpha in the paper, but we divide it out each Elements
            # in the covariance matrix for stability
            # The uncertainty and chi squared value will need to be scaled back later.
            @views scale .= (rates[:, c_ind] .* (rates[:, c_ind] .> 0)) .+ 2 .* read_var
            @views beta .= -1 .* read_var # the explicit -1 makes it fully broadcast
            beta ./= scale

            # Mask resultant differences that should be ignored.  This is half
            # of what we need to do to mask these resultant differences; the
            # rest comes later.
            @views beta .*= (diffs2use[:, 2:end] .* diffs2use[:, 1:(end - 1)])

            @views begin # recurence relations
                # All definitions and formulas here are in the paper.
                theta[:, 1:2] .= 1
                for i in 3:(ndiffs + 1)
                    theta[:, i] .= theta[:, i - 1] .- beta[:, i - 2] .^ 2 .* theta[:, i - 2]
                end

                phi[:, ndiffs] .= 1
                for i in (ndiffs - 1):-1:1
                    phi[:, i] .= phi[:, i + 1] .- beta[:, i] .^ 2 .* phi[:, i + 2]
                end

                for i in (ndiffs - 1):-1:1
                    Phi[:, i] .= Phi[:, i + 1] .* beta[:, i] .+
                                 sgn[:, i + 1] .* beta[:, i] .* phi[:, i + 2]
                end

                # This one is defined later in the paper and is used for jump
                # detection and pedestal fitting.
                #for i in (ndiffs - 1):-1:1
                #    PhiD[:, i] .= (PhiD[:, i + 1] .+
                #                   sgn[:, i + 1] .* d[:, i + 1] .* phi[:, i + 2]) .*
                #                  beta[:, i]
                #end

                Theta[:, 1] .= -1 .* theta[:, 1] # the explicit -1 makes it fully broadcast
                for i in 2:(ndiffs - 1)
                    Theta[:, i] .= Theta[:, i - 1] .* beta[:, i - 1] .+
                                   sgn[:, i] .* theta[:, i]
                end

                # the explicit -1 makes it fully broadcast
                ThetaD[:, 2] .= -1 .* d[:, 1] .* theta[:, 1]
                for i in 2:(ndiffs - 1)
                    ThetaD[:, i + 1] .= beta[:, i - 1] .* ThetaD[:, i] .+
                                        sgn[:, i] .* d[:, i] .* theta[:, i]
                end
            end

            beta_extended[:, 2:end] .= beta

            # C' in the paper
            @views dC .= sgn ./ theta[:, ndiffs + 1] .*
                         (phi[:, 2:end] .* Theta .+ theta[:, 1:(end - 1)] .* Phi) .*
                         diffs2use

            # TODO?
            #dB = sgn ./ theta_ndiffs .*
            #     (phi[2:end, :] .* ThetaD[2:end, :] .+ theta[1:(end - 1), :] .* PhiD)

            @views begin # {\cal A}, {\cal B}, {\cal C} in the paper
                A .= 0
                B .= 0
                C .= 0
                for i in 1:ndiffs
                    A .+= 2 .*
                          (d[:, i] .* sgn[:, i] .* beta_extended[:, i] .* phi[:, i + 1] .*
                           ThetaD[:, i] .+
                           d[:, i] .^ 2 .* theta[:, i] .* phi[:, i + 1]) ./
                          theta[:, ndiffs + 1]
                    B .+= d[:, i] .* dC[:, i]
                    C .+= dC[:, i]
                end
            end

            #use first countrate measurement to improve alpha,beta definitions
            #and extract better, unbiased countrate
            rates[:, c_ind] .= B ./ C
            final_vars[:, c_ind] .= scale ./ C
            final_chisqs[:, c_ind] .= (A .- B .^ 2 ./ C) ./ scale
        end
    end

    return rates, final_vars .^ (-1), final_chisqs
end

"""
    sutr_wood!(datacube, gainMat, readVarMat; firstind = 1, n_repeat = 2)

# Arguments
- `datacube` has shape (npix_x,npix_y,n_reads)
- `gainMat`: The gain for each pixel (npix_x,npix_y)
- `read_var_mat`: the read noise (as a variance) for each pixel (npix_x,npix_y)

# Keyword Arguments
- firstind 
- n_repeat: number of iterations to run, default is 2

!!! warning
    This mutates datacube. The difference images are written to datacute[:, :, firstindex+1:end]

!!! note
    This assumes that all images are sequential (ie separated by the same read time).

Written by Andrew Saydjari, based on work by Kevin McKinnon and Adam Wheeler. 
Based on [Tim Brandt's SUTR python code](https://github.com/t-brandt/fitramp).
"""
function sutr_wood!(datacube, gainMat, readVarMat; firstind = 1, n_repeat = 2)
    # Woodbury version of SUTR by Andrew Saydjari on October 17, 2024 
    # based on Tim Brandt SUTR python code (https://github.com/t-brandt/fitramp)
    # datacube has shape (npix_x,npix_y,n_reads)
    # read_var_mat has shape (npix_x,npix_y)
    # assumes all images are sequential (ie separated by one read time)

    # construct the differences images in place, overwriting datacube
    for i in size(datacube, 3):-1:(firstind + 1)
        @views datacube[:, :, i] .= gainMat .* (datacube[:, :, i] .- datacube[:, :, i - 1])
    end
    # this view is to minimize indexing headaches
    dimages = view(datacube, :, :, (firstind + 1):size(datacube, 3))

    rates = dropdims(mean(dimages; dims = 3), dims=3)
    ivars = zeros(Float64, size(datacube, 1), size(datacube, 2))
    chi2s = zeros(Float64, size(datacube, 1), size(datacube, 2))

    npix, ndiffs = size(dimages)[[1, 3]]
    # working arrays
    ones_vec = ones(ndiffs)
    KinvQones = zeros(ndiffs)
    KinvQdata = zeros(ndiffs)
    Kinv = zeros(ndiffs)

    Q = @. sin((1:ndiffs) * (1:ndiffs)' * π / (ndiffs + 1))
    Q .*= sqrt(2 / (ndiffs + 1))
    Qones = Q * ones_vec
    D = (1 .- 4 * cos.((1:ndiffs) * π / (ndiffs + 1)))

    for pixel_ind in CartesianIndices(dimages[:, :, 1])
        n = readVarMat[pixel_ind]
        @views Qdata = Q * dimages[pixel_ind, :]
        @views d1 = sum(dimages[pixel_ind, :])
        @views d2 = sum(abs2, dimages[pixel_ind, :])
        for _ in 1:n_repeat
            a = rates[pixel_ind] > 0 ? rates[pixel_ind] : 0
            x = (a + 1.5n)
            y = 2 * x / n
            Kinv .= D ./ (D .+ y)
            # is sum(.*) faster that xTx? if so, switch
            KinvQones .= Kinv .* Qones
            KinvQdata .= Kinv .* Qdata
            ivars[pixel_ind] = (ndiffs - Qones' * KinvQones) / x
            rates[pixel_ind] = (d1 - Qones' * KinvQdata) / x / ivars[pixel_ind]
            chi2s[pixel_ind] = (d2 - Qdata' * KinvQdata) / x -
                               rates[pixel_ind]^2 * ivars[pixel_ind]
        end
    end
    return rates, ivars, chi2s
end
