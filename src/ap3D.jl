# Handling the 3D data cube
using LinearAlgebra: SymTridiagonal, Diagonal, mul!
using Statistics: mean
using TimerOutputs
using FITSIO

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
    dimage = (dcubedat[:, :, end] .- dcubedat[:, :, firstind])
    # bad to use measured flux as the photon noise
    ivarimage = 1 ./ (2 .* readVarMat .+ gainMat./dimage)
    # return dimage ./ ndiffs .* gainMat, (ndiffs .^ 2) ./ (gainMat.^2) .* ivarimage, zero(dimage) #output in electrons/read
    return dimage ./ ndiffs, (ndiffs .^ 2) .* ivarimage, zero(dimage) #output in DN/read
end

function outlier_mask(dimages; clip_threshold = 20)
    mask = ones(Bool, size(dimages))
    for i in axes(dimages, 1), j in axes(dimages, 2)
        @views μ = mean(dimages[i, j, :])
        @views σ = iqr(dimages[i, j, :]) / 1.34896
        @views @. mask[i, j, :] &= (dimages[i, j, :] - μ) < (clip_threshold * σ)
    end
    mask
end

"""
    sutr_wood!(datacube, gainMat, readVarMat; firstind = 1, n_repeat = 2)

Fit the counts for each read to compute the count rate and variance for each pixel.
This assumes that all images are sequential (ie separated by the same read time).

In principle, this function is merely taking the weighted-average of the rate diffs for each pixel, 
under the assumption that their uncertainties are described by a covariance matrix with diagonal 
elements given by (2*read_variance + photon_noise) and off-diagonal elements given by -photon_noise.
In practice solving against this covariance matrix is a bit tricky.  This function uses the Woodbury
matrix identity to solve the system of equations.

# Arguments
- `datacube` has shape (npix_x,npix_y,n_reads)
- `gainMat`: The gain for each pixel (npix_x,npix_y)
- `read_var_mat`: the read noise (as a variance) for each pixel (npix_x,npix_y)

# Keyword Arguments
- firstind: the index of the first read that should be used. 
- n_repeat: number of iterations to run, default is 2

# Returns
A tuple of `(rates, ivars, chi2s)` where:
- `rates` is the best-fit count rate for each pixel
- `ivars` is the inverse variance describing the uncertainty in the count rate for each pixel
- `chi2s` is the chi squared value for each pixel

!!! warning
    This mutates datacube. The difference images are written to datacute[:, :, firstindex+1:end]

Written by Andrew Saydjari, based on work by Kevin McKinnon and Adam Wheeler. 
Based on [Tim Brandt's SUTR python code](https://github.com/t-brandt/fitramp).
"""
function sutr_wood!(datacube, gainMat, readVarMat; firstind = 1, n_repeat = 2)
    # Woodbury version of SUTR by Andrew Saydjari on October 17, 2024 
    # based on Tim Brandt SUTR python code (https://github.com/t-brandt/fitramp)

    # construct the differences images in place, overwriting datacube
    for i in size(datacube, 3):-1:(firstind + 1)
        @views datacube[:, :, i] .= gainMat .*
                                    (datacube[:, :, i] .- datacube[:, :, i - 1])
    end
    # this view is to minimize indexing headaches
    dimages = view(datacube, :, :, (firstind + 1):size(datacube, 3))

    all_good_diffs = outlier_mask(dimages)
    CRimage = sum(.!all_good_diffs, dims = 3)[:, :, 1]

    rates = dropdims(mean(dimages; dims = 3), dims = 3)
    ivars = zeros(Float64, size(datacube, 1), size(datacube, 2))
    chi2s = zeros(Float64, size(datacube, 1), size(datacube, 2))

    ndiffs = size(dimages, 3)
    # working arrays
    ones_vec = ones(ndiffs)
    KinvQones = zeros(ndiffs)
    KinvQdata = zeros(ndiffs)
    Kinv = zeros(ndiffs)

    # eigenvectors of a matrix with 1 on the diagonal, and -2 on the off-diagonals
    Q = @. sin((1:ndiffs) * (1:ndiffs)' * π / (ndiffs + 1))
    Q .*= sqrt(2 / (ndiffs + 1))
    Qones = Q * ones_vec
    # eigenvalues of that matrix
    D = (1 .- 4 * cos.((1:ndiffs) * π / (ndiffs + 1)))
    Qdata = similar(Qones)
    #TODO delete
    d1s = sum(dimages, dims = 3)
    d2s = sum(abs2, dimages, dims = 3)

    for pixel_ind in CartesianIndices(view(dimages, :, :, 1))
        good_diffs = all_good_diffs[pixel_ind, :]
        n_good_diffs = sum(good_diffs)
        if n_good_diffs == ndiffs
            read_var = readVarMat[pixel_ind] * (gainMat[pixel_ind].^2)
            if isnan(read_var)
                println((pixel_ind,read_var))
            end
            @views mul!(Qdata, Q, view(dimages, pixel_ind, :))
            d1 = d1s[pixel_ind, 1]
            d2 = d2s[pixel_ind, 1]

            Qones = Q * ones_vec

            for _ in 1:n_repeat
                rate_guess = rates[pixel_ind] > 0 ? rates[pixel_ind] : 0
                x = (rate_guess + 1.5read_var)
                y = 2 * x / read_var

                Kinv .= D ./ (D .+ y)
                KinvQones .= Kinv .* Qones
                KinvQdata .= Kinv .* Qdata

                ivars[pixel_ind] = (ndiffs - Qones' * KinvQones) / x
                rates[pixel_ind] = (d1 - Qones' * KinvQdata) / x / ivars[pixel_ind]
                chi2s[pixel_ind] = (d2 - Qdata' * KinvQdata) / x -
                                   rates[pixel_ind]^2 * ivars[pixel_ind]
                if isnan(rates[pixel_ind])
                    println((pixel_ind,x,y,rate_guess,read_var))
                end
            end
        else
            # the appendix of https://github.com/dfink/gspice/blob/main/paper/gspice.pdf
            # has a nice derivation of how to invert a matrix with masked elements
            # this implementation is super naive

            for _ in 1:n_repeat
                rate_guess = rates[pixel_ind] > 0 ? rates[pixel_ind] : 0
                read_var = readVarMat[pixel_ind] * (gainMat[pixel_ind].^2)
                @views C = SymTridiagonal(
                    (rate_guess + 2read_var) .* ones_vec, -read_var .* ones_vec[1:(end - 1)])[good_diffs, good_diffs]

                @views ivars[pixel_ind] = ones_vec[1:n_good_diffs]' * (C \
                                                                       ones_vec[1:n_good_diffs])
                @views rates[pixel_ind] = (ones_vec[1:n_good_diffs]' * (C \
                                            dimages[pixel_ind, good_diffs])) /
                                          ivars[pixel_ind]
                @views chi2s[pixel_ind] = (dimages[pixel_ind, good_diffs]' * (C \
                                            dimages[pixel_ind, good_diffs])) /
                                          ivars[pixel_ind] -
                                          rates[pixel_ind]^2 * ivars[pixel_ind]
            end
        end
    end
    # return rates ./ ndiffs, (ndiffs .^ 2) .* ivars, chi2s, CRimage # outputs in electrons/read
    return rates ./ ndiffs ./ gainMat, (ndiffs .^ 2) .* (gainMat.^2) .* ivars, chi2s, CRimage # outputs in DN/read
end

function load_gain_maps(gainReadCalDir,tele,chips)
    gainMatDict = Dict{String, Array{Float64, 2}}()
    for chip in string.(collect(chips))
        gainMatPath = gainReadCalDir*"gain_"*tele*"_"*chip*".fits"
        if isfile(gainMatPath)
            f = FITS(gainMatPath)
            dat = read(f[1])
            close(f)
            gainView = nanzeromedian(dat) .*ones(Float64, 2560, 2048)
            view(gainView, 5:2044, 5:2044) .= dat
            gainMatDict[chip] = gainView
        else
            #once we have the LCO calibrations, we should make this warning a flag that propagates and a harder error 
            warn("Gain calibration file not found for chip $chip")
            gainMatDict[chip] = 1.9 * ones(Float64, 2560, 2048) # electrons/DN
        end
    end
    return gainMatDict
end

function load_read_var_maps(gainReadCalDir,tele,chips)
    readVarMatDict = Dict{String, Array{Float64, 2}}()
    for chip in string.(collect(chips))
        readVarMatPath = gainReadCalDir*"rdnoise_"*tele*"_"*chip*".fits"
        if isfile(readVarMatPath)
            f = FITS(readVarMatPath)
            dat = read(f[1]).^2
            close(f)
            readVarView = nanzeromedian(dat) .*ones(Float64, 2560, 2048)
            view(readVarView, 5:2044, 5:2044) .= dat
            readVarView[isnanorzero.(readVarView)] .= 25
            readVarMatDict[chip] = readVarView
        else
            warn("Read noise calibration file not found for chip $chip")
            readVarMatDict[chip] = 25 * ones(Float64, 2560, 2048) # DN/read
        end
    end
    return readVarMatDict
end