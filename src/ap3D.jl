# Handling the 3D data cube
using LinearAlgebra: SymTridiagonal, Diagonal, mul!
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
Use dumb sigma clipping to remove really obvious cosmic rays.
"""
function outlier_mask(dimages; clip_threshold = 20)
    mask = ones(Bool, size(dimages))
    for i in axes(dimages, 1), j in axes(dimages, 2)
        μ = mean(dimages[i, j, :])
        σ = nanzeroiqr(dimages[i, j, :])
        @. mask[i, j, :] &= (dimages[i, j, :] - μ) < (clip_threshold * σ)
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

Based partially on Tim Brandt's SUTR python code (https://github.com/t-brandt/fitramp)

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
    # construct the differences images in place, overwriting datacube
    for i in size(datacube, 3):-1:(firstind + 1)
        @views datacube[:, :, i] .= gainMat .*
                                    (datacube[:, :, i] .- datacube[:, :, i - 1])
    end
    # this view is to minimize indexing headaches
    dimages = view(datacube, :, :, (firstind + 1):size(datacube, 3))

    all_good_diffs = outlier_mask(dimages)
    CRimage = sum(.!all_good_diffs, dims = 3)[:, :, 1]
    let n_masked_pix = sum(CRimage .> 0), n_masked_reads = sum(.!all_good_diffs)
        @info "Masked $n_masked_reads reads ($n_masked_pix pixels)"
    end

    rates = dropdims(mean(dimages; dims = 3), dims = 3)
    ivars = zeros(Float64, size(datacube, 1), size(datacube, 2))
    chi2s = zeros(Float64, size(datacube, 1), size(datacube, 2))

    ndiffs = size(dimages, 3)

    # working arrays
    n_good_diffs = ndiffs .- CRimage
    ones_vec = ones(ndiffs)
    KinvQones = zeros(ndiffs)
    KinvQdata = zeros(ndiffs)
    Kinv = zeros(ndiffs)

    # eigenvectors of a matrix with 1 on the diagonal, and -2 on the off-diagonals
    Q = @. sin((1:ndiffs) * (1:ndiffs)' * π / (ndiffs + 1))
    Q .*= sqrt(2 / (ndiffs + 1))
    Qones_unmasked = Q * ones_vec
    # eigenvalues of that matrix
    D = (1 .- 4 * cos.((1:ndiffs) * π / (ndiffs + 1)))
    Qdata_unmasked = similar(Qones_unmasked)

    for pixel_ind in CartesianIndices(view(dimages, :, :, 1))
        read_var = readVarMat[pixel_ind]
        @views good_diffs = all_good_diffs[pixel_ind, :]

        Qones = view(Qones_unmasked, good_diffs)
        Qdata = view(Qdata_unmasked, good_diffs)
        mul!(Qdata, view(Q, good_diffs, good_diffs), view(dimages, pixel_ind, good_diffs))
        # these are allocating
        @views d1 = sum(dimages[pixel_ind, good_diffs])
        @views d2 = sum(abs2, dimages[pixel_ind, good_diffs])

        for _ in 1:n_repeat
            a = rates[pixel_ind] > 0 ? rates[pixel_ind] : 0
            x = (a + 1.5read_var)
            y = 2 * x / read_var
            @views Kinv[good_diffs] .= D[good_diffs] ./ (D[good_diffs] .+ y)
            view(KinvQones, good_diffs) .= view(Kinv, good_diffs) .* Qones
            view(KinvQdata, good_diffs) .= view(Kinv, good_diffs) .* Qdata
            ivars[pixel_ind] = (n_good_diffs[pixel_ind] -
                                Qones' * view(KinvQones, good_diffs)) / x
            rates[pixel_ind] = (d1 - Qones' * view(KinvQdata, good_diffs)) / x /
                               ivars[pixel_ind]
            chi2s[pixel_ind] = (d2 - Qdata' * view(KinvQdata, good_diffs)) / x -
                               rates[pixel_ind]^2 * ivars[pixel_ind]
        end
    end
    return rates, ivars, chi2s, CRimage
end
