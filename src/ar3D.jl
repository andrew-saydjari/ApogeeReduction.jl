# Handling the 3D data cube
using LinearAlgebra: SymTridiagonal, Diagonal, mul!
using Statistics: mean
using TimerOutputs
using FITSIO, EllipsisNotation

function apz2cube(fname)
    f = FITS(fname)
    #    hdr_dict = Dict("ave" => read_header(f[2]))
    hdr_dict = Dict()
    hdr_dict[1] = read_header(f[3])
    avg_dcounts = read(f[2])
    cubedat = zeros(Float32, size(avg_dcounts)..., length(f) - 2) #use float bc of NaN
    cubedat[:, :, 1] .= read(f[3])
    for i in 2:(length(f) - 2)
        hdr_dict[i] = read_header(f[i + 2])
        cubedat[:, :, i] .= read(f[i + 2]) .+ avg_dcounts .+ cubedat[:, :, i - 1]
    end
    close(f)
    return cubedat, hdr_dict
end

function zeropoint_read_dcube!(dcube)
    ref_zpt_vec = mean(dcube[2049:end, :, :], dims = (1, 2))
    sci_zpt_vec = mean(dcube[1:4, :, :], dims = (1, 2)) +
                  mean(dcube[(end - 3):end, :, :], dims = (1, 2))

    ref_zpt_out = dropdims(ref_zpt_vec, dims = (1, 2))
    sci_zpt_out = dropdims(sci_zpt_vec, dims = (1, 2))

    dcube[1:2048, :, :] .-= sci_zpt_vec
    dcube[2049:end, :, :] .-= ref_zpt_vec

    amp_ranges = [1:512, 513:1024, 1025:1536, 1537:2048]
    ref_bot = zeros(4, length(ref_zpt_out))
    for (aindx, amp) in enumerate(amp_ranges)
        ref_bot[aindx, :] .= dropdims(mean(dcube[amp, 1:4, :], dims = (1, 2)), dims = (1, 2))
    end
    ref_top = zeros(4, length(ref_zpt_out))
    for (aindx, amp) in enumerate(amp_ranges)
        ref_top[aindx, :] .= dropdims(mean(dcube[amp, 2045:2048, :], dims = (1, 2)), dims = (1, 2))
    end

    amp_off_vec = ref_bot .- reshape((ref_bot[1, :] .+ ref_bot[4, :]) ./ 2, 1, :)
    amp_off_vec .+= ref_top .- .-reshape((ref_top[1, :] .+ ref_bot[4, :]) ./ 2, 1, :)
    amp_off_vec ./= 2
    amp_off_vec .-= reshape((amp_off_vec[1, :] .+ amp_off_vec[4, :]) ./ 2, 1, :)

    dcube[1:512, :, :] .-= reshape(amp_off_vec[1, :], 1, 1, :)
    dcube[513:1024, :, :] .-= reshape(amp_off_vec[2, :], 1, 1, :)
    dcube[1025:1536, :, :] .-= reshape(amp_off_vec[3, :], 1, 1, :)
    dcube[1537:2048, :, :] .-= reshape(amp_off_vec[4, :], 1, 1, :)

    return ref_zpt_out, sci_zpt_out, amp_off_vec
end

function dcs(dcubedat, gainMat, readVarMat; firstind = 1)
    ndiffs = size(dcubedat, 3) - firstind
    dimage = (dcubedat[:, :, end] .- dcubedat[:, :, firstind])
    # bad to use measured flux as the photon noise
    ivarimage = 1 ./ (2 .* readVarMat .+ gainMat ./ dimage)
    ivarimage[ivarimage .< 0] .= 1e-9 # set weight to zero for pixels with negative DCS IVAR

    # return dimage ./ ndiffs .* gainMat, (ndiffs .^ 2) ./ (gainMat.^2) .* ivarimage, zero(dimage) #output in electrons/read
    return dimage ./ ndiffs, (ndiffs .^ 2) .* ivarimage, zero(dimage), zeros(Int, size(dimage)) #output in DN/read
end

"""
    saturation_readmask(datacube, saturation_map; fudge_factor = 0.9, n_cutoff = 2)

Find pixels that are saturated in the datacube and return a mask of which reads to keep for each
pixel.  This uses a saturation map that specifies the saturation level for each pixel.

It finds the first read that is saturated, meaning `>= fudge_factor * saturation_map`, then
flags all reads starting from the `n_cutoff` before that.
"""
function saturation_readmask(datacube, saturation_map; fudge_factor = 0.9, n_cutoff = 2)
    readmask = ones(Bool, size(datacube))
    for i in axes(datacube, 1), j in axes(datacube, 2)
        @views first_saturated_read = findfirst(
            x -> x > saturation_map[i, j] * fudge_factor,
            datacube[i, j, :])
        if !isnothing(first_saturated_read)
            first_to_drop = max(1, first_saturated_read - n_cutoff)
            readmask[i, j, first_to_drop:end] .= false
        end
    end
    readmask
end

function outlier_mask(dimages, sat_readmask; clip_threshold = 20)
    @timeit "alloc" mask=copy(sat_readmask)

    read_buffer = zeros(eltype(dimages), size(dimages, 3))

    @timeit "loop" for i in axes(dimages, 1), j in axes(dimages, 2)
        # it's easier to write a nonallocating loop that some combo of views and broadcasting
        # this loop puts all the non-saturated reads into read_buffer
        # and counts them in n
        n = 1
        for k in axes(dimages, 3)
            if sat_readmask[i, j, k]
                read_buffer[k] = dimages[i, j, k]
                n += 1
            end
        end

        # if there's 3 of fewer unsaturated reads, mask all b/c we can't assess outliers
        if n <= 4
            mask[i, j, :] .= false
            continue
        end

        # TODO
        #smallest = minimum(read_buffer)
        #largest = maximum(read_buffer)
        ## efficiently find the second smallest and second largest values
        ## these are allocating-ish
        #second_smallest = partialsortperm(read_buffer, 2)
        #second_largest = partialsortperm(read_buffer, n - 1)

        @timeit "masked diffs view" diffs=@views read_buffer[1:(n - 1)]
        @timeit "mean" @views μ = mean(diffs)
        @timeit "iqr" @views σ = iqr(diffs) / 1.34896

        # same here, it's easier to write performant code as a loop
        @timeit "mask" for k in axes(dimages, 3)
            mask[i, j, k] = (read_buffer[k] - μ) < (clip_threshold * σ)
            mask[i, j, k] &= sat_readmask[i, j, k]
        end
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
- `gain_mat`: The gain for each pixel (npix_x,npix_y)
- `read_var_mat`: the read noise (as a variance) for each pixel (npix_x,npix_y)
- `sat_mat`: the saturation level for each pixel (npix_x,npix_y)

# Keyword Arguments
- firstind: the index of the first read that should be used.
- n_repeat: number of iterations to run, default is 2

# Returns
A tuple of:
- `rates` is the best-fit count rate for each pixel
- `ivars` is the inverse variance describing the uncertainty in the count rate for each pixel
- `chi2s` is the chi squared value for each pixel
- `CRimage` is the cosmic ray count image
- `sat_mask` is the mask of saturated pixels

!!! warning
    This mutates datacube. The difference images are written to datacute[:, :, firstindex+1:end]

Written by Andrew Saydjari, based on work by Kevin McKinnon and Adam Wheeler.
Based on [Tim Brandt's SUTR python code](https://github.com/t-brandt/fitramp).
"""
function sutr_wood!(datacube, gain_mat, read_var_mat, sat_mat; firstind = 1, n_repeat = 2)
    # Woodbury version of SUTR by Andrew Saydjari on October 17, 2024
    # based on Tim Brandt SUTR python code (https://github.com/t-brandt/fitramp)

    # identify saturated pixels
    @timeit "new sat mask" sat_readmask=saturation_readmask(datacube, sat_mat) # TODO use
    @show mean(sat_readmask)
    sat_mask = ones(Bool, size(datacube)[1:2]) # TODO remove

    # construct the differences images in place, overwriting datacube
    for i in size(datacube, 3):-1:(firstind + 1)
        @views datacube[:, :, i] .= (datacube[:, :, i] .- datacube[:, :, i - 1])
    end
    @timeit "setup views" begin
        # this view is to minimize indexing headaches
        dimages = @views datacube[:, :, (firstind + 1):end]
        sat_readmask = @views sat_readmask[:, :, (firstind + 1):end]
    end

    @timeit "CR mask" not_cosmic_ray=outlier_mask(dimages, sat_readmask)
    # don't try to do CR rejection on saturated pixels
    not_cosmic_ray[sat_mask, :] .= true
    CRimage = sum(.!not_cosmic_ray, dims = 3)[:, :, 1]

    rates = dropdims(mean(dimages; dims = 3), dims = 3)
    ivars = zeros(Float64, size(datacube, 1), size(datacube, 2))
    chi2s = zeros(Float64, size(datacube, 1), size(datacube, 2))

    @timeit "allocs" begin
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
        d1s = sum(dimages, dims = 3)
        d2s = sum(abs2, dimages, dims = 3)
    end

    @timeit "fit" for pixel_ind in CartesianIndices(view(dimages, :, :, 1))
        good_diffs = not_cosmic_ray[pixel_ind, :]
        n_good_diffs = sum(good_diffs)
        if n_good_diffs == ndiffs
            read_var = read_var_mat[pixel_ind]
            @views mul!(Qdata, Q, view(dimages, pixel_ind, :))
            d1 = d1s[pixel_ind, 1]
            d2 = d2s[pixel_ind, 1]

            Qones = Q * ones_vec

            for _ in 1:n_repeat
                rate_guess = rates[pixel_ind] > 0 ? rates[pixel_ind] / gain_mat[pixel_ind] : 0
                x = (rate_guess + 1.5read_var)
                y = 2 * x / read_var

                Kinv .= D ./ (D .+ y)
                KinvQones .= Kinv .* Qones
                KinvQdata .= Kinv .* Qdata

                ivars[pixel_ind] = (ndiffs - Qones' * KinvQones) / x
                rates[pixel_ind] = (d1 - Qones' * KinvQdata) / x / ivars[pixel_ind]
                chi2s[pixel_ind] = (d2 - Qdata' * KinvQdata) / x -
                                   rates[pixel_ind]^2 * ivars[pixel_ind]
            end
        else
            # the appendix of https://github.com/dfink/gspice/blob/main/paper/gspice.pdf
            # has a nice derivation of how to invert a matrix with masked elements
            # this implementation is super naive

            for _ in 1:n_repeat
                rate_guess = rates[pixel_ind] > 0 ? rates[pixel_ind] / gain_mat[pixel_ind] : 0
                read_var = read_var_mat[pixel_ind]
                @views C = SymTridiagonal(
                    (rate_guess + 2read_var) .* ones_vec, -read_var .* ones_vec[1:(end - 1)]
                )[good_diffs, good_diffs]

                @views ivars[pixel_ind] = ones_vec[1:n_good_diffs]' * (C \ ones_vec[1:n_good_diffs])
                @views rates[pixel_ind] = (ones_vec[1:n_good_diffs]' *
                                           (C \ dimages[pixel_ind, good_diffs])) / ivars[pixel_ind]
                @views chi2s[pixel_ind] = (dimages[pixel_ind, good_diffs]' *
                                           (C \ dimages[pixel_ind, good_diffs])) /
                                          ivars[pixel_ind] - rates[pixel_ind]^2 * ivars[pixel_ind]
            end
        end
    end
    # return rates .* gainMat, ivars ./ (gainMat .^ 2), chi2s, CRimage # outputs in electrons/read
    return rates, ivars, chi2s, CRimage, sat_readmask # outputs in DN/read
end

function load_gain_maps(gainReadCalDir, tele, chips)
    gainMatDict = Dict{String, Array{Float64, 2}}()
    for chip in string.(collect(chips))
        gainMatPath = gainReadCalDir * "gain_" * tele * "_" * chip * ".fits"
        if isfile(gainMatPath)
            f = FITS(gainMatPath)
            dat = read(f[1])
            close(f)
            gainView = nanzeromedian(dat) .* ones(Float64, 2560, 2048)
            view(gainView, 5:2044, 5:2044) .= dat
            gainMatDict[chip] = gainView
        else
            #once we have the LCO calibrations, we should make this warning a flag that propagates and a harder error
            @warn "Gain calibration file not found for chip $chip at $gainMatPath"
            gainMatDict[chip] = 1.9 * ones(Float64, 2560, 2048) # electrons/DN
        end
    end
    return gainMatDict
end

function load_read_var_maps(gainReadCalDir, tele, chips)
    readVarMatDict = Dict{String, Array{Float64, 2}}()
    for chip in string.(collect(chips))
        readVarMatPath = gainReadCalDir * "rdnoise_" * tele * "_" * chip * ".fits"
        if isfile(readVarMatPath)
            f = FITS(readVarMatPath)
            # TODO: Tim and I need to sort out this factor of 2
            dat = (read(f[1]) .^ 2) ./ 2
            close(f)
            refval = nanzeromedian(dat)
            readVarView = refval .* ones(Float64, 2560, 2048)
            view(readVarView, 5:2044, 5:2044) .= dat
            readVarView[isnanorzero.(readVarView)] .= refval
            readVarView[readVarView .== 1] .= refval
            readVarMatDict[chip] = readVarView
        else
            @warn "Read noise calibration file not found for chip $chip"
            readVarMatDict[chip] = 25 * ones(Float64, 2560, 2048) # DN/read
        end
    end
    return readVarMatDict
end
