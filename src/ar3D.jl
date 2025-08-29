# Handling the 3D data cube
using LinearAlgebra: SymTridiagonal, Diagonal, mul!
using Statistics: mean, median
using StatsBase: iqr
using TimerOutputs
using FITSIO, EllipsisNotation, AstroTime, Glob

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
    ref_zpt_vec = median(dcube[2049:end, :, :], dims = (1, 2))
    sci_zpt_vec = median(cat(dcube[1:4, :, :], dcube[2045:2048, :, :], dims = 1), dims = (1, 2))

    ref_zpt_out = dropdims(ref_zpt_vec, dims = (1, 2))
    sci_zpt_out = dropdims(sci_zpt_vec, dims = (1, 2))

    dcube[1:2048, :, :] .-= sci_zpt_vec
    dcube[2049:end, :, :] .-= ref_zpt_vec

    amp_ranges = [1:512, 513:1024, 1025:1536, 1537:2048]
    ref_bot_zpt = dropdims(
        median(cat(dcube[amp_ranges[1], 1:4, :], dcube[amp_ranges[4], 1:4, :], dims = 1),
            dims = (1, 2)),
        dims = (1, 2))
    ref_top_zpt = dropdims(
        median(
            cat(dcube[amp_ranges[1], 2045:2048, :], dcube[amp_ranges[4], 2045:2048, :], dims = 1),
            dims = (1, 2)),
        dims = (1, 2))
    amp_off_vec = zeros(4, length(ref_zpt_out))
    for (aindx, amp) in enumerate(amp_ranges)
        amp_off_vec[aindx, :] .= dropdims(
            median(
                cat(dcube[amp, 1:4, :] .- reshape(ref_bot_zpt, 1, 1, :),
                    dcube[amp, 2045:2048, :] .- reshape(ref_top_zpt, 1, 1, :),
                    dims = 2), dims = (1, 2)), dims = (1, 2))
    end
    amp_off_vec .-= reshape((amp_off_vec[1, :] .+ amp_off_vec[4, :]) ./ 2, 1, :)

    for i in eachindex(amp_ranges)
        dcube[amp_ranges[i], :, :] .-= reshape(amp_off_vec[i, :], 1, 1, :)
    end

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
    get_last_unsaturated_read(datacube, saturation_map; fudge_factor = 0.9, n_cutoff = 2)

Find pixels that are saturated in the datacube.  This uses a saturation map that specifies the 
saturation level for each pixel.

It finds the last read that is unsaturated, meaning `<= fudge_factor * saturation_map`, then
returns the index of the read before it.

For the reference array (rows 2049:end), we set all pixels to use the full number of reads
since we don't do saturation detection on them.
"""
function get_last_unsaturated_read(datacube, saturation_map; fudge_factor = 0.9, n_cutoff = 1)
    read_nums = zeros(Int, size(datacube, 1), size(datacube, 2))

    # do saturation detection on the real pixels only
    for i in 1:2048, j in 1:2048
        read_num = findlast(
            x -> x <= saturation_map[i, j] * fudge_factor,
            datacube[i, j, :])

        read_nums[i, j] = if isnothing(read_num)
            0
        else
            read_num
        end
    end

    # for the reference array, set all pixels to use the full number of reads
    read_nums[2049:end, :] .= size(datacube, 3)

    read_nums
end

function outlier_mask(dimages, last_unsaturated; clip_threshold = 20)
    @timeit "alloc" mask=ones(Bool, size(dimages))

    read_buffer = zeros(eltype(dimages), size(dimages, 3))

    # skip the reference array, just do the 2048x2048 real pixels
    @timeit "loop" for i in 1:2048, j in 1:2048
        # TODO make sure this is nonallocating. Put things in a buffer for faster calculations.
        @views read_buffer[1:last_unsaturated[i, j]] .= dimages[
            i, j, 1:last_unsaturated[i, j]]

        # if there's 2 or fewer unsaturated reads, we can't assess outliers
        if last_unsaturated[i, j] <= 2
            continue
        end

        @timeit "masked diffs view" diffs=@views read_buffer[1:last_unsaturated[i, j]]
        @timeit "mean" @views μ = mean(diffs)
        @timeit "iqr" @views σ = iqr(diffs) / 1.34896

        if σ == 0
            continue
        end

        # same here, it's easier to write performant code as a loop
        @timeit "mask" for k in 1:last_unsaturated[i, j]
            mask[i, j, k] = (read_buffer[k] - μ) < (clip_threshold * σ)
        end
    end
    mask
end

"""
Construct the difference images in place, overwriting the datacube. Returns a view of the
differences with shape (npix_x, npix_y, n_diffs) where n_diffs = n_reads - 1.
"""
function diffify_datacube!(datacube, last_unsaturated)
    # construct the differences images in place, overwriting datacube
    for i in size(datacube, 3):-1:2
        @views datacube[:, :, i] .= (datacube[:, :, i] .- datacube[:, :, i - 1])
    end
    # this view is to minimize indexing headaches
    dimages = @views datacube[:, :, 2:end]
    # modify the last_unsaturated_read accordingly.  It now refers to diffs, not reads
    last_unsaturated .-= 1
    return dimages
end

"""
    sutr_wood(dimages, gainMat, readVarMat; n_repeat = 2)

Fit the counts for each read to compute the count rate and variance for each pixel.
This assumes that all images are sequential (ie separated by the same read time).

In principle, this function is merely taking the weighted-average of the rate diffs for each pixel,
under the assumption that their uncertainties are described by a covariance matrix with diagonal
elements given by (2*read_variance + photon_noise) and off-diagonal elements given by -photon_noise.
In practice solving against this covariance matrix is a bit tricky.  This function uses the Woodbury
matrix identity to solve the system of equations.

# Arguments
- `dimages` has shape (npix_x,npix_y,n_diffs), in units of DN
- `gainMat`: The gain for each pixel (npix_x,npix_y), in units of e-/DN
- `readVarMat`: the read noise (as a variance) for each pixel (npix_x,npix_y), in units of DN/read
- `last_unsaturated`: the index of the last unsaturated read for each pixel. Has shape (npix_x,npix_y)
- `not_cosmic_ray`: a mask of which reads are not cosmic rays for each pixel. Has the same shape as dimages.

# Keyword Arguments
- n_repeat: number of iterations to run, default is 2

# Returns
A tuple of:
- `rates` is the best-fit count rate for each pixel
- `ivars` is the inverse variance describing the uncertainty in the count rate for each pixel
- `chi2s` is the chi squared value for each pixel

Written by Andrew Saydjari, based on work by Kevin McKinnon and Adam Wheeler.
Based on [Tim Brandt's SUTR python code](https://github.com/t-brandt/fitramp).
"""
function sutr_wood(dimages, gain_mat, read_var_mat, last_unsaturated, not_cosmic_ray; n_repeat = 2)
    # core Woodbury version of SUTR by Andrew Saydjari on October 17, 2024
    # based on Tim Brandt SUTR python code (https://github.com/t-brandt/fitramp)

    # initial guess for iterative flux calculation.  Median works better than mean.
    rates = dropdims(median(dimages; dims = 3), dims = 3)
    ivars = zeros(Float64, size(dimages, 1), size(dimages, 2))
    chi2s = zeros(Float64, size(dimages, 1), size(dimages, 2))

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
        # ignore cosmic rays and saturated pixels
        good_diffs = not_cosmic_ray[pixel_ind, :]
        good_diffs[max(last_unsaturated[pixel_ind] + 1, 1):end] .= false
        n_good_diffs = sum(good_diffs)

        if n_good_diffs == ndiffs # if all diffs are good, we can use a fast implementation
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
        elseif n_good_diffs == 0
            # if there are no good diffs, we can't fit anything
            rates[pixel_ind] = 0
            ivars[pixel_ind] = 0
            chi2s[pixel_ind] = 0
        else # otherwise, we need to do the slow implementation.
            # This _should_ only be needed for a small fraction of pixels.
            # Several things in this branch allocate

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
                @views r = dimages[pixel_ind, good_diffs] .- rates[pixel_ind]
                @views chi2s[pixel_ind] = r' * (C \ r)
            end
        end
    end
    return rates, ivars, chi2s # outputs in DN/read (not multiplied by gain)
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
            #            dat = (read(f[1]) .^ 2) ./ 2
            dat = (read(f[1]) .^ 2)
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

function load_saturation_maps(tel, chips; datadir = "data/saturation_maps")
    chips = replace(chips, "a" => "R", "b" => "G", "c" => "B")
    tel = lowercase(tel)
    saturationMatDict = Dict{String, Array{Float64, 2}}()
    for chip in string.(collect(chips))
        saturationMatPath = joinpath(datadir, "saturation_map_$(tel)_chip$(chip).h5")
        saturationMatDict[chip] = JLD2.load(saturationMatPath, "saturation_values")
    end
    return saturationMatDict
end

"""
Compute the mid-exposure time for an exposure a few different ways.
"""
function exposure_times(hdr_dict, firstind, lastind_loc, nread_total, tdat)
    first_image_start_time = TAIEpoch(hdr_dict[firstind]["DATE-OBS"])
    last_image_start_time = TAIEpoch(hdr_dict[lastind_loc]["DATE-OBS"])
    dtime_read = (hdr_dict[firstind]["INTOFF"] / 1000 / 3600 / 24)days #dt_read, orig in ms, convert to days
    dtime_delay = (hdr_dict[firstind]["INTDELAY"] / 3600 / 24)days #orig in seconds, convert to days
    exptime_est = (hdr_dict[firstind]["EXPTIME"] / 3600 / 24)days

    # NOT using dtime_delay because we start directly at first_image_start_time
    # (though we might need to think about this more: not sure what dtime_delay does)
    # REMEMBER to add half of a dtime_read to shift to center of exposure

    # Like DRP outputs (we think)
    mjd_mid_exposure_old = modified_julian(first_image_start_time
                                           +
                                           0.5 * exptime_est * size(tdat, 3) / nread_total)
    # Using dread_time*ndiff_USED
    mjd_mid_exposure_rough = modified_julian(first_image_start_time
                                             +
                                             dtime_read *
                                             (0.5 + 0.5 * (size(tdat, 3) - 1)))
    # Using times directly from header
    mjd_mid_exposure_precise = modified_julian(first_image_start_time + 0.5 * dtime_read
                                               +
                                               0.5 *
                                               (last_image_start_time - first_image_start_time))
    mjd_mid_exposure = mjd_mid_exposure_precise

    return mjd_mid_exposure_old, mjd_mid_exposure_rough, mjd_mid_exposure_precise, mjd_mid_exposure
end

# firstind overriden for APO dome flats
function process_3D(outdir, runname, tel, mjd, expid, chip,
        gainMatDict, readVarMatDict, saturationMatDict;
        firstind = 3, cor1fnoise = true, extractMethod = "sutr_wood",
        save3dcal = false, cluster = "sdss", suppress_warning = false, checkpoint_mode = "commit_same")
    dirName = joinpath(outdir, "apred/$(mjd)/")
    if !ispath(dirName)
        mkpath(dirName)
    end
    plotdirName = joinpath(outdir, "plots/$(mjd)/")
    if !ispath(plotdirName)
        mkpath(plotdirName)
    end

    df = read_almanac_exp_df(joinpath(outdir, "almanac/$(runname).h5"), tel, mjd)
    #        println(expid,chip,size(df.observatory),size(df.mjd),size(df.exposure_int))
    # check if chip is in the llist of chips in df.something[expid] (waiting on Andy Casey to update alamanc)
    rawpath = build_raw_path(
        df.observatory[expid], chip, df.mjd[expid], lpad(df.exposure_int[expid], 8, "0"),
        cluster = cluster, suppress_warning = suppress_warning)
    cartid = df.cartidInt[expid]

    outfname = join(
        ["ar2D", df.observatory[expid], df.mjd[expid],
            last(df.exposure_str[expid], 4), chip, df.imagetyp[expid]],
        "_")
    fname = joinpath(outdir, "apred/$(mjd)/" * outfname * ".h5")

    outfname3d = join(
        ["ar3Dcal", df.observatory[expid], df.mjd[expid],
            last(df.exposure_str[expid], 4), chip, df.imagetyp[expid]],
        "_")
    fname3d = joinpath(outdir, "apred/$(mjd)/" * outfname3d * ".h5")

    if save3dcal && check_file(fname3d, mode=checkpoint_mode) && check_file(fname, mode=checkpoint_mode)
        return fname
    end

    if check_file(fname, mode=checkpoint_mode)
        return fname
    end

    # decompress and convert apz data format to a standard 3D cube of reads
    cubedat, hdr_dict = try
        apz2cube(rawpath)
    catch
        @warn "Failed to read $(rawpath)"
        apz2cube(rawpath)
        return
    end

    nread_total = size(cubedat, 3)

    # ADD? reset anomaly fix (currently just drop first ind or two as our "fix")

    # override the firstind and extractMethod in special cases
    # we drop the first read, but might need to adjust for the few read cases (2,3,4,5)
    firstind,
    extractMethod = if ((df.imagetyp[expid] == "DOMEFLAT") &
                        (df.observatory[expid] == "apo")) # NREAD 5, and lamp gets shutoff too soon (needs to be DCS)
        2, "dcs"
    elseif ((df.imagetyp[expid] == "QUARTZFLAT") & (nread_total == 3))
        2, "dcs"
    elseif (nread_total == 3)
        #catch some weird cases (like nreads=3 with Darks)
        #but still reduce them to prevent errors later in pipeline_2d_1d
        #ULTIMATELY want to make it so these exposures are removed from runlist earlier
        2, "dcs"
    else
        firstind, extractMethod
    end

    # Might want lastind_loc as well to truncate when all reads are
    # saturated (important for calculating the exposure mid time)
    lastind_loc = size(cubedat, 3)

    tdat = @view cubedat[:, :, firstind:lastind_loc]
    ndiff_used = size(tdat, 3) - 1 # nread_used-1

    mjd_mid_exposure_old, mjd_mid_exposure_rough, mjd_mid_exposure_precise, mjd_mid_exposure = exposure_times(
        hdr_dict, firstind, lastind_loc, nread_total, tdat)

    # compute the last unsaturated read for each pixel.
    last_unsaturated = get_last_unsaturated_read(tdat, saturationMatDict[chip])

    ## zero pointing per read and for ref, sci, and relative for sci amp
    ## defer 1/f correction to 2D stage
    outdat = Float64.(tdat)
    ref_zpt_out, sci_zpt_out, amp_off_vec = zeropoint_read_dcube!(outdat)

    # ADD? reference array-based masking/correction

    # ADD? nonlinearity correction

    if save3dcal
        #make copy before it is adjusted
        outdat_cal = copy(outdat)
    end

    # extraction 3D -> 2D
    if extractMethod == "dcs"
        # TODO use last_unsaturated
        dimage, ivarimage, chisqimage = dcs(outdat, gainMatDict[chip], readVarMatDict[chip])
        CRimage = zeros(Int, size(dimage)[1:2])
    elseif extractMethod == "sutr_wood"
        dimages = diffify_datacube!(outdat, last_unsaturated)

        # try to identify any cosmic rays
        not_cosmic_ray = outlier_mask(dimages, last_unsaturated)
        CRimage = sum(.!not_cosmic_ray, dims = 3)[:, :, 1]

        # n.b. this will mutate outdat
        dimage, ivarimage, chisqimage = sutr_wood(
            dimages, gainMatDict[chip], readVarMatDict[chip], last_unsaturated, not_cosmic_ray)
    else
        error("Extraction method not recognized")
    end

    # write ar2D file
    metadata = Dict(
        "cartid" => cartid,
        "ndiff_used" => ndiff_used,
        "nread_total" => nread_total,
        "extraction_method" => extractMethod,
        "mjd_mid_exposure_old" => value(mjd_mid_exposure_old),
        "mjd_mid_exposure_rough" => value(mjd_mid_exposure_rough),
        "mjd_mid_exposure_precise" => value(mjd_mid_exposure_precise),
        "mjd_mid_exposure" => value(mjd_mid_exposure)
    )
    # need to clean up imagetyp to account for FPI versus ARCLAMP
    safe_jldsave(fname, metadata; dimage, ivarimage, chisqimage, CRimage, last_unsaturated)

    if save3dcal
        safe_jldsave(
            fname3d, metadata; dimage, ivarimage, chisqimage, CRimage, last_unsaturated,
            outdat = outdat_cal, gainMat = gainMatDict[chip], readVarMat = readVarMatDict[chip],
            ref_zpt_out = ref_zpt_out, sci_zpt_out = sci_zpt_out, amp_off_vec = amp_off_vec)
    end

    return fname
end

# come back to tuning the chi2perdofcut once more rigorously establish noise model
function process_2Dcal(fname; chi2perdofcut = 100, checkpoint_mode = "commit_same")
    sname = split(split(split(fname, "/")[end], ".h5")[1], "_")
    fnameType, tele, mjd, expnum, chip, imagetyp = sname[(end - 5):end]
    outfname = replace(fname, "ar2D" => "ar2Dcal")
    if check_file(outfname, mode = checkpoint_mode)
        return
    end

    dimage = load(fname, "dimage")
    ivarimage = load(fname, "ivarimage")
    CRimage = load(fname, "CRimage")
    chisqimage = load(fname, "chisqimage")
    last_unsaturated = load(fname, "last_unsaturated")

    metadata = read_metadata(fname)
    # n.b. every pixel does't use this many diffs. Some are truncated due to saturation
    ndiff_used = metadata["ndiff_used"]

    ### dark current subtraction
    darkRateflst = sort(glob("darkRate_$(tele)_$(chip)_*.h5", dirname(fname)))
    if length(darkRateflst) != 1
        error("I didn't just find one darkRate file for mjd $mjd, I found $(length(darkRateflst))")
    end
    darkRate = load(darkRateflst[1], "dark_rate")
    pix_bitmask = load(darkRateflst[1], "dark_pix_bitmask")

    dimage .-= darkRate
    # should I be modifying ivarimage? (uncertainty on dark rate in quad... but dark subtraction has bigger sys)

    ### flat fielding
    flatFractionflst = sort(glob("flatFraction_$(tele)_$(chip)_*.h5", dirname(fname)))
    if length(flatFractionflst) != 1
        error("I didn't just find one flatFraction file for mjd $mjd, I found $(length(flatFractionflst))")
    end
    flat_im = load(flatFractionflst[1], "flat_im")
    flat_pix_bitmask = load(flatFractionflst[1], "flat_pix_bitmask")
    dimage[5:2044, 5:2044] ./= flat_im
    ivarimage[5:2044, 5:2044] .*= flat_im .^ 2
    pix_bitmask[5:2044, 5:2044] .|= flat_pix_bitmask

    pix_bitmask .|= (CRimage .== 1) * 2^7
    pix_bitmask .|= (CRimage .> 1) * 2^8
    pix_bitmask .|= ((chisqimage ./ (last_unsaturated .- CRimage)) .> chi2perdofcut) * 2^9
    # these values are now defined for all pixels including the reference array

    pix_bitmask .|= (last_unsaturated .< ndiff_used) * 2^13
    pix_bitmask .|= (last_unsaturated .<= 0) * 2^14

    safe_jldsave(outfname, metadata; dimage, ivarimage, pix_bitmask)
    return
end
