using FastRunningMedian: running_median
using Distributions: cdf, Normal
using Interpolations: linear_interpolation, Line
using DataFrames

# this file contains the code needed to extract a 1D spectrum from a 2D images.
# trace_params is of size (n_x_pix, n_fibers, 3)
# the elements correspodn to flux, y-pixel, and gaussian witdth (sigma)

# hold off on prop ivar through until we switch to sutr_wood, also could implement a chi2 cut here
# add a condition that we should drop any x pixel where a bad bit in any of the pixels being summed is bad

function get_relFlux(fname; sig_cut = 4.5, rel_val_cut = 0.07)
    f = jldopen(fname)
    flux_1d = f["flux_1d"]
    mask_1d = f["mask_1d"]
    mask_1d_good = (mask_1d .& bad_pix_bits) .== 0
    close(f)
    metadata = read_metadata(fname)

    absthrpt = dropdims(nanzeromedian(flux_1d, 1), dims = 1)
    bitmsk_relthrpt = zeros(Int, length(absthrpt))
    relthrpt = copy(absthrpt)
    relthrpt ./= nanzeromedian(relthrpt)

    thresh = (1 .- sig_cut * nanzeroiqr(relthrpt))
    bitmsk_relthrpt[relthrpt .< thresh] .|= 2^0
    bitmsk_relthrpt[relthrpt .< rel_val_cut] .|= 2^1
    return absthrpt, relthrpt, bitmsk_relthrpt, metadata
end

"""
Regularize the trace by applying a running median filter to each param in each fiber.
Could be denoised further by fitting a low-order polynomial or similar.
"""
function regularize_trace(trace_params; window_size = 101)
    @assert isodd(window_size) # otherwise the length of the regularized array is incorrect
    n_xpix, n_fibers = size(trace_params)[1:2]
    regularized_trace = similar(trace_params)
    for fiber in 1:n_fibers, param in 1:3
        regularized_trace[:,
        fiber,
        param] = running_median(
            trace_params[:, fiber, param], window_size, :asym_trunc; nan = :ignore)
    end
    regularized_trace
end

function extract_boxcar(dimage, ivarimage, pix_bitmask, trace_params;
        boxcar_halfwidth = 2, return_resids = false)
    if return_resids
        resid_fluxes_2d = zeros(Float64, size(dimage))
        resid_ivars_2d = zeros(Float64, size(dimage))
    end
    flux_1d = extract_boxcar_core(dimage, trace_params, boxcar_halfwidth)
    var_1d = extract_boxcar_core(1 ./ ivarimage, trace_params, boxcar_halfwidth)
    ivar_1d = 1.0 ./ var_1d
    mask_1d = extract_boxcar_bitmask(pix_bitmask, trace_params, boxcar_halfwidth)
    dropped_pixel_mask_1d = zeros(Int64, N_XPIX, N_FIBERS)

    if return_resids
        return flux_1d, ivar_1d, mask_1d, dropped_pixel_mask_1d, resid_fluxes_2d, resid_ivars_2d
    else
        return flux_1d, ivar_1d, mask_1d, dropped_pixel_mask_1d
    end
end

"""
Extract a 1D spectrum using a boxcar kernel with width from trace_params. This is used twice.
Once for the flux and once for the variance.
"""
function extract_boxcar_core(dimage_in, trace_params, boxcar_halfwidth)
    out = zeros(Float64, N_XPIX, N_FIBERS)
    n_xpix, n_fibers = size(trace_params)[1:2]
    for xpix in 1:n_xpix, fib in 1:n_fibers
        _, ypixf, _ = trace_params[xpix, fib, :]
        ypix = round(Int, ypixf)
        out[xpix, fib] = sum(dimage_in[xpix, (ypix - boxcar_halfwidth):(ypix + boxcar_halfwidth)])
    end
    out
end

function extract_boxcar_bitmask(dimage_in, trace_params, boxcar_halfwidth)
    mask = zeros(Int64, N_XPIX, N_FIBERS)
    n_xpix, n_fibers = size(trace_params)[1:2]
    for xpix in 1:n_xpix, fib in 1:n_fibers
        _, ypixf, _ = trace_params[xpix, fib, :]
        ypix = round(Int, ypixf)
        mask[xpix,
        fib] = reduce(
            |, dimage_in[xpix, (ypix - boxcar_halfwidth):(ypix + boxcar_halfwidth)])
    end
    mask
end

"""
Extract a 1D spectrum using the a gaussian kernel with center and width from trace_params.

# Keyword arguments:
- `window_size` is the number of pixels to sum over.

Spacing between traces is ~ 2000 / 300 = 6.67 pixels.
"""
function extract_optimal(dimage, ivarimage, pix_bitmask, trace_params; window_half_size = 4)
    n_xpix = size(trace_params, 1)
    n_fibers = size(trace_params, 2)

    # return values to be filled
    flux_1d = Matrix{Float64}(undef, n_xpix, n_fibers)
    ivar_1d = Matrix{Float64}(undef, n_xpix, n_fibers)
    mask_1d = Matrix{Int64}(undef, n_xpix, n_fibers)

    for xpix in 1:n_xpix, fib in 1:n_fibers
        _, y_peak, y_sigma = trace_params[xpix, fib, :]

        ypixels = floor(Int, y_peak - window_half_size):ceil(
            Int, y_peak + window_half_size)
        ypix_boundaries = [ypixels .- 0.5; ypixels[end] + 0.5]
        weights = diff(cdf.(Normal(y_peak, y_sigma), ypix_boundaries))

        flux_1d[xpix,
        fib] = sum(weights .* dimage[xpix, ypixels] .* ivarimage[xpix, ypixels]) /
               sum(weights .^ 2 .* ivarimage[xpix, ypixels])
        ivar_1d[xpix, fib] = sum(weights .^ 2 .* ivarimage[xpix, ypixels])

        # bitmask
        mask_1d[xpix, fib] = reduce(|, pix_bitmask[xpix, ypixels])
    end
    flux_1d, ivar_1d, mask_1d
end

function extract_optimal_iter(dimage, ivarimage, pix_bitmask, trace_params,
        med_center_to_fiber_func, x_prof_min, x_prof_max_ind,
        n_sub, min_prof_fib, max_prof_fib, all_y_prof, all_y_prof_deriv;
        small_window_half_size = 2, fit_window_half_size = 4,
        large_window_half_size = 12, n_max_repeat = 5, flag_thresh = 0.001,
        return_resids = false, neff_thresh = 10.0)
    n_xpix = size(trace_params, 1)
    n_ypix = size(dimage, 2)
    n_fibers = size(trace_params, 2)

    # return values to be filled
    flux_1d = Matrix{Float64}(undef, n_xpix, n_fibers)
    ivar_1d = Matrix{Float64}(undef, n_xpix, n_fibers)
    mask_1d = Matrix{Int64}(undef, n_xpix, n_fibers)
    dropped_pixel_mask_1d = Matrix{Int64}(undef, n_xpix, n_fibers)

    good_pixels = ((pix_bitmask .& bad_pix_bits) .== 0) .& (ivarimage .> 0)

    model_flux_indv = zeros(Float64, n_fibers, n_ypix)
    comb_model_flux = zeros(Float64, n_ypix)
    new_comb_model_flux = zeros(Float64, n_ypix)

    model_var_indv = zeros(Float64, n_fibers, n_ypix)
    comb_model_var = zeros(Float64, n_ypix)
    new_comb_model_var = zeros(Float64, n_ypix)

    new_flux_1d = zeros(Float64, n_fibers)

    if return_resids
        resid_fluxes_2d = zeros(Float64, size(dimage))
        resid_ivars_2d = zeros(Float64, size(dimage))
    end

    for xpix in 1:n_xpix
        #iterate on best-fit fluxes
        #to remove crosstalk/contribution
        #from neighbours

        #reset model fluxes
        model_flux_indv .= 0
        comb_model_flux .= 0
        model_var_indv .= 0
        comb_model_var .= 0
        new_flux_1d .= 0

        for repeat_ind in 1:n_max_repeat
            new_comb_model_flux .= 0
            new_comb_model_var .= 0

            if repeat_ind == 1
                window_half_size = small_window_half_size
            else
                window_half_size = fit_window_half_size
            end

            for fib in 1:n_fibers
                _, y_peak, y_sigma = trace_params[xpix, fib, :]
		y_peak_round = ceil(Int,round(y_peak))

		if ((y_peak_round + window_half_size) < 1) | ((y_peak_round - window_half_size) > N_XPIX) 
	            #then the peak is so far off the edge that
		    #we can't measure it's flux from the wings
		    #so give 0 flux and skip
                    new_flux_1d[fib] = 0.0
                    ivar_1d[xpix, fib] = 0.0
                    mask_1d[xpix, fib] = reduce(|, pix_bitmask[xpix, ypixels])
                    mask_1d[xpix, fib] |= bad_1d_no_good_pix
                    mask_1d[xpix, fib] |= bad_1d_failed_extract
		    continue
		end

                #use large window to get model fluxes
		full_ypixels = max(1,y_peak_round - large_window_half_size):min(
			              N_XPIX,y_peak_round + large_window_half_size)
		y_peak_ind = y_peak_round - full_ypixels[begin] + 1
                full_ypix_boundaries = [full_ypixels .- 0.5; full_ypixels[end] + 0.5]
                #                full_model_vals = diff(cdf.(Normal(y_peak, y_sigma), full_ypix_boundaries))
                prof_fib_ind = clamp(fib, min_prof_fib, max_prof_fib)
                full_model_vals = diff(cdf_func_indv(full_ypix_boundaries, y_peak, y_sigma,
                    prof_fib_ind, x_prof_min, x_prof_max_ind,
                    n_sub, min_prof_fib, all_y_prof, all_y_prof_deriv))

		ypixels = full_ypixels[max(1,y_peak_ind - window_half_size):min(size(full_ypixels,1),y_peak_ind + window_half_size)]
		model_vals = full_model_vals[max(1,y_peak_ind - window_half_size):min(size(full_ypixels,1),y_peak_ind + window_half_size)]
                if any(good_pixels[xpix, ypixels])
                    #then mask the bad pixels, same as setting ivar=0 there
                    model_vals[.!good_pixels[xpix, ypixels]] .= 0
                    good_flux_1d = true
                else
                    good_flux_1d = false
                end

                curr_pix_ivars = 1 ./ (1 ./ ivarimage[xpix, ypixels] .+ comb_model_var[ypixels] .-
                                  model_var_indv[fib, ypixels])
                ivar_1d[xpix, fib] = sum(model_vals .^ 2 .* curr_pix_ivars)
                flux_weights = model_vals .* curr_pix_ivars ./ ivar_1d[xpix, fib]
                new_flux_1d[fib] = sum(flux_weights .*
                                       (dimage[xpix, ypixels] .- comb_model_flux[ypixels] .+
                                        model_flux_indv[fib, ypixels]))

                # bitmask
                #		curr_good_fluxes = flux_weights .>= flag_thresh
                curr_good_fluxes = model_vals .>= flag_thresh
                if !good_flux_1d
                    curr_neff = sqrt(1 / sum(model_vals .^ 2))
                    mask_1d[xpix, fib] = reduce(|, pix_bitmask[xpix, ypixels])
                    mask_1d[xpix, fib] |= bad_1d_no_good_pix
                elseif any(curr_good_fluxes)
                    curr_neff = sqrt(1 / sum(model_vals[curr_good_fluxes] .^ 2))
                    mask_1d[xpix, fib] = reduce(|, pix_bitmask[xpix, ypixels[curr_good_fluxes]])
                    if any(.!curr_good_fluxes)
                        dropped_pixel_mask_1d[xpix, fib] = reduce(
                            |, pix_bitmask[xpix, ypixels[.!curr_good_fluxes]])
                    end
                else
                    curr_neff = sqrt(1 / sum(model_vals .^ 2))
                    mask_1d[xpix, fib] = reduce(|, pix_bitmask[xpix, ypixels])
                end

                if (!isfinite(new_flux_1d[fib])) | (ivar_1d[xpix, fib] == 0.0)
                    new_flux_1d[fib] = 0.0
                    ivar_1d[xpix, fib] = 0.0
                    mask_1d[xpix, fib] |= bad_1d_failed_extract
                end

                if curr_neff > neff_thresh
                    #                    new_flux_1d[fib] = 0.0
                    #                    ivar_1d[xpix, fib] = 0.0
                    mask_1d[xpix, fib] |= bad_1d_neff
                end

                if good_flux_1d
                    model_flux_indv[fib, full_ypixels] .= max(0, new_flux_1d[fib]) * full_model_vals
                    new_comb_model_flux[full_ypixels] .+= model_flux_indv[fib, full_ypixels]

                    if ivar_1d[xpix, fib] > 0
                        model_var_indv[fib,
                        full_ypixels] .= max(0, 1 / ivar_1d[xpix, fib]) *
                                         (full_model_vals .^ 2)
                        new_comb_model_var[full_ypixels] .+= model_var_indv[fib, full_ypixels]
                    end
                end
            end

            flux_1d[xpix, :] .= new_flux_1d
            if return_resids
                resid_fluxes_2d[xpix, :] .= dimage[xpix, :] .- new_comb_model_flux
                resid_ivars_2d[xpix, :] .= 1 ./ (1 ./ ivarimage[xpix, :] .+ new_comb_model_var)
            end

            if all(abs.(new_flux_1d .- flux_1d[xpix, :]) .< 0.01) & (repeat_ind > 1)
                break
            end

            comb_model_flux .= new_comb_model_flux
            comb_model_var .= new_comb_model_var
        end
    end

    if return_resids
        return flux_1d, ivar_1d, mask_1d, dropped_pixel_mask_1d, resid_fluxes_2d, resid_ivars_2d
    else
        return flux_1d, ivar_1d, mask_1d, dropped_pixel_mask_1d
    end
end

"""
Given an open HDF.file, `f`, and the telescope, mjd, and expnum, return a dictionary
mapping fiber index (1:300 laid out on the chip) to fiber type.
"""
function get_fibTargDict(f, tele, mjd, dfindx)
    # translate confSummary/almanac terminology to AR.jl terminology
    fiber_type_names = Dict(
        # fps era
        "science" => "sci",
        "sky_boss" => "skyB",
        "standard_apogee" => "tel",
        "sky_apogee" => "sky",

        # Plate era
        "science" => "sci",
        "standard" => "tel",
        "sky" => "sky"
    )
    # TODO other fiber types:
    # "blank"s from plate era
    # FPI era "serendipitous" APOGEE fibers are those which "accidentally" point at a bright
    # star (for BOSS reasons).
    # TODO Andrew thinks the fibers with category "" might be serendipitous targets

    mjdfps2plate = get_fps_plate_divide(tele)
    configIdCol = if parse(Int, mjd) > mjdfps2plate
        "config_id"
    else
        "plate_id"
    end

    df_exp = read_almanac_exp_df(f, tele, mjd)

    if !(dfindx in df_exp.exposure)
        @warn "Exposure $(dfindx) not found in $(tele)/$(mjd)/exposures"
        return Dict(1:300 .=> "fiberTypeFail"), Dict(1:300 .=> -2)
    end
    exposure_info = df_exp[dfindx, :]
    config_id = exposure_info[configIdCol]

    fibtargDict, fiber_sdss_id_Dict = if exposure_info.image_type != "object"
        Dict(1:300 .=> "cal"), Dict(1:300 .=> -2)
    else
        # Check if config_id is -1, which should not exist in the HDF5 file
        if config_id == -1
            @warn "config_id is -1 for exposure $(dfindx) in $(tele)/$(mjd). This should have been filtered out as flagged_bad=1. Returning fiberTypeFail for all fibers."
            Dict(1:300 .=> "fiberTypeFail"), Dict(1:300 .=> -2)
        else
            try
                df_fib = DataFrame(read(f["$(tele)/$(mjd)/fibers/$(config_id)"]))
                # normalizes all column names to lowercase
                rename!(df_fib, lowercase.(names(df_fib)))

                # limit to only the APOGEE fiber/hole information
                df_fib = if configIdCol == "config_id"
                    df_fib[df_fib[!, "fiber_type"].=="APOGEE",:]
                else
                    df_fib
                end

                fiber_types = map(df_fib[!, "category"]) do t
                    if t in keys(fiber_type_names)
                        fiber_type_names[t]
                    else
                        # @warn "Unknown fiber type for $(tele)/$(mjd)/fibers/$(config_id): $(repr(t))"
                        "fiberTypeFail"
                    end
                end
                fibernumvec = df_fib[!, "fiber_id"]
                fiber_sdss_id = df_fib[!, "sdss_id"]

                #this is a Hack and Andy Casey will replace this very very soon
                msknofiberdefaults = (fibernumvec .!= -1)
                fiber_types_full = repeat(["fiberTypeFail"], N_FIBERS)
                fiber_sdss_id_full = repeat([-2], N_FIBERS)
                try
                    fiber_types_full[fiberID2fiberIndx.(fibernumvec[msknofiberdefaults])] .= fiber_types[msknofiberdefaults]
                    fiber_sdss_id_full[fiberID2fiberIndx.(fibernumvec[msknofiberdefaults])] .= fiber_sdss_id[msknofiberdefaults]
                catch e
                    @warn "Problem with getting fiber type information for $(tele)/$(mjd)/fibers/$(config_id) (exposure $(dfindx)). Returning fiberTypeFail for all fibers."
                    show(e)
                    fiber_types_full .= "fiberTypeFail"
                end
                Dict(1:N_FIBERS .=> fiber_types_full), Dict(1:N_FIBERS .=> fiber_sdss_id_full)
            catch e
                @warn "Failed to get fiber type information for $(tele)/$(mjd)/fibers/$(config_id) (exposure $(dfindx)). Returning fiberTypeFail for all fibers."
                show(e)
                Dict(1:300 .=> "fiberTypeFail"), Dict(1:300 .=> -2)
            end
        end
    end

    if parse(Int, mjd) > mjdfps2plate
        fpifib1, fpifib2 = fiberID2fiberIndx.(get_fpi_guide_fiberID(tele))
        fibtargDict[fpifib1] = "fpiguide"
        fibtargDict[fpifib2] = "fpiguide"
    end
    return fibtargDict, fiber_sdss_id_Dict
end

# hardcoded to use chip c only for now
# must use dome flats, not quartz flats (need fiber runs to telescope)
# use full exposure_id
function get_fluxing_file(dfalmanac, parent_dir, tele, mjd, dfindx, runname; fluxing_chip = "B")
    df_mjd = sort(
        dfalmanac[(dfalmanac.mjd .== parse(Int, mjd)) .& (dfalmanac.observatory .== tele), :],
        :exposure)
    expIndex = dfindx
    cartId = df_mjd.cart_id[expIndex]
    image_type = df_mjd.image_type[expIndex]

    valid_flats4fluxing_fname = joinpath(parent_dir, "almanac/valid_domeflats4fluxing_$(runname).h5")
    if !isfile(valid_flats4fluxing_fname)
        @warn "Could not find any useful relfluxing files after looking for file $(valid_flats4fluxing_fname)"
	    return 2^2, nothing
    end
    f = h5open(valid_flats4fluxing_fname, "r")
    found_tele_mjd = false
    if tele in keys(f)
        if "$(mjd)" in keys(f[tele])
            found_tele_mjd = true
	    end
    end

    if !found_tele_mjd
        close(f)
        if ((image_type == "object") | (image_type == "domeflat"))
            @warn "Could not find any useful relfluxing files in file $(valid_flats4fluxing_fname) for tele $(tele) mjd $(mjd)"
        end
        return 2^2,nothing
    end

    cal_expid_list = read(f["$(tele)/$(mjd)"])
    close(f)

    if dfindx in cal_expid_list
        #the current files is one of the dome flats that has a relfluxing file
        return 2^0,get_fluxing_file_name(
            parent_dir, tele, mjd, df_mjd.exposure[expIndex], fluxing_chip, cartId)
    end

    expIndex_before = findlast(cal_expid_list .< dfindx)
    if !isnothing(expIndex_before)
        expIndex_before = cal_expid_list[expIndex_before]
    end
    expIndex_after = findfirst(cal_expid_list .> dfindx)
    if !isnothing(expIndex_after)
        expIndex_after = cal_expid_list[expIndex_after]
    end

    valid_before = if isnothing(expIndex_before)
        0
    elseif all(df_mjd.cart_id[expIndex_before:expIndex] .== cartId)
	    1
    elseif !isnothing(expIndex_before) & (df_mjd.cart_id[expIndex_before] .== cartId)
        2
    else
        0
    end
    valid_after = if isnothing(expIndex_after)
        0
    elseif all(df_mjd.cart_id[expIndex:expIndex_after] .== cartId)
        1
    elseif !isnothing(expIndex_after) & (df_mjd.cart_id[expIndex_after] .== cartId)
        2
    else
        0
    end

    if valid_before == 1
        return 2^0, get_fluxing_file_name(
            parent_dir, tele, mjd, df_mjd.exposure[expIndex_before], fluxing_chip, cartId)
    elseif valid_after == 1
        return 2^0, get_fluxing_file_name(
            parent_dir, tele, mjd, df_mjd.exposure[expIndex_after], fluxing_chip, cartId)
        # any of the cases below here we could consider using a global file
    elseif valid_before == 2
        return 2^1, get_fluxing_file_name(
            parent_dir, tele, mjd, df_mjd.exposure[expIndex_before], fluxing_chip, cartId)
    elseif valid_after == 2
        return 2^1, get_fluxing_file_name(
            parent_dir, tele, mjd, df_mjd.exposure[expIndex_after], fluxing_chip, cartId)
    else
        return 2^2,nothing
    end
end

# TODO: switch to meta data dict and then save wavecal flags etc.
function reinterp_spectra(fname, roughwave_dict; checkpoint_mode = "commit_same", outdir = "../outdir")
    # might need to add in telluric div functionality here?
    outname = replace(replace(fname, "ar1D" => "ar1Duni"), "_$(FIRST_CHIP)_" => "_")
    if check_file(outname, mode = checkpoint_mode)
        return
    end

    sname = split(split(split(fname, "/")[end], ".h5")[1], "_")
    fnameType, tele, mjd, expnum, chip, image_type = sname[(end - 5):end]
    mjd_int = parse(Int, mjd)

    backupWave_fname = joinpath(
        outdir, "wavecal", "wavecalNightAve_$(tele)_$(mjd).h5")

    # could shift this to a preallocation step
    outflux = zeros(length(logUniWaveAPOGEE), N_FIBERS)
    outvar = zeros(length(logUniWaveAPOGEE), N_FIBERS)
    outmsk = zeros(Int, length(logUniWaveAPOGEE), N_FIBERS)
    outDropMsk = zeros(Int, length(logUniWaveAPOGEE), N_FIBERS)
    outTraceCoords = zeros(length(logUniWaveAPOGEE), N_FIBERS, 3) #(x,y,chipInt)
    cntvec = zeros(Int, length(logUniWaveAPOGEE), N_FIBERS)

    pixvec = 1:(N_CHIPS * N_XPIX)
    xpix_stack = ((pixvec .- 1) .% N_XPIX) .+ 1
    flux_stack = zeros(N_CHIPS * N_XPIX, N_FIBERS)
    ivar_stack = zeros(N_CHIPS * N_XPIX, N_FIBERS)
    mask_stack = zeros(Int, N_CHIPS * N_XPIX, N_FIBERS)
    dropped_mask_stack = zeros(Int, N_CHIPS * N_XPIX, N_FIBERS)
    wave_stack = zeros(N_CHIPS * N_XPIX, N_FIBERS)
    trace_center_stack = zeros(N_CHIPS * N_XPIX, N_FIBERS)
    chipInt_stack = zeros(N_CHIPS * N_XPIX, N_FIBERS)
    chipBit_stack = zeros(Int, N_CHIPS * N_XPIX, N_FIBERS)

    ingestBit = zeros(Int, N_FIBERS)

    # add a for loop over the exposures (stop thinking about "visits" for now)
    # probably just generate ap1D file names from the alamanc files

    # this was used for looping over exposures in the visit
    # outdir = "/uufs/chpc.utah.edu/common/home/u6039752/scratch1/working/2024_12_05/outdir/"
    # fname = outdir * "apred/$(mjd)/" * get_1d_name(parse(Int, last(expid,4)), df) * ".h5"

    wavetype_order = ["fpi", "sky"]
    found_soln = false
    wavecal_type = ""
    for wavetype in wavetype_order
        wavecal_type = "waveCalNight$(wavetype)Dither"
        wavefname = replace(replace(fname, fnameType => wavecal_type), "_$(FIRST_CHIP)_" => "_")
        if isfile(wavefname)
            f = jldopen(wavefname)
            chipWaveSoln = f["chipWaveSoln"]
            close(f)
	    if all(isfinite.(chipWaveSoln))
                found_soln = true
                break
	    end
        end
    end

    if !found_soln
        if !isfile(backupWave_fname)
	    curr_best_wave_type = "rough"
            backupWaveSoln = nothing
	    backupWave_fname = nothing
        else
            f = h5open(backupWave_fname, "r")
            curr_best_wave_type = read(f["best_wave_type"])
	    backupWaveSoln = read(f["$(curr_best_wave_type)/nightAve_wave_soln"])
	    close(f)
        end

        #this is not a great global fallback, but it works so we get something to look at
        if isnothing(backupWave_fname) || (!all(isfinite.(backupWaveSoln)))
            chipWaveSoln = zeros(N_XPIX, N_FIBERS, N_CHIPS)
            for (chipind, chip) in enumerate(CHIP_LIST)
                chipWaveSoln[:, :, chipind] .= rough_linear_wave.(
                    1:N_XPIX, a = roughwave_dict[tele][chip][1], b = roughwave_dict[tele][chip][2])
            end
            if !(image_type in ["dark", "internalflat", "quartzflat", "domeflat"])
                println("No wavecal found for $(fname), using rough linear fallback as fallback")
                flush(stdout)
            end
            wavecal_type = "error_fixed_fallback"
        else
            chipWaveSoln = backupWaveSoln
            if !(image_type in ["dark", "internalflat", "quartzflat", "domeflat"])
                println("No wavecal found for $(fname), using $(curr_best_wave_type) nightly average as fallback")
            end
            flush(stdout)
	        wavecal_type = "error_night_ave_$(curr_best_wave_type)"
        end
    end

    metadata_lst = []
    trace_lst = []
    for (chipind, chip) in enumerate(CHIP_LIST) # This needs to be the in abc RGB order, changing that will break this section
        fnameloc = replace(fname, "_$(FIRST_CHIP)_" => "_$(chip)_")
        f = jldopen(fnameloc)
        flux_1d = f["flux_1d"]
        ivar_1d = f["ivar_1d"]
        mask_1d = f["mask_1d"]
        dropped_pixels_mask_1d = f["dropped_pixels_mask_1d"]
        extract_trace_centers = f["extract_trace_centers"]
        close(f)
        push!(metadata_lst, read_metadata(fnameloc))

        flux_stack[(1:N_XPIX) .+ (3 - chipind) * N_XPIX, :] .= flux_1d[end:-1:1, :]
        ivar_stack[(1:N_XPIX) .+ (3 - chipind) * N_XPIX, :] .= ivar_1d[end:-1:1, :]
        mask_stack[(1:N_XPIX) .+ (3 - chipind) * N_XPIX, :] .= mask_1d[end:-1:1, :]
        dropped_mask_stack[(1:N_XPIX) .+ (3 - chipind) * N_XPIX, :] .= dropped_pixels_mask_1d[
            end:-1:1, :]
        try
            wave_stack[(1:N_XPIX) .+ (3 - chipind) * N_XPIX, :] .= chipWaveSoln[end:-1:1, :, chipind]
        catch
            println((typeof(wave_stack), typeof(N_XPIX), typeof(chipind), typeof(chipWaveSoln),fname, mjd_int, typeof(backupWaveSoln), typeof(backupWaveSoln)))
        end
        wave_stack[(1:N_XPIX) .+ (3 - chipind) * N_XPIX, :] .= chipWaveSoln[end:-1:1, :, chipind]
        trace_center_stack[(1:N_XPIX) .+ (3 - chipind) * N_XPIX, :] .= extract_trace_centers[
            end:-1:1, :]
        chipBit_stack[(1:N_XPIX) .+ (3 - chipind) * N_XPIX, :] .+= 2^(chipind)
        chipInt_stack[(1:N_XPIX) .+ (3 - chipind) * N_XPIX, :] .= chipind
    end

    # should add a check all entries of metadata_lst to be equal
    metadata = metadata_lst[1]
    metadata["wavecal_type"] = wavecal_type # add wavecal type to metadata

    noBadBits = (mask_stack .& bad_pix_bits .== 0)
    chipBit_stack[.!(noBadBits)] .+= 2^4 # call pixels with bad bits thrown bad
    chipBit_stack[chipBit_stack .== 0] .+= 2^4 # call missing chips bad

    # think about adding warnings for the last two cases
    good_pix = ((noBadBits) .& (chipBit_stack .& 2^4 .== 0)) .& (.!isnan.(ivar_stack)) .&
               (ivar_stack .> (10^-20))

    ## need to propagate the bit mask
    for fiberindx in 1:N_FIBERS
        good_pix_fiber = good_pix[:, fiberindx]
        flux_fiber = flux_stack[good_pix_fiber, fiberindx]
        ivar_fiber = ivar_stack[good_pix_fiber, fiberindx]
        wave_fiber = wave_stack[good_pix_fiber, fiberindx]
        trace_center_fiber = trace_center_stack[good_pix_fiber, fiberindx]
        chipBit_fiber = chipBit_stack[good_pix_fiber, fiberindx]
        chipInt_fiber = chipInt_stack[good_pix_fiber, fiberindx]
        pixindx_fiber = pixvec[good_pix_fiber]
        xpix_fiber = xpix_stack[good_pix_fiber]

        Rinv = generateInterpMatrix_sparse_inv(
            wave_fiber, chipBit_fiber, logUniWaveAPOGEE, pixindx_fiber)
        normvec = dropdims(sum(Rinv, dims = 2), dims = 2)
        msk_inter = (normvec .!= 0)

        outflux[msk_inter, fiberindx] .+= (Rinv * flux_fiber)[msk_inter]
        outvar[msk_inter, fiberindx] .+= ((Rinv .^ 2) * (1 ./ ivar_fiber))[msk_inter]
        cntvec[:, fiberindx] .+= msk_inter

        #right now, only works for a single exposure
        if length(wave_fiber) > 1
            # Check if wave_fiber has unique, sorted values
            if length(unique(wave_fiber)) != length(wave_fiber) || !issorted(wave_fiber)
                @warn "Non-unique or unsorted wavelengths for fiber $fiberindx in $fnameType $tele $mjd $expnum $chip $image_type. Cannot interpolate, filling trace coordinates with NaN."
                outTraceCoords[:, fiberindx, 1] .= NaN
                outTraceCoords[:, fiberindx, 2] .= NaN
                outTraceCoords[:, fiberindx, 3] .= NaN
            else
                outTraceCoords[:, fiberindx, 1] .= linear_interpolation(
                    wave_fiber, xpix_fiber, extrapolation_bc = Line()).(logUniWaveAPOGEE)
                outTraceCoords[:, fiberindx, 2] .= linear_interpolation(
                    wave_fiber, trace_center_fiber, extrapolation_bc = Line()).(logUniWaveAPOGEE)
                outTraceCoords[:, fiberindx, 3] .= linear_interpolation(
                    wave_fiber, chipInt_fiber, extrapolation_bc = Line()).(logUniWaveAPOGEE)
            end
        else
            # If no good pixels or only 1 good pixel for this fiber, fill with NaN
            if length(wave_fiber) == 0
                @warn "No good pixels found for fiber $fiberindx in $fnameType $tele $mjd $expnum $chip $image_type. Filling trace coordinates with NaN."
            else
                @warn "Only 1 good pixel found for fiber $fiberindx in $fnameType $tele $mjd $expnum $chip $image_type. Cannot interpolate, filling trace coordinates with NaN."
            end
            outTraceCoords[:, fiberindx, 1] .= NaN
            outTraceCoords[:, fiberindx, 2] .= NaN
            outTraceCoords[:, fiberindx, 3] .= NaN
        end

        if all(isnanorzero.(flux_fiber)) && ((ingestBit[fiberindx] & 2^1) == 0)
            ingestBit[fiberindx] += 2^1 # ap1D exposure flux are all NaNs (for at least one of the exposures)
        elseif all(.!((chipBit_fiber .& 2^4) .== 0)) && ((ingestBit[fiberindx] & 2^2) == 0)
            ingestBit[fiberindx] += 2^2 # ap1D exposure pixels are all bad by bit mask  (for at least one of the exposures)
        elseif all(isnanorzero.(ivar_fiber)) && ((ingestBit[fiberindx] & 2^3) == 0)
            ingestBit[fiberindx] += 2^3 # ap1D exposure ivars are all NaNs or zeros (for at least one of the exposures)
        end
    end

    framecnts = maximum(cntvec, dims = 1) #     framecnts = maximum(cntvec) # a little shocked that I throw it away if it is bad in even one frame
    outflux ./= framecnts
    outvar ./= (framecnts .^ 2)
    # need to update this to a bit mask that is all or any for the pixels contributing to the reinterpolation
    outmsk = (cntvec .== framecnts)

    # Write reinterpolated data
    safe_jldsave(
        outname, metadata; flux_1d = outflux, ivar_1d = 1 ./ outvar, mask_1d = outmsk,
        extract_trace_coords = outTraceCoords)
    return
end

const logUniWaveAPOGEE = 10 .^ range((start = 4.17825), step = 6.0e-6, length = 8700);

#should add a check_file call for this one
function process_1D(fname;
        outdir::String,
        runname::String,
        extraction::String,
        relFlux::Bool,
        chip_list::Vector{String} = CHIP_LIST,
        profile_path = "./data/",
        plot_path = "../outdir/$(sjd)/plots/",
        checkpoint_mode = "commit_same")
    sname = split(split(split(fname, "/")[end], ".h5")[1], "_")
    fnameType, tele, mjd, expnum, chip, image_type = sname[(end - 5):end]
    dfindx = parse(Int, expnum)

    outfname = if relFlux
        replace(fname, "ar2D" => "ar1D")
    else
        replace(replace(fname, "ar2D" => "ar1D"), "apred" => "apredrelflux")
    end

    if check_file(outfname, mode = checkpoint_mode)
        return true
    end

    # this seems annoying to load so often if we know we are doing a daily... need to ponder
    traceFname = outdir * "apred/$(mjd)/traceMain_$(tele)_$(mjd)_$(chip).h5"

    # how worried should I be about loading this every time?
    falm = h5open(joinpath(outdir, "almanac/$(runname).h5"))
    dfalmanac = read_almanac_exp_df(falm, tele, mjd)

    (med_center_to_fiber_func, x_prof_min, x_prof_max_ind, n_sub, min_prof_fib, max_prof_fib,
    all_y_prof, all_y_prof_deriv) = get_default_trace_hyperparams(tele, chip, profile_path = profile_path, plot_path = plot_path)

    fnamecal = if (fnameType == "ar2D")
        replace(fname, "ar2D" => "ar2Dcal")
    else
        fname
    end

    dimage = load(fname, "dimage")
    ivarimage = load(fname, "ivarimage")
    pix_bitmask = load(fnamecal, "pix_bitmask")
    metadata = read_metadata(fname)

    regularized_trace_params = try
        load(traceFname, "regularized_trace_params")
    catch
        @warn "No regularized trace params found for $(traceFname)"
        return false
    end
    trace_metadata = read_metadata(traceFname)
    metadata = merge(metadata, trace_metadata)

    flux_1d, ivar_1d,
    mask_1d,
    dropped_pixels_mask_1d,
    resid_flux,
    resid_ivar = if extraction == "boxcar"
        extract_boxcar(
            dimage, ivarimage, pix_bitmask, regularized_trace_params, return_resids = true)
    elseif extraction == "optimal"
        #            extract_optimal(dimage, ivarimage, pix_bitmask, regularized_trace_params)
        extract_optimal_iter(dimage, ivarimage, pix_bitmask, regularized_trace_params,
            med_center_to_fiber_func, x_prof_min, x_prof_max_ind, n_sub,
            min_prof_fib, max_prof_fib, all_y_prof, all_y_prof_deriv, return_resids = true)
    else
        error("Extraction method $(extraction) not recognized")
    end

    resid_outfname = replace(fname, "ar2D" => "ar2Dresiduals")
    safe_jldsave(resid_outfname, metadata; resid_flux, resid_ivar, trace_used_param_fname = traceFname)
    if relFlux
        # relative fluxing (using B (last chip) only for now)
        # this is the path to the underlying fluxing file.
        # it is symlinked below to an exposure-specific file (linkPath).
        relflux_bit,calPath = get_fluxing_file(
            dfalmanac, outdir, tele, mjd, dfindx, runname, fluxing_chip = chip_list[end])
        fibtargDict, fiber_sdss_id_Dict = get_fibTargDict(falm, tele, mjd, dfindx)
        fiberTypeList = map(x -> fibtargDict[x], 1:300)

        if isnothing(calPath)
            # TODO uncomment this
            if (image_type == "object") | (image_type == "domeflat")
                @warn "No fluxing file available for $(tele) $(mjd) $(dfindx) $(chip)"
            end
            relthrpt = ones(size(flux_1d, 2))
            bitmsk_relthrpt = 2^2 * ones(Int, size(flux_1d, 2))
        elseif !isfile(calPath)
            error("Fluxing file $(calPath) for $(tele) $(mjd) $(dfindx) $(chip) does not exist")
        else
            linkPath = abspath(joinpath(
                dirname(fname), "relFlux_$(tele)_$(mjd)_$(dfindx_fname_format(dfindx))_$(chip).h5"))
            if !islink(linkPath)
                symlink(abspath(calPath), linkPath)
            end
            relthrpt = load(linkPath, "relthrpt")
            relthrptr = reshape(relthrpt, (1, length(relthrpt)))
            bitmsk_relthrpt = load(linkPath, "bitmsk_relthrpt")
        end

        # don't flux broken fibers (don't use warn fibers for sky)
        msk_goodwarn = (bitmsk_relthrpt .== 0) .| (bitmsk_relthrpt .& 2^0) .== 2^0
        if any(msk_goodwarn)
            flux_1d[:, msk_goodwarn] ./= relthrptr[:, msk_goodwarn]
            ivar_1d[:, msk_goodwarn] .*= relthrptr[:, msk_goodwarn] .^ 2
        end

	metadata["bitmsk_relFluxFile"] = relflux_bit
        # we probably want to append info from the fiber dictionary from alamanac into the file name
        safe_jldsave(outfname, metadata; flux_1d, ivar_1d, mask_1d, dropped_pixels_mask_1d,
            extract_trace_centers = regularized_trace_params[:, :, 2],
            relthrpt, bitmsk_relthrpt, fiberTypeList,
            trace_used_param_fname = traceFname)
    else
        dirName = dirname(outfname)
        if !ispath(dirName)
            mkpath(dirName)
        end
        safe_jldsave(
            outfname, metadata; flux_1d, ivar_1d, mask_1d, dropped_pixels_mask_1d,
            extract_trace_centers = regularized_trace_params[:, :, 2],
            trace_used_param_fname = traceFname)
    end
    close(falm)
    return true
end
