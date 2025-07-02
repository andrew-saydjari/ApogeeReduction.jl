using FastRunningMedian: running_median
using Distributions: cdf, Normal
using Interpolations: linear_interpolation, Line
using DataFrames

include("./wavecal.jl")
include("./skyline_peaks.jl")
include("./fileNameHandling.jl")
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

    if return_resids
        return flux_1d, ivar_1d, mask_1d, resid_fluxes_2d, resid_ivars_2d
    else
        return flux_1d, ivar_1d, mask_1d
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

                #use large window to get model fluxes
                full_ypixels = floor(
                    Int, y_peak - large_window_half_size):ceil(
                    Int, y_peak + large_window_half_size)
                full_ypix_boundaries = [full_ypixels .- 0.5; full_ypixels[end] + 0.5]
                #                full_model_vals = diff(cdf.(Normal(y_peak, y_sigma), full_ypix_boundaries))
                prof_fib_ind = clamp(fib, min_prof_fib, max_prof_fib)
                full_model_vals = diff(cdf_func_indv(full_ypix_boundaries, y_peak, y_sigma,
                    prof_fib_ind, x_prof_min, x_prof_max_ind,
                    n_sub, min_prof_fib, all_y_prof, all_y_prof_deriv))

                ypixels = full_ypixels[(begin + (large_window_half_size - window_half_size)):(end - (large_window_half_size - window_half_size))]
                model_vals = full_model_vals[(begin + (large_window_half_size - window_half_size)):(end - (large_window_half_size - window_half_size))]
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
        return flux_1d, ivar_1d, mask_1d, resid_fluxes_2d, resid_ivars_2d
    else
        return flux_1d, ivar_1d, mask_1d
    end
end

"""
Given an open HDF.file, `f`, and the telescope, mjd, and expnum, return a dictionary
mapping fiber id to fiber type.
"""
function get_fibTargDict(f, tele, mjd, expnum)
    exposure_id = short_expid_to_long(mjd, expnum)

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
    configName, configIdCol, target_type_col = if parse(Int, mjd) > mjdfps2plate
        "fps", "configid", "category"
    else
        "plates", "plateid", "target_type" # TODO should this be source_type?
    end

    df_exp = read_almanac_exp_df(f, tele, mjd)

    if !(exposure_id in df_exp.exposure_str)
        @warn "Exposure $(exposure_id) not found in $(tele)/$(mjd)/exposures"
        return Dict(1:300 .=> "fiberTypeFail")
    end
    exposure_info = df_exp[findfirst(df_exp[!, "exposure_str"] .== exposure_id), :]
    configid = exposure_info[configIdCol]

    fibtargDict = if exposure_info.exptype == "OBJECT"
        try
            df_fib = DataFrame(read(f["$(tele)/$(mjd)/fibers/$(configName)/$(configid)"]))
            # normalizes all column names to lowercase
            rename!(df_fib, lowercase.(names(df_fib)))

            # sometimes the fiber id column is called "fiber_id" and sometimes it is called "fiberid"
            fiberid_col = if "fiber_id" in names(df_fib)
                "fiber_id"
            elseif "fiberid" in names(df_fib)
                "fiberid"
            else
                error("fiber_id or fiberid column not found in $(configName)/$(configid). Available columns: $(names(df_fib))")
            end

            fiber_types = map(df_fib[!, target_type_col]) do t
                if t in keys(fiber_type_names)
                    fiber_type_names[t]
                else
                    # @warn "Unknown fiber type for $(tele)/$(mjd)/fibers/$(configName)/$(configid): $(repr(t))"
                    "fiberTypeFail"
                end
            end
            fibernum_col = df_fib[!, fiberid_col]
            # println(typeof(fibernum_col))
            fibernumvec = if fibernum_col isa AbstractVector{<:Integer}
                fibernum_col
            elseif fibernum_col isa AbstractVector{<:String}
                parse.(Int, fibernum_col)
            else
                @warn "Fiber numbers are neither integers or strings"
                fibernum_col
            end
            Dict(fibernumvec .=> fiber_types)
        catch e
            rethrow(e)
            @warn "Failed to get any fiber type information for $(tele)/$(mjd)/fibers/$(configName)/$(configid) (exposure $(exposure_id)). Returning fiberTypeFail for all fibers."
            show(e)
            Dict(1:300 .=> "fiberTypeFail")
        end
    else
        Dict(1:300 .=> "cal")
    end

    if parse(Int, mjd) > mjdfps2plate
        fpifib1, fpifib2 = get_fpi_guide_fiberID(tele)
        fibtargDict[fpifib1] = "fpiguide"
        fibtargDict[fpifib2] = "fpiguide"
    end
    return fibtargDict
end

# hardcoded to use chip c only for now
# must use dome flats, not quartz flats (need fiber runs to telescope)
# use full exposure_id
function get_fluxing_file(dfalmanac, parent_dir, tele, mjd, expnum; fluxing_chip = "B")
    exposure_id = parse(Int, short_expid_to_long(mjd, expnum))
    df_mjd = sort(
        dfalmanac[(dfalmanac.mjd .== parse(Int, mjd)) .& (dfalmanac.observatory .== tele), :],
        :exposure)
    expIndex = findfirst(df_mjd.exposure_int .== exposure_id)
    cartId = df_mjd.cartidInt[expIndex]
    # this needs to have cuts that match those in make_runlist_dome_flats.jl
    expIndex_before = findlast((df_mjd.imagetyp .== "DomeFlat") .&
                               (df_mjd.exposure_int .< exposure_id) .& (df_mjd.nreadInt .> 3))
    expIndex_after = findfirst((df_mjd.imagetyp .== "DomeFlat") .&
                               (df_mjd.exposure_int .> exposure_id) .& (df_mjd.nreadInt .> 3))
    valid_before = if !isnothing(expIndex_before)
        all(df_mjd.cartidInt[expIndex_before:expIndex] .== cartId) * 1
    elseif !isnothing(expIndex_before)
        (df_mjd.cartidInt[expIndex_after] .== cartId) * 2
    else
        0
    end
    valid_after = if !isnothing(expIndex_after)
        all(df_mjd.cartidInt[expIndex:expIndex_after] .== cartId) * 1
    elseif !isnothing(expIndex_after)
        (df_mjd.cartidInt[expIndex_after] .== cartId) * 2
    else
        0
    end

    if valid_before == 1
        return get_fluxing_file_name(
            parent_dir, tele, mjd, last(df_mjd.exposure_str[expIndex_before], 4), fluxing_chip, cartId)
    elseif valid_after == 1
        return get_fluxing_file_name(
            parent_dir, tele, mjd, last(df_mjd.exposure_str[expIndex_after], 4), fluxing_chip, cartId)
        # any of the cases below here we could consider using a global file
    elseif valid_before == 2
        return get_fluxing_file_name(
            parent_dir, tele, mjd, last(df_mjd.exposure_str[expIndex_before], 4), fluxing_chip, cartId)
    elseif valid_after == 2
        return get_fluxing_file_name(
            parent_dir, tele, mjd, last(df_mjd.exposure_str[expIndex_after], 4), fluxing_chip, cartId)
    else
        return nothing
    end
end

# TODO: switch to meta data dict and then save wavecal flags etc.
function reinterp_spectra(fname; backupWaveSoln = nothing)
    # might need to add in telluric div functionality here?

    sname = split(split(split(fname, "/")[end], ".h5")[1], "_")
    fnameType, tele, mjd, expnum, chip, exptype = sname[(end - 5):end]

    # could shift this to a preallocation step
    outflux = zeros(length(logUniWaveAPOGEE), N_FIBERS)
    outvar = zeros(length(logUniWaveAPOGEE), N_FIBERS)
    outmsk = zeros(Int, length(logUniWaveAPOGEE), N_FIBERS)
    outTraceCoords = zeros(length(logUniWaveAPOGEE), N_FIBERS, 3) #(x,y,chipInt)
    cntvec = zeros(Int, length(logUniWaveAPOGEE), N_FIBERS)

    pixvec = 1:(N_CHIPS * N_XPIX)
    xpix_stack = ((pixvec .- 1) .% N_XPIX) .+ 1
    flux_stack = zeros(N_CHIPS * N_XPIX, N_FIBERS)
    ivar_stack = zeros(N_CHIPS * N_XPIX, N_FIBERS)
    mask_stack = zeros(Int, N_CHIPS * N_XPIX, N_FIBERS)
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
            found_soln = true
            break
        end
    end
    if !found_soln
        #this is not a great global fallback, but it works so we get something to look at
        if isnothing(backupWaveSoln)
            chipWaveSoln = zeros(N_XPIX, N_FIBERS, N_CHIPS)
            for (chipind, chip) in enumerate(CHIP_LIST)
                chipWaveSoln[:, :, chipind] .= rough_linear_wave.(
                    1:N_XPIX, a = roughwave_dict[tele][chip][1], b = roughwave_dict[tele][chip][2])
            end
            println("No wavecal found for $(fname), using rough linear fallback")
            flush(stdout)
            wavecal_type = "error_fixed_fallback"
        else
            chipWaveSoln = backupWaveSoln
            println("No wavecal found for $(fname), using nightly average as fallback")
            flush(stdout)
            wavecal_type = "error_night_ave_fallback"
        end
    end

    metadata_lst = []
    for (chipind, chip) in enumerate(CHIP_LIST) # This needs to be the in abc RGB order, changing that will break this section
        fnameloc = replace(fname, "_$(FIRST_CHIP)_" => "_$(chip)_")
        f = jldopen(fnameloc)
        flux_1d = f["flux_1d"]
        ivar_1d = f["ivar_1d"]
        mask_1d = f["mask_1d"]
	    extract_trace_centers = f["extract_trace_centers"]
        close(f)
        push!(metadata_lst,read_metadata(fnameloc))

        flux_stack[(1:N_XPIX) .+ (3 - chipind) * N_XPIX, :] .= flux_1d[end:-1:1, :]
        ivar_stack[(1:N_XPIX) .+ (3 - chipind) * N_XPIX, :] .= ivar_1d[end:-1:1, :]
        mask_stack[(1:N_XPIX) .+ (3 - chipind) * N_XPIX, :] .= mask_1d[end:-1:1, :]
        wave_stack[(1:N_XPIX) .+ (3 - chipind) * N_XPIX, :] .= chipWaveSoln[end:-1:1, :, chipind]
        trace_center_stack[(1:N_XPIX) .+ (3 - chipind) * N_XPIX, :] .= extract_trace_centers[end:-1:1, :]
        chipBit_stack[(1:N_XPIX) .+ (3 - chipind) * N_XPIX, :] .+= 2^(chipind)
        chipInt_stack[(1:N_XPIX) .+ (3 - chipind) * N_XPIX, :] .= chipind
    end

    # should add a check all entries of metadata_lst to be equal
    metadata = metadata_lst[1]

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
	outTraceCoords[:, fiberindx, 1] .= linear_interpolation(wave_fiber, xpix_fiber, extrapolation_bc = Line()).(logUniWaveAPOGEE)
	outTraceCoords[:, fiberindx, 2] .= linear_interpolation(wave_fiber, trace_center_fiber, extrapolation_bc = Line()).(logUniWaveAPOGEE)
	outTraceCoords[:, fiberindx, 3] .= linear_interpolation(wave_fiber, chipInt_fiber, extrapolation_bc = Line()).(logUniWaveAPOGEE)

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
    outname = replace(replace(fname, "ar1D" => "ar1Duni"), "_$(FIRST_CHIP)_" => "_")
    safe_jldsave(
        outname, metadata; flux_1d = outflux, ivar_1d = 1 ./ outvar, mask_1d = outmsk, extract_trace_coords = outTraceCoords, wavecal_type)
    return
end

logUniWaveAPOGEE = 10 .^ range((start = 4.17825), step = 6.0e-6, length = 8700);
