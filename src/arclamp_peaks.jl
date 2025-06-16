#using StatsBase

using Polynomials: fit
using SpecialFunctions: erf
using Interpolations: linear_interpolation, Line

function read_fpiPeakLoc_coeffs(tele, chip;
        data_path = "./data/")

    fname = joinpath(data_path,"fpiPeakSpacingParams_$(tele)_$(chip).h5")

    # opening the file with a "do" closure guarantees that the file is closed
    # (analogous to "with open() as f:" in Python)
    coeffs_peak_ind_to_x,
    coeffs_x_to_peak_ind = jldopen(fname) do params
        (params["coeffs_peak_ind_to_x"], params["coeffs_x_to_peak_ind"])
    end

    return coeffs_peak_ind_to_x,coeffs_x_to_peak_ind
end

function gaussian_int_pdf(x_vals, params; return_deriv = false)
    #calculate the integrated gaussian probabilities

    #params shape = (n_peaks,n_params)
    #x_vals shape = (n_peaks,n_pix)

    #params = (height,mean,width)
    #x_vals = integers (ie pixels)

    z_above = ((x_vals .+ 0.5) .- params[:, 2]) ./ params[:, 3]
    z_below = ((x_vals .- 0.5) .- params[:, 2]) ./ params[:, 3]

    cdf = (0.5 * erf.(z_above ./ (2^0.5))) .- (0.5 * erf.(z_below ./ (2^0.5)))

    if return_deriv
        dcdf_dmu = (2 * π)^(-0.5) .*
                   (exp.(-1 .* z_above .^ 2 ./ 2) .- exp.(-1 .* z_below .^ 2 ./ 2)) .*
                   (-1 ./ params[:, 3])
        dcdf_dsig = (2 * π)^(-0.5) .*
                    (exp.(-1 .* z_above .^ 2 ./ 2) .* (-z_above ./ params[:, 3]) .-
                     exp.(-1 .* z_below .^ 2 ./ 2) .* (-z_below ./ params[:, 3]))
        #
        #        dmu = 0.0001
        #        z_above_above = ((x_vals .+ 0.5) .- (params[:, 2] .+ dmu)) ./ params[:, 3]
        #        z_below_above = ((x_vals .- 0.5) .- (params[:, 2] .+ dmu)) ./ params[:, 3]
        #        cdf_above = (0.5 * erf.(z_above_above ./ (2 ^ 0.5)))            .- (0.5 * erf.(z_below_above ./ (2 ^ 0.5))) 
        #        z_above_below = ((x_vals .+ 0.5) .- (params[:, 2] .- dmu)) ./ params[:, 3]
        #        z_below_below = ((x_vals .- 0.5) .- (params[:, 2] .- dmu)) ./ params[:, 3]
        #        cdf_below = (0.5 * erf.(z_above_below ./ (2 ^ 0.5)))            .- (0.5 * erf.(z_below_below ./ (2 ^ 0.5))) 
        #        dcdf_dmu_emp = (cdf_above .- cdf_below) ./ (2*dmu)
        #    
        #        dsig = 0.001
        #        z_above_above = ((x_vals .+ 0.5) .- (params[:, 2] .+ 0.0)) ./ (params[:, 3] .+ dsig)
        #        z_below_above = ((x_vals .- 0.5) .- (params[:, 2] .+ 0.0)) ./ (params[:, 3] .+ dsig)
        #        cdf_above = (0.5 * erf.(z_above_above ./ (2 ^ 0.5)))            .- (0.5 * erf.(z_below_above ./ (2 ^ 0.5))) 
        #        z_above_below = ((x_vals .+ 0.5) .- (params[:, 2] .+ 0.0)) ./ (params[:, 3] .- dsig)
        #        z_below_below = ((x_vals .- 0.5) .- (params[:, 2] .+ 0.0)) ./ (params[:, 3] .- dsig)
        #        cdf_below = (0.5 * erf.(z_above_below ./ (2 ^ 0.5))) .- (0.5 * erf.(z_below_below ./ (2 ^ 0.5))) 
        #        dcdf_dsig_emp = (cdf_above .- cdf_below) ./ (2*dsig)
        #
        #	fib_ind = 10
        #	println("diff ",(dcdf_dmu .- dcdf_dmu_emp)[fib_ind,:])
        #	println("orig ",dcdf_dmu[fib_ind,:])
        #	println("emp ",dcdf_dmu_emp[fib_ind,:])
        #
        #	println("diff ",(dcdf_dsig .- dcdf_dsig_emp)[fib_ind,:])
        #	println("orig ",dcdf_dsig[fib_ind,:])
        #	println("emp ",dcdf_dsig_emp[fib_ind,:])

        return cdf, dcdf_dmu, dcdf_dsig
    else
        return cdf
    end
end

function fit_gauss_and_bias(all_rel_fluxes, all_rel_ivars, first_guess_params,
        fit_inds, best_model_fit_inds, offset_inds;
        n_iter = 10, dmu = 0.001, dsig = 0.001, return_cov = false,
        use_first_guess_heights = false, max_center_move = 3,
        min_widths = 0.1, max_widths = 3.0, n_pad = 0, n_offset = 10)
    n_pixels = size(all_rel_fluxes, 1)
    fit_inds .+= n_pad
    best_model_fit_inds .+= n_pad

    fit_fluxes = copy(all_rel_fluxes[fit_inds])'
    fit_ivars = copy(all_rel_ivars[fit_inds])'
    good_ivars = fit_ivars .> 0
    if any(good_ivars)
        fit_ivars[fit_ivars .== 0] .= minimum(fit_ivars[fit_ivars .> 0]) / (10.0^2)
    else
        fit_ivars[fit_ivars .== 0] .= 1000^-2
    end

    n_params = size(first_guess_params, 2)

    curr_guess = copy(first_guess_params)
    new_params = copy(first_guess_params)
    curr_n_peaks = size(curr_guess, 1)
    v_hat = zeros(Float64, (n_params, curr_n_peaks))
    v_hat_cov = zeros(Float64, (n_params, n_params, curr_n_peaks))

    #(n_pix,n_params,n_peaks)
    M_vects = zeros(Float64,
        (size(fit_fluxes, 1), n_params, curr_n_peaks))

    comb_model_fluxes = zeros(size(all_rel_fluxes))
    param_offsets = zeros(Float64, size(curr_guess))

    flux_diffs = copy(fit_fluxes)

    for r_ind in 1:n_iter

        # CDF version
        model_fluxes_unit_height, dmodel_dmu, dmodel_dsig = gaussian_int_pdf(
            fit_inds, curr_guess, return_deriv = true)

        if r_ind == 1
            #then fit only for the heights using the first guesses on centers and widths
            if !use_first_guess_heights
                #don't use the first guess heights to define a best model
                #just use the first guess mean and widths to get the best
                #guess for height
                flux_diffs[:, :] .= fit_fluxes[:, :]
            else
                #use the initial provided heights to define best models to
                #subtract off, then get height updates
                model_fluxes = (curr_guess[:, 1] .* model_fluxes_unit_height)

                best_model_fluxes_unit_height = gaussian_int_pdf(
                    best_model_fit_inds, curr_guess)
                best_model_fluxes = (curr_guess[:, 1] .* best_model_fluxes_unit_height)
                best_model_fluxes[.!isfinite.(curr_guess[:, 1]), :] .= 0

                comb_model_fluxes .*= 0
                for j in 1:size(best_model_fit_inds, 1)
                    comb_model_fluxes[best_model_fit_inds[j, :]] .+= best_model_fluxes[j, :]
                end

                flux_diffs .= fit_fluxes .- model_fluxes'
                flux_diffs .-= comb_model_fluxes[fit_inds]'
                flux_diffs .+= best_model_fluxes[
                    :, offset_inds .+ (size(best_model_fluxes, 2) .÷ 2) .+ 1]'
            end

            #get best fit flux heights
            ivar_heights = (nansum((model_fluxes_unit_height .^ 2)' .* fit_ivars, 1))[1, :]
            mean_heights = nansum(model_fluxes_unit_height' .* flux_diffs .* fit_ivars, 1)[1, :] ./
                           ivar_heights
            #             @show mean_heights

            #     	    println(0," ",first_guess_params[:,1])

            if !use_first_guess_heights
                #then mean_heights gives the first guess heights
                first_guess_params[:, 1] .= max.(0.01, mean_heights)
            else
                #then mean_heights gives the offset heights to update by
                first_guess_params[:, 1] .= max.(0.01, first_guess_params[:, 1] .+ mean_heights)
            end
            curr_guess .= first_guess_params
            #     	    println(0," ",first_guess_params[:,1])
        end

        model_fluxes = (curr_guess[:, 1] .* model_fluxes_unit_height) .+ curr_guess[:, 4]

        best_model_fluxes_unit_height = gaussian_int_pdf(
            best_model_fit_inds, curr_guess)
        best_model_fluxes = (curr_guess[:, 1] .* best_model_fluxes_unit_height)

        best_model_fluxes[.!isfinite.(curr_guess[:, 1]), :] .= 0

        comb_model_fluxes .*= 0
        for j in 1:size(best_model_fit_inds, 1)
            comb_model_fluxes[best_model_fit_inds[j, :]] .+= best_model_fluxes[j, :]
        end

        # remove combined gaussians
        flux_diffs .= fit_fluxes .- comb_model_fluxes[fit_inds]'
        # add back in the individual gaussian for each peak so only neighbours are removed
        flux_diffs .+= best_model_fluxes[
            :, offset_inds .+ (size(best_model_fluxes, 2) .÷ 2) .+ 1]'
        flux_diffs .-= model_fluxes'

        dmodel_dmu .*= curr_guess[:, 1]
        dmodel_dsig .*= curr_guess[:, 1]

        # derivative matrix
        M_vects[:, 1, :] .= model_fluxes_unit_height'  # d flux/d height
        M_vects[:, 2, :] .= dmodel_dmu'      # d flux/d mu
        M_vects[:, 3, :] .= dmodel_dsig'     # d flux/d sigma
        M_vects[:, 4, :] .= 1.0              # d flux/d bias

        for j in 1:curr_n_peaks
            #preallocate curr_M_T_dot_V_inv
            curr_M_T_dot_V_inv = M_vects[:, :, j]' .* fit_ivars[:, j]'

            if (return_cov) & (r_ind == n_iter)
                # only invert matrix on the last iteration
                try
                    curr_V = inv(curr_M_T_dot_V_inv * M_vects[:, :, j])
                    v_hat[:, j] .= curr_V * (curr_M_T_dot_V_inv * flux_diffs[:, j])
                    v_hat_cov[:, :, j] .= curr_V
                catch
                    v_hat[:, j] .= 0.0
                    v_hat_cov[:, :, j] .= NaN
                end
            else
                try
                    v_hat[:, j] .= (curr_M_T_dot_V_inv * M_vects[:, :, j]) \
                                   (curr_M_T_dot_V_inv * flux_diffs[:, j])
                catch
                    v_hat[:, j] .= 0.0
                end
            end
        end

        #	println(r_ind," ",v_hat[1,:])
        #	println(r_ind," ",v_hat[2,:])
        #	println(r_ind," ",v_hat[3,:])
        #	println(r_ind," ",v_hat[4,:])

        new_params .= curr_guess .+ v_hat'
        new_params[:, 1] .= max.(new_params[:, 1], 0.01)
        new_params[:, 2] .= clamp.(new_params[:, 2],
            max.((n_offset + 1), first_guess_params[:, 2] .- max_center_move),
            min.(n_pixels - (n_offset + 1), first_guess_params[:, 2] .+ max_center_move))
        new_params[:, 3] .= clamp.(new_params[:, 3], min_widths, max_widths)

        curr_guess .= new_params
    end

    if !return_cov
        return curr_guess
    else
        return curr_guess, permutedims(v_hat_cov, (3, 1, 2))
    end
end

"""
Returns frequency values for input list of pixel positions

Arguments:
- `x`: pixel coordinates to be changed to frequency, size (n_pixels)
- `freqParams`: coefficients of frequency solution, size (order+1)
- `chipPolyParams`: coefficients of chip transformation, size (n_chips,order+1)
- `chipIndx`: chip index [1,2,3], type Int

Returns:
- `freq`: frequency values of x given current frequecny solution, size (n_pixels)
"""
function apply_freq_soln(x, freqParams, chipPolyParams, chipIndx)
    porder = size(freqParams, 1) - 1
    ximport = (x .- 1024) ./ 2048
    xt = transform_x_chips(ximport, chipPolyParams[chipIndx, :])
    Ax = positional_poly_mat(xt, porder = porder)
    freq = Ax * freqParams
    return freq
end

"""
Uses the FPI data, skyline wavelength solution & chip polynomial solution
to return estimates for the FPI parameters by fitting Gaussians to FPI peaks.

Returns:
- ``:
"""
function get_firstGuessParams(
        fluxfpi, ivarfpi, goodfpi, fiber_waveParams, fiber_chipPolyParams; bias_porder = 10)
    n_pixels = size(fluxfpi, 1)
    n_fibers = size(fluxfpi, 2)
    n_chips = size(fluxfpi, 3)

    x = collect(1:n_pixels) #same for all fibers
    x_boundaries = [x .- 0.5; x[end] + 0.5] #pixel edges

    fib_inds = 1:n_fibers
    chip_inds = 1:n_chips

    Ax_bias = positional_poly_mat(x, porder = bias_porder)

    #don't use bad pixels in fitting
    #    fluxfpi[.!goodfpi] = 0.0
    ivarfpi[.!goodfpi] = 0.0

    #output = (flux,mean,width,bias)
    n_params = 4

    #fit Gaussians and bias to FPI fluxes
    for chipIndx in chip_inds
        for fibIndx in fib_inds
            fit_initial_fpi_peaks(x, fluxfpi[fibIndx, :, chipIndx],
                ivarfpi[fibIndx, :, chipIndx])
            curr_guess_params = zeros(Float64, ())
            fit_gaussians(fluxfpi[fibIndx, :, chipIndx],
                ivarfpi[fibIndx, :, chipIndx] .^ -0.5,
                curr_guess_params,
                fit_inds, best_model_fit_inds, offset_inds;
                n_iter = 10, dmu = 0.001, dsig = 0.001, return_cov = false,
                use_first_guess_heights = false, max_center_move = 3,
                min_widths = 0.1, max_widths = 5.0)

            waves = apply_freq_soln(x,
                fiber_waveParams[fibIndx, :],
                fiber_chipPolyParams[fibIndx, :],
                chipIndx)
            freqs = 3.0e5 ./ waves
        end
    end

    #intrinsic width of FPI summed-Lorentzian profile
    gamma = 0.02 * peak_freq_spacing

    return gamma
end

"""
Returns the FPI fluxes, ivars, good_pixels for a given expid,sjd,tele.

Returns:
- `fluxOut`: fluxes, size (n_pixels,n_fibers,n_chips):
- `ivarOut`: inverse-variances, size (n_pixels,n_fibers,n_chips):
- `goodPixOut`: boolean of whether flux is good/usedul, type Boolean, size (n_pixels,n_fibers,n_chips):
"""
function get_fpi_data(expid, sjd, tele; ftype = "ar1Dcal", datapath = "")
    fname = datapath * "$(ftype)_$(tele)_$(sjd)_$(expid)_ARCLAMP.jld2"
    chips = ["a", "b", "c"]
    fluxOut = zeros(Float64, (2048, 300, 3))
    ivarOut = zeros(Float64, (2048, 300, 3))
    maskOut = zeros(Int, (2048, 300, 3))

    for chip in chips
        chipIndx = getChipIndx(chip)
        fluxOut[:, :, chipIndx] = load(fname, "flux_1d")
        ivarOut[:, :, chipIndx] = load(fname, "ivar_1d")
        maskOut[:, :, chipIndx] = load(fname, "mask_1d")
    end

    goodPixOut = ((maskOut .& bad_pix_bits) .== 0)

    return fluxOut, ivarOut, goodPixOut
end

"""
Fits Gaussians and Bias function to flux,ivar data

Arguments:
- `flux`: size (n_pixels)
- `ivar`: size (n_pixels)
- `coeffs_peak_ind_to_x`: size (n_params), coefficients to estimate peak locations from integers
- `coeffs_x_to_peak_ind`: size (n_params), coefficients to estimate peak integers from locations
"""
function get_initial_fpi_peaks(flux, ivar,
		coeffs_peak_ind_to_x, coeffs_x_to_peak_ind)
    #smooth the fluxes, then find local maxima
    #fit Gaussians to those local maxima
    #use the positions of those Gaussians to infer the positions of missing peaks
    #do a second pass on fitting 
    #return bias polynomial coefficients and peak (width,height,center)

    peak_to_x_func = Polynomial(coeffs_peak_ind_to_x)
    x_to_peak_func = Polynomial(coeffs_x_to_peak_ind)

    n_pixels = size(flux, 1)
    x = collect(1:n_pixels)

    good_ivars = (ivar .> 0)
#    if !any(good_ivars)
#        return [], zeros((0, 4)), zeros((0, 4, 4))
#    end

    sigma = 1.0 # smoothing length, pixels
    n_smooth_pix = max(5, round(Int, 3 * sigma) + 1)
    smooth_inds = range(start = -n_smooth_pix, stop = n_smooth_pix, step = 1)

    smooth_weights = exp.(-0.5 * (smooth_inds ./ sigma) .^ 2)
    smooth_weights ./= sum(smooth_weights)

    all_smooth_inds = (((smooth_inds' .+ x) .- 1 .+ 2048) .% 2048) .+ 1
    smoothed_fluxes = nansum(((flux .* good_ivars)[all_smooth_inds])' .* smooth_weights, 1)'
    smoothed_fluxes ./= nansum(good_ivars[all_smooth_inds]' .* smooth_weights, 1)'

    #    flux_thresh = nanzeropercentile(smoothed_fluxes,percent_vec=[70])[1] 
    flux_thresh = nanzeropercentile(smoothed_fluxes, percent_vec = [60])[1]

    # find local maxima and minima
    slopes = smoothed_fluxes[(begin + 1):end] .- smoothed_fluxes[begin:(end - 1)]
    local_max_inds = findall((slopes[1:(end - 1)] .>= 0) .& (slopes[2:end] .<= 0) .&
                             (.!((slopes[1:(end - 1)] .== 0) .& (slopes[2:end] .== 0))) .&
                             (smoothed_fluxes[2:(end - 1)] .>= flux_thresh)) .+ 1
    local_min_inds = findall((slopes[1:(end - 1)] .<= 0) .& (slopes[2:end] .>= 0) .&
                             (smoothed_fluxes[2:(end - 1)] .<= flux_thresh)) .+ 1

    poss_local_max_fluxes = smoothed_fluxes[local_max_inds]
    poss_local_min_fluxes = smoothed_fluxes[local_min_inds]
    med_min_fluxes = nanmedian(poss_local_min_fluxes)
    med_max_fluxes = nanmedian(poss_local_max_fluxes)

    local_min_waves = x[local_min_inds]
    local_max_waves = x[local_max_inds]

    # fit polynomial functions to min(Y) and max(Y) to define good peaks
    keep_min = (.!(isnan.(poss_local_min_fluxes))) .&
               (abs.(poss_local_min_fluxes ./ med_min_fluxes .- 1.0) .< 0.5)
    keep_max = (.!(isnan.(poss_local_max_fluxes))) .&
               (abs.(poss_local_max_fluxes ./ med_max_fluxes .- 1.0) .< 0.5)

    for r_ind in 1:2
        x_min, y_min = local_min_waves[keep_min], poss_local_min_fluxes[keep_min]
        x_max, y_max = local_max_waves[keep_max], poss_local_max_fluxes[keep_max]

        min_func = fit(x_min, y_min, 3)
        max_func = fit(x_max, y_max, 3)

        resids_min = poss_local_min_fluxes .- min_func.(local_min_waves)
        resids_max = poss_local_max_fluxes .- max_func.(local_max_waves)

        resid_summary_min = nanzeropercentile(resids_min[keep_min], percent_vec = [16, 50, 64])
        resid_summary_min = [resid_summary_min[2],
            resid_summary_min[2] - resid_summary_min[1],
            resid_summary_min[3] - resid_summary_min[2]]
        resid_summary_max = nanzeropercentile(resids_max[keep_max], percent_vec = [16, 50, 64])
        resid_summary_max = [resid_summary_max[2],
            resid_summary_max[2] - resid_summary_max[1],
            resid_summary_max[3] - resid_summary_max[2]]
        n_sigma = 5
        keep_min .= (resids_min .>=
                     resid_summary_min[1] - n_sigma * resid_summary_min[2]) .&
                    (resids_min .<= resid_summary_min[1] + n_sigma * resid_summary_min[3])
        keep_max .= (resids_max .>=
                     resid_summary_max[1] - n_sigma * resid_summary_max[2]) .&
                    (resids_max .<= resid_summary_max[1] + n_sigma * resid_summary_max[3])
    end

    x_min, y_min = local_min_waves[keep_min], poss_local_min_fluxes[keep_min]
    x_max, y_max = local_max_waves[keep_max], poss_local_max_fluxes[keep_max]

    min_func = fit(x_min, y_min, 3)
    max_func = fit(x_max, y_max, 3)

    resids_min = poss_local_min_fluxes .- min_func.(local_min_waves)
    resids_max = poss_local_max_fluxes .- max_func.(local_max_waves)

    med_min_fluxes = min_func.(local_max_waves)
    med_max_fluxes = max_func.(local_max_waves)

    #only keep peaks that are substantially above the minimum function
    relative_fluxes = (poss_local_max_fluxes .- med_min_fluxes) ./
                      (med_max_fluxes .- med_min_fluxes)

    n_offset = 15
    n_pad = 20

    good_max_inds = (abs.(relative_fluxes .- 1.0) .< 0.2) .& (local_max_waves .>= (n_offset + 1)) .&
                    (local_max_waves .<= n_pixels - (n_offset + 1))

    good_y_vals = local_max_waves[good_max_inds]

    peak_int_guess = ceil.(Int,round.(x_to_peak_func.(good_y_vals)))
    peak_x_guess = peak_to_x_func.(peak_int_guess)
    curr_offset = nanmedian(peak_x_guess .- good_y_vals)
    min_max_peak_ints = ceil.(Int,round.(x_to_peak_func.([1,2048])))
    peak_ints = collect(min_max_peak_ints[1]:min_max_peak_ints[2])
    good_y_vals = ceil.(Int,round.(peak_to_x_func.(peak_ints) .- curr_offset)) 

    good_max_inds = (good_y_vals .>= (n_offset + 1)) .&
                    (good_y_vals .<= n_pixels - (n_offset + 1))
    good_y_vals = good_y_vals[good_max_inds]

    #fit 1D gaussians to each identified peak, using offset_inds around each peak
    offset_inds = range(start = -4, stop = 4, step = 1)
    #use a larger number of pixels for removing the contribution from neighbouring peaks
    best_model_offset_inds = range(start = -n_offset, stop = n_offset, step = 1)

    fit_inds = good_y_vals .+ offset_inds'
    best_model_fit_inds = good_y_vals .+ best_model_offset_inds'

    fit_fluxes = flux[fit_inds]

    first_guess_params = zeros(Float64, size(fit_fluxes, 1), 4)
    first_guess_params[:, 1] .= max.(0.01, maximum(fit_fluxes, dims = 2)) # height
    first_guess_params[:, 2] .= good_y_vals # center position index for mu
    first_guess_params[:, 3] .= 1.3 # sigma
    #    first_guess_params[:, 3] .= 1.0 # sigma
    first_guess_params[:, 4] .= 0 # bias

    first_guess_params[:, 1] .*= (2 * π)^0.5 * first_guess_params[:, 3] # change to integrated height

    #pad the fluxes
    flux_pad = zeros(size(flux, 1) + 2 * n_pad)
    ivar_pad = zeros(size(flux, 1) + 2 * n_pad)
    flux_pad[(begin + n_pad):(end - n_pad)] .= flux
    ivar_pad[(begin + n_pad):(end - n_pad)] .= ivar
    #add on n_pad to the measured peak centers
    first_guess_params[:, 2] .+= n_pad

    #get first guess peak parameters
    new_params = fit_gauss_and_bias(flux_pad, ivar_pad, first_guess_params,
        fit_inds, best_model_fit_inds, offset_inds;
        n_iter = 10, dmu = 0.001, dsig = 0.001, return_cov = false,
        use_first_guess_heights = true, max_center_move = 3,
        min_widths = 0.1, max_widths = 3.0, n_pad = n_pad, n_offset = n_offset)
    #subtract off n_pad from the center positions
    new_params[:, 2] .-= n_pad

    # use the location of the current peaks to fill in any gaps for low-throughput fibers
    peak_locs = new_params[1:(end - 1), 2]
    peak_spacing = new_params[2:end, 2] .- new_params[1:(end - 1), 2]
    peak_space_spacing = peak_spacing[2:end] .- peak_spacing[1:(end - 1)]
    peak_space_spacing_locs = peak_locs[2:end]
    med_peak_space_spacing = nanmedian(peak_space_spacing ./
                                       (peak_locs[2:end] .- peak_locs[1:(end - 1)]))
    med_peak_spacing = nanmedian(peak_spacing .- peak_locs .* med_peak_space_spacing)
    peak_space_func = Polynomial([med_peak_spacing, med_peak_space_spacing])
    expect_peak_spacing = peak_space_func.(peak_locs)

    keep_space = abs.(peak_spacing ./ expect_peak_spacing .- 1.0) .< 0.5

    sigma = 5.0
    n_smooth_pix = round(3 * sigma) + 1
    smooth_inds = range(start = -n_smooth_pix, stop = n_smooth_pix, step = 1)

    smooth_weights = exp.(-0.5 * (smooth_inds ./ sigma) .^ 2)
    smooth_weights ./= sum(smooth_weights)

    peak_counts = range(start = 0,
        stop = size(new_params, 1) - 1,
        step = 1)
    all_smooth_inds = abs.(smooth_inds' .+ peak_counts)
    all_smooth_inds[all_smooth_inds .> size(new_params,
    1) .- 1] .= size(new_params, 1) .- 1 .+
                (size(new_params, 1) .- 1 .-
                 all_smooth_inds[all_smooth_inds .> size(new_params, 1) .- 1])
    all_smooth_inds = round.(Int, all_smooth_inds) .+ 1

    smoothed_biases = nansum(new_params[all_smooth_inds, 4]' .* smooth_weights, 1)'
    smoothed_widths = nansum(new_params[all_smooth_inds, 3]' .* smooth_weights, 1)'
    smoothed_heights = nansum(new_params[all_smooth_inds, 1]' .* smooth_weights, 1)'
    smoothed_width_locs = copy(new_params[:, 2])

    keep_widths = ones(Bool, size(new_params, 1))
    keep_heights = ones(Bool, size(new_params, 1))
    keep_biases = ones(Bool, size(new_params, 1))

    for r_ind in 1:2
        space_func = fit(peak_locs[keep_space], peak_spacing[keep_space], 3)
        resids = peak_spacing .- space_func.(peak_locs)
        resid_summary = nanzeropercentile(resids[keep_space], percent_vec = [16, 50, 64])
        resid_summary = [resid_summary[2],
            resid_summary[2] - resid_summary[1],
            resid_summary[3] - resid_summary[2]]
        n_sigma = 5
        keep_space .= (resids .>= resid_summary[1] - n_sigma * resid_summary[2]) .&
                      (resids .<= resid_summary[1] + n_sigma * resid_summary[3])

        resids = new_params[:, 3] .- smoothed_widths
        resid_summary = nanzeropercentile(resids[keep_widths], percent_vec = [16, 50, 64])
        resid_summary = [resid_summary[2],
            resid_summary[2] - resid_summary[1],
            resid_summary[3] - resid_summary[2]]
        n_sigma = 3
        keep_widths .= (resids .>= resid_summary[1] - n_sigma * resid_summary[2]) .&
                       (resids .<= resid_summary[1] + n_sigma * resid_summary[3])
        smoothed_widths .= nansum(
            ((new_params[:, 3] .* keep_widths)[all_smooth_inds])' .* smooth_weights, 1)'
        smoothed_widths ./= nansum(keep_widths[all_smooth_inds]' .* smooth_weights, 1)'

        resids = new_params[:, 1] .- smoothed_heights
        resid_summary = nanzeropercentile(resids[keep_heights], percent_vec = [16, 50, 64])
        resid_summary = [resid_summary[2],
            resid_summary[2] - resid_summary[1],
            resid_summary[3] - resid_summary[2]]
        n_sigma = 3
        keep_widths .= (resids .>= resid_summary[1] - n_sigma * resid_summary[2]) .&
                       (resids .<= resid_summary[1] + n_sigma * resid_summary[3])
        smoothed_heights .= nansum(
            ((new_params[:, 1] .* keep_heights)[all_smooth_inds])' .* smooth_weights, 1)'
        smoothed_heights ./= nansum(keep_heights[all_smooth_inds]' .* smooth_weights, 1)'

        resids = new_params[:, 4] .- smoothed_biases
        resid_summary = nanzeropercentile(resids[keep_biases], percent_vec = [16, 50, 64])
        resid_summary = [resid_summary[2],
            resid_summary[2] - resid_summary[1],
            resid_summary[3] - resid_summary[2]]
        n_sigma = 3
        keep_biases .= (resids .>= resid_summary[1] - n_sigma * resid_summary[2]) .&
                       (resids .<= resid_summary[1] + n_sigma * resid_summary[3])
        smoothed_biases .= nansum(
            ((new_params[:, 4] .* keep_biases)[all_smooth_inds])' .* smooth_weights, 1)'
        smoothed_biases ./= nansum(keep_biases[all_smooth_inds]' .* smooth_weights, 1)'
    end

    space_func = fit(peak_locs[keep_space], peak_spacing[keep_space], 3)

    space_func_eval = space_func.(peak_locs)
    int_spacing = floor.(Int, round.(peak_spacing ./ space_func_eval))

    peak_ints = zeros(Int, size(int_spacing, 1) + 1)
    peak_ints[2:end] .= cumsum(int_spacing, dims = 1)

    keep_peaks = ones(Bool, size(peak_ints, 1))
    for r_ind in 1:2
        peak_func = fit(peak_ints[keep_peaks], new_params[keep_peaks, 2], 3)

        resids = new_params[:, 2] .- peak_func.(peak_ints)
        resid_summary = nanzeropercentile(resids[keep_peaks], percent_vec = [16, 50, 64])
        resid_summary = [resid_summary[2],
            resid_summary[2] - resid_summary[1],
            resid_summary[3] - resid_summary[2]]
        n_sigma = 3
        keep_peaks .= (resids .>= resid_summary[1] - n_sigma * resid_summary[2]) .&
                      (resids .<= resid_summary[1] + n_sigma * resid_summary[3])
    end

#    peak_func = fit(peak_ints[keep_peaks], new_params[keep_peaks, 2], 3)
#
#    all_peak_ints = range(
#        start = max(peak_ints[1], 0) - 50, stop = min(1000, peak_ints[end]) + 50, step = 1)
#    all_peak_locs = peak_func.(all_peak_ints)
#    keep_peak_ints = (all_peak_locs .>= 1) .& (all_peak_locs .<= 2048)
#    all_peak_ints = all_peak_ints[keep_peak_ints]
#    all_peak_locs = all_peak_locs[keep_peak_ints]

    peak_ints = collect(min_max_peak_ints[1]:min_max_peak_ints[2])
    peak_int_guess = ceil.(Int,round.(x_to_peak_func.(new_params[:, 2])))
    peak_x_guess = peak_to_x_func.(peak_int_guess)
    curr_offset = nanmedian(peak_x_guess .- new_params[:, 2])
    all_peak_locs = peak_to_x_func.(peak_ints) .- curr_offset 
    keep_peak_ints = (all_peak_locs .>= 1) .& (all_peak_locs .<= 2048)
    all_peak_ints = peak_ints[keep_peak_ints]
    all_peak_locs = all_peak_locs[keep_peak_ints]

    offset_inds = range(start = -7, stop = 7, step = 1)
    fit_inds = floor.(Int, round.(all_peak_locs)) .+ offset_inds'
    best_model_fit_inds = floor.(Int, round.(all_peak_locs)) .+ best_model_offset_inds'

    first_guess_params = zeros(Float64, length(all_peak_locs), 4)
    first_guess_params[:, 2] .= all_peak_locs # center position index for mu

    for j in 1:size(first_guess_params, 1)
        match_ind = argmin(abs.(all_peak_locs[j] .- smoothed_width_locs))
        first_guess_params[j, 1] = max(0.01, smoothed_heights[match_ind])
        first_guess_params[j, 3] = smoothed_widths[match_ind]
        first_guess_params[j, 4] = smoothed_biases[match_ind]
    end

    #add on n_pad to the measured peak centers
    first_guess_params[:, 2] .+= n_pad

    #repeat fit to get measurements of all peaks
    new_params, v_hat_cov = fit_gauss_and_bias(flux_pad, ivar_pad, first_guess_params,
        fit_inds, best_model_fit_inds, offset_inds;
        n_iter = 10, dmu = 0.001, dsig = 0.001, return_cov = true,
        use_first_guess_heights = true, max_center_move = 3,
        min_widths = 0.1, max_widths = 3.0, n_pad = n_pad, n_offset = n_offset)

    #subtract off n_pad from the center positions
    new_params[:, 2] .-= n_pad

    return all_peak_ints, new_params, v_hat_cov
end

"""
Fits Gaussians and Bias function to flux,ivar data

Arguments:
- `flux`: size (n_pixels)
- `ivar`: size (n_pixels)
"""
function get_initial_arclamp_peaks(flux, ivar)
    #smooth the fluxes, then find local maxima
    #fit Gaussians to those local maxima

    n_pixels = size(flux, 1)
    x = collect(1:n_pixels)

    good_ivars = (ivar .> 0)
    if !any(good_ivars)
        return [], zeros((0, 4)), zeros((0, 4, 4))
    end

    sigma = 1.0 # smoothing length, pixels
    n_smooth_pix = max(5, round(Int, 3 * sigma) + 1)
    smooth_inds = range(start = -n_smooth_pix, stop = n_smooth_pix, step = 1)

    smooth_weights = exp.(-0.5 * (smooth_inds ./ sigma) .^ 2)
    smooth_weights ./= sum(smooth_weights)

    all_smooth_inds = (((smooth_inds' .+ x) .- 1 .+ 2048) .% 2048) .+ 1
    smoothed_fluxes = nansum(((flux .* good_ivars)[all_smooth_inds])' .* smooth_weights, 1)'
    smoothed_fluxes ./= nansum(good_ivars[all_smooth_inds]' .* smooth_weights, 1)'

    flux_thresh = nanzeropercentile(smoothed_fluxes, percent_vec = [95])[1]

    # find local maxima and minima
    slopes = smoothed_fluxes[(begin + 1):end] .- smoothed_fluxes[begin:(end - 1)]
    local_max_inds = findall((slopes[1:(end - 1)] .>= 0) .& (slopes[2:end] .<= 0) .&
                             (.!((slopes[1:(end - 1)] .== 0) .& (slopes[2:end] .== 0))) .&
                             (smoothed_fluxes[2:(end - 1)] .>= flux_thresh)) .+ 1
    local_min_inds = findall((slopes[1:(end - 1)] .<= 0) .& (slopes[2:end] .>= 0) .&
                             (smoothed_fluxes[2:(end - 1)] .<= flux_thresh)) .+ 1

    poss_local_max_fluxes = smoothed_fluxes[local_max_inds]
    poss_local_min_fluxes = smoothed_fluxes[local_min_inds]
    med_min_fluxes = nanmedian(poss_local_min_fluxes)
    med_max_fluxes = nanmedian(poss_local_max_fluxes)

    local_min_waves = x[local_min_inds]
    local_max_waves = x[local_max_inds]

    n_offset = 15
    n_pad = 20

    good_max_inds = ((local_max_waves .>= (n_offset + 1)) .&
                     (local_max_waves .<= n_pixels - (n_offset + 1)))
    good_y_vals = local_max_waves[good_max_inds]

    #fit 1D gaussians to each identified peak, using offset_inds around each peak
    offset_inds = range(start = -4, stop = 4, step = 1)
    #use a larger number of pixels for removing the contribution from neighbouring peaks
    best_model_offset_inds = range(start = -n_offset, stop = n_offset, step = 1)

    fit_inds = good_y_vals .+ offset_inds'
    best_model_fit_inds = good_y_vals .+ best_model_offset_inds'

    fit_fluxes = flux[fit_inds]

    first_guess_params = zeros(Float64, size(fit_fluxes, 1), 4)
    first_guess_params[:, 1] .= max.(0.01, maximum(fit_fluxes, dims = 2)) # height
    first_guess_params[:, 2] .= good_y_vals # center position index for mu
    first_guess_params[:, 3] .= 1.3 # sigma
    #    first_guess_params[:, 3] .= 1.0 # sigma
    first_guess_params[:, 4] .= 0 # bias

    first_guess_params[:, 1] .*= (2 * π)^0.5 * first_guess_params[:, 3] # change to integrated height

    #pad the fluxes
    flux_pad = zeros(size(flux, 1) + 2 * n_pad)
    ivar_pad = zeros(size(flux, 1) + 2 * n_pad)
    flux_pad[(begin + n_pad):(end - n_pad)] .= flux
    ivar_pad[(begin + n_pad):(end - n_pad)] .= ivar
    #add on n_pad to the measured peak centers
    first_guess_params[:, 2] .+= n_pad

    #get first guess peak parameters
    new_params, v_hat_cov = fit_gauss_and_bias(flux_pad, ivar_pad, first_guess_params,
        fit_inds, best_model_fit_inds, offset_inds;
        n_iter = 10, dmu = 0.001, dsig = 0.001, return_cov = true,
        use_first_guess_heights = true, max_center_move = 3,
        min_widths = 0.1, max_widths = 3.0, n_pad = n_pad, n_offset = n_offset)
    #subtract off n_pad from the center positions
    new_params[:, 2] .-= n_pad

    return collect(1:size(new_params, 1)), new_params, v_hat_cov
end

function get_and_save_fpi_peaks(fname)
    sname = split(fname, "_")
    tele, mjd, expid, chip = sname[(end - 4):(end - 1)]
    f = jldopen(fname, "r+")
    flux_1d = f["flux_1d"]
    ivar_1d = f["ivar_1d"]
    mask_1d = f["mask_1d"]
    extract_trace_centers = f["extract_trace_centers"]
    close(f)

    n_pixels = size(flux_1d, 1)
    n_fibers = size(flux_1d, 2)

    fib_inds = 1:n_fibers

    good_pix = ((mask_1d .& bad_pix_bits) .== 0)
    ivar_1d[.!good_pix] .= 0.0
    flux_1d[.!good_pix] .= 0.0
    x = collect(1:n_pixels)

    coeffs_peak_ind_to_x,
    coeffs_x_to_peak_ind = read_fpiPeakLoc_coeffs(
	tele, chip; data_path = "./data/")

    function get_peaks_partial(intup)
        flux_1d, ivar_1d, 
	coeffs_peak_ind_to_x, coeffs_x_to_peak_ind = intup
        get_initial_fpi_peaks(flux_1d, ivar_1d, 
	    coeffs_peak_ind_to_x, coeffs_x_to_peak_ind)
    end
    in2do = Iterators.zip(eachcol(flux_1d), eachcol(ivar_1d), 
		eachcol(coeffs_peak_ind_to_x), eachcol(coeffs_x_to_peak_ind))
    pout = map(get_peaks_partial, in2do)

    max_peaks = 0
    for j in 1:size(pout, 1)
        max_peaks = max(max_peaks, length(pout[j][1]))
    end
    n_params = size(pout[1][2], 2)
    fpi_trace_centers = zeros(
        Float64, max_peaks, n_fibers)
    fpi_line_mat = zeros(
        Float64, max_peaks, n_params, n_fibers)
    fpi_line_cov_mat = zeros(
        Float64, max_peaks, n_params, n_params, n_fibers)
    fill!(fpi_line_mat, NaN)
    fill!(fpi_line_cov_mat, NaN)

    x_pixels = collect(1:N_XPIX)
    for i in fib_inds
        if size(pout[i][1], 1) == 0
            continue
        end
        for peak_ind in 1:size(pout[i][1], 1)
            fpi_line_mat[peak_ind, :, i] .= pout[i][2][peak_ind, :]
            fpi_line_cov_mat[peak_ind, :, :, i] .= pout[i][3][peak_ind, :, :]
        end
	fpi_trace_centers[:, i] .= linear_interpolation(x_pixels, extract_trace_centers[:, i], extrapolation_bc = Line()).(fpi_line_mat[:,2,i])
    end

    outname = replace(replace(fname, "ar1Dcal" => "fpiPeaks"), "ar1D" => "fpiPeaks")
    f = h5open(outname, "w")

    # Write original data
    write(f, "fpi_line_mat", fpi_line_mat)
    attrs(f["fpi_line_mat"])["axis_1"] = "peak_index"
    attrs(f["fpi_line_mat"])["axis_2"] = "fit_info"
    attrs(f["fpi_line_mat"])["axis_3"] = "fibers"

    write(f, "fpi_line_cov_mat", fpi_line_cov_mat)
    attrs(f["fpi_line_cov_mat"])["axis_1"] = "peak_index"
    attrs(f["fpi_line_cov_mat"])["axis_2"] = "fit_info"
    attrs(f["fpi_line_cov_mat"])["axis_3"] = "fit_info"
    attrs(f["fpi_line_cov_mat"])["axis_4"] = "fibers"

    write(f, "fpi_line_trace_centers", fpi_trace_centers)
    attrs(f["fpi_line_trace_centers"])["axis_1"] = "peak_index"
    attrs(f["fpi_line_trace_centers"])["axis_2"] = "fibers"
    close(f)

    return outname
end

function get_and_save_arclamp_peaks(fname)
    sname = split(fname, "_")
    tele, mjd, chip, expid = sname[(end - 4):(end - 1)]
    f = jldopen(fname, "r+")
    flux_1d = f["flux_1d"]
    ivar_1d = f["ivar_1d"]
    mask_1d = f["mask_1d"]
    extract_trace_centers = f["extract_trace_centers"]
    close(f)

    n_pixels = size(flux_1d, 1)
    n_fibers = size(flux_1d, 2)

    fib_inds = 1:n_fibers

    good_pix = ((mask_1d .& bad_pix_bits) .== 0)
    ivar_1d[.!good_pix] .= 0.0
    flux_1d[.!good_pix] .= 0.0
    x = collect(1:n_pixels)

    function get_peaks_partial(intup)
        flux_1d, ivar_1d = intup
        get_initial_arclamp_peaks(flux_1d, ivar_1d)
    end
    in2do = Iterators.zip(eachcol(flux_1d), eachcol(ivar_1d))
    pout = map(get_peaks_partial, in2do)

    max_peaks = 0
    for i in 1:size(pout, 1)
        max_peaks = max(max_peaks, length(pout[i][1]))
    end

    n_params = size(pout[1][2], 2)
    fpi_trace_centers = zeros(
        Float64, max_peaks, n_fibers)
    fpi_line_mat = zeros(
        Float64, max_peaks, n_params, n_fibers)
    fpi_line_cov_mat = zeros(
        Float64, max_peaks, n_params, n_params, n_fibers)
    fill!(fpi_line_mat, NaN)
    fill!(fpi_line_cov_mat, NaN)

    x_pixels = collect(1:N_XPIX)
    for i in fib_inds
        if size(pout[i][1], 1) == 0
            continue
        end
        for peak_ind in 1:size(pout[i][1], 1)
            fpi_line_mat[peak_ind, :, i] .= pout[i][2][peak_ind, :]
            fpi_line_cov_mat[peak_ind, :, :, i] .= pout[i][3][peak_ind, :, :]
        end
	fpi_trace_centers[:, i] .= linear_interpolation(x_pixels, extract_trace_centers[:, i], extrapolation_bc = Line()).(fpi_line_mat[:,2,i])
    end

    outname = replace(replace(fname, "ar1Dcal" => "arclamp_peaks"), "ar1D" => "arclamp_peaks")
    f = h5open(outname, "w")

    # Write original data
    write(f, "arclamp_line_mat", fpi_line_mat)
    attrs(f["arclamp_line_mat"])["axis_1"] = "peak_index"
    attrs(f["arclamp_line_mat"])["axis_2"] = "fit_info"
    attrs(f["arclamp_line_mat"])["axis_3"] = "fibers"

    write(f, "arclamp_line_cov_mat", fpi_line_cov_mat)
    attrs(f["arclamp_line_cov_mat"])["axis_1"] = "peak_index"
    attrs(f["arclamp_line_cov_mat"])["axis_2"] = "fit_info"
    attrs(f["arclamp_line_cov_mat"])["axis_3"] = "fit_info"
    attrs(f["arclamp_line_cov_mat"])["axis_4"] = "fibers"

    write(f, "arclamp_line_trace_centers", fpi_trace_centers)
    attrs(f["arclamp_line_trace_centers"])["axis_1"] = "peak_index"
    attrs(f["arclamp_line_trace_centers"])["axis_2"] = "fibers"
    close(f)

    return outname
end

#using FITSIO, HDF5, FileIO, JLD2, Glob, CSV
#using DataFrames, EllipsisNotation, StatsBase
#using ParallelDataTransfer, ProgressMeter

#src_dir = "../"
#include(src_dir * "src/ar1D.jl")
#include(src_dir * "src/spectraInterpolation.jl")
#include(src_dir * "src/fileNameHandling.jl")
#include(src_dir * "src/utils.jl")
#fname = "../outdir//apred/60816/ar1Dcal_lco_60816_0008_a_ARCLAMP.h5"
#fname = "../outdir//apred/60816/ar1Dcal_apo_60816_0009_a_ARCLAMP.h5"

#get_and_save_fpi_peaks(fname)

#mjd = 60816
#tele = "apo"
#expids = ["0008","0009","0073","0074"]
#for expid in expids
#    for chip in ["a","b","c"]
#        fname = "../outdir//apred/60816/ar1Dcal_$(tele)_$(mjd)_$(expid)_$(chip)_ARCLAMP.h5"
#	println(fname)
#        get_and_save_fpi_peaks(fname)
#    end
#end


