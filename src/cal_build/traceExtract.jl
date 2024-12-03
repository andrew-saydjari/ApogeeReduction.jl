using Polynomials: fit
using SpecialFunctions: erf
using Interpolations: linear_interpolation, Line 


function trace_extract(image_data, ivar_image, tele, mjd, chip, expid; image_mask = nothing)
    noise_image = 1 ./ sqrt.(ivar_image)
    if isnothing(image_mask)
        image_mask = ones(Bool, size(image_data))
    end
    image_mask .= image_mask .& (ivar_image .> 0)
    n_center_cols = 100 # +/- n_cols to use from middle to sum to find peaks
    x_center = 1024 + 1

    x_inds = [-n_center_cols, n_center_cols + 1] .+ x_center

    # Cutout image in x direction
    cutout_fluxes = copy(image_data[
        (x_center - n_center_cols):(x_center + n_center_cols), begin:end])'
    cutout_errs = copy(noise_image[
        (x_center - n_center_cols):(x_center + n_center_cols), begin:end])'
    cutout_masks = copy(image_mask[
        (x_center - n_center_cols):(x_center + n_center_cols), begin:end])'

    # Mask bad pixels
    cutout_fluxes[.!cutout_masks] .= 0
    cutout_ivars = cutout_errs .^ -2
    cutout_ivars[.!cutout_masks] .= 0
    cutout_errs[.!cutout_masks] .= Inf

    #use median to mask large outlier fluxes from average
    med_flux = nanzeromedian(cutout_fluxes, 2)
    sigma_dists = (cutout_fluxes .- med_flux) .* (cutout_ivars .^ (0.5))
    good_sigmas = (abs.(sigma_dists) .< 5) .& (cutout_ivars .!= 0)

    cutout_fluxes[.!good_sigmas] .= 0
    cutout_ivars[.!good_sigmas] .= 0
    cutout_errs[.!good_sigmas] .= Inf

    y_vals = range(start = 1, stop = size(cutout_fluxes, 1), step = 1)

    comb_ivars = nansum(cutout_ivars, 2)
    good_ivars = comb_ivars .> 0

    # ivar-weight combine fluxes
    comb_errs = zeros(size(comb_ivars)) .+ Inf
    comb_fluxes = zeros(size(comb_ivars))
    comb_fluxes[good_ivars] = nansum((cutout_fluxes .* cutout_ivars), 2)[good_ivars] ./
                              comb_ivars[good_ivars]

    comb_errs[good_ivars] = comb_ivars[good_ivars] .^ (-0.5)

    sigma = 0.5 # smoothing length, pixels

    n_smooth_pix = max(5,round(Int, 3 * sigma) + 1)
    smooth_inds = range(start = -n_smooth_pix, stop = n_smooth_pix, step = 1)

    smooth_weights = exp.(-0.5 * (smooth_inds ./ sigma) .^ 2)
    smooth_weights ./= sum(smooth_weights)

    all_smooth_inds = (((smooth_inds' .+ y_vals) .- 1 .+ 2048) .% 2048) .+ 1
#    smoothed_fluxes = nansum(comb_fluxes[all_smooth_inds]' .* smooth_weights, 1)'
    smoothed_fluxes = nansum(((comb_fluxes .* good_ivars)[all_smooth_inds])' .* smooth_weights,1)'
    smoothed_fluxes ./= nansum(good_ivars[all_smooth_inds]' .* smooth_weights,1)'

    # find local maxima and minima
    slopes = smoothed_fluxes[(begin + 1):end] - smoothed_fluxes[begin:(end - 1)]
    local_max_inds = findall((slopes[1:(end - 1)] .>= 0) .& (slopes[2:end] .<= 0) .&
                             .!((slopes[1:(end - 1)] .== 0) .& (slopes[2:end] .== 0))) .+ 1
    local_min_inds = findall((slopes[1:(end - 1)] .<= 0) .& (slopes[2:end] .>= 0)) .+ 1

    poss_local_max_fluxes = smoothed_fluxes[local_max_inds]
    poss_local_min_fluxes = smoothed_fluxes[local_min_inds]
    med_min_fluxes = nanmedian(poss_local_min_fluxes)
    med_max_fluxes = nanmedian(poss_local_max_fluxes)

    local_min_waves = y_vals[local_min_inds]
    local_max_waves = y_vals[local_max_inds]

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

        resid_summary_min = quantile(resids_min[keep_min], [0.16, 0.5, 0.84])
        resid_summary_min = [resid_summary_min[2],
            resid_summary_min[2] - resid_summary_min[1],
            resid_summary_min[3] - resid_summary_min[2]]
        resid_summary_max = quantile(resids_max[keep_max], [0.16, 0.5, 0.84])
        resid_summary_max = [resid_summary_max[2],
            resid_summary_max[2] - resid_summary_max[1],
            resid_summary_max[3] - resid_summary_max[2]]
        n_sigma = 2
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
    good_max_inds = abs.(relative_fluxes .- 1.0) .< 0.5

    good_y_vals = local_max_waves[good_max_inds]

    #fit 1D gaussians to each identified peak, using offset_inds around each peak
    offset_inds = range(start = -4, stop = 4, step = 1)
    #use a larger number of pixels for removing the contribution from neighbouring peaks
    best_model_offset_inds = range(start = -10, stop = 10, step = 1)

    all_min_fluxes = min_func.(y_vals)
    all_max_fluxes = max_func.(y_vals)

    #fit scaled fluxes
#    all_rel_fluxes = (comb_fluxes) ./ (all_max_fluxes)
#    all_rel_errs = comb_errs ./ (all_max_fluxes)

    all_rel_fluxes = comb_fluxes
    all_rel_errs = comb_errs

    # first guess parameters
    fit_inds = good_y_vals .+ offset_inds'
    best_model_fit_inds = good_y_vals .+ best_model_offset_inds'
    fit_fluxes = copy(all_rel_fluxes[fit_inds])
    fit_errs = copy(all_rel_errs[fit_inds])
    fit_ivars = fit_errs .^ -2

    first_guess_params = zeros(Float64, size(fit_fluxes, 1), 3)
    first_guess_params[:, 1] .= max.(0.01, maximum(fit_fluxes, dims = 2)) # height
    first_guess_params[:, 2] .= good_y_vals # center position index for mu
    first_guess_params[:, 3] .= 0.7 # sigma

    first_guess_params[:, 1] .*= (2 * π)^0.5 * first_guess_params[:, 3] # change to integrated height

    curr_guess = copy(first_guess_params)

    dmu = 0.001
    dsig = 0.01
    group_inv(x, y) = mapslices(inv, x, dims = y)

    n_iter = 10
    for r_ind in 1:n_iter
        if r_ind > 1
            # use previous iteration parameters
#            fit_inds = floor.(Int, round.(curr_guess[:, 2])) .+ offset_inds'
#            best_model_fit_inds = floor.(Int, round.(curr_guess[:, 2])) .+
#                                  best_model_offset_inds'
            fit_fluxes = copy(all_rel_fluxes[fit_inds])
            fit_errs = copy(all_rel_errs[fit_inds])
            fit_ivars = fit_errs .^ -2
        end

        # CDF version
        model_fluxes_unit_height = 0.5 *
                                   erf.(((fit_inds .+ 0.5) .- curr_guess[:, 2]) ./
                                        (2^0.5 * curr_guess[:, 3])) .-
                                   0.5 *
                                   erf.(((fit_inds .- 0.5) .- curr_guess[:, 2]) ./
                                        (2^0.5 * curr_guess[:, 3]))

	if r_ind == 1
            model_fluxes = (curr_guess[:, 1] .* model_fluxes_unit_height)
            #get best fit flux height	    
	    flux_diffs = fit_fluxes 

            # derivative matrix
            M_vects = zeros(Float64,
                (size(model_fluxes, 1), size(model_fluxes, 2), 1))
	    M_vects[:,:,1] .= model_fluxes_unit_height #d flux/d height

            M_T_dot_V_inv = M_vects .* fit_ivars

            @einsum M_T_dot_V_inv_dot_M[n, j, k] := M_T_dot_V_inv[n, i, j] * M_vects[n, i, k]
            @einsum M_T_dot_V_inv_dot_y[n, j] := M_T_dot_V_inv[n, i, j] * flux_diffs[n, i]

            v_hat_cov = group_inv(M_T_dot_V_inv_dot_M, (2, 3))

            # update vector
            @einsum v_hat[n, i] := v_hat_cov[n, i, j] * M_T_dot_V_inv_dot_y[n, j]
	    first_guess_params[:,1] .= max.(0.01,v_hat[:,1])
	    curr_guess = copy(first_guess_params)

	end

        model_fluxes = (curr_guess[:, 1] .* model_fluxes_unit_height)

        model_fluxes_unit_height_above_mu = 0.5 *
                                            erf.(((fit_inds .+ 0.5) .-
                                                  (curr_guess[:, 2] .+ dmu)) ./
                                                 (2^0.5 * curr_guess[:, 3])) .-
                                            0.5 *
                                            erf.(((fit_inds .- 0.5) .-
                                                  (curr_guess[:, 2] .+ dmu)) ./
                                                 (2^0.5 * curr_guess[:, 3]))
        model_fluxes_unit_height_below_mu = 0.5 *
                                            erf.(((fit_inds .+ 0.5) .-
                                                  (curr_guess[:, 2] .- dmu)) ./
                                                 (2^0.5 * curr_guess[:, 3])) .-
                                            0.5 *
                                            erf.(((fit_inds .- 0.5) .-
                                                  (curr_guess[:, 2] .- dmu)) ./
                                                 (2^0.5 * curr_guess[:, 3]))
        model_fluxes_unit_height_above_sig = 0.5 *
                                             erf.(((fit_inds .+ 0.5) .- curr_guess[:, 2]) ./
                                                  (2^0.5 * (curr_guess[:, 3] .+ dsig))) .-
                                             0.5 *
                                             erf.(((fit_inds .- 0.5) .- curr_guess[:, 2]) ./
                                                  (2^0.5 * (curr_guess[:, 3] .+ dsig)))
        model_fluxes_unit_height_below_sig = 0.5 *
                                             erf.(((fit_inds .+ 0.5) .- curr_guess[:, 2]) ./
                                                  (2^0.5 * (curr_guess[:, 3] .- dsig))) .-
                                             0.5 *
                                             erf.(((fit_inds .- 0.5) .- curr_guess[:, 2]) ./
                                                  (2^0.5 * (curr_guess[:, 3] .- dsig)))

        best_model_fluxes_unit_height = 0.5 *
                                        erf.(((best_model_fit_inds .+ 0.5) .-
                                              curr_guess[:, 2]) ./
                                             (2^0.5 * curr_guess[:, 3])) .-
                                        0.5 *
                                        erf.(((best_model_fit_inds .- 0.5) .-
                                              curr_guess[:, 2]) ./
                                             (2^0.5 * curr_guess[:, 3]))
        best_model_fluxes = (curr_guess[:, 1] .* best_model_fluxes_unit_height)
        best_model_fluxes[.!isfinite.(curr_guess[:, 1]), :] .= 0

        comb_model_fluxes = zeros(size(all_rel_fluxes))
        for j in 1:size(best_model_fit_inds, 1)
            comb_model_fluxes[best_model_fit_inds[j, :]] .+= best_model_fluxes[j, :]
        end

        # remove combined gaussians
        fit_fluxes .-= comb_model_fluxes[fit_inds]
        # add back in the individual gaussian for each peak so only neighbours are removed
        fit_fluxes .+= best_model_fluxes[
            :, offset_inds .+ (size(best_model_fluxes, 2) .÷ 2) .+ 1]

        dmodel_dheight = model_fluxes_unit_height
        dmodel_dmu = (curr_guess[:, 1] .* (model_fluxes_unit_height_above_mu .-
                       model_fluxes_unit_height_below_mu)) ./ dmu
        dmodel_dsig = (curr_guess[:, 1] .* (model_fluxes_unit_height_above_sig .-
                        model_fluxes_unit_height_below_sig)) ./ dsig

        flux_diffs = fit_fluxes .- model_fluxes

        # derivative matrix
        M_vects = zeros(Float64,
            (size(model_fluxes, 1), size(model_fluxes, 2), size(first_guess_params, 2)))

        M_vects[:, :, 1] .= dmodel_dheight  # d flux/d height
        M_vects[:, :, 2] .= dmodel_dmu      # d flux/d mu
        M_vects[:, :, 3] .= dmodel_dsig     # d flux/d sigma

        M_T_dot_V_inv = M_vects .* fit_ivars

        @einsum M_T_dot_V_inv_dot_M[n, j, k] := M_T_dot_V_inv[n, i, j] * M_vects[n, i, k]
        @einsum M_T_dot_V_inv_dot_y[n, j] := M_T_dot_V_inv[n, i, j] * flux_diffs[n, i]

#        scales = M_T_dot_V_inv_dot_M[1:size(M_T_dot_V_inv_dot_M, 1), 1, 1]
#        v_hat_cov = group_inv(M_T_dot_V_inv_dot_M ./ scales, (2, 3)) ./ scales
        v_hat_cov = group_inv(M_T_dot_V_inv_dot_M, (2, 3))


        # update vector
        @einsum v_hat[n, i] := v_hat_cov[n, i, j] * M_T_dot_V_inv_dot_y[n, j]

        # restrict how far the update can move from the starting point
#        v_hat[:, 1] .= clamp.(v_hat[:, 1], -curr_guess[:, 1] * 0.9, curr_guess[:, 1] * 5)
#        v_hat[:, 2] .= clamp.(v_hat[:, 2], -2, 2)
#        v_hat[:, 3] .= clamp.(v_hat[:, 3], -curr_guess[:, 3] * 0.9, curr_guess[:, 3] * 5)

        new_params = curr_guess .+ v_hat
        new_params[:, 1] .= max.(new_params[:, 1], 0.01)
        new_params[:, 2] .= clamp.(new_params[:, 2], first_guess_params[:,2] .- 2, first_guess_params[:,2] .+ 2)
        new_params[:, 3] .= clamp.(new_params[:, 3], 0.5, 2.0)

        curr_guess .= new_params

    end

    new_params = curr_guess

    curr_best_params = copy(new_params)

    # use the location of the current peaks to fill in any gaps for low-throughput fibers
    peak_locs = new_params[1:(end - 1), 2]
    peak_spacing = new_params[2:end, 2] .- new_params[1:(end - 1), 2]
    med_peak_spacing = nanmedian(peak_spacing)

    good_peak_spacing = abs.(peak_spacing ./ med_peak_spacing .- 1.0) .< 0.5
    keep_space = copy(good_peak_spacing)

    sigma = 10.0
    n_smooth_pix = round(3*sigma)+1
    smooth_inds = range(start = -n_smooth_pix, stop = n_smooth_pix, step = 1)

    smooth_weights = exp.(-0.5 * (smooth_inds ./ sigma) .^ 2)
    smooth_weights ./= sum(smooth_weights)

    peak_counts = range(start = 0,
		        stop = size(new_params,1) - 1,
		        step = 1)
    all_smooth_inds = abs.(smooth_inds' .+ peak_counts)
    all_smooth_inds[all_smooth_inds .> size(new_params,1) .- 1] .= size(new_params,1) .- 1 .+ (size(new_params,1) .- 1 .- all_smooth_inds[all_smooth_inds .> size(new_params,1) .- 1])
    all_smooth_inds = round.(Int,all_smooth_inds) .+ 1

    smoothed_widths = nansum(new_params[all_smooth_inds,3]' .* smooth_weights, 1)'
    smoothed_heights = nansum(new_params[all_smooth_inds,1]' .* smooth_weights, 1)'
    smoothed_width_locs = copy(new_params[:,2])

    keep_widths = ones(Bool, size(new_params, 1))
    keep_heights = ones(Bool, size(new_params, 1))

    for r_ind in 1:2
        space_func = fit(peak_locs[keep_space], peak_spacing[keep_space], 3)
        resids = peak_spacing .- space_func.(peak_locs)
        resid_summary = quantile(resids[keep_space], [0.16, 0.5, 0.84])
        resid_summary = [resid_summary[2],
            resid_summary[2] - resid_summary[1],
            resid_summary[3] - resid_summary[2]]
        n_sigma = 5
        keep_space .= (resids .>= resid_summary[1] - n_sigma * resid_summary[2]) .&
                      (resids .<= resid_summary[1] + n_sigma * resid_summary[3])


	resids = new_params[:,3] .- smoothed_widths
	resid_summary = quantile(resids[keep_widths], [0.16, 0.5, 0.84])
        resid_summary = [resid_summary[2],
            resid_summary[2] - resid_summary[1],
            resid_summary[3] - resid_summary[2]]
        n_sigma = 3
        keep_widths .= (resids .>= resid_summary[1] - n_sigma * resid_summary[2]) .&
                       (resids .<= resid_summary[1] + n_sigma * resid_summary[3])
	smoothed_widths .= nansum(((new_params[:,3] .* keep_widths)[all_smooth_inds])' .* smooth_weights,1)'
	smoothed_widths ./= nansum(keep_widths[all_smooth_inds]' .* smooth_weights,1)'

	resids = new_params[:,2] .- smoothed_heights
	resid_summary = quantile(resids[keep_heights], [0.16, 0.5, 0.84])
        resid_summary = [resid_summary[2],
            resid_summary[2] - resid_summary[1],
            resid_summary[3] - resid_summary[2]]
        n_sigma = 3
        keep_widths .= (resids .>= resid_summary[1] - n_sigma * resid_summary[2]) .&
                       (resids .<= resid_summary[1] + n_sigma * resid_summary[3])
	smoothed_heights .= nansum(((new_params[:,1] .* keep_heights)[all_smooth_inds])' .* smooth_weights,1)'
	smoothed_heights ./= nansum(keep_heights[all_smooth_inds]' .* smooth_weights,1)'
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
        resid_summary = quantile(resids[keep_peaks], [0.16, 0.5, 0.84])
        resid_summary = [resid_summary[2],
            resid_summary[2] - resid_summary[1],
            resid_summary[3] - resid_summary[2]]
        n_sigma = 3
        keep_peaks .= (resids .>= resid_summary[1] - n_sigma * resid_summary[2]) .&
                      (resids .<= resid_summary[1] + n_sigma * resid_summary[3])
    end

    peak_func = fit(peak_ints[keep_peaks], new_params[keep_peaks, 2], 3)

    all_peak_ints = range(
        start = max(peak_ints[1], 0) - 50, stop = min(1000, peak_ints[end]) + 50, step = 1)
    all_peak_locs = peak_func.(all_peak_ints)
    keep_peak_ints = (all_peak_locs .> 1 + 10) .& (all_peak_locs .< 2048 - 10)
    all_peak_ints = all_peak_ints[keep_peak_ints]
    all_peak_locs = all_peak_locs[keep_peak_ints]

    missing_ints = zeros(Bool, size(all_peak_ints, 1))
    for (ind, peak_int) in enumerate(all_peak_ints)
        if !(peak_int in peak_ints)
            missing_ints[ind] = true
        end
    end

    missing_ints = all_peak_ints[missing_ints]
    
    peak_offsets = new_params[keep_peaks,2] .- peak_func.(peak_ints[keep_peaks])
    extrap = linear_interpolation(new_params[keep_peaks,2], peak_offsets, extrapolation_bc = Line())
    all_peak_locs .+= extrap.(all_peak_locs)

    #use more pixels around each peak, otherwise large-width peak fitting fails
    #(see weird discontinuities in height and width vs X)
    #This ONLY works because we are fitting all the peaks simultaneously 
    #and removing contamination from neighbouring fibers
    offset_inds = range(-7, 7, step = 1)

    curr_best_heights = zeros(Float64, size(all_peak_locs)) .+ 0.00001
    for j in 1:size(new_params,1)
        match_ind = argmin(abs.(all_peak_locs .- new_params[j,2]))
	curr_best_heights[match_ind] = new_params[j,1]
    end
    curr_best_widths = zeros(Float64, size(all_peak_locs))
    for j in 1:size(curr_best_widths,1)
        match_ind = argmin(abs.(all_peak_locs[j] .- smoothed_width_locs))
	curr_best_widths[j] = smoothed_widths[match_ind]
    end

    #using the full list of peaks, get the best-fit parameters
    all_rel_fluxes = comb_fluxes
    all_rel_errs = comb_errs

    # first guess parameters
    fit_inds = floor.(Int, round.(all_peak_locs)) .+ offset_inds'
    best_model_fit_inds = floor.(Int, round.(all_peak_locs)) .+ best_model_offset_inds'
    fit_fluxes = copy(all_rel_fluxes[fit_inds])
    fit_errs = copy(all_rel_errs[fit_inds])
    fit_ivars = fit_errs .^ -2

    first_guess_params = zeros(Float64, size(fit_fluxes, 1), 3)
    first_guess_params[:, 1] .= max.(0.01, curr_best_heights) # height
    first_guess_params[:, 2] .= all_peak_locs #center position index for mu
    first_guess_params[:, 3] .= curr_best_widths # sigma

    curr_guess = copy(first_guess_params)

    dmu = 0.001
    dsig = 0.001

    n_iter = 20
    for r_ind in 1:n_iter
        if r_ind > 1
            # use previous iteration parameters
#            fit_inds = floor.(Int, round.(curr_guess[:, 2])) .+ offset_inds'
#            best_model_fit_inds = floor.(Int, round.(curr_guess[:, 2])) .+
#                                  best_model_offset_inds'
            fit_fluxes = copy(all_rel_fluxes[fit_inds])
            fit_errs = copy(all_rel_errs[fit_inds])
            fit_ivars = fit_errs .^ -2
        end
        # CDF version
        model_fluxes_unit_height = 0.5 *
                                   erf.(((fit_inds .+ 0.5) .- curr_guess[:, 2]) ./
                                        (2^0.5 * curr_guess[:, 3])) .-
                                   0.5 *
                                   erf.(((fit_inds .- 0.5) .- curr_guess[:, 2]) ./
                                        (2^0.5 * curr_guess[:, 3]))

        if r_ind == 1
            model_fluxes = (curr_guess[:, 1] .* model_fluxes_unit_height)
            best_model_fluxes_unit_height = 0.5 *
                                            erf.(((best_model_fit_inds .+ 0.5) .-
                                                  curr_guess[:, 2]) ./
                                                 (2^0.5 * curr_guess[:, 3])) .-
                                            0.5 *
                                            erf.(((best_model_fit_inds .- 0.5) .-
                                                  curr_guess[:, 2]) ./
                                                 (2^0.5 * curr_guess[:, 3]))
            best_model_fluxes = (curr_guess[:, 1] .* best_model_fluxes_unit_height)
            best_model_fluxes[.!isfinite.(curr_guess[:, 1]), :] .= 0
            comb_model_fluxes = zeros(size(all_rel_fluxes))
            for j in 1:size(best_model_fit_inds, 1)
                comb_model_fluxes[best_model_fit_inds[j, :]] .+= best_model_fluxes[j, :]
            end

            flux_diffs = fit_fluxes .- model_fluxes
 	    flux_diffs .-= comb_model_fluxes[fit_inds]
            flux_diffs .+= best_model_fluxes[
                :, offset_inds .+ (size(best_model_fluxes, 2) .÷ 2) .+ 1]

            # derivative matrix
            M_vects = zeros(Float64,
                (size(model_fluxes, 1), size(model_fluxes, 2), 1))
	    M_vects[:,:,1] .= model_fluxes_unit_height #d flux/d height
            M_T_dot_V_inv = M_vects .* fit_ivars

            @einsum M_T_dot_V_inv_dot_M[n, j, k] := M_T_dot_V_inv[n, i, j] * M_vects[n, i, k]
            @einsum M_T_dot_V_inv_dot_y[n, j] := M_T_dot_V_inv[n, i, j] * flux_diffs[n, i]

            v_hat_cov = group_inv(M_T_dot_V_inv_dot_M, (2, 3))

            # update vector
            @einsum v_hat[n, i] := v_hat_cov[n, i, j] * M_T_dot_V_inv_dot_y[n, j]
       	    first_guess_params[:,1] .= max.(0.01,first_guess_params[:,1] .+ v_hat[:,1])
	    curr_guess = copy(first_guess_params)
	end

        model_fluxes = (curr_guess[:, 1] .* model_fluxes_unit_height)

        model_fluxes_unit_height_above_mu = 0.5 *
                                            erf.(((fit_inds .+ 0.5) .-
                                                  (curr_guess[:, 2] .+ dmu)) ./
                                                 (2^0.5 * curr_guess[:, 3])) .-
                                            0.5 *
                                            erf.(((fit_inds .- 0.5) .-
                                                  (curr_guess[:, 2] .+ dmu)) ./
                                                 (2^0.5 * curr_guess[:, 3]))
        model_fluxes_unit_height_below_mu = 0.5 *
                                            erf.(((fit_inds .+ 0.5) .-
                                                  (curr_guess[:, 2] .- dmu)) ./
                                                 (2^0.5 * curr_guess[:, 3])) .-
                                            0.5 *
                                            erf.(((fit_inds .- 0.5) .-
                                                  (curr_guess[:, 2] .- dmu)) ./
                                                 (2^0.5 * curr_guess[:, 3]))
        model_fluxes_unit_height_above_sig = 0.5 *
                                             erf.(((fit_inds .+ 0.5) .- curr_guess[:, 2]) ./
                                                  (2^0.5 * (curr_guess[:, 3] .+ dsig))) .-
                                             0.5 *
                                             erf.(((fit_inds .- 0.5) .- curr_guess[:, 2]) ./
                                                  (2^0.5 * (curr_guess[:, 3] .+ dsig)))
        model_fluxes_unit_height_below_sig = 0.5 *
                                             erf.(((fit_inds .+ 0.5) .- curr_guess[:, 2]) ./
                                                  (2^0.5 * (curr_guess[:, 3] .- dsig))) .-
                                             0.5 *
                                             erf.(((fit_inds .- 0.5) .- curr_guess[:, 2]) ./
                                                  (2^0.5 * (curr_guess[:, 3] .- dsig)))

        best_model_fluxes_unit_height = 0.5 *
                                        erf.(((best_model_fit_inds .+ 0.5) .-
                                              curr_guess[:, 2]) ./
                                             (2^0.5 * curr_guess[:, 3])) .-
                                        0.5 *
                                        erf.(((best_model_fit_inds .- 0.5) .-
                                              curr_guess[:, 2]) ./
                                             (2^0.5 * curr_guess[:, 3]))
        best_model_fluxes = (curr_guess[:, 1] .* best_model_fluxes_unit_height)
        best_model_fluxes[.!isfinite.(curr_guess[:, 1]), :] .= 0

        comb_model_fluxes = zeros(size(all_rel_fluxes))
        for j in 1:size(best_model_fit_inds, 1)
            comb_model_fluxes[best_model_fit_inds[j, :]] .+= best_model_fluxes[j, :]
        end

        # remove combined gaussians
        fit_fluxes .-= comb_model_fluxes[fit_inds]
        # add back in the individual gaussian for each peak so only neighbours are removed
        fit_fluxes .+= best_model_fluxes[
            :, offset_inds .+ (size(best_model_fluxes, 2) .÷ 2) .+ 1]

        dmodel_dheight = model_fluxes_unit_height
        dmodel_dmu = (curr_guess[:, 1] .* (model_fluxes_unit_height_above_mu .-
                       model_fluxes_unit_height_below_mu)) ./ dmu
        dmodel_dsig = (curr_guess[:, 1] .* (model_fluxes_unit_height_above_sig .-
                        model_fluxes_unit_height_below_sig)) ./ dsig

        flux_diffs = fit_fluxes .- model_fluxes

        # derivative matrix
        M_vects = zeros(Float64,
            (size(model_fluxes, 1), size(model_fluxes, 2), size(first_guess_params, 2)))

        M_vects[:, :, 1] .= dmodel_dheight  # d flux/d height
        M_vects[:, :, 2] .= dmodel_dmu      # d flux/d mu
        M_vects[:, :, 3] .= dmodel_dsig     # d flux/d sigma

        M_T_dot_V_inv = M_vects .* fit_ivars
        @einsum M_T_dot_V_inv_dot_M[n, j, k] := M_T_dot_V_inv[n, i, j] * M_vects[n, i, k]
        @einsum M_T_dot_V_inv_dot_y[n, j] := M_T_dot_V_inv[n, i, j] * flux_diffs[n, i]

#        scales = M_T_dot_V_inv_dot_M[1:size(M_T_dot_V_inv_dot_M, 1), 1, 1]
#        v_hat_cov = group_inv(M_T_dot_V_inv_dot_M ./ scales, (2, 3)) ./ scales
        v_hat_cov = group_inv(M_T_dot_V_inv_dot_M, (2, 3))

        # update vector
        @einsum v_hat[n, i] := v_hat_cov[n, i, j] * M_T_dot_V_inv_dot_y[n, j]

        # restrict how far the update can move from the starting point
#        v_hat[:, 1] .= clamp.(v_hat[:, 1], -curr_guess[:, 1] * 0.9, curr_guess[:, 1] * 5)
#        v_hat[:, 2] .= clamp.(v_hat[:, 2], -2, 2)
#        v_hat[:, 3] .= clamp.(v_hat[:, 3], -curr_guess[:, 3] * 0.9, curr_guess[:, 3] * 5)

        new_params = curr_guess .+ v_hat
        new_params[:, 1] .= max.(new_params[:, 1], 0.01)
        new_params[:, 2] .= clamp.(
            clamp.(new_params[:, 2], first_guess_params[:, 2] .- 3, first_guess_params[:, 2] .+ 3), 11, 2048 - 11)
        new_params[:, 3] .= clamp.(new_params[:, 3], 0.8 .* first_guess_params[:, 3], 1.2 .* first_guess_params[:, 2])

        curr_guess .= new_params
	println(size(v_hat_cov),v_hat_cov[1,:,:])
	println(flux_diffs[1,:])
	println(v_hat[1,:])
	println(M_T_dot_V_inv_dot_y[1,:])
    end

    new_params = curr_guess
    println("new_params 1",new_params[:,1])
    println("new_params 2",new_params[:,2])
    println("new_params 3",new_params[:,3])
    aksjhskjh

    #save the middle-of-detector best fit parameters for first guess
    best_fit_ave_params = copy(new_params)

    #work outwards from the middle of the detector
    #using the previous analyses to give good first guesses
    #when iterating over all the X pixels

    peak_locs = new_params[1:(end - 1), 2]
    peak_spacing = new_params[2:end, 2] - new_params[1:(end - 1), 2]
    med_peak_spacing = nanmedian(peak_spacing, 1)

    good_peak_spacing = abs.(peak_spacing ./ med_peak_spacing .- 1.0) .< 0.5
    med_flux = nanmedian(new_params[:, 1], 1)

    sigma = 10.0
    n_smooth_pix = round(3*sigma)+1
    smooth_inds = range(start = -n_smooth_pix, stop = n_smooth_pix, step = 1)

    smooth_weights = exp.(-0.5 * (smooth_inds ./ sigma) .^ 2)
    smooth_weights ./= sum(smooth_weights)

    peak_counts = range(start = 0,
		        stop = size(new_params,1) - 1,
		        step = 1)
    all_smooth_inds = abs.(smooth_inds' .+ peak_counts)
    all_smooth_inds[all_smooth_inds .> size(new_params,1) .- 1] .= size(new_params,1) .- 1 .+ (size(new_params,1) .- 1 .- all_smooth_inds[all_smooth_inds .> size(new_params,1) .- 1])
    all_smooth_inds = round.(Int,all_smooth_inds) .+ 1

    smoothed_widths = nansum(new_params[all_smooth_inds,3]' .* smooth_weights, 1)'
    smoothed_heights = nansum(new_params[all_smooth_inds,1]' .* smooth_weights, 1)'
    smoothed_width_locs = copy(new_params[:,2])
    println("width",smoothed_widths)
    println("height",smoothed_heights)

    keep_widths = ones(Bool, size(new_params, 1))
    keep_heights = ones(Bool, size(new_params, 1))

    for r_ind in 1:2

	resids = new_params[:,3] .- smoothed_widths
	resid_summary = quantile(resids[keep_widths], [0.16, 0.5, 0.84])
        resid_summary = [resid_summary[2],
            resid_summary[2] - resid_summary[1],
            resid_summary[3] - resid_summary[2]]
        n_sigma = 3
        keep_widths .= (resids .>= resid_summary[1] - n_sigma * resid_summary[2]) .&
                       (resids .<= resid_summary[1] + n_sigma * resid_summary[3])
	smoothed_widths .= nansum(((new_params[:,3] .* keep_widths)[all_smooth_inds])' .* smooth_weights,1)'
	smoothed_widths ./= nansum(keep_widths[all_smooth_inds]' .* smooth_weights,1)'

	resids = new_params[:,2] .- smoothed_heights
	resid_summary = quantile(resids[keep_heights], [0.16, 0.5, 0.84])
        resid_summary = [resid_summary[2],
            resid_summary[2] - resid_summary[1],
            resid_summary[3] - resid_summary[2]]
        n_sigma = 3
        keep_widths .= (resids .>= resid_summary[1] - n_sigma * resid_summary[2]) .&
                       (resids .<= resid_summary[1] + n_sigma * resid_summary[3])
	smoothed_heights .= nansum(((new_params[:,1] .* keep_heights)[all_smooth_inds])' .* smooth_weights,1)'
	smoothed_heights ./= nansum(keep_heights[all_smooth_inds]' .* smooth_weights,1)'
    end

    println("width",smoothed_widths)
    println("height",smoothed_heights)

    #remove the edge possible peaks if they have no throughput, because they likely don't exist
    good_throughput_fibers = (best_fit_ave_params[:, 1] ./ med_flux) .> 0.1
    low_throughput_fibers = findall(.!good_throughput_fibers)
    if size(low_throughput_fibers, 1) > 0
        if low_throughput_fibers[1] == 1
            #then there is a truncation applied to the left
            curr_cut_int = low_throughput_fibers[1]
            j = 2
            while j < size(low_throughput_fibers, 1) + 1
                if curr_cut_int + 1 == low_throughput_fibers[j]
                    curr_cut_int = low_throughput_fibers[j]
                else
                    break
                end
                j += 1
            end
            # println("Left Cut Integer ",curr_cut_int)
            left_cut_ind = curr_cut_int + 1
        else
            left_cut_ind = 1
        end

        if low_throughput_fibers[end] == size(new_params, 1)
            #then there is a truncation applied to the right
            curr_cut_int = low_throughput_fibers[end]
            j = size(low_throughput_fibers, 1) - 1
            while j > 1
                if curr_cut_int - 1 == low_throughput_fibers[j]
                    curr_cut_int = low_throughput_fibers[j]
                else
                    break
                end
                j -= 1
            end
            # println("Right Cut Integer ",curr_cut_int)
            right_cut_ind = curr_cut_int - 1
        else
            right_cut_ind = size(new_params, 1)
        end
    end

    best_fit_ave_params = best_fit_ave_params[left_cut_ind:right_cut_ind, :]
        
    med_flux = nanmedian(best_fit_ave_params[: ,1], 1)
    good_throughput_fibers = (best_fit_ave_params[:, 1] ./ med_flux) .> 0.1
    low_throughput_fibers = findall(.!good_throughput_fibers)

    x_inds = range(1, 2048, step = 1)
    final_param_outputs = zeros((
        size(x_inds, 1), size(best_fit_ave_params, 1), size(best_fit_ave_params, 2)))
    final_param_output_covs = zeros((size(x_inds, 1), size(best_fit_ave_params, 1),
        size(best_fit_ave_params, 2), size(best_fit_ave_params, 2)))

    best_fit_ave_params = best_fit_ave_params[good_throughput_fibers,:]

    println("Final number of peaks:", size(best_fit_ave_params))

    curr_best_widths = zeros(Float64, size(best_fit_ave_params,1))
    for j in 1:size(curr_best_widths,1)
        match_ind = argmin(abs.(best_fit_ave_params[j,2] .- smoothed_width_locs))
	curr_best_widths[j] = smoothed_widths[match_ind]
    end
    best_fit_ave_params[:,3] .= curr_best_widths

    #work from the center outward to use constraints from previous analysis 
    sorted_x_inds = x_inds[sortperm(abs.(x_inds .- 1025))]

    param_outputs = zeros((
        size(x_inds, 1), size(best_fit_ave_params, 1), size(best_fit_ave_params, 2)))
    param_output_covs = zeros((size(x_inds, 1), size(best_fit_ave_params, 1),
        size(best_fit_ave_params, 2), size(best_fit_ave_params, 2)))

    all_rel_fluxes_mat = copy(image_data)
    all_rel_errs_mat = copy(noise_image)
    all_rel_masks_mat = copy(image_mask)

    all_rel_fluxes_mat[.!all_rel_masks_mat] .= 0
    all_rel_errs_mat[.!all_rel_masks_mat] .= Inf

    all_rel_ivars_mat = 1 ./ (all_rel_errs_mat .^ 2)
    all_rel_ivars_mat[.!all_rel_masks_mat] .= 0
    all_rel_ivars_mat[all_rel_ivars_mat .== 0] .= 1e-5

    n_use = 200
    @showprogress for (ind, x_ind) in enumerate(sorted_x_inds)
        #use previous analyses to constrain first guesses
	
        if abs(x_ind - 1025) < 100
	    n_use = 100
        elseif abs(x_ind - 1025) < 200
	    n_use = 200
        elseif abs(x_ind - 1025) < 300
	    n_use = 300
	else
	    n_use = 500
	end

        if abs(x_ind - 1025) < 0.5 * n_use
            first_guess_params = copy(best_fit_ave_params)
        elseif x_ind > 1025
            first_guess_params = nanmedian(param_outputs[(x_ind - (n_use + 1)):(x_ind - 1), :, :], 1)[
                1, :, :]
	    param_slopes = nanmedian(diff(param_outputs[(x_ind - (n_use + 1)):(x_ind - 1), :, :],1),1)
	    first_guess_params .+= 0.5 * n_use .* param_slopes
        elseif x_ind < 1025
            first_guess_params = nanmedian(param_outputs[(x_ind + 1):(x_ind + (n_use + 1)), :, :], 1)[
                1, :, :]
	    param_slopes = nanmedian(diff(param_outputs[(x_ind + 1):(x_ind + (n_use + 1)), :, :], 1), 1)
	    first_guess_params .+= -0.5 * n_use .* param_slopes
        end

	first_guess_params[:,3] = best_fit_ave_params[:,3]*nanmedian(first_guess_params[:,3]/best_fit_ave_params[:,3])
	first_guess_params[:,1] = best_fit_ave_params[:,1]*nanmedian(first_guess_params[:,1]/best_fit_ave_params[:,1])

        all_rel_fluxes = all_rel_fluxes_mat[x_ind, :]
        all_rel_errs = all_rel_ivars_mat[x_ind, :]
        all_rel_ivars = all_rel_ivars_mat[x_ind, :]
        all_rel_masks = all_rel_errs_mat[x_ind, :]

        # first guess parameters
        fit_inds = floor.(Int, round.(first_guess_params[:, 2])) .+ offset_inds'
        best_model_fit_inds = floor.(Int, round.(first_guess_params[:, 2])) .+
                              best_model_offset_inds'
        fit_fluxes = copy(all_rel_fluxes[fit_inds])
        fit_errs = copy(all_rel_errs[fit_inds])
        fit_ivars = fit_errs .^ -2

        curr_guess = copy(first_guess_params)

        dmu = 0.001
        dsig = 0.001

        n_iter = 30
        n_iter = 10
        for r_ind in 1:n_iter
            if r_ind > 1
                # use previous iteration parameters
 #               fit_inds = floor.(Int, round.(curr_guess[:, 2])) .+ offset_inds'
 #               best_model_fit_inds = floor.(Int, round.(curr_guess[:, 2])) .+
 #                                     best_model_offset_inds'
                fit_fluxes = copy(all_rel_fluxes[fit_inds])
                fit_errs = copy(all_rel_errs[fit_inds])
                fit_ivars = fit_errs .^ -2
            end

            fit_ivars[fit_ivars .== 0] .= 1e-5

            # CDF version
            model_fluxes_unit_height = 0.5 *
                                       erf.(((fit_inds .+ 0.5) .- curr_guess[:, 2]) ./
                                            (2^0.5 * curr_guess[:, 3])) .-
                                       0.5 *
                                       erf.(((fit_inds .- 0.5) .- curr_guess[:, 2]) ./
                                            (2^0.5 * curr_guess[:, 3]))


            if r_ind == 1
                model_fluxes = (curr_guess[:, 1] .* model_fluxes_unit_height)
                best_model_fluxes_unit_height = 0.5 *
                                                erf.(((best_model_fit_inds .+ 0.5) .-
                                                      curr_guess[:, 2]) ./
                                                     (2^0.5 * curr_guess[:, 3])) .-
                                                0.5 *
                                                erf.(((best_model_fit_inds .- 0.5) .-
                                                      curr_guess[:, 2]) ./
                                                     (2^0.5 * curr_guess[:, 3]))
                best_model_fluxes = (curr_guess[:, 1] .* best_model_fluxes_unit_height)
                best_model_fluxes[.!isfinite.(curr_guess[:, 1]), :] .= 0

                comb_model_fluxes = zeros(size(all_rel_fluxes))
                for j in 1:size(best_model_fit_inds, 1)
                    comb_model_fluxes[best_model_fit_inds[j, :]] .+= best_model_fluxes[j, :]
                end

                flux_diffs = fit_fluxes .- model_fluxes
		flux_diffs .-= comb_model_fluxes[fit_inds]
                flux_diffs .+= best_model_fluxes[
                    :, offset_inds .+ (size(best_model_fluxes, 2) .÷ 2) .+ 1]

                # derivative matrix
                M_vects = zeros(Float64,
                    (size(model_fluxes, 1), size(model_fluxes, 2), 1))
     	        M_vects[:,:,1] .= model_fluxes_unit_height #d flux/d height
                M_T_dot_V_inv = M_vects .* fit_ivars

                @einsum M_T_dot_V_inv_dot_M[n, j, k] := M_T_dot_V_inv[n, i, j] * M_vects[n, i, k]
                @einsum M_T_dot_V_inv_dot_y[n, j] := M_T_dot_V_inv[n, i, j] * flux_diffs[n, i]

                v_hat_cov = group_inv(M_T_dot_V_inv_dot_M, (2, 3))

                # update vector
                @einsum v_hat[n, i] := v_hat_cov[n, i, j] * M_T_dot_V_inv_dot_y[n, j]
       	        first_guess_params[:,1] .= max.(0.01,first_guess_params[:,1] .+ v_hat[:,1])
     	        curr_guess = copy(first_guess_params)
            end

            model_fluxes = (curr_guess[:, 1] .* model_fluxes_unit_height)
            model_fluxes_unit_height_above_mu = 0.5 *
                                                erf.(((fit_inds .+ 0.5) .-
                                                      (curr_guess[:, 2] .+ dmu)) ./
                                                     (2^0.5 * curr_guess[:, 3])) .-
                                                0.5 *
                                                erf.(((fit_inds .- 0.5) .-
                                                      (curr_guess[:, 2] .+ dmu)) ./
                                                     (2^0.5 * curr_guess[:, 3]))
            model_fluxes_unit_height_below_mu = 0.5 *
                                                erf.(((fit_inds .+ 0.5) .-
                                                      (curr_guess[:, 2] .- dmu)) ./
                                                     (2^0.5 * curr_guess[:, 3])) .-
                                                0.5 *
                                                erf.(((fit_inds .- 0.5) .-
                                                      (curr_guess[:, 2] .- dmu)) ./
                                                     (2^0.5 * curr_guess[:, 3]))
            model_fluxes_unit_height_above_sig = 0.5 *
                                                 erf.(((fit_inds .+ 0.5) .-
                                                       curr_guess[:, 2]) ./ (2^0.5 *
                                                       (curr_guess[:, 3] .+ dsig))) .-
                                                 0.5 *
                                                 erf.(((fit_inds .- 0.5) .-
                                                       curr_guess[:, 2]) ./
                                                      (2^0.5 * (curr_guess[:, 3] .+ dsig)))
            model_fluxes_unit_height_below_sig = 0.5 *
                                                 erf.(((fit_inds .+ 0.5) .-
                                                       curr_guess[:, 2]) ./ (2^0.5 *
                                                       (curr_guess[:, 3] .- dsig))) .-
                                                 0.5 *
                                                 erf.(((fit_inds .- 0.5) .-
                                                       curr_guess[:, 2]) ./
                                                      (2^0.5 * (curr_guess[:, 3] .- dsig)))

            best_model_fluxes_unit_height = 0.5 *
                                            erf.(((best_model_fit_inds .+ 0.5) .-
                                                  curr_guess[:, 2]) ./
                                                 (2^0.5 * curr_guess[:, 3])) .-
                                            0.5 *
                                            erf.(((best_model_fit_inds .- 0.5) .-
                                                  curr_guess[:, 2]) ./
                                                 (2^0.5 * curr_guess[:, 3]))
            best_model_fluxes = (curr_guess[:, 1] .* best_model_fluxes_unit_height)
            best_model_fluxes[.!isfinite.(curr_guess[:, 1]), :] .= 0

            comb_model_fluxes = zeros(size(all_rel_fluxes))
            for j in 1:size(best_model_fit_inds, 1)
                comb_model_fluxes[best_model_fit_inds[j, :]] .+= best_model_fluxes[j, :]
            end

            # remove combined gaussians
            fit_fluxes .-= comb_model_fluxes[fit_inds]
            # add back in the individual gaussian for each peak so only neighbours are removed
            fit_fluxes .+= best_model_fluxes[
                :, offset_inds .+ (size(best_model_fluxes, 2) .÷ 2) .+ 1]

            dmodel_dheight = model_fluxes_unit_height
            dmodel_dmu = (curr_guess[:, 1] .* (model_fluxes_unit_height_above_mu .-
                           model_fluxes_unit_height_below_mu)) ./ dmu
            dmodel_dsig = (curr_guess[:, 1] .* (model_fluxes_unit_height_above_sig .-
                            model_fluxes_unit_height_below_sig)) ./ dsig

            flux_diffs = fit_fluxes .- model_fluxes

            # derivative matrix
            M_vects = zeros(Float64,
                (size(model_fluxes, 1), size(model_fluxes, 2), size(first_guess_params, 2)))

            M_vects[:, :, 1] .= dmodel_dheight  # d flux/d height
            M_vects[:, :, 2] .= dmodel_dmu      # d flux/d mu
            M_vects[:, :, 3] .= dmodel_dsig     # d flux/d sigma

            M_T_dot_V_inv = M_vects .* fit_ivars
            @einsum M_T_dot_V_inv_dot_M[n, j, k] := M_T_dot_V_inv[n, i, j] *
                                                    M_vects[n, i, k]
            @einsum M_T_dot_V_inv_dot_y[n, j] := M_T_dot_V_inv[n, i, j] * flux_diffs[n, i]

#            scales = M_T_dot_V_inv_dot_M[1:size(M_T_dot_V_inv_dot_M, 1), 1, 1]
#            v_hat_cov = group_inv(M_T_dot_V_inv_dot_M ./ scales, (2, 3)) ./ scales
            v_hat_cov = group_inv(M_T_dot_V_inv_dot_M, (2, 3))

            # update vector
            @einsum v_hat[n, i] := v_hat_cov[n, i, j] * M_T_dot_V_inv_dot_y[n, j]

            # restrict how far the update can move from the starting point
#            v_hat[:, 1] .= clamp.(
#                v_hat[:, 1], -curr_guess[:, 1] * 0.9, curr_guess[:, 1] * 5)
#            v_hat[:, 2] .= clamp.(v_hat[:, 2], -2, 2)
#            v_hat[:, 3] .= clamp.(
#                v_hat[:, 3], -curr_guess[:, 3] * 0.9, curr_guess[:, 3] * 5)

            new_params = curr_guess .+ v_hat
            new_params[:, 1] .= max.(
                new_params[:, 1], 0.01)
            new_params[:, 2] .= clamp.(
                clamp.(new_params[:, 2], first_guess_params[:, 2] .- 1,
                    first_guess_params[:, 2] .+ 1),
                11,
                2048 - 11)
            new_params[:, 3] .= clamp.(new_params[:, 3], first_guess_params[:, 3] .* 0.5, first_guess_params[:, 3] .* 2.0)

            curr_guess .= new_params
            param_outputs[x_ind, :, :] .= new_params
            param_output_covs[x_ind, :, :, :] .= v_hat_cov
        end
    end

    final_param_outputs[:,good_throughput_fibers] .= param_outputs
    final_param_output_covs[:,good_throughput_fibers] .= param_output_covs 
	    
    #extrapolate from nearby good fibers for low throughput ones
    good_throughput_fibers = findall(good_throughput_fibers)
    final_param_outputs[:,low_throughput_fibers,1] .= 0
    for j in 1:size(low_throughput_fibers)
        nearest_inds = argsort(abs.(low_throughput_fibers[j] .- good_throughput_fibers))[[1,2]]
        nearest_inds = good_throughput_fibers[nearest_inds]
        param_slopes = (diff(final_param_outputs[:,nearest_inds],2)[:,1]) ./ (nearest_inds[2] - nearest_inds[1])
        final_param_outputs[:,low_throughput_fibers[j],[2,3]] .= final_param_outputs[:,nearest_inds[1],[2,3]] \
                                                                 .+ param_slopes[:,[2,3]] .* (low_throughput_fibers[j]-nearest_inds[1])
    end

    param_outputs = final_param_outputs
    param_output_covs = final_param_output_covs

    # serialize(dirNamePlots * "param_outputs_$(tele)_$(mjd)_$(expid)_$(chip).jld2", param_outputs)
    return param_outputs
end
