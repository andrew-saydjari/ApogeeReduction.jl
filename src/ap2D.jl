# Handling the 2D image

import Pkg
Pkg.add("Polynomials")
Pkg.add("SpecialFunctions")
Pkg.add("Einsum")
#Pkg.add(url="https://github.com/jishnub/PolyFit.jl")
#using PolyFit
#Pkg.instantiate();
#Pkg.precompile();

using Polynomials, JLD2, StatsBase, LinearAlgebra
using ProgressMeter, ArgParse, SlackThreads
using SpecialFunctions
using Einsum
using Polynomials: fit

src_dir = "./"
include(src_dir * "/utils.jl")
include(src_dir * "/plotutils.jl")

#Pkg.add(url="https://github.com/jishnub/PolyFit.jl")
#using PolyFit


dirNamePlots = "../plots/"
if !ispath(dirNamePlots)
    mkpath(dirNamePlots)
end


function gauss_fit_linear(fluxes,flux_errs,first_guess_params)
    
    return nothing
end

function trace_extract(image_data,ivar_image,tele,mjd,chip,expid;image_mask=nothing)
    # Last editted by Kevin McKinnon on Sept 10, 2024 
    # image_data has shape (npix_x,npix_y)
    # ivar_image has shape (npix_x,npix_y)
    # where X is the wavelength dimension
	        
    plotter = true
    noise_image = ivar_image .^(-0.5)

    if isnothing(image_mask)
        image_mask = ones(Bool, size(image_data))
    end

    n_center_cols = 100 # +/- n_cols to use from middle to sum to find peaks
    x_center = 1024 + 1
   
    x_inds = [-n_center_cols, n_center_cols + 1] .+ x_center
   
    # Cutout image in x direction
    cutout_fluxes = copy(image_data[x_center-n_center_cols:x_center+n_center_cols,begin:end])'
    cutout_errs = copy(noise_image[x_center-n_center_cols:x_center+n_center_cols,begin:end])'
    cutout_masks = copy(image_mask[x_center-n_center_cols:x_center+n_center_cols,begin:end])'

    # Mask bad pixels
    cutout_fluxes[.!cutout_masks] .= 0
    cutout_ivars = cutout_errs .^ -2
    cutout_ivars[.!cutout_masks] .= 0
    cutout_errs[.!cutout_masks] .= Inf

    y_vals = range(start = 1, stop = size(cutout_fluxes, 1), step = 1)

    comb_ivars = nansum(cutout_ivars, 2)
    good_ivars = comb_ivars .> 0

    # ivar-weight combine fluxes
    comb_errs = zeros(size(comb_ivars)) .+ Inf
    comb_fluxes = zeros(size(comb_ivars))
    comb_fluxes[good_ivars] = nansum((cutout_fluxes .* cutout_ivars), 2)[good_ivars] ./ comb_ivars[good_ivars]

    comb_errs[good_ivars] = comb_ivars[good_ivars] .^ (-0.5)

    sigma = 1.0 # smoothing length, pixels

    n_smooth_pix = round(Int, 3 * sigma) + 1
    smooth_inds = range(start = -n_smooth_pix, stop = n_smooth_pix, step = 1) .+ 1

    smooth_weights = exp.(-0.5 * (smooth_inds ./ sigma) .^ 2)
    smooth_weights ./= sum(smooth_weights)

    all_smooth_inds = (((smooth_inds' .+ y_vals) .- 1 .+ 2048) .% 2048) .+ 1
    smoothed_fluxes = nansum(comb_fluxes[all_smooth_inds]' .* smooth_weights, 1)'

    if plotter

        keep = (abs.(y_vals .- 0)) .< 100

        # Figures for QA
        fig = plt.figure(figsize = (10, 5), dpi = 300)
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(y_vals[keep],smoothed_fluxes[keep],lw=1,color="C0")
        ax.errorbar(y_vals[keep],comb_fluxes[keep],yerr=comb_errs[keep],lw=1,color="C1")
        ttl = plt.text(
            0.5, 1.01, "Tele: $(tele), MJD: $(mjd), Chip: $(chip), Expid: $(expid)",
            ha = "center", va = "bottom", transform = ax.transAxes)
	imagePath = dirNamePlots * 
			"smoothed_fluxes_$(tele)_$(mjd)_$(expid)_$(chip).png"
        fig.savefig(imagePath, bbox_inches = "tight", pad_inches = 0.1)
	plt.close()
    end

    # find local maxima and minima
    slopes = smoothed_fluxes[begin+1:end]-smoothed_fluxes[begin:end-1]
    local_max_inds = findall((slopes[1:end-1] .>= 0) .& (slopes[2:end] .<= 0) .& .!((slopes[1:end-1] .== 0) .& (slopes[2:end] .== 0))) .+ 1
    local_min_inds = findall((slopes[1:end-1] .<= 0) .& (slopes[2:end] .>= 0)) .+ 1

    poss_local_max_fluxes = smoothed_fluxes[local_max_inds]
    poss_local_min_fluxes = smoothed_fluxes[local_min_inds]
    med_min_fluxes = nanmedian(poss_local_min_fluxes)
    med_max_fluxes = nanmedian(poss_local_max_fluxes)

    local_min_waves = y_vals[local_min_inds]
    local_max_waves = y_vals[local_max_inds]

    # fit polynomial functions to min(Y) and max(Y) to define good peaks

    keep_min = (.!(isnan.(poss_local_min_fluxes))) .&
               (abs.(poss_local_min_fluxes./med_min_fluxes .- 1.0) .< 0.5)
    keep_max = (.!(isnan.(poss_local_max_fluxes))) .&
               (abs.(poss_local_max_fluxes./med_max_fluxes .- 1.0) .< 0.5)
#    keep_min[begin:begin+5] .= false
#    keep_min[end-5:end] .= false
#    keep_max[begin:begin+5] .= false
#    keep_max[end-5:end] .= false
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
        keep_min .= (resids_min .>= resid_summary_min[1] - n_sigma * resid_summary_min[2]) .&
                    (resids_min .<= resid_summary_min[1] + n_sigma * resid_summary_min[3])
        keep_max .= (resids_max .>= resid_summary_max[1] - n_sigma * resid_summary_max[2]) .&
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

    if plotter
        # Figures for QA
        fig = plt.figure(figsize = (10, 5), dpi = 300)
        ax = fig.add_subplot(1, 1, 1)
	ax.plot(local_min_waves,poss_local_min_fluxes)
	ax.plot(local_max_waves,poss_local_max_fluxes)
	ax.plot(local_max_waves,med_min_fluxes)
	ax.plot(local_max_waves,med_max_fluxes)
        ttl = plt.text(
            0.5, 1.01, "Tele: $(tele), MJD: $(mjd), Chip: $(chip), Expid: $(expid)",
            ha = "center", va = "bottom", transform = ax.transAxes)
	imagePath = dirNamePlots * 
			"local_max_min_polynomials_$(tele)_$(mjd)_$(expid)_$(chip).png"
        fig.savefig(imagePath, bbox_inches = "tight", pad_inches = 0.1)
	plt.close()
    end

    #only keep peaks that are substantially above the minimum function
    relative_fluxes = (poss_local_max_fluxes .- med_min_fluxes) ./ (med_max_fluxes .- med_min_fluxes)
    good_max_inds = abs.(relative_fluxes .- 1.0) .< 0.5

    good_y_vals = local_max_waves[good_max_inds]

    #fit 1D gaussians to each identified peak, using offset_inds around each peak
    offset_inds = [-2,-1,0,1,2] .+ 1
    #use a larger number of pixels for removing the contribution from neighbouring peaks
    best_model_offset_inds = range(start=-10,stop=10,step=1) .+ 1 

    all_min_fluxes = min_func.(y_vals)
    all_max_fluxes = max_func.(y_vals)

    #fit scaled fluxes
    all_rel_fluxes = (comb_fluxes)./(all_max_fluxes)
    all_rel_errs = comb_errs./(all_max_fluxes)

    # first guess parameters
    fit_inds = good_y_vals .+ offset_inds'
    best_model_fit_inds = good_y_vals .+ best_model_offset_inds'
    fit_fluxes = all_rel_fluxes[fit_inds]
    fit_errs = all_rel_errs[fit_inds]
    fit_ivars = fit_errs .^ -2

    first_guess_params = zeros(Float64, size(fit_fluxes, 1), 3)
    first_guess_params[:, 1] .= max.(0.1, maximum(fit_fluxes, dims=2)) # height
#    first_guess_params[:, 2] .= (fit_inds[range(1,size(fit_inds,1),step=1),argmax(fit_fluxes, dims=2)])' # center position index for mu
    first_guess_params[:, 2] .= good_y_vals # center position index for mu
    first_guess_params[:, 3] .= 0.7 # sigma

    first_guess_params[:, 1] .*= (2 * π) ^ 0.5 * first_guess_params[:, 3] # change to integrated height

    curr_guess = copy(first_guess_params)

    dmu = 0.001
    dsig = 0.01
    group_inv(x, y) = mapslices(inv, x, dims = y)

    n_iter = 10
    for r_ind in 1:n_iter
        if r_ind > 1
        # use previous iteration parameters
	    fit_inds = floor.(Int,round.(curr_guess[:, 2])) .+ offset_inds'
	    best_model_fit_inds = floor.(Int,round.(curr_guess[:, 2])) .+ best_model_offset_inds'
	    fit_fluxes = copy(all_rel_fluxes[fit_inds])
	    fit_errs = all_rel_errs[fit_inds]
	    fit_ivars = fit_errs .^ -2
        end
        # CDF version
        model_fluxes_unit_height = 0.5 * erf.(((fit_inds .+ 0.5) .- curr_guess[:, 2]) ./ (2 ^ 0.5 * curr_guess[:, 3])) .- 
                                   0.5 * erf.(((fit_inds .- 0.5) .- curr_guess[:, 2]) ./ (2 ^ 0.5 * curr_guess[:, 3]))
    
        model_fluxes_unit_height_above_mu = 0.5 * erf.(((fit_inds .+ 0.5) .- (curr_guess[:, 2] .+ dmu)) ./ (2 ^ 0.5 * curr_guess[:, 3])) .- 
                                            0.5 * erf.(((fit_inds .- 0.5) .- (curr_guess[:, 2] .+ dmu)) ./ (2 ^ 0.5 * curr_guess[:, 3]))
        model_fluxes_unit_height_below_mu = 0.5 * erf.(((fit_inds .+ 0.5) .- (curr_guess[:, 2] .- dmu)) ./ (2 ^ 0.5 * curr_guess[:, 3])) .- 
                                            0.5 * erf.(((fit_inds .- 0.5) .- (curr_guess[:, 2] .- dmu)) ./ (2 ^ 0.5 * curr_guess[:, 3]))
        model_fluxes_unit_height_above_sig = 0.5 * erf.(((fit_inds .+ 0.5) .- curr_guess[:, 2]) ./ (2 ^ 0.5 * (curr_guess[:, 3] .+ dsig))) .- 
                                             0.5 * erf.(((fit_inds .- 0.5) .- curr_guess[:, 2]) ./ (2 ^ 0.5 * (curr_guess[:, 3] .+ dsig)))
        model_fluxes_unit_height_below_sig = 0.5 * erf.(((fit_inds .+ 0.5) .- curr_guess[:, 2]) ./ (2 ^ 0.5 * (curr_guess[:, 3] .- dsig))) .- 
                                             0.5 * erf.(((fit_inds .- 0.5) .- curr_guess[:, 2]) ./ (2 ^ 0.5 * (curr_guess[:, 3] .- dsig)))    

        model_fluxes = (curr_guess[:, 1] .* model_fluxes_unit_height)
        best_model_fluxes_unit_height = 0.5 * erf.(((best_model_fit_inds .+ 0.5) .- curr_guess[:, 2]) ./ (2 ^ 0.5 * curr_guess[:, 3])) .- 
                                        0.5 * erf.(((best_model_fit_inds .- 0.5) .- curr_guess[:, 2]) ./ (2 ^ 0.5 * curr_guess[:, 3]))
        best_model_fluxes = (curr_guess[:, 1] .* best_model_fluxes_unit_height)
        best_model_fluxes[.!isfinite.(curr_guess[:, 1]),:] .= 0

        comb_model_fluxes = zeros(size(all_rel_fluxes))
        for j in 1:size(best_model_fit_inds, 1)
            comb_model_fluxes[best_model_fit_inds[j,:]] .+= best_model_fluxes[j,:]
        end

        # remove combined gaussians
        fit_fluxes .-= comb_model_fluxes[fit_inds]
        # add back in the individual gaussian for each peak so only neighbours are removed
        fit_fluxes .+= best_model_fluxes[:, offset_inds .+ (size(best_model_fluxes, 2) .÷ 2)]

        dmodel_dheight = model_fluxes_unit_height
        dmodel_dmu = (curr_guess[:, 1] .* (model_fluxes_unit_height_above_mu .- model_fluxes_unit_height_below_mu)) ./ dmu
        dmodel_dsig = (curr_guess[:, 1] .* (model_fluxes_unit_height_above_sig .- model_fluxes_unit_height_below_sig)) ./ dsig

        flux_diffs = fit_fluxes .- model_fluxes

        # derivative matrix
        M_vects = zeros(Float64, (size(model_fluxes, 1), size(model_fluxes, 2), size(first_guess_params,2)))

        M_vects[:, :, 1] .= dmodel_dheight  # d flux/d height
        M_vects[:, :, 2] .= dmodel_dmu      # d flux/d mu
        M_vects[:, :, 3] .= dmodel_dsig     # d flux/d sigma

        M_T_dot_V_inv = M_vects .* fit_ivars
	@einsum M_T_dot_V_inv_dot_M[n,j,k] := M_T_dot_V_inv[n,i,j] * M_vects[n,i,k]
	@einsum M_T_dot_V_inv_dot_y[n,j] := M_T_dot_V_inv[n,i,j] * flux_diffs[n,i]

        scales = M_T_dot_V_inv_dot_M[1:size(M_T_dot_V_inv_dot_M, 1), 1, 1]
        v_hat_cov = group_inv(M_T_dot_V_inv_dot_M ./ scales,(2,3)) ./ scales

        # update vector
        @einsum v_hat[n,i] := v_hat_cov[n,i,j] * M_T_dot_V_inv_dot_y[n,j]

        # restrict how far the update can move from the starting point
        v_hat[:, 1] .= clamp.(v_hat[:, 1], -curr_guess[:, 1] * 0.9, curr_guess[:, 1] * 5)
        v_hat[:, 2] .= clamp.(v_hat[:, 2], -2, 2)
        v_hat[:, 3] .= clamp.(v_hat[:, 3], -curr_guess[:, 3] * 0.9, curr_guess[:, 3] * 5)

        new_params = curr_guess .+ v_hat
        curr_guess .= new_params
    end

    new_params = curr_guess

    # use the location of the current peaks to fill in any gaps for low-throughput fibers
    peak_locs = new_params[1:end-1, 2]
    peak_spacing = new_params[2:end, 2] .- new_params[1:end-1, 2]
    med_peak_spacing = nanmedian(peak_spacing)

    good_peak_spacing = abs.(peak_spacing ./ med_peak_spacing .- 1.0) .< 0.5

    if plotter
        # Figures for QA
        fig = plt.figure(figsize = (10, 5), dpi = 300)
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(new_params[:, 2],new_params[:, 1],lw=1,color="C0")
        ttl = plt.text(
            0.5, 1.01, "Tele: $(tele), MJD: $(mjd), Chip: $(chip), Expid: $(expid)",
            ha = "center", va = "bottom", transform = ax.transAxes)
	imagePath = dirNamePlots * 
			"ave_middle_heights_$(tele)_$(mjd)_$(expid)_$(chip).png"
        fig.savefig(imagePath, bbox_inches = "tight", pad_inches = 0.1)
	plt.close()

        fig = plt.figure(figsize = (10, 5), dpi = 300)
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(new_params[:, 2],new_params[:, 3],lw=1,color="C0")
        ttl = plt.text(
            0.5, 1.01, "Tele: $(tele), MJD: $(mjd), Chip: $(chip), Expid: $(expid)",
            ha = "center", va = "bottom", transform = ax.transAxes)
	imagePath = dirNamePlots * 
			"ave_middle_widths_$(tele)_$(mjd)_$(expid)_$(chip).png"
        fig.savefig(imagePath, bbox_inches = "tight", pad_inches = 0.1)
	plt.close()

        fig = plt.figure(figsize = (10, 5), dpi = 300)
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(peak_locs,peak_spacing, marker="o")
        ax.plot(peak_locs[good_peak_spacing],peak_spacing[good_peak_spacing], marker="o")
        ttl = plt.text(
            0.5, 1.01, "Tele: $(tele), MJD: $(mjd), Chip: $(chip), Expid: $(expid)",
            ha = "center", va = "bottom", transform = ax.transAxes)
	imagePath = dirNamePlots * 
			"ave_middle_peak_spacing_$(tele)_$(mjd)_$(expid)_$(chip).png"
        fig.savefig(imagePath, bbox_inches = "tight", pad_inches = 0.1)
	plt.close()
    end

    keep_space = copy(good_peak_spacing)

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
    end
    space_func = fit(peak_locs[keep_space], peak_spacing[keep_space], 3)

    space_func_eval = space_func.(peak_locs)
    int_spacing = floor.(Int,round.(peak_spacing ./ space_func_eval))

    if plotter
        fig = plt.figure(figsize = (10, 5), dpi = 300)
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(peak_locs[good_peak_spacing],peak_spacing[good_peak_spacing], marker="o")
        ax.plot(peak_locs[keep_space],peak_spacing[keep_space], marker="o")
        ax.plot(peak_locs,space_func_eval)
        ttl = plt.text(
            0.5, 1.01, "Tele: $(tele), MJD: $(mjd), Chip: $(chip), Expid: $(expid)",
            ha = "center", va = "bottom", transform = ax.transAxes)
	imagePath = dirNamePlots * 
			"ave_middle_peak_spacing_func_$(tele)_$(mjd)_$(expid)_$(chip).png"
        fig.savefig(imagePath, bbox_inches = "tight", pad_inches = 0.1)
	plt.close()

        fig = plt.figure(figsize = (10, 5), dpi = 300)
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(peak_locs,peak_spacing ./ space_func_eval, marker="o")
        ax.plot(peak_locs,int_spacing)
        ttl = plt.text(
            0.5, 1.01, "Tele: $(tele), MJD: $(mjd), Chip: $(chip), Expid: $(expid)",
            ha = "center", va = "bottom", transform = ax.transAxes)
	imagePath = dirNamePlots * 
			"ave_middle_int_spacing_$(tele)_$(mjd)_$(expid)_$(chip).png"
        fig.savefig(imagePath, bbox_inches = "tight", pad_inches = 0.1)
	plt.close()
    end

    peak_ints = zeros(Int, size(int_spacing,1) + 1)
    peak_ints[2:end] .= cumsum(int_spacing, dims=1)

    keep_peaks = ones(Bool, size(peak_ints,1))
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

    all_peak_ints = range(start=max(peak_ints[1], 0) - 50,stop=min(1000, peak_ints[end]) + 50,step=1)
    all_peak_locs = peak_func.(all_peak_ints)
    keep_peak_ints = (all_peak_locs .> 1 + 10) .& (all_peak_locs .< 2048 - 10)
    all_peak_ints = all_peak_ints[keep_peak_ints]
    all_peak_locs = all_peak_locs[keep_peak_ints]

    missing_ints = zeros(Bool,size(all_peak_ints,1))
    for (ind, peak_int) in enumerate(all_peak_ints)
        if !(peak_int in peak_ints)
            missing_ints[ind] = true
        end
    end

    missing_ints = all_peak_ints[missing_ints]

    if plotter
        fig = plt.figure(figsize = (10, 5), dpi = 300)
        ax = fig.add_subplot(1, 1, 1)
        ax.scatter(peak_ints,new_params[:,2])
        ax.scatter(peak_ints[keep_peaks],new_params[keep_peaks,2])
        ax.plot(peak_ints,peak_func.(peak_ints),c="r")
	for ind in missing_ints
	    ax.axvline(ind,c="w",lw=1)
	end
        ttl = plt.text(
            0.5, 1.01, "Tele: $(tele), MJD: $(mjd), Chip: $(chip), Expid: $(expid)",
            ha = "center", va = "bottom", transform = ax.transAxes)
	imagePath = dirNamePlots * 
			"ave_middle_missing_peaks_$(tele)_$(mjd)_$(expid)_$(chip).png"
        fig.savefig(imagePath, bbox_inches = "tight", pad_inches = 0.1)
	plt.close()

        fig = plt.figure(figsize = (10, 5), dpi = 300)
        ax = fig.add_subplot(1, 1, 1)
        ax.scatter(peak_ints,new_params[:,2]-peak_func.(peak_ints))
        ax.scatter(peak_ints[keep_peaks],new_params[keep_peaks,2]-peak_func.(peak_ints)[keep_peaks])
	ax.axhline(c="r")
	for ind in missing_ints
	    ax.axvline(ind,c="w",lw=1)
	end
        ttl = plt.text(
            0.5, 1.01, "Tele: $(tele), MJD: $(mjd), Chip: $(chip), Expid: $(expid)",
            ha = "center", va = "bottom", transform = ax.transAxes)
	imagePath = dirNamePlots * 
			"ave_middle_missing_peak_resids_$(tele)_$(mjd)_$(expid)_$(chip).png"
        fig.savefig(imagePath, bbox_inches = "tight", pad_inches = 0.1)
	plt.close()

    end

    #using the full list of peaks, get the best-fit parameters
    all_rel_fluxes = comb_fluxes
    all_rel_errs = comb_errs

    # first guess parameters
    fit_inds = floor.(Int,round.(all_peak_locs)) .+ offset_inds'
    best_model_fit_inds = floor.(Int,round.(all_peak_locs)) .+ best_model_offset_inds'
    fit_fluxes = all_rel_fluxes[fit_inds]
    fit_errs = all_rel_errs[fit_inds]
    fit_ivars = fit_errs .^ -2

    first_guess_params = zeros(Float64, size(fit_fluxes, 1), 3)
    first_guess_params[:, 1] .= max.(0.1, maximum(fit_fluxes, dims=2)) # height
    first_guess_params[:, 2] .= all_peak_locs #center position index for mu
    first_guess_params[:, 3] .= nanmedian(new_params[:,3]) # sigma

    first_guess_params[:, 1] .*= (2 * π) ^ 0.5 * first_guess_params[:, 3] # change to integrated height

    curr_guess = copy(first_guess_params)

    dmu = 0.001
    dsig = 0.01

    n_iter = 20
    for r_ind in 1:n_iter
        if r_ind > 1
        # use previous iteration parameters
	    fit_inds = floor.(Int,round.(curr_guess[:, 2])) .+ offset_inds'
	    best_model_fit_inds = floor.(Int,round.(curr_guess[:, 2])) .+ best_model_offset_inds'
	    fit_fluxes = copy(all_rel_fluxes[fit_inds])
	    fit_errs = all_rel_errs[fit_inds]
	    fit_ivars = fit_errs .^ -2
        end
        # CDF version
        model_fluxes_unit_height = 0.5 * erf.(((fit_inds .+ 0.5) .- curr_guess[:, 2]) ./ (2 ^ 0.5 * curr_guess[:, 3])) .- 
                                   0.5 * erf.(((fit_inds .- 0.5) .- curr_guess[:, 2]) ./ (2 ^ 0.5 * curr_guess[:, 3]))
    
        model_fluxes_unit_height_above_mu = 0.5 * erf.(((fit_inds .+ 0.5) .- (curr_guess[:, 2] .+ dmu)) ./ (2 ^ 0.5 * curr_guess[:, 3])) .- 
                                            0.5 * erf.(((fit_inds .- 0.5) .- (curr_guess[:, 2] .+ dmu)) ./ (2 ^ 0.5 * curr_guess[:, 3]))
        model_fluxes_unit_height_below_mu = 0.5 * erf.(((fit_inds .+ 0.5) .- (curr_guess[:, 2] .- dmu)) ./ (2 ^ 0.5 * curr_guess[:, 3])) .- 
                                            0.5 * erf.(((fit_inds .- 0.5) .- (curr_guess[:, 2] .- dmu)) ./ (2 ^ 0.5 * curr_guess[:, 3]))
        model_fluxes_unit_height_above_sig = 0.5 * erf.(((fit_inds .+ 0.5) .- curr_guess[:, 2]) ./ (2 ^ 0.5 * (curr_guess[:, 3] .+ dsig))) .- 
                                             0.5 * erf.(((fit_inds .- 0.5) .- curr_guess[:, 2]) ./ (2 ^ 0.5 * (curr_guess[:, 3] .+ dsig)))
        model_fluxes_unit_height_below_sig = 0.5 * erf.(((fit_inds .+ 0.5) .- curr_guess[:, 2]) ./ (2 ^ 0.5 * (curr_guess[:, 3] .- dsig))) .- 
                                             0.5 * erf.(((fit_inds .- 0.5) .- curr_guess[:, 2]) ./ (2 ^ 0.5 * (curr_guess[:, 3] .- dsig)))    

        model_fluxes = (curr_guess[:, 1] .* model_fluxes_unit_height)
        best_model_fluxes_unit_height = 0.5 * erf.(((best_model_fit_inds .+ 0.5) .- curr_guess[:, 2]) ./ (2 ^ 0.5 * curr_guess[:, 3])) .- 
                                        0.5 * erf.(((best_model_fit_inds .- 0.5) .- curr_guess[:, 2]) ./ (2 ^ 0.5 * curr_guess[:, 3]))
        best_model_fluxes = (curr_guess[:, 1] .* best_model_fluxes_unit_height)
        best_model_fluxes[.!isfinite.(curr_guess[:, 1]),:] .= 0

        comb_model_fluxes = zeros(size(all_rel_fluxes))
        for j in 1:size(best_model_fit_inds, 1)
            comb_model_fluxes[best_model_fit_inds[j,:]] .+= best_model_fluxes[j,:]
        end

        # remove combined gaussians
        fit_fluxes .-= comb_model_fluxes[fit_inds]
        # add back in the individual gaussian for each peak so only neighbours are removed
        fit_fluxes .+= best_model_fluxes[:, offset_inds .+ (size(best_model_fluxes, 2) .÷ 2)]

        dmodel_dheight = model_fluxes_unit_height
        dmodel_dmu = (curr_guess[:, 1] .* (model_fluxes_unit_height_above_mu .- model_fluxes_unit_height_below_mu)) ./ dmu
        dmodel_dsig = (curr_guess[:, 1] .* (model_fluxes_unit_height_above_sig .- model_fluxes_unit_height_below_sig)) ./ dsig

        flux_diffs = fit_fluxes .- model_fluxes

        # derivative matrix
        M_vects = zeros(Float64, (size(model_fluxes, 1), size(model_fluxes, 2), size(first_guess_params,2)))

        M_vects[:, :, 1] .= dmodel_dheight  # d flux/d height
        M_vects[:, :, 2] .= dmodel_dmu      # d flux/d mu
        M_vects[:, :, 3] .= dmodel_dsig     # d flux/d sigma

        M_T_dot_V_inv = M_vects .* fit_ivars
	@einsum M_T_dot_V_inv_dot_M[n,j,k] := M_T_dot_V_inv[n,i,j] * M_vects[n,i,k]
	@einsum M_T_dot_V_inv_dot_y[n,j] := M_T_dot_V_inv[n,i,j] * flux_diffs[n,i]

        scales = M_T_dot_V_inv_dot_M[1:size(M_T_dot_V_inv_dot_M, 1), 1, 1]
        v_hat_cov = group_inv(M_T_dot_V_inv_dot_M ./ scales,(2,3)) ./ scales

        # update vector
        @einsum v_hat[n,i] := v_hat_cov[n,i,j] * M_T_dot_V_inv_dot_y[n,j]


        # restrict how far the update can move from the starting point
        v_hat[:, 1] .= clamp.(v_hat[:, 1], -curr_guess[:, 1] * 0.9, curr_guess[:, 1] * 5)
        v_hat[:, 2] .= clamp.(v_hat[:, 2], -2, 2)
        v_hat[:, 3] .= clamp.(v_hat[:, 3], -curr_guess[:, 3] * 0.9, curr_guess[:, 3] * 5)

        new_params = curr_guess .+ v_hat
	new_params[:, 1] .= clamp.(new_params[:, 1], 0.01, first_guess_params[:, 1]*10)
	new_params[:, 2] .= clamp.(new_params[:, 2], all_peak_locs .- 3, all_peak_locs .+ 3)
	new_params[:, 3] .= clamp.(new_params[:, 3], 0.5, 2.0)

        curr_guess .= new_params
        if plotter & (r_ind == n_iter)
            fig = plt.figure(figsize = (10, 5), dpi = 300)
            ax = fig.add_subplot(1, 1, 1)
            ax.plot(all_rel_fluxes)
            ax.plot(comb_model_fluxes)
            ttl = plt.text(
                0.5, 1.01, "Tele: $(tele), MJD: $(mjd), Chip: $(chip), Expid: $(expid)",
                ha = "center", va = "bottom", transform = ax.transAxes)
    	    imagePath = dirNamePlots * 
			"ave_middle_flux_comp_$(tele)_$(mjd)_$(expid)_$(chip).png"
            fig.savefig(imagePath, bbox_inches = "tight", pad_inches = 0.1)
   	    plt.close()

            fig = plt.figure(figsize = (10, 5), dpi = 300)
            ax = fig.add_subplot(1, 1, 1)
            ax.plot(all_rel_fluxes-comb_model_fluxes)
            ttl = plt.text(
                0.5, 1.01, "Tele: $(tele), MJD: $(mjd), Chip: $(chip), Expid: $(expid)",
                ha = "center", va = "bottom", transform = ax.transAxes)
    	    imagePath = dirNamePlots * 
			"ave_middle_flux_comp_resids_$(tele)_$(mjd)_$(expid)_$(chip).png"
            fig.savefig(imagePath, bbox_inches = "tight", pad_inches = 0.1)
	    plt.close()
        end
    end

    new_params = curr_guess

    #save the middle-of-detector best fit parameters for first guess
    best_fit_ave_params = copy(new_params)

    #work outwards from the middle of the detector
    #using the previous analyses to give good first guesses
    #when iterating over all the X pixels

    peak_locs = new_params[1:end-1,2]
    peak_spacing = new_params[2:end,2]-new_params[1:end-1,2]
    med_peak_spacing = nanmedian(peak_spacing,1)

    good_peak_spacing = abs.(peak_spacing ./ med_peak_spacing .- 1.0) .< 0.5
    med_flux = nanmedian(new_params[:,1],1)

    if plotter
        fig = plt.figure(figsize = (10, 5), dpi = 300)
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(new_params[:,2],new_params[:,1],marker="o")
        ttl = plt.text(
            0.5, 1.01, "Tele: $(tele), MJD: $(mjd), Chip: $(chip), Expid: $(expid)",
            ha = "center", va = "bottom", transform = ax.transAxes)
        imagePath = dirNamePlots * 
			"ave_middle_heights_ALL_$(tele)_$(mjd)_$(expid)_$(chip).png"
        fig.savefig(imagePath, bbox_inches = "tight", pad_inches = 0.1)
	plt.close()

        fig = plt.figure(figsize = (10, 5), dpi = 300)
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(new_params[:,2],new_params[:,3],marker="o")
        ttl = plt.text(
            0.5, 1.01, "Tele: $(tele), MJD: $(mjd), Chip: $(chip), Expid: $(expid)",
            ha = "center", va = "bottom", transform = ax.transAxes)
        imagePath = dirNamePlots * 
			"ave_middle_widths_ALL_$(tele)_$(mjd)_$(expid)_$(chip).png"
        fig.savefig(imagePath, bbox_inches = "tight", pad_inches = 0.1)
	plt.close()

        fig = plt.figure(figsize = (10, 5), dpi = 300)
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(peak_locs,peak_spacing,marker="o")
        ax.plot(peak_locs[good_peak_spacing],peak_spacing[good_peak_spacing],marker="o")
        ttl = plt.text(
            0.5, 1.01, "Tele: $(tele), MJD: $(mjd), Chip: $(chip), Expid: $(expid)",
            ha = "center", va = "bottom", transform = ax.transAxes)
        imagePath = dirNamePlots * 
			"ave_middle_peak_spacing_ALL_$(tele)_$(mjd)_$(expid)_$(chip).png"
        fig.savefig(imagePath, bbox_inches = "tight", pad_inches = 0.1)
	plt.close()

        fig = plt.figure(figsize = (10, 5), dpi = 300)
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(new_params[:,2],new_params[:,1]./med_flux,marker="o")
        ttl = plt.text(
            0.5, 1.01, "Tele: $(tele), MJD: $(mjd), Chip: $(chip), Expid: $(expid)",
            ha = "center", va = "bottom", transform = ax.transAxes)
        imagePath = dirNamePlots * 
			"ave_middle_heightsScaled_ALL_$(tele)_$(mjd)_$(expid)_$(chip).png"
        fig.savefig(imagePath, bbox_inches = "tight", pad_inches = 0.1)
	plt.close()

    end

    #remove the edge possible peaks if they have no throughput, because they likely don't exist
    good_throughput_fibers = (best_fit_ave_params[:,1] ./ med_flux) .> 0.1
    low_throughput_fibers = findall(.!good_throughput_fibers)
    if size(low_throughput_fibers,1) > 0
        if low_throughput_fibers[1] == 1
            #then there is a truncation applied to the left
            curr_cut_int = low_throughput_fibers[1]
            j = 2
            while j < size(low_throughput_fibers,1)+1
                if curr_cut_int+1 == low_throughput_fibers[j]
                    curr_cut_int = low_throughput_fibers[j]
	        else
	            break
		end
	        j += 1
	    end
	    println("Left Cut Integer ",curr_cut_int)
            left_cut_ind = curr_cut_int+1
        else
            left_cut_ind = 1
	end

        if low_throughput_fibers[end] == size(new_params,1)
            #then there is a truncation applied to the right
	    curr_cut_int = low_throughput_fibers[end]
	    j = size(low_throughput_fibers,1)-1
	    while j > 1
	        if curr_cut_int-1 == low_throughput_fibers[j]
	            curr_cut_int = low_throughput_fibers[j]
	        else
	            break
		end
	        j -= 1
	    end
	    println("Right Cut Integer ",curr_cut_int)
	    right_cut_ind = curr_cut_int-1
	else
	    right_cut_ind = size(new_params,1)
	end

        best_fit_ave_params = best_fit_ave_params[left_cut_ind:right_cut_ind,:]
    end

    println("Final number of peaks:",size(best_fit_ave_params))

    #use more pixels around each peak, otherwise large-width peak fitting fails
    #(see weird discontinuities in height and width vs X)
    #This ONLY works because we are fitting all the peaks simultaneously 
    #and removing contamination from neighbouring fibers
    offset_inds = range(-7,7,step=1) .+ 1

    x_inds = range(1,2048,step=1)
    #work from the center outward to use constraints from previous analysis 
    sorted_x_inds = x_inds[sortperm(abs.(x_inds .- 1025))]

    param_outputs = zeros((size(x_inds,1),size(best_fit_ave_params,1),size(best_fit_ave_params,2)))
    param_output_covs = zeros((size(x_inds,1),size(best_fit_ave_params,1),size(best_fit_ave_params,2),size(best_fit_ave_params,2))) 

    @showprogress for (ind,x_ind) in enumerate(sorted_x_inds)
        #use previous analyses to constrain first guesses
        if abs(x_ind - 1025) <= 100
            first_guess_params = copy(best_fit_ave_params)
        elseif x_ind > 1025
	    first_guess_params = nanmedian(param_outputs[x_ind-51:x_ind-1,:,:],1)[1,:,:]
	elseif x_ind < 1025
	    first_guess_params = nanmedian(param_outputs[x_ind+1:x_ind+51,:,:],1)[1,:,:]
	end
					    
        all_rel_fluxes = copy(image_data[x_ind,:])
        all_rel_errs = copy(noise_image[x_ind,:])
        all_rel_masks = copy(image_mask[x_ind,:])

        all_rel_fluxes[.!all_rel_masks] .= 0
        all_rel_ivars = all_rel_errs .^ (-2)
        all_rel_ivars[.!all_rel_masks] .= 0
        all_rel_errs[.!all_rel_masks] .= Inf

        # first guess parameters
        fit_inds = floor.(Int,round.(first_guess_params[:,2])) .+ offset_inds'
        best_model_fit_inds = floor.(Int,round.(first_guess_params[:,2])) .+ best_model_offset_inds'
        fit_fluxes = all_rel_fluxes[fit_inds]
        fit_errs = all_rel_errs[fit_inds]
        fit_ivars = fit_errs .^ -2

        curr_guess = copy(first_guess_params)

        dmu = 0.001
        dsig = 0.01

        n_iter = 30
        n_iter = 20
        for r_ind in 1:n_iter
            if r_ind > 1
                # use previous iteration parameters
                fit_inds = floor.(Int,round.(curr_guess[:, 2])) .+ offset_inds'
                best_model_fit_inds = floor.(Int,round.(curr_guess[:, 2])) .+ best_model_offset_inds'
    	        fit_fluxes = copy(all_rel_fluxes[fit_inds])
    	        fit_errs = all_rel_errs[fit_inds]
	        fit_ivars = fit_errs .^ -2
            end

            fit_ivars[fit_ivars .== 0] .= 1e-5

            # CDF version
            model_fluxes_unit_height = 0.5 * erf.(((fit_inds .+ 0.5) .- curr_guess[:, 2]) ./ (2 ^ 0.5 * curr_guess[:, 3])) .- 
                                       0.5 * erf.(((fit_inds .- 0.5) .- curr_guess[:, 2]) ./ (2 ^ 0.5 * curr_guess[:, 3]))
    
            model_fluxes_unit_height_above_mu = 0.5 * erf.(((fit_inds .+ 0.5) .- (curr_guess[:, 2] .+ dmu)) ./ (2 ^ 0.5 * curr_guess[:, 3])) .- 
                                                0.5 * erf.(((fit_inds .- 0.5) .- (curr_guess[:, 2] .+ dmu)) ./ (2 ^ 0.5 * curr_guess[:, 3]))
            model_fluxes_unit_height_below_mu = 0.5 * erf.(((fit_inds .+ 0.5) .- (curr_guess[:, 2] .- dmu)) ./ (2 ^ 0.5 * curr_guess[:, 3])) .- 
                                                0.5 * erf.(((fit_inds .- 0.5) .- (curr_guess[:, 2] .- dmu)) ./ (2 ^ 0.5 * curr_guess[:, 3]))
            model_fluxes_unit_height_above_sig = 0.5 * erf.(((fit_inds .+ 0.5) .- curr_guess[:, 2]) ./ (2 ^ 0.5 * (curr_guess[:, 3] .+ dsig))) .- 
                                                 0.5 * erf.(((fit_inds .- 0.5) .- curr_guess[:, 2]) ./ (2 ^ 0.5 * (curr_guess[:, 3] .+ dsig)))
            model_fluxes_unit_height_below_sig = 0.5 * erf.(((fit_inds .+ 0.5) .- curr_guess[:, 2]) ./ (2 ^ 0.5 * (curr_guess[:, 3] .- dsig))) .- 
                                                 0.5 * erf.(((fit_inds .- 0.5) .- curr_guess[:, 2]) ./ (2 ^ 0.5 * (curr_guess[:, 3] .- dsig)))    

            model_fluxes = (curr_guess[:, 1] .* model_fluxes_unit_height)
            best_model_fluxes_unit_height = 0.5 * erf.(((best_model_fit_inds .+ 0.5) .- curr_guess[:, 2]) ./ (2 ^ 0.5 * curr_guess[:, 3])) .- 
                                            0.5 * erf.(((best_model_fit_inds .- 0.5) .- curr_guess[:, 2]) ./ (2 ^ 0.5 * curr_guess[:, 3]))
            best_model_fluxes = (curr_guess[:, 1] .* best_model_fluxes_unit_height)
            best_model_fluxes[.!isfinite.(curr_guess[:, 1]),:] .= 0

            comb_model_fluxes = zeros(size(all_rel_fluxes))
            for j in 1:size(best_model_fit_inds, 1)
                 comb_model_fluxes[best_model_fit_inds[j,:]] .+= best_model_fluxes[j,:]
            end

            # remove combined gaussians
            fit_fluxes .-= comb_model_fluxes[fit_inds]
            # add back in the individual gaussian for each peak so only neighbours are removed
            fit_fluxes .+= best_model_fluxes[:, offset_inds .+ (size(best_model_fluxes, 2) .÷ 2)]

            dmodel_dheight = model_fluxes_unit_height
            dmodel_dmu = (curr_guess[:, 1] .* (model_fluxes_unit_height_above_mu .- model_fluxes_unit_height_below_mu)) ./ dmu
            dmodel_dsig = (curr_guess[:, 1] .* (model_fluxes_unit_height_above_sig .- model_fluxes_unit_height_below_sig)) ./ dsig

            flux_diffs = fit_fluxes .- model_fluxes

            # derivative matrix
            M_vects = zeros(Float64, (size(model_fluxes, 1), size(model_fluxes, 2), size(first_guess_params,2)))

            M_vects[:, :, 1] .= dmodel_dheight  # d flux/d height
            M_vects[:, :, 2] .= dmodel_dmu      # d flux/d mu
            M_vects[:, :, 3] .= dmodel_dsig     # d flux/d sigma

            M_T_dot_V_inv = M_vects .* fit_ivars
    	    @einsum M_T_dot_V_inv_dot_M[n,j,k] := M_T_dot_V_inv[n,i,j] * M_vects[n,i,k]
	    @einsum M_T_dot_V_inv_dot_y[n,j] := M_T_dot_V_inv[n,i,j] * flux_diffs[n,i]

            scales = M_T_dot_V_inv_dot_M[1:size(M_T_dot_V_inv_dot_M, 1), 1, 1]
            v_hat_cov = group_inv(M_T_dot_V_inv_dot_M ./ scales,(2,3)) ./ scales

            # update vector
            @einsum v_hat[n,i] := v_hat_cov[n,i,j] * M_T_dot_V_inv_dot_y[n,j]


            # restrict how far the update can move from the starting point
            v_hat[:, 1] .= clamp.(v_hat[:, 1], -curr_guess[:, 1] * 0.9, curr_guess[:, 1] * 5)
            v_hat[:, 2] .= clamp.(v_hat[:, 2], -2, 2)
            v_hat[:, 3] .= clamp.(v_hat[:, 3], -curr_guess[:, 3] * 0.9, curr_guess[:, 3] * 5)

#	    println(r_ind," ",x_ind," ",v_hat[begin:begin+3,2])

            new_params = curr_guess .+ v_hat
   	    new_params[:, 1] .= clamp.(new_params[:, 1], 0.01, first_guess_params[:, 1]*10)
       	    new_params[:, 2] .= clamp.(new_params[:, 2], first_guess_params[:, 2] .- 3, first_guess_params[:, 2] .+ 3)
	    new_params[:, 3] .= clamp.(new_params[:, 3], 0.5, 2.0)

#	    println(r_ind," ",x_ind," ",new_params[begin:begin+3,2])

            curr_guess .= new_params

            if plotter & (r_ind == n_iter) & (x_ind == 1025)
                fig = plt.figure(figsize = (10, 5), dpi = 300)
                ax = fig.add_subplot(1, 1, 1)
                ax.plot(all_rel_fluxes)
                ax.plot(comb_model_fluxes)
                ttl = plt.text(
                    0.5, 1.01, "Tele: $(tele), MJD: $(mjd), Chip: $(chip), Expid: $(expid)",
                    ha = "center", va = "bottom", transform = ax.transAxes)
                imagePath = dirNamePlots * 
   			"x1025_flux_comp_$(tele)_$(mjd)_$(expid)_$(chip).png"
                fig.savefig(imagePath, bbox_inches = "tight", pad_inches = 0.1)
	        plt.close()

                fig = plt.figure(figsize = (10, 5), dpi = 300)
                ax = fig.add_subplot(1, 1, 1)
                ax.plot(all_rel_fluxes-comb_model_fluxes)
                ttl = plt.text(
                    0.5, 1.01, "Tele: $(tele), MJD: $(mjd), Chip: $(chip), Expid: $(expid)",
                    ha = "center", va = "bottom", transform = ax.transAxes)
                imagePath = dirNamePlots * 
   			"x1025_flux_comp_resids_$(tele)_$(mjd)_$(expid)_$(chip).png"
                fig.savefig(imagePath, bbox_inches = "tight", pad_inches = 0.1)
	        plt.close()
            end

            param_outputs[x_ind,:,:] .= new_params
            param_output_covs[x_ind,:,:,:] .= v_hat_cov
        end
    end

    if plotter
        fig = plt.figure(figsize = (10, 5), dpi = 300)
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(nanmedian(param_outputs[:,:,1],1)[1,:])
        ttl = plt.text(
            0.5, 1.01, "Tele: $(tele), MJD: $(mjd), Chip: $(chip), Expid: $(expid)",
            ha = "center", va = "bottom", transform = ax.transAxes)
        imagePath = dirNamePlots * 
			"final_median_flux_$(tele)_$(mjd)_$(expid)_$(chip).png"
        fig.savefig(imagePath, bbox_inches = "tight", pad_inches = 0.1)
	plt.close()

	fiber_inds = collect(0:50:300) .+ 1
        fiber_inds[1] = 1
        fiber_inds[end] = size(param_outputs,2)

        median_heights = nanmedian(param_outputs[:,:,1],1)
        keep_params = (param_outputs[:,:,3] .< 1.7) .& (param_outputs[:,:,3] .> 0.5) .& 
                      (param_outputs[:,:,1] ./ median_heights .> 0.1) .&
                      (param_outputs[:,:,1] ./ median_heights .< 2)

        fig = plt.figure(figsize = (10, 5), dpi = 300)
        ax = fig.add_subplot(1, 1, 1)
	ax.set_xlabel("X")
#	ax.set_ylabel("Y$-\langle$Y$\rangle$")
	ax.set_ylabel("Y")
	for fiber_ind in fiber_inds
            curr_keep = keep_params[:,fiber_ind]
	    ax.plot(x_inds[curr_keep],
		    param_outputs[curr_keep,fiber_ind,2] .- nanmedian(param_outputs[curr_keep,fiber_ind,2]),
		    marker="o",lw=1,label=fiber_ind,alpha=0.7,ms=1)
	end
	leg = ax.legend(loc=6,bbox_to_anchor=(1.05,0.5))
        ttl = plt.text(
            0.5, 1.01, "Tele: $(tele), MJD: $(mjd), Chip: $(chip), Expid: $(expid)",
            ha = "center", va = "bottom", transform = ax.transAxes)
        imagePath = dirNamePlots * 
			"final_fiber_trace_mean_$(tele)_$(mjd)_$(expid)_$(chip).png"
        fig.savefig(imagePath, bbox_inches = "tight", pad_inches = 0.1)
	plt.close()

        fig = plt.figure(figsize = (10, 5), dpi = 300)
        ax = fig.add_subplot(1, 1, 1)
	ax.set_xlabel("X")
	ax.set_ylabel("Flux")
	for fiber_ind in fiber_inds
            curr_keep = keep_params[:,fiber_ind]
	    ax.plot(x_inds[curr_keep],
		    param_outputs[curr_keep,fiber_ind,1],
		    marker="o",lw=1,label=fiber_ind,alpha=0.7,ms=1)
	end
	leg = ax.legend(loc=6,bbox_to_anchor=(1.05,0.5))
        ttl = plt.text(
            0.5, 1.01, "Tele: $(tele), MJD: $(mjd), Chip: $(chip), Expid: $(expid)",
            ha = "center", va = "bottom", transform = ax.transAxes)
        imagePath = dirNamePlots * 
			"final_fiber_trace_flux_$(tele)_$(mjd)_$(expid)_$(chip).png"
        fig.savefig(imagePath, bbox_inches = "tight", pad_inches = 0.1)
	plt.close()

        fig = plt.figure(figsize = (10, 5), dpi = 300)
        ax = fig.add_subplot(1, 1, 1)
	ax.set_xlabel("X")
	ax.set_ylabel("Y Width")
	for fiber_ind in fiber_inds
            curr_keep = keep_params[:,fiber_ind]
	    ax.plot(x_inds[curr_keep],
		    param_outputs[curr_keep,fiber_ind,3],
		    marker="o",lw=1,label=fiber_ind,alpha=0.7,ms=1)
	end
	leg = ax.legend(loc=6,bbox_to_anchor=(1.05,0.5))
        ttl = plt.text(
            0.5, 1.01, "Tele: $(tele), MJD: $(mjd), Chip: $(chip), Expid: $(expid)",
            ha = "center", va = "bottom", transform = ax.transAxes)
        imagePath = dirNamePlots * 
			"final_fiber_trace_width_$(tele)_$(mjd)_$(expid)_$(chip).png"
        fig.savefig(imagePath, bbox_inches = "tight", pad_inches = 0.1)
	plt.close()

        all_x_inds = zeros(Float64,size(param_outputs,1),size(param_outputs,2))
        all_x_inds[:,:] .= range(1,2048,step=1)

        ravel_x_inds = collect(Iterators.flatten(all_x_inds))
        ravel_y_centers = collect(Iterators.flatten(param_outputs[:,:,2]))
 
        curr_params = collect(Iterators.flatten(param_outputs[:,:,3]))
        param_summary = quantile(filter(!isnan, curr_params),[0.16, 0.5, 0.84])
        param_summary = [param_summary[2],
     		 param_summary[2]-param_summary[1],
		 param_summary[3]-param_summary[2]]
        vmin = param_summary[1] - 2*param_summary[2]
        vmax = param_summary[1] + 2*param_summary[3]
        fig = plt.figure(figsize = (15, 13), dpi = 300)
        ax = fig.add_subplot(1, 1, 1)
        img = ax.scatter(ravel_x_inds,ravel_y_centers,c=curr_params,
            	    alpha=0.7,s=5,vmin=vmin,vmax=vmax)
	ax.set_xlabel("X")
	ax.set_xlabel("Y")
	ax.set_xlim(1,2048)
	ax.set_ylim(1,2048)
        divider = mpltk.make_axes_locatable(ax)
        cax = divider.append_axes("right", size = "5%", pad = 0.05)
        cbar = plt.colorbar(img, cax = cax, orientation = "vertical",label="Y Width")
        ttl = plt.text(
            0.5, 1.01, "Tele: $(tele), MJD: $(mjd), Chip: $(chip), Expid: $(expid)",
            ha = "center", va = "bottom", transform = ax.transAxes)
        imagePath = dirNamePlots * 
			"final_trace_width_map_$(tele)_$(mjd)_$(expid)_$(chip).png"
        fig.savefig(imagePath, bbox_inches = "tight", pad_inches = 0.1)
	plt.close()

        curr_params = collect(Iterators.flatten(param_outputs[:,:,1]))
        param_summary = quantile(filter(!isnan, curr_params),[0.16, 0.5, 0.84])
        param_summary = [param_summary[2],
     		 param_summary[2]-param_summary[1],
		 param_summary[3]-param_summary[2]]
        vmin = param_summary[1] - 2*param_summary[2]
        vmax = param_summary[1] + 2*param_summary[3]
        fig = plt.figure(figsize = (15, 13), dpi = 300)
        ax = fig.add_subplot(1, 1, 1)
        img = ax.scatter(ravel_x_inds,ravel_y_centers,c=curr_params,
            	    alpha=0.7,s=5,vmin=vmin,vmax=vmax)
	ax.set_xlabel("X")
	ax.set_xlabel("Y")
	ax.set_xlim(1,2048)
	ax.set_ylim(1,2048)
        divider = mpltk.make_axes_locatable(ax)
        cax = divider.append_axes("right", size = "5%", pad = 0.05)
        cbar = plt.colorbar(img, cax = cax, orientation = "vertical",label="Flux")
        ttl = plt.text(
            0.5, 1.01, "Tele: $(tele), MJD: $(mjd), Chip: $(chip), Expid: $(expid)",
            ha = "center", va = "bottom", transform = ax.transAxes)
        imagePath = dirNamePlots * 
			"final_trace_flux_map_$(tele)_$(mjd)_$(expid)_$(chip).png"
        fig.savefig(imagePath, bbox_inches = "tight", pad_inches = 0.1)
	plt.close()
    end

    return nothing 
end

dataPath = "/uufs/chpc.utah.edu/common/home/u6057633/projects/outdir/ap2D/"

tele = "apo"
mjd = "60546"
chip = "c"
expid = "49840013"
ftype = "DOMEFLAT"
fname = dataPath * "ap2D_$(tele)_$(mjd)_$(chip)_$(expid)_$(ftype).jld2"


f = jldopen(fname)
image_data = f["dimage"][1 : 2048, 1 : 2048]
ivar_image = f["ivarimage"][1 : 2048, 1 : 2048]
close(f)

trace_funcs = trace_extract(image_data,ivar_image,tele,mjd,chip,expid;image_mask=nothing)


