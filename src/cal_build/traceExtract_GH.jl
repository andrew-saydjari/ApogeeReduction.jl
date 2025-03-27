using Polynomials: fit, Polynomial
using SpecialFunctions: erf
using Interpolations: linear_interpolation, Line

#profile_path = "/uufs/chpc.utah.edu/common/home/u6057633/scratch/20250226/outdir/trace_profile/"
profile_path = "../../data/"

function int_gauss_hermite_term(x_bins, n; mean = 0.0, width = 1.0, return_deriv = false)
    #integrated version of Gaussian-Hermite terms
    #rule is cdf H(n) = -1*H(n-1)

    zscore = (x_bins .- mean) ./ width

    if !return_deriv
        if n == 0
            #integrate the Gaussian
            return 0.5 .* (1 .+ erf.(zscore ./ (2^0.5)))
        elseif n == 1
            factor = ones(size(zscore))
        elseif n == 2
            factor = zscore
        elseif n == 3
            factor = (zscore .^ 2) .- 1
        elseif n == 4
            factor = (zscore .^ 3) .- 3 .* zscore
        elseif n == 5
            factor = (zscore .^ 4) .- 6 .* (zscore .^ 2) .+ 3
        elseif n == 6
            factor = (zscore .^ 5) .- 10 .* (zscore .^ 3) .+ 15 .* zscore
        elseif n == 7
            factor = (zscore .^ 6) .- 15 .* (zscore .^ 4) .+ 45 .* (zscore .^ 2) .- 15
        elseif n == 8
            factor = (zscore .^ 7) .- 21 .* (zscore .^ 5) .+ 105 .* (zscore .^ 3) .- 105 .* zscore
        elseif n == 9
            factor = (zscore .^ 8) .- 28 .* (zscore .^ 6) .+ 210 .* (zscore .^ 4) .-
                     420 .* (zscore .^ 2) .+ 105
        elseif n == 10
            factor = (zscore .^ 9) .- 36 .* (zscore .^ 7) .+ 378 .* (zscore .^ 5) .-
                     1260 .* (zscore .^ 3) .+ 945 .* zscore
        elseif n == 11
            factor = (zscore .^ 10) .- 45 .* (zscore .^ 8) .+ 630 .* (zscore .^ 6) .-
                     3150 .* (zscore .^ 4) .+ 4725 .* (zscore .^ 2) .- 945
            #    else
            #        raise ValueError(f'ERROR: order n={n} for Hermite polynomial is not an option. Must be leq 11.')
        end

        gauss = exp.(-0.5 .* (zscore .^ 2)) ./ ((2 * π)^0.5) ./ width

        return -1 .* gauss .* factor ./ (width .^ (n - 1))
    else
        if n == 0
            dgauss_dz = (1 / ((2 * π)^0.5)) .* exp.(-0.5 .* (zscore .^ 2))
            #integrate the Gaussian
            return 0.5 .* (1 .+ erf.(zscore ./ (2^0.5))), dgauss_dz
        elseif n == 1
            factor = ones(size(zscore))
            dfactor_dz = zeros(size(zscore))
        elseif n == 2
            factor = zscore
            dfactor_dz = ones(size(zscore))
        elseif n == 3
            factor = (zscore .^ 2) .- 1
            dfactor_dz = 2 .* (zscore)
        elseif n == 4
            factor = (zscore .^ 3) .- 3 .* zscore
            dfactor_dz = 3 .* (zscore .^ 2) .- 3
        elseif n == 5
            factor = (zscore .^ 4) .- 6 .* (zscore .^ 2) .+ 3
            dfactor_dz = 4 .* (zscore .^ 3) .- (6 * 2) .* zscore
        elseif n == 6
            factor = (zscore .^ 5) .- 10 .* (zscore .^ 3) .+ 15 .* zscore
            dfactor_dz = 5 .* (zscore .^ 4) .- (10 * 3) .* (zscore .^ 2) .+ 15
        elseif n == 7
            factor = (zscore .^ 6) .- 15 .* (zscore .^ 4) .+ 45 .* (zscore .^ 2) .- 15
            dfactor_dz = 6 .* (zscore .^ 5) .- (15 * 4) .* (zscore .^ 3) .+ (45 * 2) .* (zscore)
        elseif n == 8
            factor = (zscore .^ 7) .- 21 .* (zscore .^ 5) .+ 105 .* (zscore .^ 3) .- 105 .* zscore
            dfactor_dz = 7 .* (zscore .^ 6) .- (21 * 5) .* (zscore .^ 4) .+
                         (105 * 3) .* (zscore .^ 2) .- 105
        elseif n == 9
            factor = (zscore .^ 8) .- 28 .* (zscore .^ 6) .+ 210 .* (zscore .^ 4) .-
                     420 .* (zscore .^ 2) .+ 105
            dfactor_dz = 8 .* (zscore .^ 7) .- (28 * 6) .* (zscore .^ 5) .+
                         (210 * 4) .* (zscore .^ 3) .- (420 * 2) .* zscore
        elseif n == 10
            factor = (zscore .^ 9) .- 36 .* (zscore .^ 7) .+ 378 .* (zscore .^ 5) .-
                     1260 .* (zscore .^ 3) .+ 945 .* zscore
            dfactor_dz = 9 .* (zscore .^ 8) .- (36 * 7) .* (zscore .^ 6) .+
                         (378 * 5) .* (zscore .^ 4) .- (1260 * 3) .* (zscore .^ 2) .+ 945
        elseif n == 11
            factor = (zscore .^ 10) .- 45 .* (zscore .^ 8) .+ 630 .* (zscore .^ 6) .-
                     3150 .* (zscore .^ 4) .+ 4725 .* (zscore .^ 2) .- 945
            dfactor_dz = 10 .* (zscore .^ 9) .- (45 * 8) .* (zscore .^ 7) .+
                         (630 * 6) .* (zscore .^ 5) .- (3150 * 4) .* (zscore .^ 3) .+
                         (4725 * 2) .* zscore
            #    else
            #        raise ValueError(f'ERROR: order n={n} for Hermite polynomial is not an option. Must be leq 11.')
        end

        gauss = exp.(-0.5 .* (zscore .^ 2)) ./ ((2 * π)^0.5) ./ width
        dgauss_dz = gauss .* (-1 .* zscore)

        return -1 .* gauss .* factor ./ (width .^ (n - 1)),
        -1 .* (gauss .* dfactor_dz .+ dgauss_dz .* factor) ./ (width .^ (n - 1))
    end
end

function gh_profiles(tele, mjd, chip, expid; n_sub = 100, make_plots = false)
    if tele == "apo"
        profile_mjd = "59549"
        profile_expid = "39870035"
    elseif tele == "lco"
        profile_mjd = "60044"
        profile_expid = "44820015"
    end

    profile_fname = profile_path *
                    "quartzTraceProfileParams_$(tele)_$(profile_mjd)_$(profile_expid)_$(chip).jld2"

    prof_params = jldopen(profile_fname)

    prof_fiber_inds = prof_params["fiber_index"] .+ 1
    prof_fiber_centers = prof_params["fiber_median_y_center"] .+ 1
    fiber_inds = collect(minimum(prof_fiber_inds):maximum(prof_fiber_inds))

    prof_width_coeffs = prof_params["y_width_coeffs"]
    gh_order = prof_params["gh_order"][1]
    n_gauss = gh_order + 1
    poly_order = prof_params["poly_order"][1]
    prof_height_coeffs = zeros((n_gauss, poly_order + 1))
    smooth_new_indv_widths = ones(size(fiber_inds, 1))
    smooth_new_indv_heights = zeros((size(fiber_inds, 1), n_gauss))
    for j in 1:n_gauss
        prof_height_coeffs[j, :] .= reverse(prof_params["gh_$(j-1)_height_coeffs"]) #numpy to Julia Polynomials 
        smooth_new_indv_heights[:, j] .= Polynomial(prof_height_coeffs[j, :]).(fiber_inds)
    end

    close(prof_params)

    if make_plots
        dirNamePlots = "../outdir/plots/"

        for j in 1:n_gauss
            fig = Figure(size = (800, 800))
            ax = Axis(fig[1, 1],
                xlabel = "Fiber Index",
                ylabel = "GH Height $(j-1)",
                title = "Trace Profile GH Heights\nTele: $(tele), MJD: $(profile_mjd), Chip: $(chip), Expid: $(profile_expid)")

            scatter!(ax, fiber_inds, smooth_new_indv_heights[:, j])

            tracePlot_heights_Path = dirNamePlots *
                                     "GH$(j-1)_heights_$(tele)_$(profile_mjd)_$(profile_expid)_$(chip).png"
            save(tracePlot_heights_Path, fig)
        end
    end

    med_center_to_fiber_func = fit(prof_fiber_centers, prof_fiber_inds, 3)

    min_prof_fib = minimum(prof_fiber_inds)
    max_prof_fib = maximum(prof_fiber_inds)

    smoothed_cdf_zeros = zeros(size(fiber_inds, 1))
    smoothed_cdf_scales = zeros(size(fiber_inds, 1))
    n_offset_pix = 15

    x = range(start = -n_offset_pix, stop = n_offset_pix, step = 1)
    x_bins = range(start = x[1] - 0.5, stop = x[end] + 0.5, step = 1 / n_sub)
    cdf = zeros(size(x_bins, 1))
    dcdf_dz = zeros(size(x_bins, 1))

    all_y_prof = zeros((size(fiber_inds, 1), size(x_bins, 1)))
    all_y_prof_deriv = zeros((size(fiber_inds, 1), size(x_bins, 1)))

    for ind in 1:size(fiber_inds, 1)
        fiber_ind = fiber_inds[ind]
        cdf[:] .= 0
        dcdf_dz[:] .= 0

        for j in 1:n_gauss
            vals, deriv_vals = int_gauss_hermite_term(x_bins, j - 1, return_deriv = true)
            cdf .+= smooth_new_indv_heights[ind, j] * vals
            dcdf_dz .+= smooth_new_indv_heights[ind, j] * deriv_vals
        end

        smoothed_cdf_zeros[ind] = cdf[1]
        cdf .-= smoothed_cdf_zeros[ind]
        smoothed_cdf_scales[ind] = cdf[end]
        cdf ./= smoothed_cdf_scales[ind]
        dcdf_dz ./= smoothed_cdf_scales[ind]

        all_y_prof[ind, :] .= cdf
        all_y_prof_deriv[ind, :] .= dcdf_dz
    end

    x_prof_min = x_bins[1]
    x_prof_max_ind = size(x_bins, 1)

    if make_plots
        x_centers = 0.5 .* (x_bins[(begin + 1):end] .+ x_bins[begin:(end - 1)])

        fig = Figure(size = (800, 800))
        ax = Axis(fig[1, 1],
            xlabel = "dy (pix)",
            ylabel = "Profile Height",
            limits = ((-5, 5), nothing),
            title = "Trace Profiles\nTele: $(tele), MJD: $(profile_mjd), Chip: $(chip), Expid: $(profile_expid)")

        for j in 1:size(prof_fiber_inds, 1)
            fiber_ind = prof_fiber_inds[j]
            ind = findall(fiber_inds .== fiber_ind)[1]
            cdf[:] .= all_y_prof[ind, :]
            pdf = diff(cdf) ./ diff(x_bins)

            lines!(ax, x_centers, pdf .+ 0.02 * (j - 1))
        end

        tracePlot_heights_Path = dirNamePlots *
                                 "GH_profiles_$(tele)_$(profile_mjd)_$(profile_expid)_$(chip).png"
        save(tracePlot_heights_Path, fig)
    end

    return med_center_to_fiber_func,
    x_prof_min, x_prof_max_ind, n_sub, min_prof_fib, max_prof_fib, all_y_prof, all_y_prof_deriv
end

function cdf_func_indv(x_bins, mean, width, fiber_ind, x_prof_min, x_prof_max_ind,
        n_sub, min_prof_fib, all_y_prof, all_y_prof_deriv)
    z_bins = (x_bins .- mean) ./ width
    z_ints = ceil.(Int, round.((z_bins .- x_prof_min) .* n_sub)) .+ 1
    z_ints = clamp.(z_ints, 1, x_prof_max_ind)

    cdf = all_y_prof[fiber_ind .- min_prof_fib .+ 1, z_ints]
    return cdf
end

function pdf_func_mult(x, mean, width, fiber_inds, x_prof_min, x_prof_max_ind, n_sub,
        min_prof_fib, all_y_prof, all_y_prof_deriv; return_deriv = false)
    x_bins = x[:, 1] .+ range(-0.5, size(x, 2))'
    z_bins = (x_bins .- mean) ./ width
    z_ints = ceil.(Int, round.((z_bins .- x_prof_min) .* n_sub)) .+ 1
    z_ints = clamp.(z_ints, 1, x_prof_max_ind)

    fib_ints = (fiber_inds .- min_prof_fib) .+ 1 .+ zeros(Int, size(z_ints))
    pdfs = diff(getindex.(Ref(all_y_prof), fib_ints, z_ints), dims = 2)
    if !return_deriv
        return pdfs
    else
        dpdf_dmus = diff(getindex.(Ref(all_y_prof_deriv), fib_ints, z_ints), dims = 2) .*
                    (-1 ./ width)
        dpdf_dsigs = diff(
            getindex.(Ref(all_y_prof_deriv), fib_ints, z_ints) .* (-1 .* z_bins ./ width), dims = 2)
        return pdfs, dpdf_dmus, dpdf_dsigs
    end
end

function fit_gaussians(all_rel_fluxes, all_rel_errs, first_guess_params,
        fit_inds, best_model_fit_inds, offset_inds,
        fiber_inds, x_prof_min, x_prof_max_ind,
        n_sub, min_prof_fib, all_y_prof, all_y_prof_deriv;
        n_iter = 10, dmu = 0.01, dsig = 0.01, return_cov = false,
        use_first_guess_heights = false, max_center_move = 3,
        min_widths = 0.5, max_widths = 2.0)
    fit_fluxes = copy(all_rel_fluxes[fit_inds])'
    fit_errs = copy(all_rel_errs[fit_inds])'
    fit_ivars = fit_errs .^ -2
    #    fit_ivars[fit_ivars .== 0] .= 1e-10
    fit_ivars[fit_ivars .== 0] .= minimum(fit_ivars[fit_ivars .!= 0]) / 1000

    curr_guess = copy(first_guess_params)
    new_params = copy(first_guess_params)
    curr_n_peaks = size(curr_guess, 1)
    v_hat = zeros(Float64, (3, curr_n_peaks))
    v_hat_cov = zeros(Float64, (3, 3, curr_n_peaks))

    #(n_pix,n_params,n_peaks)
    M_vects = zeros(Float64,
        (size(fit_fluxes, 1), size(first_guess_params, 2), curr_n_peaks))
    curr_M_T_dot_V_inv = zeros(Float64,
        (size(first_guess_params, 2), size(fit_fluxes, 1)))

    comb_model_fluxes = zeros(size(all_rel_fluxes))
    param_offsets = zeros(Float64, size(curr_guess))

    for r_ind in 1:n_iter
        fit_fluxes .= all_rel_fluxes[fit_inds]'
        fit_errs .= all_rel_errs[fit_inds]'
        fit_ivars = fit_errs .^ -2
        #        fit_ivars[fit_ivars .== 0] .= 1e-10
        fit_ivars[fit_ivars .== 0] .= minimum(fit_ivars[fit_ivars .!= 0]) / 1000

        # CDF version
        model_fluxes_unit_height, dmodel_dmu, dmodel_dsig = pdf_func_mult(
            fit_inds, curr_guess[:, 2], curr_guess[:, 3],
            fiber_inds, x_prof_min, x_prof_max_ind,
            n_sub, min_prof_fib, all_y_prof, all_y_prof_deriv, return_deriv = true)

        best_model_fluxes_unit_height = pdf_func_mult(
            best_model_fit_inds, curr_guess[:, 2], curr_guess[:, 3],
            fiber_inds, x_prof_min, x_prof_max_ind,
            n_sub, min_prof_fib, all_y_prof, all_y_prof_deriv)

        if r_ind == 1
            #then fit only for the heights using the first guesses on centers and widths
            if !use_first_guess_heights
                #don't use the first guess heights to define a best model
                #just use the first guess mean and widths to get the best
                #guess for height
                flux_diffs = copy(fit_fluxes)
            else
                #use the initial provided heights to define best models to
                #subtract off, then get height updates
                model_fluxes = (curr_guess[:, 1] .* model_fluxes_unit_height)

                best_model_fluxes = (curr_guess[:, 1] .* best_model_fluxes_unit_height)
                best_model_fluxes[.!isfinite.(curr_guess[:, 1]), :] .= 0

                comb_model_fluxes .*= 0
                for j in 1:size(best_model_fit_inds, 1)
                    comb_model_fluxes[best_model_fit_inds[j, :]] .+= best_model_fluxes[j, :]
                end

                flux_diffs = fit_fluxes .- model_fluxes'
                flux_diffs .-= comb_model_fluxes[fit_inds]'
                flux_diffs .+= best_model_fluxes[
                    :, offset_inds .+ (size(best_model_fluxes, 2) .÷ 2) .+ 1]'
            end

            #get best fit flux heights
            ivar_heights = (nansum((model_fluxes_unit_height .^ 2)' .* fit_ivars, 1))[1, :]
            mean_heights = nansum(model_fluxes_unit_height' .* flux_diffs .* fit_ivars, 1)[1, :] ./
                           ivar_heights
            #             @show mean_heights

            if !use_first_guess_heights
                #then mean_heights gives the first guess heights
                first_guess_params[:, 1] .= max.(0.01, mean_heights)
            else
                #then mean_heights gives the offset heights to update by
                first_guess_params[:, 1] .= max.(0.01, first_guess_params[:, 1] .+ mean_heights)
            end
            curr_guess .= first_guess_params
        end

        model_fluxes = (curr_guess[:, 1] .* model_fluxes_unit_height)

        dmodel_dmu .*= curr_guess[:, 1]
        dmodel_dsig .*= curr_guess[:, 1]

        best_model_fluxes = (curr_guess[:, 1] .* best_model_fluxes_unit_height)

        best_model_fluxes[.!isfinite.(curr_guess[:, 1]), :] .= 0

        comb_model_fluxes .*= 0
        for j in 1:size(best_model_fit_inds, 1)
            comb_model_fluxes[best_model_fit_inds[j, :]] .+= best_model_fluxes[j, :]
        end

        # remove combined gaussians
        fit_fluxes .-= comb_model_fluxes[fit_inds]'
        # add back in the individual gaussian for each peak so only neighbours are removed
        fit_fluxes .+= best_model_fluxes[
            :, offset_inds .+ (size(best_model_fluxes, 2) .÷ 2) .+ 1]'

        dmodel_dheight = model_fluxes_unit_height

        flux_diffs = fit_fluxes .- model_fluxes'

        # derivative matrix
        M_vects[:, 1, :] .= dmodel_dheight'  # d flux/d height
        M_vects[:, 2, :] .= dmodel_dmu'      # d flux/d mu
        M_vects[:, 3, :] .= dmodel_dsig'     # d flux/d sigma

        for j in 1:curr_n_peaks
            #preallocate curr_M_T_dot_V_inv
            curr_M_T_dot_V_inv[:] = M_vects[:, :, j]' .* fit_ivars[:, j]'

            try
                v_hat[:, j] = (curr_M_T_dot_V_inv * M_vects[:, :, j]) \
                              (curr_M_T_dot_V_inv * flux_diffs[:, j])
                if (return_cov) & (r_ind == n_iter)
                    # only invert matrix on the last iteration
                    v_hat_cov[:, :, j] = inv(curr_M_T_dot_V_inv * M_vects[:, :, j])
                end
            catch
                #put back to the first guess
                v_hat[:, j] .= first_guess_params[j, :] .- curr_guess[j, :]
                v_hat_cov[:, :, j] .= NaN
            end
        end

        new_params[:] = curr_guess .+ v_hat'
        new_params[:, 1] .= max.(new_params[:, 1], 0.01)
        new_params[:, 2] .= clamp.(new_params[:, 2],
            max.(11, first_guess_params[:, 2] .- max_center_move),
            min.(2048 - 11, first_guess_params[:, 2] .+ max_center_move))
        new_params[:, 3] .= clamp.(new_params[:, 3], min_widths, max_widths)

        bad_params = (new_params[:, 3] .<= min_widths .* 1.01) .|
                     (new_params[:, 3] .>= max_widths .* 0.99) .|
                     (new_params[:, 2] .<= first_guess_params[:, 2] .- max_center_move .+ 0.01) .|
                     (new_params[:, 2] .>= first_guess_params[:, 2] .+ max_center_move .- 0.01)
        new_params[bad_params, :] .= first_guess_params[bad_params, :]
        v_hat_cov[:, :, bad_params] .= NaN

        curr_guess .= new_params
    end

    if !return_cov
        return curr_guess
    else
        return curr_guess, permutedims(v_hat_cov, (3, 1, 2))
    end
end

function trace_extract(image_data, ivar_image, tele, mjd, chip, expid,
        med_center_to_fiber_func, x_prof_min, x_prof_max_ind, n_sub, min_prof_fib,
        max_prof_fib, all_y_prof, all_y_prof_deriv;
        good_pixels = nothing, mid = 1025, n_center_cols = 100, verbose = false)
    noise_image = 1 ./ sqrt.(ivar_image)
    if isnothing(good_pixels)
        good_pixels = ones(Bool, size(image_data))
    end
    good_pixels .= good_pixels .& (ivar_image .> 0)
    #     n_center_cols = 100 # +/- n_cols to use from middle to sum to find peaks
    x_center = mid

    x_inds = [-n_center_cols, n_center_cols + 1] .+ x_center

    # Cutout image in x direction
    cutout_fluxes = copy(image_data[
        (x_center - n_center_cols):(x_center + n_center_cols), begin:end])'
    cutout_errs = copy(noise_image[
        (x_center - n_center_cols):(x_center + n_center_cols), begin:end])'
    cutout_masks = copy(good_pixels[
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

    n_smooth_pix = max(5, round(Int, 3 * sigma) + 1)
    smooth_inds = range(start = -n_smooth_pix, stop = n_smooth_pix, step = 1)

    smooth_weights = exp.(-0.5 * (smooth_inds ./ sigma) .^ 2)
    smooth_weights ./= sum(smooth_weights)

    all_smooth_inds = (((smooth_inds' .+ y_vals) .- 1 .+ 2048) .% 2048) .+ 1
    #    smoothed_fluxes = nansum(comb_fluxes[all_smooth_inds]' .* smooth_weights, 1)'
    smoothed_fluxes = nansum(((comb_fluxes .* good_ivars)[all_smooth_inds])' .* smooth_weights, 1)'
    smoothed_fluxes ./= nansum(good_ivars[all_smooth_inds]' .* smooth_weights, 1)'

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
    good_max_inds = (abs.(relative_fluxes .- 1.0) .< 0.5) .& (local_max_waves .>= 11) .&
                    (local_max_waves .<= 2048 - 11)
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

    fit_inds = good_y_vals .+ offset_inds'
    best_model_fit_inds = good_y_vals .+ best_model_offset_inds'

    fit_fluxes = copy(all_rel_fluxes[fit_inds])
    fit_errs = copy(all_rel_errs[fit_inds])
    fit_ivars = fit_errs .^ -2

    curr_fiber_inds = clamp.(ceil.(Int, round.(med_center_to_fiber_func.(float.(good_y_vals)))),
        min_prof_fib, max_prof_fib)

    # first guess parameters
    first_guess_params = zeros(Float64, size(fit_fluxes, 1), 3)
    first_guess_params[:, 1] .= max.(0.01, maximum(fit_fluxes, dims = 2)) # height
    first_guess_params[:, 2] .= good_y_vals # center position index for mu
    first_guess_params[:, 3] .= 0.7 # sigma

    first_guess_params[:, 1] .*= (2 * π)^0.5 * first_guess_params[:, 3] # change to integrated height

    dmu = 0.01
    dsig = 0.01

    dmu = 0.001
    dsig = 0.001

    new_params = fit_gaussians(all_rel_fluxes, all_rel_errs, first_guess_params,
        fit_inds, best_model_fit_inds, offset_inds,
        curr_fiber_inds, x_prof_min, x_prof_max_ind,
        n_sub, min_prof_fib, all_y_prof, all_y_prof_deriv,
        n_iter = 10, dmu = dmu, dsig = dsig, return_cov = false,
        use_first_guess_heights = false, max_center_move = 2,
        min_widths = 0.5, max_widths = 2.0)

    curr_best_params = copy(new_params)

    # use the location of the current peaks to fill in any gaps for low-throughput fibers
    peak_locs = new_params[1:(end - 1), 2]
    peak_spacing = new_params[2:end, 2] .- new_params[1:(end - 1), 2]
    med_peak_spacing = nanmedian(peak_spacing)

    good_peak_spacing = abs.(peak_spacing ./ med_peak_spacing .- 1.0) .< 0.5
    keep_space = copy(good_peak_spacing)

    sigma = 10.0
    n_smooth_pix = round(3 * sigma) + 1
    smooth_inds = range(start = -n_smooth_pix, stop = n_smooth_pix, step = 1)

    smooth_weights = exp.(-0.5 * (smooth_inds ./ sigma) .^ 2)
    smooth_weights ./= sum(smooth_weights)

    peak_counts = range(start = 0,
        stop = size(new_params, 1) - 1,
        step = 1)
    all_smooth_inds = abs.(smooth_inds' .+ peak_counts)
    all_smooth_inds[all_smooth_inds .> size(new_params, 1) .- 1] .= size(new_params, 1) .- 1 .+
                                                                    (size(new_params, 1) .- 1 .-
                                                                     all_smooth_inds[all_smooth_inds .> size(new_params, 1) .- 1])
    all_smooth_inds = round.(Int, all_smooth_inds) .+ 1

    smoothed_widths = nansum(new_params[all_smooth_inds, 3]' .* smooth_weights, 1)'
    smoothed_heights = nansum(new_params[all_smooth_inds, 1]' .* smooth_weights, 1)'
    smoothed_width_locs = copy(new_params[:, 2])

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

        resids = new_params[:, 3] .- smoothed_widths
        resid_summary = quantile(resids[keep_widths], [0.16, 0.5, 0.84])
        resid_summary = [resid_summary[2],
            resid_summary[2] - resid_summary[1],
            resid_summary[3] - resid_summary[2]]
        n_sigma = 3
        keep_widths .= (resids .>= resid_summary[1] - n_sigma * resid_summary[2]) .&
                       (resids .<= resid_summary[1] + n_sigma * resid_summary[3])
        smoothed_widths .= nansum(
            ((new_params[:, 3] .* keep_widths)[all_smooth_inds])' .* smooth_weights, 1)'
        smoothed_widths ./= nansum(keep_widths[all_smooth_inds]' .* smooth_weights, 1)'

        resids = new_params[:, 2] .- smoothed_heights
        resid_summary = quantile(resids[keep_heights], [0.16, 0.5, 0.84])
        resid_summary = [resid_summary[2],
            resid_summary[2] - resid_summary[1],
            resid_summary[3] - resid_summary[2]]
        n_sigma = 3
        keep_widths .= (resids .>= resid_summary[1] - n_sigma * resid_summary[2]) .&
                       (resids .<= resid_summary[1] + n_sigma * resid_summary[3])
        smoothed_heights .= nansum(
            ((new_params[:, 1] .* keep_heights)[all_smooth_inds])' .* smooth_weights, 1)'
        smoothed_heights ./= nansum(keep_heights[all_smooth_inds]' .* smooth_weights, 1)'
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

    peak_offsets = new_params[keep_peaks, 2] .- peak_func.(peak_ints[keep_peaks])
    extrap = linear_interpolation(
        new_params[keep_peaks, 2], peak_offsets, extrapolation_bc = Line())
    all_peak_locs .+= extrap.(all_peak_locs)

    #use more pixels around each peak, otherwise large-width peak fitting fails
    #(see weird discontinuities in height and width vs X)
    #This ONLY works because we are fitting all the peaks simultaneously
    #and removing contamination from neighbouring fibers
    # offset_inds = range(-7, 7, step = 1)
    offset_inds = range(-4, 4, step = 1)

    curr_best_heights = zeros(Float64, size(all_peak_locs)) .+ 0.00001
    for j in 1:size(new_params, 1)
        match_ind = argmin(abs.(all_peak_locs .- new_params[j, 2]))
        curr_best_heights[match_ind] = new_params[j, 1]
    end
    curr_best_widths = zeros(Float64, size(all_peak_locs))
    for j in 1:size(curr_best_widths, 1)
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

    curr_fiber_inds = clamp.(
        ceil.(Int, round.(med_center_to_fiber_func.(all_peak_locs))), min_prof_fib, max_prof_fib)

    new_params = fit_gaussians(all_rel_fluxes, all_rel_errs, first_guess_params,
        fit_inds, best_model_fit_inds, offset_inds,
        curr_fiber_inds, x_prof_min, x_prof_max_ind,
        n_sub, min_prof_fib, all_y_prof, all_y_prof_deriv,
        n_iter = 20, dmu = dmu, dsig = dsig, return_cov = false,
        use_first_guess_heights = true, max_center_move = 3,
        min_widths = 0.8 .* first_guess_params[:, 3],
        max_widths = 1.2 .* first_guess_params[:, 3])

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
    n_smooth_pix = round(3 * sigma) + 1
    smooth_inds = range(start = -n_smooth_pix, stop = n_smooth_pix, step = 1)

    smooth_weights = exp.(-0.5 * (smooth_inds ./ sigma) .^ 2)
    smooth_weights ./= sum(smooth_weights)

    peak_counts = range(start = 0,
        stop = size(new_params, 1) - 1,
        step = 1)
    all_smooth_inds = abs.(smooth_inds' .+ peak_counts)
    all_smooth_inds[all_smooth_inds .> size(new_params, 1) .- 1] .= size(new_params, 1) .- 1 .+
                                                                    (size(new_params, 1) .- 1 .-
                                                                     all_smooth_inds[all_smooth_inds .> size(new_params, 1) .- 1])
    all_smooth_inds = round.(Int, all_smooth_inds) .+ 1

    smoothed_widths = nansum(new_params[all_smooth_inds, 3]' .* smooth_weights, 1)'
    smoothed_heights = nansum(new_params[all_smooth_inds, 1]' .* smooth_weights, 1)'
    smoothed_width_locs = copy(new_params[:, 2])

    keep_widths = ones(Bool, size(new_params, 1))
    keep_heights = ones(Bool, size(new_params, 1))

    for r_ind in 1:2
        resids = new_params[:, 3] .- smoothed_widths
        resid_summary = quantile(resids[keep_widths], [0.16, 0.5, 0.84])
        resid_summary = [resid_summary[2],
            resid_summary[2] - resid_summary[1],
            resid_summary[3] - resid_summary[2]]
        n_sigma = 3
        keep_widths .= (resids .>= resid_summary[1] - n_sigma * resid_summary[2]) .&
                       (resids .<= resid_summary[1] + n_sigma * resid_summary[3])
        smoothed_widths .= nansum(
            ((new_params[:, 3] .* keep_widths)[all_smooth_inds])' .* smooth_weights, 1)'
        smoothed_widths ./= nansum(keep_widths[all_smooth_inds]' .* smooth_weights, 1)'

        resids = new_params[:, 2] .- smoothed_heights
        resid_summary = quantile(resids[keep_heights], [0.16, 0.5, 0.84])
        resid_summary = [resid_summary[2],
            resid_summary[2] - resid_summary[1],
            resid_summary[3] - resid_summary[2]]
        n_sigma = 3
        keep_widths .= (resids .>= resid_summary[1] - n_sigma * resid_summary[2]) .&
                       (resids .<= resid_summary[1] + n_sigma * resid_summary[3])
        smoothed_heights .= nansum(
            ((new_params[:, 1] .* keep_heights)[all_smooth_inds])' .* smooth_weights, 1)'
        smoothed_heights ./= nansum(keep_heights[all_smooth_inds]' .* smooth_weights, 1)'
    end

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
            right_cut_ind = curr_cut_int - 1
        else
            right_cut_ind = size(new_params, 1)
        end
    end

    best_fit_ave_params = best_fit_ave_params[left_cut_ind:right_cut_ind, :]

    #    curr_fiber_inds = ceil.(Int,round.(med_center_to_fiber_func.(best_fit_ave_params[:,2])))
    curr_fiber_inds = clamp.(range(1, size(best_fit_ave_params, 1)), min_prof_fib, max_prof_fib)

    med_flux = nanmedian(best_fit_ave_params[:, 1], 1)
    good_throughput_fibers = (best_fit_ave_params[:, 1] ./ med_flux) .> 0.1
    low_throughput_fibers = findall(.!good_throughput_fibers)

    #     x_inds = range(1, 2048, step = 1)
    #     x_inds = range(1,size(image_data,1),step=1)
    x_inds = axes(image_data, 1)
    final_param_outputs = zeros((
        size(x_inds, 1), size(best_fit_ave_params, 1), size(best_fit_ave_params, 2)))
    final_param_output_covs = zeros((size(x_inds, 1), size(best_fit_ave_params, 1),
        size(best_fit_ave_params, 2), size(best_fit_ave_params, 2)))

    best_fit_ave_params = best_fit_ave_params[good_throughput_fibers, :]
    curr_fiber_inds = curr_fiber_inds[good_throughput_fibers]

    verbose && println("Final number of peaks:", size(best_fit_ave_params))

    curr_best_widths = zeros(Float64, size(best_fit_ave_params, 1))
    for j in 1:size(curr_best_widths, 1)
        match_ind = argmin(abs.(best_fit_ave_params[j, 2] .- smoothed_width_locs))
        curr_best_widths[j] = smoothed_widths[match_ind]
    end
    best_fit_ave_params[:, 3] .= curr_best_widths

    #work from the center outward to use constraints from previous analysis
    sorted_x_inds = x_inds[sortperm(abs.(x_inds .- mid))]

    param_outputs = zeros((
        size(x_inds, 1), size(best_fit_ave_params, 1), size(best_fit_ave_params, 2)))
    param_output_covs = zeros((size(x_inds, 1), size(best_fit_ave_params, 1),
        size(best_fit_ave_params, 2), size(best_fit_ave_params, 2)))

    all_rel_fluxes_mat = copy(image_data)
    all_rel_errs_mat = copy(noise_image)
    all_rel_masks_mat = copy(good_pixels)

    all_rel_fluxes_mat[.!all_rel_masks_mat] .= 0
    all_rel_errs_mat[.!all_rel_masks_mat] .= Inf

    all_rel_ivars_mat = 1 ./ (all_rel_errs_mat .^ 2)
    all_rel_ivars_mat[.!all_rel_masks_mat] .= 0
    #    all_rel_ivars_mat[all_rel_ivars_mat .== 0] .= 1e-10
    all_rel_ivars_mat[all_rel_ivars_mat .== 0] .= minimum(all_rel_ivars_mat[all_rel_ivars_mat .!= 0]) /
                                                  1000
    all_rel_errs_mat = all_rel_ivars_mat .^ -0.5

    first_guess_params = copy(best_fit_ave_params)
    fit_inds = floor.(Int, round.(first_guess_params[:, 2])) .+ offset_inds'
    best_model_fit_inds = floor.(Int, round.(first_guess_params[:, 2])) .+
                          best_model_offset_inds'

    @showprogress enabled=verbose desc="Fitting traces" for (ind, x_ind) in enumerate(sorted_x_inds)
        #use previous analyses to constrain first guesses

        if abs(x_ind - mid) < n_center_cols
            n_use = n_center_cols
        elseif abs(x_ind - mid) < 200
            n_use = 200
        elseif abs(x_ind - mid) < 300
            n_use = 300
        else
            n_use = 500
        end

        if abs(x_ind - mid) < 0.5 * n_use
            first_guess_params .= best_fit_ave_params
        elseif x_ind > mid
            first_guess_params .= nanmedian(
                param_outputs[(x_ind - (n_use + 1)):(x_ind - 1), :, :], 1)[
                1, :, :]
            param_slopes = nanmedian(
                diff(param_outputs[(x_ind - (n_use + 1)):(x_ind - 1), :, :], dims = 1), 1)[1, :, :]
            first_guess_params .+= 0.5 * n_use .* param_slopes
        elseif x_ind < mid
            first_guess_params .= nanmedian(
                param_outputs[(x_ind + 1):(x_ind + (n_use + 1)), :, :], 1)[
                1, :, :]
            param_slopes = nanmedian(
                diff(param_outputs[(x_ind + 1):(x_ind + (n_use + 1)), :, :], dims = 1), 1)[1, :, :]
            first_guess_params .+= -0.5 * n_use .* param_slopes
        end

        first_guess_params[:, 3] .= best_fit_ave_params[:, 3] .*
                                    nanmedian(first_guess_params[:, 3] ./ best_fit_ave_params[:, 3])
        first_guess_params[:, 1] .= best_fit_ave_params[:, 1] .*
                                    nanmedian(first_guess_params[:, 1] ./ best_fit_ave_params[:, 1])

        # first guess parameters
        fit_inds .= floor.(Int, round.(first_guess_params[:, 2])) .+ offset_inds'
        best_model_fit_inds .= floor.(Int, round.(first_guess_params[:, 2])) .+
                               best_model_offset_inds'

        curr_guess = copy(first_guess_params)

        new_params, v_hat_cov = fit_gaussians(all_rel_fluxes_mat[x_ind, :],
            all_rel_errs_mat[x_ind, :],
            first_guess_params,
            fit_inds, best_model_fit_inds, offset_inds,
            curr_fiber_inds, x_prof_min, x_prof_max_ind,
            n_sub, min_prof_fib, all_y_prof, all_y_prof_deriv,
            n_iter = 10, dmu = dmu, dsig = dsig, return_cov = true,
            use_first_guess_heights = true, max_center_move = 1,
            min_widths = 0.5 .* first_guess_params[:, 3],
            max_widths = 2.0 .* first_guess_params[:, 3])

        param_outputs[x_ind, :, :] .= new_params
        param_output_covs[x_ind, :, :, :] .= v_hat_cov
    end

    final_param_outputs[:, good_throughput_fibers, :] .= param_outputs[:, :, :]
    final_param_output_covs[:, good_throughput_fibers, :, :] .= param_output_covs[:, :, :, :]

    #extrapolate from nearby good fibers for low throughput ones
    good_throughput_fibers = findall(good_throughput_fibers)
    final_param_outputs[:, low_throughput_fibers, 1] .= 0
    for j in 1:size(low_throughput_fibers, 1)
        nearest_inds = sortperm(abs.(low_throughput_fibers[j] .- good_throughput_fibers))[[1, 2]]
        nearest_inds = good_throughput_fibers[nearest_inds]
        param_slopes = (diff(final_param_outputs[:, nearest_inds, :], dims = 2)[:, 1, :]) ./
                       (nearest_inds[2] - nearest_inds[1])
        final_param_outputs[:, low_throughput_fibers[j], [2, 3]] .= (final_param_outputs[
            :, nearest_inds[1], [2, 3]] .+
                                                                     param_slopes[:, [2, 3]] .*
                                                                     (low_throughput_fibers[j] -
                                                                      nearest_inds[1]))
    end

    param_outputs = final_param_outputs
    param_output_covs = final_param_output_covs

    return param_outputs, param_output_covs
end
