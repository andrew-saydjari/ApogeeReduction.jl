using HDF5
using Polynomials: Polynomial, fit
using LinearAlgebra: qr

function linear_loss_fit(x, y; wporder = 2, returnL2only = false, linparam = nothing)
    A = positional_poly_mat(x, porder = wporder)
    if isnothing(linparam)
        linparam = A \ y
    end

    if returnL2only
        return sum((y .- A * linparam) .^ 2)
    else
        return y .- A * linparam, linparam
    end
end

function positional_poly_mat(x; porder = 2)
    if porder == 0
        return [ones(length(x))]
    elseif porder == 1
        return [ones(length(x)) x]
    elseif porder == 2
        return [ones(length(x)) x x .^ 2]
    elseif porder == 3
        return [ones(length(x)) x x .^ 2 x .^ 3]
    elseif porder == 4
        return [ones(length(x)) x x .^ 2 x .^ 3 x .^ 4]
    elseif porder == 5
        return [ones(length(x)) x x .^ 2 x .^ 3 x .^ 4 x .^ 5]
    elseif porder == 6
        return [ones(length(x)) x x .^ 2 x .^ 3 x .^ 4 x .^ 5 x .^ 6]
    elseif porder == 7
        return [ones(length(x)) x x .^ 2 x .^ 3 x .^ 4 x .^ 5 x .^ 6 x .^ 7]
    elseif porder == 8
        return [ones(length(x)) x x .^ 2 x .^ 3 x .^ 4 x .^ 5 x .^ 6 x .^ 7 x .^ 8]
    elseif porder == 9
        return [ones(length(x)) x x .^ 2 x .^ 3 x .^ 4 x .^ 5 x .^ 6 x .^ 7 x .^ 8 x .^ 9]
    elseif porder == 10
        return [ones(length(x)) x x .^ 2 x .^ 3 x .^ 4 x .^ 5 x .^ 6 x .^ 7 x .^ 8 x .^ 9 x .^ 10]
    elseif porder == 11
        return [ones(length(x)) x x .^ 2 x .^ 3 x .^ 4 x .^ 5 x .^ 6 x .^ 7 x .^ 8 x .^ 9 x .^ 10 x .^
                                                                                                  11]
    elseif porder == 12
        return [ones(length(x)) x x .^ 2 x .^ 3 x .^ 4 x .^ 5 x .^ 6 x .^ 7 x .^ 8 x .^ 9 x .^ 10 x .^
                                                                                                  11 x .^
                                                                                                     12]
    elseif porder == 13
        return [ones(length(x)) x x .^ 2 x .^ 3 x .^ 4 x .^ 5 x .^ 6 x .^ 7 x .^ 8 x .^ 9 x .^ 10 x .^
                                                                                                  11 x .^
                                                                                                     12 x .^
                                                                                                        13]
    elseif porder == 14
        return [ones(length(x)) x x .^ 2 x .^ 3 x .^ 4 x .^ 5 x .^ 6 x .^ 7 x .^ 8 x .^ 9 x .^ 10 x .^
                                                                                                  11 x .^
                                                                                                     12 x .^
                                                                                                        13 x .^
                                                                                                           14]
    elseif porder == 15
        return [ones(length(x)) x x .^ 2 x .^ 3 x .^ 4 x .^ 5 x .^ 6 x .^ 7 x .^ 8 x .^ 9 x .^ 10 x .^
                                                                                                  11 x .^
                                                                                                     12 x .^
                                                                                                        13 x .^
                                                                                                           14 x .^
                                                                                                              15]
    elseif porder == 16
        return [ones(length(x)) x x .^ 2 x .^ 3 x .^ 4 x .^ 5 x .^ 6 x .^ 7 x .^ 8 x .^ 9 x .^ 10 x .^
                                                                                                  11 x .^
                                                                                                     12 x .^
                                                                                                        13 x .^
                                                                                                           14 x .^
                                                                                                              15 x .^
                                                                                                                 16]
    elseif porder == 17
        return [ones(length(x)) x x .^ 2 x .^ 3 x .^ 4 x .^ 5 x .^ 6 x .^ 7 x .^ 8 x .^ 9 x .^ 10 x .^
                                                                                                  11 x .^
                                                                                                     12 x .^
                                                                                                        13 x .^
                                                                                                           14 x .^
                                                                                                              15 x .^
                                                                                                                 16 x .^
                                                                                                                    17]
    elseif porder == 18
        return [ones(length(x)) x x .^ 2 x .^ 3 x .^ 4 x .^ 5 x .^ 6 x .^ 7 x .^ 8 x .^ 9 x .^ 10 x .^
                                                                                                  11 x .^
                                                                                                     12 x .^
                                                                                                        13 x .^
                                                                                                           14 x .^
                                                                                                              15 x .^
                                                                                                                 16 x .^
                                                                                                                    17 x .^
                                                                                                                       18]
    elseif porder == 19
        return [ones(length(x)) x x .^ 2 x .^ 3 x .^ 4 x .^ 5 x .^ 6 x .^ 7 x .^ 8 x .^ 9 x .^ 10 x .^
                                                                                                  11 x .^
                                                                                                     12 x .^
                                                                                                        13 x .^
                                                                                                           14 x .^
                                                                                                              15 x .^
                                                                                                                 16 x .^
                                                                                                                    17 x .^
                                                                                                                       18 x .^
                                                                                                                          19]
    elseif porder == 20
        return [ones(length(x)) x x .^ 2 x .^ 3 x .^ 4 x .^ 5 x .^ 6 x .^ 7 x .^ 8 x .^ 9 x .^ 10 x .^
                                                                                                  11 x .^
                                                                                                     12 x .^
                                                                                                        13 x .^
                                                                                                           14 x .^
                                                                                                              15 x .^
                                                                                                                 16 x .^
                                                                                                                    17 x .^
                                                                                                                       18 x .^
                                                                                                                          19 x .^
                                                                                                                             20]
    end
end

function dither_nonlinear_loss_fit(
        inparams, x, y, linparam; wporder = 2, returnL2only = false)
    xt = transform_x_chips(x, inparams)
    At = positional_poly_mat(xt, porder = wporder)
    if returnL2only
        return sum((y .- At * linparam) .^ 2)
    else
        return y .- At * linparam
    end
end

# make sure to pass this a copy so you can modify in place
function comb_exp_nonlinear_loss_fit!(
        chipPolyParams, ditherPolyParams,
        inparams, x, y, chipInt, expInt;
        wporder = 2, cporder = 1,
        returnL2only = false, linparam = nothing, dporder = 1,
        dither_only = false, chip_only = false)
    n_fnames = maximum(expInt)
    n_dither_coeffs = (dporder + 1) * (n_fnames - 1)
    if (!dither_only) & (!chip_only)
        params2ChipPolyParams!(chipPolyParams, inparams[(begin + n_dither_coeffs):end], cporder)
        comb_exp_params2DitherPolyParams!(
            ditherPolyParams, inparams[begin:(begin + n_dither_coeffs - 1)], dporder)
    elseif dither_only
        comb_exp_params2DitherPolyParams!(ditherPolyParams, inparams, dporder)
    elseif chip_only
        params2ChipPolyParams!(chipPolyParams, inparams, cporder)
    else
        throw("Cannot have both dither_only and chip_only be true.")
    end
    xt = zeros(Float64, length(x))
    for i in 1:N_CHIPS
        msk = (chipInt .== i)
        xt[msk] .= transform_x_chips(x[msk], chipPolyParams[i, :])
    end
    for e_ind in 1:n_fnames
        msk = (expInt .== e_ind)
        xt[msk] .= transform_x_chips(xt[msk], ditherPolyParams[e_ind, :])
    end
    return linear_loss_fit(
        xt, y, wporder = wporder, returnL2only = returnL2only, linparam = linparam)
end

# make sure to pass this a copy so you can modify in place
function nonlinear_loss_fit!(
        chipPolyParams, inparams, x, y, chipInt; wporder = 2, cporder = 1,
        returnL2only = false, linparam = nothing)
    params2ChipPolyParams!(chipPolyParams, inparams, cporder)
    xt = zeros(Float64, length(x))
    for i in 1:N_CHIPS
        msk = chipInt .== i
        xt[msk] .= transform_x_chips(x[msk], chipPolyParams[i, :])
    end
    return linear_loss_fit(
        xt, y, wporder = wporder, returnL2only = returnL2only, linparam = linparam)
end

function transform_x_chips(x, chipPolyParams)
    porder = length(chipPolyParams) - 1
    if porder == 0
        #assumes 1.0 multiplies x
        return chipPolyParams[1] .+ x
    elseif porder == 1
        return chipPolyParams[1] .+ chipPolyParams[2] .* x
    elseif porder == 2
        return chipPolyParams[1] .+ chipPolyParams[2] .* x .+ chipPolyParams[3] .* x .^ 2
    elseif porder == 3
        return chipPolyParams[1] .+ chipPolyParams[2] .* x .+ chipPolyParams[3] .* x .^ 2 .+
               chipPolyParams[4] .* x .^ 3
    end
end

function comb_exp_params2ChipPolyParams!(chipPolyParams, inparams, cporder)
    n_fnames = size(chipPolyParams, 1)
    n_params = size(chipPolyParams, 3)

    for fname_ind in 1:n_fnames
        for chip_ind in 1:N_CHIPS
            start = 1 + (fname_ind - 1) * N_CHIPS * n_params + (chip_ind - 1) * n_params
            stop = start + n_params - 1
            chipPolyParams[fname_ind, chip_ind, :] .= inparams[start:stop]
        end
    end

    return nothing
end

function comb_exp_ChipPolyParams2Params(chipPolyParams)
    n_fnames = size(chipPolyParams, 1)
    n_params = size(chipPolyParams, 3)

    inparams = zeros(Float64, (N_CHIPS * n_fnames * n_params))
    for fname_ind in 1:n_fnames
        for chip_ind in 1:N_CHIPS
            start = 1 + (fname_ind - 1) * N_CHIPS * n_params + (chip_ind - 1) * n_params
            stop = start + n_params - 1
            inparams[start:stop] .= chipPolyParams[fname_ind, chip_ind, :]
        end
    end

    return inparams
end

function comb_exp_params2DitherPolyParams!(ditherPolyParams, inparams, dporder; skip_first = true)
    n_fnames = size(ditherPolyParams, 1)
    n_params = size(ditherPolyParams, 2)

    if skip_first
        ditherPolyParams[1, :] .= 0.0
        if dporder > 0
            ditherPolyParams[1, 2] = 1.0
        end
        first_ind = 2
    else
        first_ind = 1
    end

    for fname_ind in first_ind:n_fnames
        start = 1 + (fname_ind - first_ind) * n_params
        stop = start + n_params - 1
        ditherPolyParams[fname_ind, :] .= inparams[start:stop]
    end

    return nothing
end

function comb_exp_DitherPolyParams2Params(ditherPolyParams; skip_first = true)
    n_fnames = size(ditherPolyParams, 1)
    n_params = size(ditherPolyParams, 2)

    if skip_first
        first_ind = 2
    else
        first_ind = 1
    end
    inparams = zeros(Float64, ((n_fnames - (first_ind - 1)) * n_params))
    for fname_ind in first_ind:n_fnames
        start = 1 + (fname_ind - first_ind) * n_params
        stop = start + n_params - 1
        inparams[start:stop] .= ditherPolyParams[fname_ind, :]
    end

    return inparams
end

function params2ChipPolyParams!(chipPolyParams, inparams, cporder; fibIndx = nothing)
    if isnothing(fibIndx)
        chipPolyParams[1, :] .= inparams[1:(cporder + 1), :]
        chipPolyParams[3, :] .= inparams[(cporder + 2):(2 * (cporder + 1)), :]
    else
        chipPolyParams[fibIndx, 1, :] .= inparams[fibIndx, 1:(cporder + 1), :]
        chipPolyParams[fibIndx, 3, :] .= inparams[fibIndx, (cporder + 2):(2 * (cporder + 1)), :]
    end
    return nothing
end

function ChipPolyParams2Params(chipPolyParams)
    return vcat(chipPolyParams[1, :], chipPolyParams[3, :])
end

# Sky line wavecal
"""
    sky_wavecal_chipPolyParams0(fname; cporder = 0, fibInds = 1:N_FIBERS)

Initial guesses for the (low-order) chip-transform polynomial parameters of
the per-fiber sky wavelength solution, chosen per telescope from the file
name. Factored out of `get_and_save_sky_wavecal` so offline validation can
run the identical fit.
"""
function sky_wavecal_chipPolyParams0(fname; cporder = 0, fibInds = 1:N_FIBERS)
    # initial guess for the (low-order)chip polynomial parameters
    if occursin("_apo_", fname)
        # chipPolyParams0 = [-1.0716 1.00111
        #                    0 1
        #                    1.07009 0.98803]
        offset_func_chip1 = Polynomial([-1.07221500e+00, 4.52450367e-06])
        scale_func_chip1 = Polynomial([1.00093875e+00, -4.41886670e-07])
        offset_func_chip3 = Polynomial([1.06972294e+00, 2.77444782e-06])
        scale_func_chip3 = Polynomial([9.87857338e-01, 1.09350510e-06])

        #don't use the polynomials above for now...
        offset_func_chip1 = Polynomial([-1.070, 0.0])
        scale_func_chip1 = Polynomial([1.0, 0.0])
        offset_func_chip3 = Polynomial([1.0755, 0.0])
        scale_func_chip3 = Polynomial([1.0, 0.0])
    elseif occursin("_lco_", fname)
        # chipPolyParams0 = [-1.0748 1.00168
        #                    0 1
        #                    1.07089 0.98763]
        offset_func_chip1 = Polynomial([-1.07456222e+00, -1.52100076e-07])
        scale_func_chip1 = Polynomial([1.00123795e+00, 5.30751281e-07])
        offset_func_chip3 = Polynomial([1.07199520e+00, -7.11920517e-06])
        scale_func_chip3 = Polynomial([9.87968936e-01, -2.76150881e-07])

        #don't use the polynomials above for now...
        offset_func_chip1 = Polynomial([-1.07295, 0.0])
        scale_func_chip1 = Polynomial([1.0, 0.0])
        offset_func_chip3 = Polynomial([1.076, 0.0])
        scale_func_chip3 = Polynomial([1.0, 0.0])
    else
        # chipPolyParams0 = [-1.070 1
        #                    0 1
        #                    1.076 1]
        offset_func_chip1 = Polynomial([-1.070, 0.0])
        scale_func_chip1 = Polynomial([1.0, 0.0])
        offset_func_chip3 = Polynomial([1.076, 0.0])
        scale_func_chip3 = Polynomial([1.0, 0.0])
    end

    chipPolyParams0 = zeros(Float64, (size(fibInds, 1), N_CHIPS, cporder + 1))
    chipPolyParams0[:, 1, 1] .= offset_func_chip1.(fibInds)
    chipPolyParams0[:, 2, 1] .= 0.0
    chipPolyParams0[:, 3, 1] .= offset_func_chip3.(fibInds)
    if cporder > 0
        chipPolyParams0[:, 1, 2] .= scale_func_chip1.(fibInds)
        chipPolyParams0[:, 2, 2] .= 1.0
        chipPolyParams0[:, 3, 2] .= scale_func_chip3.(fibInds)
    end
    return chipPolyParams0
end

#function get_and_save_sky_wavecal(fname; cporder = 1, wporder = 2)
function get_and_save_sky_wavecal(fname; cporder = 0, wporder = 4, checkpoint_mode = "commit_same")
    outname = replace(replace(fname, "skyLinePeaks" => "waveCalSkyLine"), "_$(FIRST_CHIP)_" => "_")
    if check_file(outname, mode = checkpoint_mode)
        return outname
    end

    fibInds = 1:N_FIBERS
    chipPolyParams0 = sky_wavecal_chipPolyParams0(fname; cporder = cporder, fibInds = fibInds)

    sky_line_uxlst, sky_line_fwlst, sky_line_chipInt = ingest_skyLines_exp(fname)
    if size(sky_line_uxlst, 1) < 1
        return nothing
    end
    linParams, nlParams, resid_vec,
    loo_resid, loo_dropped, loo_status = get_sky_wavecal(
        sky_line_uxlst, sky_line_fwlst, sky_line_chipInt,
        chipPolyParams0; cporder = cporder, wporder = wporder)
    # remove LOO-rejected assignments from the ingest arrays so the
    # interpolated-parameter constant-offset refit below cannot see them
    for fibIndx in fibInds
        bad = findall(view(loo_dropped, fibIndx, :))
        sky_line_uxlst[bad, fibIndx] .= NaN
        sky_line_fwlst[bad, fibIndx] .= NaN
    end
    n_wholesale = count(loo_status .== LOO_WHOLESALE)
    if n_wholesale > 0
        println("LOO WHOLESALE flag on $(n_wholesale) fibers in $(fname): " *
                "association-scale failure survived to the wavelength solution")
    end
    resid_xt = zeros(Float64, (size(fibInds, 1), size(sky_line_uxlst, 1)))
    fill!(resid_xt, NaN)

    # iterpolate between fibers to
    # regularize the nlParams and linParams
    interp_nlParams, interp_linParams = interpolate_wave_params(
        fibInds, nlParams, linParams)

    # do final pass with interp params to
    # get a new constant-term for wave soln
    # (which we do not think should be interpolated)

    interp_resid_vec = zeros(Float64, (size(fibInds, 1), size(sky_line_uxlst, 1)))
    fill!(interp_resid_vec, NaN)
    interp_resid_xt = zeros(Float64, (size(fibInds, 1), size(sky_line_uxlst, 1)))
    fill!(interp_resid_xt, NaN)

    for fibIndx in fibInds
        xv = sky_line_uxlst[:, fibIndx]
        yv = sky_line_fwlst[:, fibIndx]
        chipIntv = sky_line_chipInt[:, fibIndx]
        sky_msk = .!isnan.(xv)
        if sum(sky_msk) == 0
            continue
        end

        chipPolyParams = copy(chipPolyParams0[fibIndx, :, :])
        inparams = nlParams[fibIndx, :]

        params2ChipPolyParams!(chipPolyParams, inparams, cporder)
        xt = zeros(Float64, length(xv[sky_msk]))
        for i in 1:N_CHIPS
            msk = chipIntv[sky_msk] .== i
            xt[msk] .= transform_x_chips(xv[sky_msk][msk], chipPolyParams[i, :])
        end
        resid_xt[fibIndx, sky_msk] .= xt

        inparams = interp_nlParams[fibIndx, :]

        params2ChipPolyParams!(chipPolyParams, inparams, cporder)
        xt = zeros(Float64, length(xv[sky_msk]))
        for i in 1:N_CHIPS
            msk = chipIntv[sky_msk] .== i
            xt[msk] .= transform_x_chips(xv[sky_msk][msk], chipPolyParams[i, :])
        end
        interp_resid_xt[fibIndx, sky_msk] .= xt

        A = positional_poly_mat(xt, porder = wporder)
        curr_resids = yv[sky_msk] .- A * interp_linParams[fibIndx, :]
        const_offset = mean(curr_resids) #could also be median
        interp_linParams[fibIndx, 1] += const_offset
        interp_resid_vec[fibIndx, sky_msk] .= curr_resids .- const_offset

        # NOTE: using the interp params and doing a final FULL update
        # for all params (e.g. fit full linParams, then get best nlParams)
        # yielded the same outliers in wave soln that we want to avoid
        # so probably shouldn't do it
    end

    chipWaveSoln = zeros(Float64, N_XPIX, N_FIBERS, N_CHIPS)
    interp_chipWaveSoln = zeros(Float64, N_XPIX, N_FIBERS, N_CHIPS)
    x = 1:N_XPIX
    ximport = (x .- (N_XPIX ÷ 2)) ./ N_XPIX
    for chip in CHIP_LIST
        chipIndx = getChipIndx(chip)
        for fibIndx in 1:N_FIBERS
            params2ChipPolyParams!(chipPolyParams0, nlParams, cporder, fibIndx = fibIndx)
            xt = transform_x_chips(ximport, chipPolyParams0[fibIndx, chipIndx, :])
            Ax = positional_poly_mat(xt, porder = wporder)
            yt = Ax * linParams[fibIndx, :]
            chipWaveSoln[:, fibIndx, chipIndx] .= yt

            params2ChipPolyParams!(chipPolyParams0, interp_nlParams, cporder, fibIndx = fibIndx)
            xt = transform_x_chips(ximport, chipPolyParams0[fibIndx, chipIndx, :])
            Ax = positional_poly_mat(xt, porder = wporder)
            yt = Ax * interp_linParams[fibIndx, :]
            interp_chipWaveSoln[:, fibIndx, chipIndx] .= yt
        end
    end

    # the best parameters are likely the interpolated ones
    # but return the raw measurements as well
    safe_jldsave(outname; linParams = interp_linParams, nlParams = interp_nlParams,
        resid_vec = interp_resid_vec, chipWaveSoln = interp_chipWaveSoln,
        resid_xt = interp_resid_xt,
        raw_linParams = linParams, raw_nlParams = nlParams,
        raw_resid_vec = resid_vec, raw_chipWaveSoln = chipWaveSoln,
        raw_resid_xt = resid_xt,
        loo_resid = loo_resid, loo_dropped = collect(loo_dropped),
        loo_status = collect(loo_status), no_metadata = true) # KEVIN FIX THIS!!!!
    return outname
end

#read in a list of fnames, then take the 
#average/median of the wavelength solution
#parameters (removing dither differences)
#to have a stable average for the night
function get_ave_night_wave_soln(
        fname_list; fit_dither = false, wavetype = "sky", dporder = 1, dither_fname_list = nothing)
    #open first to get order of linear and non-linear parameters

    #guard against repetition
    for chip in CHIP_LIST
        fname_list .= replace.(fname_list, "_$(chip)_" => "_")
    end
    unique_fname_list = unique(fname_list)

    f = h5open(unique_fname_list[1], "r+")

    linParams, nlParams = try
        read(f["linParams"]), read(f["nlParams"])
    catch
        println(unique_fname_list[1])
        read(f["linParams"]), read(f["nlParams"])
    end
    close(f)

    n_fnames = size(unique_fname_list, 1)
    n_lin_coeffs = size(linParams, 2)
    n_nl_coeffs = size(nlParams, 2)

    all_linParams = zeros(Float64, (N_FIBERS, n_lin_coeffs, n_fnames))
    all_nlParams = zeros(Float64, (N_FIBERS, n_nl_coeffs, n_fnames))
    all_ditherParams = zeros(Float64, (N_FIBERS, dporder + 1, n_fnames))
    if dporder > 0
        all_ditherParams[:, 2, :] .= 1.0
    end
    fill!(all_ditherParams, NaN)
    fill!(all_linParams, NaN)
    fill!(all_nlParams, NaN)

    for fname_ind in 1:n_fnames
        f = h5open(unique_fname_list[fname_ind], "r+")

        linParams, nlParams = try
            read(f["linParams"]), read(f["nlParams"])
        catch
            println(unique_fname_list[fname_ind])
            read(f["linParams"]), read(f["nlParams"])
        end
        all_linParams[:, :, fname_ind] .= linParams[:, :]
        all_nlParams[:, :, fname_ind] .= nlParams[:, :]
        close(f)
    end

    med_linParams = nanmedian(all_linParams, 3)
    med_nlParams = nanmedian(all_nlParams, 3)

    #remove difference between dithers
    linParam_offsets = nanmedian(all_linParams .- med_linParams, 2)
    nlParam_offsets = nanmedian(all_nlParams .- med_nlParams, 2)

    if fit_dither
        med_nlParams .= nanmedian((all_nlParams .- nlParam_offsets), 3)[:, :, 1]
        med_linParams .= nanmedian((all_linParams .- linParam_offsets), 3)

        #find exposure nearest to the median wave solution
        scores = sum((all_linParams .- med_linParams) .^ 2, dims=(1,2))[1, 1, :]
	scores[.!isfinite.(scores)] .= Inf
        best_exp_ind = argmin(scores)
	if all(isfinite.(all_linParams[:, :, best_exp_ind]))
            #use those parameters to measure dither offsets
            med_linParams .= copy(all_linParams[:, :, best_exp_ind])
	end

        #use the inverse of the dither transformation
        #to put all linParams on the same scale
        #to get a better median solution
        updated_linParams = zeros(Float64, size(all_linParams))
        fill!(updated_linParams, NaN)

        #repeat a second time with better median
        for r_ind in 1:2
            for fname_ind in 1:n_fnames
                curr_path = join(split(unique_fname_list[fname_ind], "/")[begin:(end - 1)], "/")
                curr_tele, curr_mjd, curr_expid, curr_imagetyp = split(
                    split(split(unique_fname_list[fname_ind], "/")[end], ".h5")[1], "_")[(end - 3):end]
                curr_peak_fname = "$(curr_path)/$(wavetype)LinePeaks_$(curr_tele)_$(curr_mjd)_$(curr_expid)_$(FIRST_CHIP)_$(curr_imagetyp).h5"

                dither_outputs = get_sky_dither_per_fiber(
                            curr_peak_fname, med_linParams, med_nlParams; dporder = dporder)
                if isnothing(dither_outputs)
                    continue
                end
                all_ditherParams[:, :, fname_ind] .= dither_outputs[1]

                for fibIndx in 1:N_FIBERS
                    med_wave_poly = Polynomial(all_linParams[fibIndx, :, fname_ind])
                    #use the inverse transformation (i.e. x = -b/a + 1/a * xt)
                    if dporder > 0
                        curr_dither_poly = Polynomial([
                            -all_ditherParams[fibIndx, 1, fname_ind], 1.0] ./
                                                      all_ditherParams[fibIndx, 2, fname_ind])
                    else
                        curr_dither_poly = Polynomial([-all_ditherParams[fibIndx, 1, fname_ind]])
                    end
		    curr_coeffs = med_wave_poly(curr_dither_poly).coeffs
		    n_curr_use_coeffs = min(size(curr_coeffs,1),n_lin_coeffs)
		    if n_curr_use_coeffs == 0
		        continue
		    end

		    updated_linParams[fibIndx, :, fname_ind] .= 0
                    updated_linParams[fibIndx, begin:n_curr_use_coeffs, fname_ind] .= curr_coeffs[begin:n_curr_use_coeffs] 
                end
            end
            med_linParams .= nanmedian(updated_linParams, 3)[:, :, 1]
        end

    else
        med_nlParams = nanmedian((all_nlParams .- nlParam_offsets), 3)[:, :, 1]
        med_linParams = nanmedian((all_linParams .- linParam_offsets), 3)[:, :, 1]
    end

    wporder = n_lin_coeffs - 1
    cporder = Int(n_nl_coeffs / 2) - 1
    chipPolyParams0 = zeros(Float64, (N_FIBERS, N_CHIPS, cporder + 1))
    if cporder > 0
        chipPolyParams0[:, 2, 2] .= 1.0 #chip b scale
    end
    chipWaveSoln = zeros(Float64, N_XPIX, N_FIBERS, N_CHIPS)
    x = 1:N_XPIX
    ximport = (x .- (N_XPIX ÷ 2)) ./ N_XPIX
    for chip in CHIP_LIST
        chipIndx = getChipIndx(chip)
        for fibIndx in 1:N_FIBERS
            params2ChipPolyParams!(chipPolyParams0, med_nlParams, cporder, fibIndx = fibIndx)
            xt = transform_x_chips(ximport, chipPolyParams0[fibIndx, chipIndx, :])
            Ax = positional_poly_mat(xt, porder = wporder)
            yt = Ax * med_linParams[fibIndx, :]
            chipWaveSoln[:, fibIndx, chipIndx] .= yt
        end
    end

    return med_linParams, med_nlParams, chipWaveSoln
end

function get_and_save_sky_dither_per_fiber(
        fname, linParams, nlParams; dporder = 1, wavetype = "sky", max_offset = 1.0, checkpoint_mode = "commit_same")
    outname = replace(replace(fname, "skyLinePeaks" => "waveCalNight$(wavetype)Dither"), "_$(FIRST_CHIP)_" => "_")
    if check_file(outname, mode = checkpoint_mode)
        return outname
    end

    dither_outputs = get_sky_dither_per_fiber(
        fname, linParams, nlParams; dporder = dporder, return_waveSoln = true, max_offset = max_offset)
    if isnothing(dither_outputs)
        return nothing
    end
    ditherParams, resid_vec, chipWaveSoln, resid_xt = dither_outputs
    
    safe_jldsave(outname; linParams, nlParams, ditherParams,
        resid_vec, resid_xt, chipWaveSoln, no_metadata = true) # KEVIN FIX THIS!!!!

    return outname
end

#using the given linParams,nlParams that define a wavelength solution,
#calculate the per-fiber dither shifts
#xt = nlParams applied to orig_x
#xt_new = dither params applied to xt
#wave = linParams applied to xt_new
function get_sky_dither_per_fiber(
        fname, linParams, nlParams; dporder = 1, return_waveSoln = false, max_offset = 1.0)
    fibInds = 1:N_FIBERS
    ditherPolyParams = zeros(Float64, (size(fibInds, 1), dporder + 1))
    if dporder > 0
        ditherPolyParams[:, 2] .= 1.0 #expect very close to 1
    end

    n_lin_coeffs = size(linParams, 2)
    n_nl_coeffs = size(nlParams, 2)
    wporder = n_lin_coeffs - 1
    cporder = Int(n_nl_coeffs / 2) - 1
    chipPolyParams = zeros(Float64, (N_CHIPS, cporder + 1))
    if cporder > 0
        chipPolyParams[2, 2] = 1.0
    end

    #read in sky peak positions and wavelengths
    sky_line_uxlst, sky_line_fwlst, sky_line_chipInt = ingest_skyLines_exp(fname)

    resid_vec = zeros(Float64, N_FIBERS, size(sky_line_uxlst, 1))
    resid_xt = zeros(Float64, N_FIBERS, size(sky_line_uxlst, 1))
    fill!(resid_vec, NaN)

    if size(sky_line_uxlst, 1) == 0
        return nothing
    end

    good_fibers = ones(Bool, N_FIBERS)

    for fibIndx in fibInds
        xv = sky_line_uxlst[:, fibIndx]
        yv = sky_line_fwlst[:, fibIndx]
        chipIntv = sky_line_chipInt[:, fibIndx]
        sky_msk = .!isnan.(xv)
        if sum(sky_msk) == 0
            good_fibers[fibIndx] = false
            continue
        end

        curr_nlParams = nlParams[fibIndx, :]
        params2ChipPolyParams!(chipPolyParams, curr_nlParams, cporder)
        xt = zeros(Float64, length(xv))
        for i in 1:N_CHIPS
            msk = chipIntv .== i
            xt[msk] .= transform_x_chips(xv[msk], chipPolyParams[i, :])
        end
        resid_xt[fibIndx, :] .= xt

        resid_vec[fibIndx, sky_msk] .= dither_nonlinear_loss_fit(ditherPolyParams[fibIndx, :],
            xt[sky_msk], yv[sky_msk], linParams[fibIndx, :];
            wporder = wporder, returnL2only = false)
        good_resids = abs.(resid_vec[fibIndx, sky_msk]) .<= max_offset
        sky_msk[sky_msk] .&= good_resids
        if sum(sky_msk) == 0
            good_fibers[fibIndx] = false
            continue
        end

        function dither_nonlinear_loss_fit_partial(inparams)
            dither_nonlinear_loss_fit(inparams, xt[sky_msk], yv[sky_msk], linParams[fibIndx, :];
                wporder = wporder, returnL2only = true)
        end
        res = optimize(
            dither_nonlinear_loss_fit_partial, ditherPolyParams[fibIndx, :], LBFGS(), Optim.Options(show_trace = false))
        ditherParamsOpt = Optim.minimizer(res)
        ditherPolyParams[fibIndx, :] .= ditherParamsOpt
        resid_vec[fibIndx, sky_msk] .= dither_nonlinear_loss_fit(
            ditherParamsOpt, xt[sky_msk], yv[sky_msk], linParams[fibIndx, :];
            wporder = wporder, returnL2only = false)
    end

    #extrapolate from nearby good fibers for bad ones
    bad_fibers = findall(.!good_fibers)
    good_fibers = findall(good_fibers)

    if size(good_fibers, 1) == 0
        return nothing
    elseif size(good_fibers, 1) == 1
        for j in 1:size(bad_fibers, 1)
            ditherPolyParams[bad_fibers[j], :] .= ditherPolyParams[good_fibers[1], :]
        end
    else
        for j in 1:size(bad_fibers, 1)
            nearest_inds = sortperm(abs.(bad_fibers[j] .- good_fibers))[[1, 2]]
            nearest_inds = good_fibers[nearest_inds]
            param_slopes = (diff(ditherPolyParams[nearest_inds, :], dims = 1))[1, :] ./
                           (nearest_inds[2] - nearest_inds[1])

            ditherPolyParams[bad_fibers[j], :] .= ditherPolyParams[nearest_inds[1], :] .+
                                                  param_slopes .* (bad_fibers[j] - nearest_inds[1])
        end
    end


    if return_waveSoln
        chipWaveSoln = zeros(Float64, N_XPIX, N_FIBERS, N_CHIPS)
        x = 1:N_XPIX
        ximport = (x .- (N_XPIX ÷ 2)) ./ N_XPIX
        for chip in CHIP_LIST
            chipIndx = getChipIndx(chip)
            for fibIndx in 1:N_FIBERS
                params2ChipPolyParams!(chipPolyParams, nlParams[fibIndx, :], cporder)
                xt = transform_x_chips(ximport, chipPolyParams[chipIndx, :])
                xt = transform_x_chips(xt, ditherPolyParams[fibIndx, :])
                Ax = positional_poly_mat(xt, porder = wporder)
                yt = Ax * linParams[fibIndx, :]
                chipWaveSoln[:, fibIndx, chipIndx] .= yt
            end
        end
        return ditherPolyParams, resid_vec, chipWaveSoln, resid_xt
    else
        return ditherPolyParams, resid_vec
    end
end

"""
    loo_linear_residuals(A, y)

Exact leave-one-out (LOO) residuals of the least-squares fit `y ~ A c`,
computed in closed form from one QR factorization: `r_i = e_i / (1 - h_ii)`,
where `e` are the in-sample residuals and `h_ii` the leverages (diagonal of
the hat matrix). `r_i` equals `y_i` minus the prediction of the fit performed
WITHOUT point i -- no per-point refits are needed. Returns `(r_loo, leverage)`;
points whose leverage is within 1e-10 of 1 (the point fully determines the
fit, so its LOO prediction is undefined) carry NaN.
"""
function loo_linear_residuals(A, y)
    F = qr(A)
    c = F \ y
    e = y .- A * c
    Qthin = Matrix(F.Q)
    h = vec(sum(abs2, Qthin, dims = 2))
    r = similar(e)
    for i in eachindex(e)
        r[i] = h[i] < 1 - 1e-10 ? e[i] / (1 - h[i]) : NaN
    end
    return r, h
end

"""
    each_combination(f, n, k)

Call `f(idx)` for every k-element combination `idx` (an Int vector, reused
in place -- copy it if you keep it) of `1:n`. Standard lexicographic
next-combination iteration; no dependencies.
"""
function each_combination(f, n, k)
    k <= n || return nothing
    idx = collect(1:k)
    while true
        f(idx)
        j = k
        while j >= 1 && idx[j] == n - k + j
            j -= 1
        end
        j == 0 && break
        idx[j] += 1
        for l in (j + 1):k
            idx[l] = idx[l - 1] + 1
        end
    end
    return nothing
end

"""
    lms_poly_classify(x, y; porder = 2, nsigma = 7.0, scatter_floor = 0.0)

Least-median-of-squares-style robust outlier classification against a
polynomial trend: every (porder+2)-element subset of the points is fit by
least squares, the winning core is the one minimizing the median absolute
residual over ALL points, the scatter scale is 1.4826 x that minimal median
(bounded below by `scatter_floor`), and a point is kept when its |residual|
against the core fit, normalized by its prediction-leverage factor
`sqrt(1 + h)` (a point far from the core -- e.g. a chip-edge line -- is
predicted by extrapolation and legitimately misses by more), is <= `nsigma`
x scale. Because the core is chosen combinatorially, the classification
survives up to ~half the points being corrupted (where sequential
leave-one-out deletion suffers masking: the corrupted fit inflates every
residual and no single point stands out). The threshold is a ratio against
the data's own scatter -- never an absolute wavelength number. Returns
`(keep::BitVector, core_coeffs, scale)`.
"""
function lms_poly_classify(x, y; porder = 2, nsigma = 7.0, scatter_floor = 0.0)
    n = length(x)
    k = porder + 2
    if n < k + 1
        # too few points to distinguish core from contamination
        A = positional_poly_mat(x, porder = porder)
        return trues(n), A \ y, NaN
    end
    best_score = Inf
    best_coeffs = zeros(porder + 1)
    best_idx = collect(1:k)
    resid = zeros(n)
    Afull = positional_poly_mat(x, porder = porder)
    each_combination(n, k) do idx
        A = positional_poly_mat(x[idx], porder = porder)
        c = try
            A \ y[idx]
        catch
            return
        end
        resid .= abs.(y .- Afull * c)
        score = median(resid)
        if score < best_score
            best_score = score
            best_coeffs .= c
            best_idx .= idx
        end
    end
    s = max(1.4826 * best_score, scatter_floor)
    Acore = positional_poly_mat(x[best_idx], porder = porder)
    keep = trues(n)
    for i in 1:n
        h = i in best_idx ? 0.0 : prediction_leverage(Acore, view(Afull, i, :))
        keep[i] = abs(y[i] - (view(Afull, i, :)' * best_coeffs)) <=
                  nsigma * s * sqrt(1 + h)
    end
    # refit the polynomial on everything kept (more stable than the tiny core)
    ki = findall(keep)
    coeffs = length(ki) >= porder + 1 ?
             positional_poly_mat(x[ki], porder = porder) \ y[ki] : best_coeffs
    return keep, coeffs, s
end

"""
    loo_poly_reject(x, y; porder = 2, nsigma = 7.0, max_drop_frac = 1/3,
                    scatter_floor = 0.0)

Robust polynomial fit `y ~ poly(x)`: an LMS combinatorial core
(`lms_poly_classify`) first protects against heavy contamination, then
sequential leave-one-out rejection refines the kept set. In the sequential
stage a point is an outlier when its |LOO residual| exceeds `nsigma` times
the robust scatter (1.4826 x median |LOO residual|) of the OTHER points --
all thresholds are ratios against the data's own scatter, never absolute
wavelength numbers. Dropping stops when clean or when `max_drop_frac` of
the points are out (`capped = true`; the caller should treat the fit as
suspect). `scatter_floor` bounds the scatter estimate from below (e.g. by a
resolution scale such as the segment-centroid pixel quantization) so a run
of unusually consistent points cannot make the threshold arbitrarily tight.
Returns `(keep::BitVector, coeffs, capped)`.
"""
function loo_poly_reject(x, y; porder = 2, nsigma = 7.0, max_drop_frac = 1 / 3,
        scatter_floor = 0.0)
    n = length(x)
    keep, _, _ = lms_poly_classify(x, y; porder = porder, nsigma = nsigma,
        scatter_floor = scatter_floor)
    maxdrop = max(floor(Int, n * max_drop_frac), n - count(keep))
    capped = n - count(keep) >= maxdrop && count(keep) < n
    while true
        ki = findall(keep)
        A = positional_poly_mat(x[ki], porder = porder)
        r, h = loo_linear_residuals(A, y[ki])
        # studentized: e/sqrt(1-h) has common variance sigma^2, so edge
        # (high-leverage) points are not over-penalized by the 1/(1-h) LOO
        # inflation
        absr = abs.(r) .* sqrt.(max.(1 .- h, 0.0))
        worst, worst_ratio = 0, 0.0
        for i in eachindex(ki)
            isnan(absr[i]) && continue
            others = [absr[j] for j in eachindex(ki) if j != i && !isnan(absr[j])]
            isempty(others) && continue
            s = max(1.4826 * median(others), scatter_floor)
            ratio = s > 0 ? absr[i] / s : Inf
            if ratio > worst_ratio
                worst_ratio, worst = ratio, i
            end
        end
        if worst_ratio > nsigma && worst > 0
            if n - count(keep) >= maxdrop
                capped = true
                break
            end
            keep[ki[worst]] = false
        else
            break
        end
    end
    ki = findall(keep)
    A = positional_poly_mat(x[ki], porder = porder)
    return keep, A \ y[ki], capped
end

"""
    prediction_leverage(A, a)

`a' (A'A)^-1 a` for a design matrix `A` and a prediction row `a` -- the
leverage of an out-of-sample point. The variance of a holdout residual is
`sigma^2 (1 + h)`, so holdout residuals must be normalized by `sqrt(1 + h)`
before comparison against the in-sample scatter (an excluded point at the
edge of the fit range is predicted by extrapolation and legitimately misses
by more than an interior one).
"""
function prediction_leverage(A, a)
    R = qr(A).R
    z = R' \ collect(a)
    return sum(abs2, z)
end

# loo_status codes recorded per fiber in the waveCalSkyLine product
const LOO_OK = Int8(0)        # every measured line passed the LOO test
const LOO_DROPPED = Int8(1)   # >=1 line failed, was dropped from the fit + flagged
const LOO_WHOLESALE = Int8(2) # drop cap hit: association-scale failure, fit suspect
const LOO_UNTESTED = Int8(3)  # too few measured lines to LOO-test

"""
    fit_fiber_sky_wavecal(xv, yv, chipIntv, chipPolyParams; wporder, cporder)

One per-fiber sky wavelength-solution fit (the loop body of
`get_sky_wavecal`): LBFGS over the chip-transform (nonlinear) parameters with
the wavelength polynomial solved linearly at each step. On return
`chipPolyParams` is left consistent with the OPTIMAL nonlinear parameters.
Returns `(nlParamsOpt, linParamsOpt, linResid)`.
"""
function fit_fiber_sky_wavecal(xv, yv, chipIntv, chipPolyParams; wporder = 2, cporder = 1)
    inparams = ChipPolyParams2Params(chipPolyParams)
    function nonlinear_loss_fit_partial(p)
        nonlinear_loss_fit!(chipPolyParams, p, xv, yv, chipIntv;
            wporder = wporder, cporder = cporder, returnL2only = true)
    end
    res = optimize(
        nonlinear_loss_fit_partial, inparams, LBFGS(), Optim.Options(show_trace = false))
    nlParamsOpt = Optim.minimizer(res)
    linResid,
    linParamsOpt = nonlinear_loss_fit!(
        chipPolyParams, nlParamsOpt, xv, yv, chipIntv;
        wporder = wporder, cporder = cporder, returnL2only = false)
    # leave chipPolyParams at the optimum (the last LBFGS evaluation need not be it)
    params2ChipPolyParams!(chipPolyParams, nlParamsOpt, cporder)
    return nlParamsOpt, linParamsOpt, linResid
end

"""
    loo_validate_fiber(xv, yv, chipIntv, chipPolyParams0;
                       wporder = 4, cporder = 0, loo_nsigma = 12.0,
                       chip_screen_nsigma = 12.0, loo_max_drop_frac = 1/3)

Leave-one-out validation of every sky line feeding one fiber's wavelength
solution, in three passes:

1. **Per-chip robust screen**: within each chip, the (x, catalog-wavelength)
   pairs are classified against a quadratic trend with an LMS combinatorial
   core (`lms_poly_classify`). Misassignments -- the measured peak actually
   belongs to a neighboring feature -- sit off the trend by the line
   separation (the known theft modes: 5.7 A for 15546->15540, ~10 A for
   16702->16692, 27.4 A for the 16442 label landing on 16414), while correct
   lines follow it to the within-chip truncation scale (~0.1 A). The
   combinatorial core survives several simultaneous misassignments, where a
   joint least-squares fit would smear the corruption into every residual
   and mask it.
2. **Joint fit + sequential LOO**: the full 3-chip solution is fit on the
   surviving lines, and each line's exact LOO residual (`loo_linear_residuals`
   at the linear stage, nonlinear chip parameters refit between drops) is
   tested against `loo_nsigma` times the robust scatter of the OTHER lines'
   LOO residuals. Correct assignments sit at the measurement-noise level
   (5-30 mA in the healthy corpus) -- 2-4 orders of magnitude below a
   misassignment.
3. **Rescue**: every screened-out line gets a holdout residual against the
   final cleaned solution (it took no part in that fit -- this IS its LOO
   residual); lines consistent with the kept lines' own LOO scatter are
   re-admitted and the solution refit once, so an over-eager per-chip screen
   cannot silently thin the fit.

Every threshold is a ratio against the fiber's own scatter -- never an
absolute wavelength number -- so no hardcoded pixel->wavelength lookup is
load-bearing for correctness. If more than `loo_max_drop_frac` of the lines
end up dropped the fiber is flagged `LOO_WHOLESALE` (association-scale
failure; the self-calibrating re-association in `get_sky_peaks` should have
prevented it).

Returns `(dropped, loo_resid, status, nlParams, linParams, resid)`; `resid`
is NaN for dropped lines (they are excluded from the returned fit), and
`loo_resid` for a dropped line is its holdout residual against the final
cleaned solution.
"""
function loo_validate_fiber(xv, yv, chipIntv, chipPolyParams0;
        wporder = 4, cporder = 0, loo_nsigma = 12.0, chip_screen_nsigma = 12.0,
        loo_max_drop_frac = 1 / 3)
    nline = length(xv)
    dropped = falses(nline)
    loo_resid = fill(NaN, nline)
    resid = fill(NaN, nline)
    msk = .!isnan.(xv)
    n0 = count(msk)
    nparam = (wporder + 1) + 2 * (cporder + 1)
    chipPolyParams = copy(chipPolyParams0)
    if n0 < nparam + 3
        # not enough lines for a meaningful LOO test; fit as before, flag untested
        nl, lin, lr = fit_fiber_sky_wavecal(
            xv[msk], yv[msk], chipIntv[msk], chipPolyParams;
            wporder = wporder, cporder = cporder)
        resid[msk] .= lr
        return dropped, loo_resid, LOO_UNTESTED, nl, lin, resid
    end

    # ---- pass 1: per-chip LMS screen --------------------------------------
    # Two passes: first collect each chip's own LMS scatter scale, then
    # classify with every chip's scale floored by the median of the chips'
    # scales. On a 5-line chip the LMS median is the 3rd-smallest residual of
    # a 4-point core -- essentially interpolation noise -- and without the
    # cross-chip floor the threshold collapses and honest edge lines get
    # flagged; the floor is still derived entirely from this fiber's own
    # measurements (never an absolute wavelength number).
    chip_ci = [findall((chipIntv .== c) .& msk) for c in 1:N_CHIPS]
    chip_s = fill(NaN, N_CHIPS)
    for c in 1:N_CHIPS
        length(chip_ci[c]) >= 5 || continue
        _, _, s = lms_poly_classify(xv[chip_ci[c]], yv[chip_ci[c]];
            porder = 2, nsigma = chip_screen_nsigma)
        chip_s[c] = s
    end
    s_floor = any(.!isnan.(chip_s)) ? median(filter(!isnan, chip_s)) : 0.0
    for c in 1:N_CHIPS
        length(chip_ci[c]) >= 5 || continue
        keep, _, _ = lms_poly_classify(xv[chip_ci[c]], yv[chip_ci[c]];
            porder = 2, nsigma = chip_screen_nsigma, scatter_floor = s_floor)
        dropped[chip_ci[c][.!keep]] .= true
    end
    active = msk .& .!dropped
    if count(active) < nparam + 1
        # screen left too little to fit -- distrust it, fall back to all lines
        dropped .= false
        active = copy(msk)
    end

    # ---- pass 2: joint fit + sequential LOO -------------------------------
    maxdrop = max(floor(Int, n0 * loo_max_drop_frac), count(dropped))
    status = LOO_OK
    local nl, lin, lr, ai, A_active, r_std
    function joint_fit_and_loo!()
        chipPolyParams .= chipPolyParams0
        ai = findall(active)
        nl, lin, lr = fit_fiber_sky_wavecal(
            xv[ai], yv[ai], chipIntv[ai], chipPolyParams;
            wporder = wporder, cporder = cporder)
        xt = zeros(Float64, length(ai))
        for c in 1:N_CHIPS
            m = chipIntv[ai] .== c
            xt[m] .= transform_x_chips(xv[ai][m], chipPolyParams[c, :])
        end
        A_active = positional_poly_mat(xt, porder = wporder)
        r, h = loo_linear_residuals(A_active, yv[ai])
        loo_resid[ai] .= r
        # studentized: e/sqrt(1-h) has common variance sigma^2, so the ratio
        # test does not over-penalize high-leverage (chip-edge) lines, whose
        # raw LOO residuals are inflated by 1/(1-h)
        r_std = r .* sqrt.(max.(1 .- h, 0.0))
        return r_std
    end
    function worst_loo(r)
        absr = abs.(r)
        worst, worst_ratio = 0, 0.0
        for i in eachindex(r)
            isnan(absr[i]) && continue
            others = [absr[j] for j in eachindex(r) if j != i && !isnan(absr[j])]
            isempty(others) && continue
            s = 1.4826 * median(others)
            ratio = s > 0 ? absr[i] / s : Inf
            if ratio > worst_ratio
                worst_ratio, worst = ratio, i
            end
        end
        return worst, worst_ratio
    end
    while true
        r = joint_fit_and_loo!()
        worst, worst_ratio = worst_loo(r)
        if worst_ratio > loo_nsigma && worst > 0
            if count(dropped) >= maxdrop
                status = LOO_WHOLESALE
                break
            end
            dropped[ai[worst]] = true
            active[ai[worst]] = false
        else
            break
        end
    end

    # ---- pass 3: rescue screened lines consistent with the clean fit ------
    di = findall(dropped)
    if !isempty(di) && status != LOO_WHOLESALE
        s_kept = 1.4826 * median(abs.(filter(!isnan, r_std)))
        rescued = false
        for i in di
            xt_i = transform_x_chips(xv[i], chipPolyParams[chipIntv[i], :])
            Ai = positional_poly_mat([xt_i], porder = wporder)
            loo_resid[i] = yv[i] - (Ai * lin)[1]
            # a holdout residual has variance sigma^2 (1 + h_out): normalize
            # before comparing to the kept lines' (studentized) scatter
            h_out = prediction_leverage(A_active, vec(Ai))
            if s_kept > 0 && abs(loo_resid[i]) / sqrt(1 + h_out) <= loo_nsigma * s_kept
                dropped[i] = false
                active[i] = true
                rescued = true
            end
        end
        if rescued
            r = joint_fit_and_loo!()
            # one more sequential guard after the refit
            while true
                worst, worst_ratio = worst_loo(r)
                if !(worst_ratio > loo_nsigma && worst > 0)
                    break
                end
                if count(dropped) >= maxdrop
                    status = LOO_WHOLESALE
                    break
                end
                dropped[ai[worst]] = true
                active[ai[worst]] = false
                r = joint_fit_and_loo!()
            end
        end
    end

    if status != LOO_WHOLESALE && any(dropped)
        status = LOO_DROPPED
    end
    if count(dropped) > n0 ÷ 2
        status = LOO_WHOLESALE
    end
    # cross-chip wholesale tripwire: under the JOINT model the LOO residuals
    # of correct lines are at the (chip-independent) measurement-noise level,
    # so a chip whose kept lines still scatter far above its siblings' has a
    # corrupted association that the within-chip screen could not see (>=
    # half a chip consistently wrong is invisible from inside the chip --
    # the other chips supply the scale anchor, still from the night's own
    # data). Flag only; recovery belongs to the association stage.
    if status != LOO_WHOLESALE
        chip_loo_scale = fill(NaN, N_CHIPS)
        for c in 1:N_CHIPS
            vals = [abs(r_std[k]) for (k, i) in enumerate(ai)
                    if chipIntv[i] == c && !isnan(r_std[k])]
            length(vals) >= 3 && (chip_loo_scale[c] = 1.4826 * median(vals))
        end
        for c in 1:N_CHIPS
            isnan(chip_loo_scale[c]) && continue
            others = [chip_loo_scale[cc]
                      for cc in 1:N_CHIPS if cc != c && !isnan(chip_loo_scale[cc])]
            isempty(others) && continue
            if chip_loo_scale[c] > loo_nsigma * max(median(others), eps())
                status = LOO_WHOLESALE
            end
        end
    end
    resid[ai] .= lr
    # final holdout residuals of everything still dropped
    for i in findall(dropped)
        xt_i = transform_x_chips(xv[i], chipPolyParams[chipIntv[i], :])
        Ai = positional_poly_mat([xt_i], porder = wporder)
        loo_resid[i] = yv[i] - (Ai * lin)[1]
        resid[i] = NaN
    end
    return dropped, loo_resid, status, nl, lin, resid
end

function get_sky_wavecal(
        sky_line_uxlst, sky_line_fwlst, sky_line_chipInt, chipPolyParams0;
        cporder = 1, wporder = 2,
        loo_reject = true, loo_nsigma = 12.0, loo_max_drop_frac = 1 / 3)
    nline = size(sky_line_uxlst, 1)
    linParams = zeros(Float64, N_FIBERS, wporder + 1)
    nlParams = zeros(Float64, N_FIBERS, 2 * (cporder + 1))
    resid_vec = zeros(Float64, N_FIBERS, nline)
    fill!(resid_vec, NaN)
    loo_resid = fill(NaN, N_FIBERS, nline)
    loo_dropped = falses(N_FIBERS, nline)
    loo_status = fill(LOO_UNTESTED, N_FIBERS)
    for i in 1:N_FIBERS
        xv = sky_line_uxlst[:, i]
        yv = sky_line_fwlst[:, i]
        chipIntv = sky_line_chipInt[:, i]
        msk = .!isnan.(xv)
        if count(msk) == 0
            continue
        end
        if loo_reject
            d, rl, st, nlp, lp, rv = loo_validate_fiber(
                xv, yv, chipIntv, chipPolyParams0[i, :, :];
                wporder = wporder, cporder = cporder,
                loo_nsigma = loo_nsigma, loo_max_drop_frac = loo_max_drop_frac)
            loo_dropped[i, :] .= d
            loo_resid[i, :] .= rl
            loo_status[i] = st
            nlParams[i, :] .= nlp
            linParams[i, :] .= lp
            resid_vec[i, :] .= rv
        else
            chipPolyParams = copy(chipPolyParams0[i, :, :])
            nlParamsOpt, linParamsOpt, linResid = fit_fiber_sky_wavecal(
                xv[msk], yv[msk], chipIntv[msk], chipPolyParams;
                wporder = wporder, cporder = cporder)
            nlParams[i, :] = nlParamsOpt
            linParams[i, :] = linParamsOpt
            resid_vec[i, msk] = linResid
        end
    end
    # cross-fiber wholesale tripwire: a fiber whose kept-line residual
    # scatter is far above the exposure's typical fiber has no
    # self-consistent solution. "Consistently bad" is invisible from inside
    # one fiber when most of its lines are corrupted together (the per-fiber
    # ratio tests are scale-free by design); the other fibers of the same
    # exposure supply the missing scale anchor -- still the night's own data,
    # never an absolute number.
    if loo_reject
        s_fib = fill(NaN, N_FIBERS)
        for i in 1:N_FIBERS
            vals = abs.(filter(!isnan, resid_vec[i, :]))
            length(vals) >= 3 && (s_fib[i] = 1.4826 * median(vals))
        end
        med_s = median(filter(!isnan, s_fib))
        if !isnan(med_s) && med_s > 0
            for i in 1:N_FIBERS
                if !isnan(s_fib[i]) && s_fib[i] > loo_nsigma * med_s &&
                   loo_status[i] != LOO_UNTESTED
                    loo_status[i] = LOO_WHOLESALE
                end
            end
        end
    end
    return linParams, nlParams, resid_vec, loo_resid, loo_dropped, loo_status
end

# FPI line wavecal
"""
    get_fpi_exp_m0s(fpi_line_expInt, peak_ints, n_fnames) -> Matrix{Float64}

For each (exposure, fiber), the minimum FPI peak integer m among that fiber's
ingested FPI peaks in that exposure. Entries where a fiber has NO ingested
peaks for an exposure are NaN, which `nanmedian` then skips when the nightly
median m0 is formed.

The matrix is deliberately Float64: it used to be an Int array and the NaN
assignment for a peakless (fiber, exposure) threw an InexactError that killed
the whole night's FPI wavecal, forcing the fallback to the sky solution (A5).
"""
function get_fpi_exp_m0s(fpi_line_expInt, peak_ints, n_fnames)
    n_fibers = size(fpi_line_expInt, 2)
    exp_m0s = fill(NaN, n_fnames, n_fibers)
    for fname_ind in 1:n_fnames
        for fibIndx in 1:n_fibers
            in_exp = findall(fpi_line_expInt[:, fibIndx] .== fname_ind)
            if size(in_exp, 1) == 0
                continue
            end
            exp_m0s[fname_ind, fibIndx] = minimum(peak_ints[in_exp, fibIndx])
        end
    end
    return exp_m0s
end

function comb_exp_get_and_save_fpi_wavecal(
        fname_list, initial_linParams, initial_nlParams;
        cporder = 1, wporder = 4, dporder = 2,
        n_sigma = 5, max_ang_sigma = 0.2, max_iter = 3, verbose = true, checkpoint_mode = "commit_same")
    n_fnames = size(fname_list, 1)
    fname = fname_list[1]
    tele, mjd, expid = split(fname, "_")[(end - 4):(end - 2)]
    outname = replace(
        replace(replace(fname, "fpiPeaks" => "waveCalFPI"), "_$(FIRST_CHIP)_" => "_"),
        "_$(expid)_" => "_")
    if check_file(outname, mode = checkpoint_mode)
        linParams, nlParams, chipWaveSoln = load(outname, "linParams", "nlParams", "chipWaveSoln")
        return outname, linParams, nlParams, chipWaveSoln
    end

    #first guess FPI cavity size, Angstroms
    if occursin("_lco_", fname)
        init_cavity_size = 3.73610e7
        init_m_offset = 0.20
    elseif occursin("_apo_", fname)
        init_cavity_size = 3.7363125e7
        init_m_offset = 0.225
    else
        println("ERROR: No initial cavity size determined for filename $(fname)")
        return nothing
    end

    cavity_size = init_cavity_size
    m_offset = init_m_offset

    initial_n_lin_coeffs = size(initial_linParams, 2)
    initial_n_nl_coeffs = size(initial_nlParams, 2)
    fiber_inds = collect(1:N_FIBERS)

    initial_wporder = initial_n_lin_coeffs - 1
    initial_cporder = Int(initial_n_nl_coeffs / 2) - 1

    #read in the peak positions
    fpi_line_uxlst, fpi_line_uxlst_errs, fpi_line_chipInt, fpi_line_expInt, fpi_line_peakInt = ingest_fpiLines(fname_list)
    n_peaks = size(fpi_line_uxlst, 1)
    println(n_peaks)

    #use the initial wavelength solution and the 
    #chip b peaks to define the m integers 
    #(m = 2*cavity_size/wavelength) of the 
    #peaks on all the chips and for all fibers
    #(they might be offset from each other if
    # some peaks were missed during fitting)

    #transform peak positions to wavelengths
    peak_waves = zeros(Float64, n_peaks, N_FIBERS)
    peak_ints = copy(fpi_line_peakInt) #integer m values
    fpi_ximport = (fpi_line_uxlst .- (N_XPIX ÷ 2)) ./ N_XPIX
    chipPolyParams0 = zeros(Float64, (N_FIBERS, N_CHIPS, initial_cporder + 1))
    if initial_cporder > 0
        chipPolyParams0[:, 2, 2] .= 1.0 #chip b scaling
    end
    for fname_ind in 1:n_fnames
        for chip in CHIP_LIST
            chipIndx = getChipIndx(chip)
            for fibIndx in 1:N_FIBERS
                on_chip = findall((fpi_line_chipInt[:, fibIndx] .== chipIndx) .&
                                  (fpi_line_expInt[:, fibIndx] .== fname_ind))
                chipPolyParams = copy(chipPolyParams0[fibIndx, :, :])
                params2ChipPolyParams!(
                    chipPolyParams, initial_nlParams[fibIndx, :], initial_cporder)
                chipPolyParams0[fibIndx, :, :] .= chipPolyParams
                if size(on_chip, 1) == 0
                    continue
                end

                xt = transform_x_chips(
                    fpi_ximport[on_chip, fibIndx], chipPolyParams[chipIndx, :])
                Ax = positional_poly_mat(xt, porder = initial_wporder)
                yt = Ax * initial_linParams[fibIndx, :]
                peak_waves[on_chip, fibIndx] .= yt

		if size(on_chip,1) == 1
		    int_offset = round(nanmedian(2 * cavity_size ./
                                                                peak_waves[on_chip, fibIndx] .-
                                                                (peak_ints[on_chip, fibIndx] .+
                                                                 m_offset)))
		    if !isfinite(int_offset)
                        int_offset = 0
		    end
                    peak_ints[on_chip, fibIndx] += int_offset
		else
		    int_offset = round(nanmedian(2 * cavity_size ./
                                                                peak_waves[on_chip, fibIndx] .-
                                                                (peak_ints[on_chip, fibIndx] .+
                                                                 m_offset)))
		    if !isfinite(int_offset)
                        int_offset = 0
		    end
                    peak_ints[on_chip, fibIndx] .+= int_offset
		end
            end
        end

    end
    exp_m0s = get_fpi_exp_m0s(fpi_line_expInt, peak_ints, n_fnames)
    med_m0 = round(nanmedian(nanmedian(exp_m0s, 2), 1)[1, 1])
    expect_peak_waves = 2 .* cavity_size ./ (peak_ints .+ m_offset)
    verbose && println("First FPI peak integer m0 = $(med_m0)")

    test_peak_ints = round.(2 .* cavity_size ./ peak_waves .- m_offset)
    #    bad_peaks = (abs.(test_peak_ints .- peak_ints) .> 5) .| (abs.(peak_waves .- expect_peak_waves) .> 2.0)
    #    bad_peaks = (abs.(peak_waves .- expect_peak_waves) .> 2.0) .| isnan.(fpi_ximport)
    bad_peaks = (abs.(test_peak_ints .- peak_ints) .> 5) .| isnan.(fpi_ximport) .| isnan.(peak_ints)
    #    bad_peaks = isnan.(fpi_ximport)
    verbose &&
        println("Found $(sum(bad_peaks)) bad peaks (out of a total $(length(bad_peaks)) peaks) that are too offset from the expected peak positions")
    good_peaks = (.!isnan.(fpi_ximport)) .& (.!bad_peaks)

    linParams = zeros(Float64, N_FIBERS, wporder + 1)
    #every image but the first will have an offset and scaling for dithers
    nlParams = zeros(Float64, N_FIBERS, 2 * (cporder + 1))
    ditherParams = zeros(Float64, N_FIBERS, (n_fnames - 1) * (dporder + 1))
    resid_vec = zeros(Float64, N_FIBERS, size(fpi_line_uxlst, 1))
    fill!(resid_vec, NaN)
    resid_xt = zeros(Float64, N_FIBERS, size(fpi_line_uxlst, 1))
    fill!(resid_xt, NaN)
    n_lin_coeffs = size(linParams, 2)
    n_nl_coeffs = size(nlParams, 2)
    n_dither_coeffs = size(ditherParams, 2)

    linParams[:, begin:size(initial_linParams, 2)] .= initial_linParams
    ditherPolyParams = zeros(Float64, (n_fnames, N_FIBERS, dporder + 1))
    if dporder > 0
        ditherPolyParams[:, :, 2] .= 1.0
    end

    first_dither_offsets = zeros(N_FIBERS)
    for i in 1:N_FIBERS
        on_first = findall((fpi_line_expInt[:, i] .== 1))
        if size(on_first, 1) == 0
            continue
        elseif .!any(good_peaks[on_first, i])
            continue
        end
        dlam = (expect_peak_waves[on_first, i] .- peak_waves[on_first, i])[(begin + 1):end]
        dlam_dx = diff(peak_waves[on_first, i]) ./ diff(fpi_ximport[on_first, i])
        dx = nanmedian(dlam ./ dlam_dx)
        first_dither_offsets[i] = dx
        #	println(i," ",dx," ",dx * 2048)
    end
    #    println(nanmedian(first_dither_offsets)," ",nanmedian(first_dither_offsets)*2048)
    first_dither_offset = nanmedian(first_dither_offsets)
    dither_poly = Polynomial([first_dither_offset, 1])

    fname_ind = 1
    for i in 1:N_FIBERS
        #change wavelength solution to account for first dither offset

	curr_coeffs = Polynomial(linParams[i, :])(dither_poly).coeffs
	n_curr_use_coeffs = min(size(curr_coeffs,1),n_lin_coeffs)
	linParams[i, :] .= 0
        linParams[i, begin:n_curr_use_coeffs] .= curr_coeffs[begin:n_curr_use_coeffs]
        for chip in CHIP_LIST
            chipIndx = getChipIndx(chip)
            on_chip = findall((fpi_line_chipInt[:, i] .== chipIndx) .&
                              (fpi_line_expInt[:, i] .== fname_ind))
            chipPolyParams = copy(chipPolyParams0[i, :, :])
            params2ChipPolyParams!(
                chipPolyParams, initial_nlParams[i, :], initial_cporder)
            chipPolyParams0[i, :, :] .= chipPolyParams
            if size(on_chip, 1) == 0
                continue
            end

            xt = transform_x_chips(
                fpi_ximport[on_chip, i], chipPolyParams[chipIndx, :])
            Ax = positional_poly_mat(xt, porder = wporder)
            yt = Ax * linParams[i, :]
            peak_waves[on_chip, i] .= yt
        end
    end

    #    first_dither_offsets = zeros(N_FIBERS)
    #    for i in 1:N_FIBERS
    #        on_first = findall((fpi_line_expInt[:, i] .== 1))
    #        if size(on_first,1) == 0
    #            continue
    #        elseif .!any(good_peaks[on_first, i])
    #            continue
    #        end
    #	dlam = (expect_peak_waves[on_first,i] .- peak_waves[on_first,i])[begin+1:end]
    #	dlam_dx = diff(peak_waves[on_first,i]) ./ diff(fpi_ximport[on_first,i])
    #	dx = nanmedian(dlam ./ dlam_dx)
    #	first_dither_offsets[i] = dx
    #	println(i," ",dx," ",dx * 2048)
    #    end
    #    println(nanmedian(first_dither_offsets)," ",nanmedian(first_dither_offsets)*2048)
    #    first_dither_offset = nanmedian(first_dither_offsets)

    #measure dither offset from first observation
    dither_offsets = zeros(Float64, (N_FIBERS, N_CHIPS, n_fnames))
    fill!(dither_offsets, NaN)
    dither_offsets[:, :, 1] .= 0.0

    for i in 1:N_FIBERS
        for chip in CHIP_LIST
            chipIndx = getChipIndx(chip)
            on_first = findall((fpi_line_chipInt[:, i] .== chipIndx) .&
                               (fpi_line_expInt[:, i] .== 1))
            if size(on_first, 1) == 0
                continue
            elseif .!any(good_peaks[on_first, i])
                continue
            end

            n_first = sum(good_peaks[on_first, i])
            for fname_ind in 2:n_fnames
                on_other = findall((fpi_line_chipInt[:, i] .== chipIndx) .&
                                   (fpi_line_expInt[:, i] .== fname_ind))
                if size(on_other, 1) == 0
                    continue
                elseif .!any(good_peaks[on_other, i])
                    continue
                end
                n_other = sum(good_peaks[on_other, i])

	        n_use = min(n_first,n_other)

		peak_offset = nanmedian(peak_ints[on_first,i][begin:n_use]-peak_ints[on_other,i][begin:n_use])
		if isfinite(peak_offset) 
    		    peak_offset = ceil(Int,peak_offset)
		end

		if !isfinite(peak_offset)
		    pix_offset = NaN
		elseif peak_offset >= 0
		    new_peak_offset = nanmedian(peak_ints[on_first,i][begin:n_use-peak_offset]-peak_ints[on_other,i][begin+peak_offset:n_use])
		    if isfinite(new_peak_offset)
                        new_peak_offset = ceil(Int,new_peak_offset)
                        pix_offset = nanmedian(fpi_ximport[on_first,i][begin:n_use-peak_offset]-fpi_ximport[on_other,i][begin+peak_offset:n_use])
		    else
                        new_peak_offset = NaN
			pix_offset = NaN
		    end
		else
		    peak_offset = abs(peak_offset)
		    new_peak_offset = nanmedian(peak_ints[on_first,i][begin+peak_offset:n_use]-peak_ints[on_other,i][begin:n_use-peak_offset])
                    if isfinite(new_peak_offset)
                        new_peak_offset = ceil(Int,new_peak_offset)
                        pix_offset = nanmedian(fpi_ximport[on_first,i][begin+peak_offset:n_use]-fpi_ximport[on_other,i][begin:n_use-peak_offset])
		        peak_offset *= -1
		    else
                        new_peak_offset = NaN
			pix_offset = NaN
		    end
		end
		# println(i," ",chipIndx," ",fname_ind," ",peak_offset," ",new_peak_offset," ",pix_offset," ",pix_offset * 2048)

                dither_offsets[i, chipIndx, fname_ind] = pix_offset
            end
        end
    end
    exp_med_dithers = nanmedian(nanmedian(dither_offsets, 2), 1)[1, 1, :]
    # println("dithers ",size(exp_med_dithers),exp_med_dithers,exp_med_dithers .* 2048)

    for i in 1:N_FIBERS
        start = 1
        if cporder > 0
            nlParams[i, start:(start + (initial_cporder + 1) - 1)] .= chipPolyParams0[i, 1, :]
            start = 1 + (cporder + 1)
            nlParams[i, start:(start + (initial_cporder + 1) - 1)] .= chipPolyParams0[i, 3, :]
            if initial_cporder == 0
                start = 1
                nlParams[i, (start + (initial_cporder + 1) - 1) + 1] = 1.0
                start = 1 + (cporder + 1)
                nlParams[i, (start + (initial_cporder + 1) - 1) + 1] = 1.0
            end
        else
            nlParams[i, 1] = chipPolyParams0[i, 1, 1]
            nlParams[i, 2] = chipPolyParams0[i, 3, 1]
        end
        ditherPolyParams[:, i, 1] .= exp_med_dithers
        ditherParams[i, :] .= comb_exp_DitherPolyParams2Params(ditherPolyParams[:, i, :])
    end

    chipPolyParams0 = zeros(Float64, (N_FIBERS, N_CHIPS, cporder + 1))
    if cporder > 0
        chipPolyParams0[:, 2, 2] .= 1.0 #chip b scaling
    end

    for fname_ind in 2:n_fnames
        for fibIndx in 1:N_FIBERS
            for chip in CHIP_LIST
                chipIndx = getChipIndx(chip)
                on_chip = findall((fpi_line_chipInt[:, fibIndx] .== chipIndx) .&
                                  (fpi_line_expInt[:, fibIndx] .== fname_ind))
                chipPolyParams = copy(chipPolyParams0[fibIndx, :, :])
                params2ChipPolyParams!(
                    chipPolyParams, nlParams[fibIndx, :], cporder)
                chipPolyParams0[fibIndx, :, :] .= chipPolyParams
                if size(on_chip, 1) == 0
                    continue
                end

                xt = transform_x_chips(
                    fpi_ximport[on_chip, fibIndx], chipPolyParams[chipIndx, :])
                #		println(fname_ind,chip,fibIndx,size(xt),ditherPolyParams[fname_ind, fibIndx, :],chipPolyParams[chipIndx,:])
                xt = transform_x_chips(xt, ditherPolyParams[fname_ind, fibIndx, :])
                Ax = positional_poly_mat(xt, porder = wporder)
                yt = Ax * linParams[fibIndx, :]
                peak_waves[on_chip, fibIndx] .= yt
            end
        end
    end

    #need to fit additional parameters
    #to the chip offset/scaling and wave soln,
    #specifically, the FPI cavity size and
    #the non-integer offset to m0 integer
    #(should be smaller than 1 in size)

    #    #estimate non-integer offset to m0
    #    m_offset = nanmedian(2.0 * cavity_size ./ peak_waves[good_peaks] .- peak_ints[good_peaks])
    #    #update cavity size guess
    #    cavity_size = nanmedian(peak_waves[good_peaks] .* (peak_ints[good_peaks] .+ m_offset) ./ 2.0)
    #    verbose && println("$(m_offset) $(cavity_size) $(cavity_size-init_cavity_size)")

    #mask outliers
    expect_peak_waves .= 2 .* cavity_size ./ (peak_ints .+ m_offset)
    resids = (peak_waves .- expect_peak_waves)

    r_ind = 0
    verbose && println("Fitting FPI peaks for best wavelength solutions")
    verbose &&
        println("ind   num_good_peaks/num_total   m_offset   cavity_size   delta_from_init_cavity_size   resid_summary_[16,50,84]")
    verbose && println(
        "$(r_ind) $(sum(good_peaks))/$(length(good_peaks)) $(m_offset) $(cavity_size) $(cavity_size-init_cavity_size) ",
        nanzeropercentile(resids[good_peaks], percent_vec = [16, 50, 84]))

    for r_ind in 1:max_iter

        #outlier masking

        for fname_ind in 1:n_fnames
            for chip in CHIP_LIST
                chipIndx = getChipIndx(chip)
                for i in 1:N_FIBERS
                    on_chip = findall((fpi_line_chipInt[:, i] .== chipIndx) .&
                                      (fpi_line_expInt[:, i] .== fname_ind))
                    if size(on_chip, 1) == 0
                        continue
                    elseif .!any(good_peaks[on_chip, i])
                        continue
                    end
                    resid_summary = nanzeropercentile(
                        resids[on_chip, i][good_peaks[on_chip, i]], percent_vec = [16, 50, 84])
                    resid_summary = [resid_summary[2],
                        min(max_ang_sigma, 0.5 * (resid_summary[3] - resid_summary[1])),
                        min(max_ang_sigma, 0.5 * (resid_summary[3] - resid_summary[1]))]
                    good_peaks[on_chip, i] .= (resids[on_chip, i] .>=
                                               resid_summary[1] - n_sigma * resid_summary[2]) .&
                                              (resids[on_chip, i] .<=
                                               resid_summary[1] + n_sigma * resid_summary[3]) .&
                                              (.!bad_peaks[on_chip, i])
                end
            end
        end

        good_fibers = any(good_peaks, dims = 1)[1, :]
        bad_fibers = .!good_fibers
        #	if r_ind == 1
        #            #best fit non-integer offset to m0
        #            m_offset = nanmedian(2.0 * cavity_size ./ peak_waves[good_peaks] .- peak_ints[good_peaks])
        #            #update cavity size guess
        #            cavity_size = nanmedian(peak_waves[good_peaks] .* (peak_ints[good_peaks] .+ m_offset) ./ 2.0)
        #            #expected wavelengths
        #            expect_peak_waves .= 2 * cavity_size ./ (peak_ints .+ m_offset)
        #	else
        if true
            #extrapolate from nearby good fibers for bad ones
            good_fibers = findall(good_fibers)
            bad_fibers = findall(bad_fibers)

            for j in 1:size(bad_fibers, 1)
                if length(good_fibers) >= 2
                    nearest_inds = sortperm(abs.(bad_fibers[j] .- good_fibers))[[1, 2]]
                    nearest_inds = good_fibers[nearest_inds]
                    param_slopes = (diff(linParams[nearest_inds, :], dims = 1))[1, :] ./
                                   (nearest_inds[2] - nearest_inds[1])
                    linParams[bad_fibers[j], :] .= linParams[nearest_inds[1], :] .+
                                                   param_slopes .* (bad_fibers[j] - nearest_inds[1])
                    
                    #for the linear 0th order term, measure change
                    #from the skyline/initial solution and use that 
                    change = diff(linParams[nearest_inds, 1] .- initial_linParams[nearest_inds, 1])[1] ./
                             (nearest_inds[2] - nearest_inds[1])

                    linParams[bad_fibers[j], 1] = initial_linParams[bad_fibers[j], 1] .+
                                                  (linParams[nearest_inds[1], 1] .-
                                                   initial_linParams[nearest_inds[1], 1]) .+
                                                  change .* (bad_fibers[j] - nearest_inds[1])

                    param_slopes = (diff(nlParams[nearest_inds, :], dims = 1))[1, :] ./
                                   (nearest_inds[2] - nearest_inds[1])
                    nlParams[bad_fibers[j], :] .= nlParams[nearest_inds[1], :] .+
                                                  param_slopes .* (bad_fibers[j] - nearest_inds[1])
                    
                    param_slopes = (diff(ditherParams[nearest_inds, :], dims = 1))[1, :] ./
                                   (nearest_inds[2] - nearest_inds[1])
                    ditherParams[bad_fibers[j], :] .= ditherParams[nearest_inds[1], :] .+
                                                      param_slopes .* (bad_fibers[j] - nearest_inds[1])
                elseif length(good_fibers) == 1
                    # If only one good fiber, just copy its parameters
                    linParams[bad_fibers[j], :] .= linParams[good_fibers[1], :]
                    nlParams[bad_fibers[j], :] .= nlParams[good_fibers[1], :]
                    ditherParams[bad_fibers[j], :] .= ditherParams[good_fibers[1], :]
                else
                    # If no good fibers, keep the initial parameters (already set)
                    continue
                end
            end
        end

        #update chip offset and wavelength solution parameters
        for i in 1:N_FIBERS
            resid_vec[i, :] .= NaN
            xv = @view fpi_ximport[:, i]
            yv = @view expect_peak_waves[:, i]
            chipIntv = @view fpi_line_chipInt[:, i]
            expIntv = @view fpi_line_expInt[:, i]
            msk = @view good_peaks[:, i]
            chipPolyParams = copy(chipPolyParams0[i, :, :])
            params2ChipPolyParams!(
                chipPolyParams, nlParams[i, :], cporder)
            curr_ditherPolyParams = copy(ditherPolyParams[:, i, :])
            comb_exp_params2DitherPolyParams!(
                curr_ditherPolyParams, ditherParams[i, :], dporder)
            #	    inparams = vcat(ditherParams[i,:], nlParams[i, :])
            inparams = nlParams[i, :]

            xt = zeros(Float64, length(xv))
            for chipIndx in 1:N_CHIPS
                chip_msk = chipIntv .== chipIndx
                xt[chip_msk] .= transform_x_chips(xv[chip_msk], chipPolyParams[chipIndx, :])
            end

            for e_ind in 1:n_fnames
                in_exp = findall(expIntv .== e_ind)

                curr_xt = @view xt[in_exp]
                curr_y = @view yv[in_exp]
                curr_msk = @view msk[in_exp]

                function dither_nonlinear_loss_fit_partial(inparams)
                    dither_nonlinear_loss_fit(
                        inparams, curr_xt[curr_msk], curr_y[curr_msk], linParams[i, :];
                        wporder = wporder, returnL2only = true)
                end
                res = optimize(
                    dither_nonlinear_loss_fit_partial, curr_ditherPolyParams[e_ind, :], LBFGS(), Optim.Options(show_trace = false))
                curr_ditherParamsOpt = Optim.minimizer(res)
                curr_ditherPolyParams[e_ind, :] .= curr_ditherParamsOpt
                if e_ind == 1
                    #push first exposure dither params onto the wavelength solution
                    #		    println(i," ",curr_ditherParamsOpt," poly ",Polynomial(curr_ditherParamsOpt))
                    #		    new_coeffs = Polynomial(linParams[i,:])(Polynomial(curr_ditherParamsOpt)).coeffs
                    #		    println(linParams[i,:],new_coeffs)
                    linParams[i, :] .= (Polynomial(linParams[i, :])(Polynomial(curr_ditherParamsOpt))).coeffs[begin:n_lin_coeffs]
                    curr_ditherPolyParams[e_ind, :] .= 0.0
                    if dporder > 0
                        curr_ditherPolyParams[e_ind, 2] = 1.0
                    end
                    curr_ditherParamsOpt .= curr_ditherPolyParams[e_ind, :]
                end
                resid_vec[i, in_exp[curr_msk]] .= dither_nonlinear_loss_fit(
                    curr_ditherParamsOpt, curr_xt[curr_msk], curr_y[curr_msk], linParams[i, :];
                    wporder = wporder, returnL2only = false)
            end
            #transfer the first dither parameters to the others
            if dporder > 1
                curr_ditherPolyParams[:, 2] ./= curr_ditherPolyParams[1, 2]
                curr_ditherPolyParams[:, 1] .-= curr_ditherPolyParams[1, 1]
            elseif dporder == 1
                curr_ditherPolyParams .= [curr_ditherPolyParams[:, 1] .-
                                          curr_ditherPolyParams[1, 1] .*
                                          curr_ditherPolyParams[:, 2] ./ curr_ditherPolyParams[1, 2] curr_ditherPolyParams[:,
                    2] ./ curr_ditherPolyParams[1,
                    2]]
            else
                curr_ditherPolyParams .-= curr_ditherPolyParams[1, 1]
            end

            ditherParams[i, :] .= comb_exp_DitherPolyParams2Params(curr_ditherPolyParams)

            function comb_exp_nonlinear_loss_fit_partial(inparams)
                comb_exp_nonlinear_loss_fit!(chipPolyParams, curr_ditherPolyParams, inparams,
                    xv[msk], yv[msk], chipIntv[msk], expIntv[msk];
                    wporder = wporder, cporder = cporder, returnL2only = true, dporder = dporder, chip_only = true)
            end

            if sum(msk) < (n_lin_coeffs + n_nl_coeffs + n_dither_coeffs)
                #then turn off fitting for this peak
                #and don't fit
                good_peaks[:, i] .= false
                linParamsOpt = @view linParams[i, :]
                nlParamsOpt = @view nlParams[i, :]
                ditherParamsOpt = ditherParams[i, :]
                #		outparamsOpt = vcat(ditherParamsOpt, nlParamsOpt)
                outparamsOpt = nlParamsOpt
            else
                res = optimize(
                    comb_exp_nonlinear_loss_fit_partial, inparams, LBFGS(), Optim.Options(show_trace = false))
                outparamsOpt = Optim.minimizer(res)
                #		ditherParamsOpt = outparamsOpt[begin:begin+n_dither_coeffs-1]
                #		nlParamsOpt = outparamsOpt[begin+n_dither_coeffs:end] 
                ditherParamsOpt = @view ditherParams[i, :]
                nlParamsOpt = outparamsOpt
                linResid,
                linParamsOpt = comb_exp_nonlinear_loss_fit!(
                    chipPolyParams, curr_ditherPolyParams, outparamsOpt, xv[msk],
                    yv[msk], chipIntv[msk], expIntv[msk];
                    wporder = wporder, cporder = cporder, returnL2only = false, dporder = dporder, chip_only = true)
            end
            #get resids for all the peaks, even masked ones
            #but don't repeat fit (by passing linParamsOpt)
            linResid,
            linParamsOpt = comb_exp_nonlinear_loss_fit!(
                chipPolyParams, curr_ditherPolyParams, outparamsOpt, xv, yv, chipIntv, expIntv;
                wporder = wporder, cporder = cporder, returnL2only = false,
                linparam = linParamsOpt, dporder = dporder, chip_only = true)
            nlParams[i, :] .= nlParamsOpt
            linParams[i, :] .= linParamsOpt
            ditherParams[i, :] .= ditherParamsOpt
            resid_vec[i, :] .= linResid
            curr_ditherPolyParams = copy(ditherPolyParams[:, i, :])
            comb_exp_params2DitherPolyParams!(
                curr_ditherPolyParams, ditherParams[i, :], dporder)
        end

        #change peak_waves to be measured peak wavelengths (not expected)
        peak_waves .= expect_peak_waves .- resid_vec'

        mults = [1.0e3, 1.0]
        expparams0 = [0.0, m_offset] ./ mults
        function nonlinear_expparams_fit(expparams)
            return nansum((peak_waves[good_peaks] .-
                           (2.0 * (cavity_size + expparams[1] * mults[1])) ./
                           (peak_ints[good_peaks] .+ expparams[2] * mults[2])) .^ 2)
        end
        res = optimize(
            nonlinear_expparams_fit, expparams0, LBFGS(), Optim.Options(show_trace = false))
        expparamsOpt = Optim.minimizer(res) .* mults
        cavity_size += expparamsOpt[1]
        m_offset = expparamsOpt[2]

        #        #EM definition
        #        m_offset = nanmedian(2.0 * cavity_size ./ peak_waves .- peak_ints)
        #        cavity_size = nanmedian(peak_waves .* (peak_ints .+ m_offset) ./ 2.0)

        #expected wavelengths
        expect_peak_waves .= 2 * cavity_size ./ (peak_ints .+ m_offset)

        #outlier rejection
        resids .= (peak_waves .- expect_peak_waves)
        resid_vec .= resids'

        verbose && println(
            "$(r_ind) $(sum(good_peaks))/$(length(good_peaks)) $(m_offset) $(cavity_size) $(cavity_size-init_cavity_size) ",
            nanzeropercentile(resids[good_peaks], percent_vec = [16, 50, 84]))
    end

    chipWaveSoln = zeros(Float64, N_XPIX, N_FIBERS, N_CHIPS)
    x = 1:N_XPIX
    ximport = (x .- (N_XPIX ÷ 2)) ./ N_XPIX
    for chip in CHIP_LIST
        chipIndx = getChipIndx(chip)
        for fibIndx in 1:N_FIBERS
            chipPolyParams = copy(chipPolyParams0[fibIndx, :, :])
            params2ChipPolyParams!(chipPolyParams, nlParams[fibIndx, :], cporder)
            curr_ditherPolyParams = copy(ditherPolyParams[:, fibIndx, :])
            comb_exp_params2DitherPolyParams!(
                curr_ditherPolyParams, ditherParams[fibIndx, :], dporder)
            ditherPolyParams[:, fibIndx, :] .= curr_ditherPolyParams
            xt = transform_x_chips(ximport, chipPolyParams[chipIndx, :])
            Ax = positional_poly_mat(xt, porder = wporder)
            yt = Ax * linParams[fibIndx, :]
            chipWaveSoln[:, fibIndx, chipIndx] .= yt

            curr_keep = findall(fpi_line_chipInt[:, fibIndx] .== chipIndx)
            resid_xt[fibIndx, curr_keep] .= transform_x_chips(
                fpi_ximport[curr_keep, fibIndx], chipPolyParams[chipIndx, :])
        end
    end

    #loop over the exposures and make the expected waveCalNightfpiDither output
    for fname_ind in 1:n_fnames
        curr_outname = replace(replace(
            replace(fname_list[fname_ind], "fpiPeaks" => "waveCalNightfpiDither"),
            "_$(FIRST_CHIP)_" => "_"))

        for chip in CHIP_LIST
            chipIndx = getChipIndx(chip)
            for fibIndx in 1:N_FIBERS
                chipPolyParams = copy(chipPolyParams0[fibIndx, :, :])
                params2ChipPolyParams!(chipPolyParams, nlParams[fibIndx, :], cporder)
                xt = transform_x_chips(ximport, chipPolyParams[chipIndx, :])
                xt = transform_x_chips(xt, ditherPolyParams[fname_ind, fibIndx, :])
                Ax = positional_poly_mat(xt, porder = wporder)
                yt = Ax * linParams[fibIndx, :]
                chipWaveSoln[:, fibIndx, chipIndx] .= yt
            end
        end

        in_exp = findall(fpi_line_expInt[1, :] .== fname_ind)
        safe_jldsave(curr_outname;
            linParams = linParams, nlParams = nlParams, ditherParams = ditherPolyParams[
                fname_ind, :, :],
            resid_vec = resid_vec[:, in_exp], resid_xt = resid_xt[:, in_exp],
            chipWaveSoln = chipWaveSoln, no_metadata = true) # KEVIN FIX THIS!!!!
    end
    safe_jldsave(outname; linParams = linParams, nlParams = nlParams, ditherParams = ditherParams,
        exposure_names = fname_list, resid_vec = resid_vec, resid_xt = resid_xt,
        resid_chipInts = fpi_line_chipInt, resid_peak_ints = peak_ints,
        resid_exp_ints = fpi_line_expInt, resid_used_in_fit = good_peaks,
        chipWaveSoln = chipWaveSoln,
        fpi_m0 = med_m0, fpi_m0_offset = m_offset, fpi_cavity_size = cavity_size, no_metadata = true) # KEVIN FIX THIS!!!!

    return outname, linParams, nlParams, chipWaveSoln
end

function get_fpi_wavecal(
        fpi_line_uxlst, fpi_line_uxlst_errs, peak_ints,
        fpi_line_chipInt, cavity_size, m_offset,
        initial_linParams, initial_nlParams,
        chipPolyParams0; cporder = 1, wporder = 4)
    linParams = zeros(Float64, N_FIBERS, wporder + 1)
    nlParams = zeros(Float64, N_FIBERS, 2 * (cporder + 1) + 2) #extra 2 for integer offset and cavity size
    resid_vec = zeros(Float64, N_FIBERS, size(fpi_line_uxlst, 1))
    fill!(resid_vec, NaN)

    linParams[:, begin:(begin + size(initial_linParams, 2))] .= initial_linParams
    nlParams[:, begin:(begin + size(initial_nlParams, 2))] .= initial_nlParams

    for chip in CHIP_LIST
        chipIndx = getChipIndx(chip)
        for fibIndx in 1:N_FIBERS
            on_chip = findall(fpi_line_chipInt[:, fibIndx] .== chipIndx)
            params2ChipPolyParams!(
                chipPolyParams0, initial_nlParams, initial_cporder, fibIndx = fibIndx)
            xt = transform_x_chips(
                fpi_ximport[on_chip, fibIndx], chipPolyParams0[fibIndx, chipIndx, :])
            Ax = positional_poly_mat(xt, porder = initial_wporder)
            yt = Ax * initial_linParams[fibIndx, :]
            peak_waves[on_chip, fibIndx] .= yt

            peak_ints[on_chip, fibIndx] .= collect(1:size(on_chip, 1)) .- 1
            fiber_m0s[fibIndx] = nanmedian(2 * cavity_size ./ peak_waves[on_chip, fibIndx] .-
                                           peak_ints[on_chip, fibIndx])
            peak_ints[on_chip, fibIndx] .+= round(fiber_m0s[fibIndx])
        end
    end

    #best fit non-integer offset to m0
    m_offset = nanmedian(2 * cavity_size ./ peak_waves .- peak_ints)
    #update cavity size guess
    cavity_size = nanmedian(peak_waves .* (peak_ints .+ m_offset) ./ 2.0)

    for i in 1:N_FIBERS
        xv = fpi_line_uxlst[:, i]
        yv = fpi_line_fwlst[:, i]
        chipIntv = fpi_line_chipInt[:, i]
        msk = .!isnan.(xv)
        chipPolyParams = copy(chipPolyParams0[i, :, :])
        inparams = ChipPolyParams2Params(chipPolyParams)
        function nonlinear_loss_fit_partial(inparams)
            nonlinear_loss_fit!(chipPolyParams, inparams, xv[msk], yv[msk], chipIntv[msk];
                wporder = wporder, cporder = cporder, returnL2only = true)
        end
        res = optimize(
            nonlinear_loss_fit_partial, inparams, LBFGS(), Optim.Options(show_trace = false))
        nlParamsOpt = Optim.minimizer(res)
        linResid,
        linParamsOpt = nonlinear_loss_fit!(
            chipPolyParams, nlParamsOpt, xv[msk], yv[msk], chipIntv[msk];
            wporder = wporder, cporder = cporder, returnL2only = false)
        nlParams[i, :] = nlParamsOpt
        linParams[i, :] = linParamsOpt
        resid_vec[i, msk] = linResid
    end

    #append the 
    shared_params = [cavity_size, m_offset]

    for i in 1:N_FIBERS
        xv = fpi_line_uxlst[:, i]
        yv = fpi_line_fwlst[:, i]
        chipIntv = fpi_line_chipInt[:, i]
        msk = .!isnan.(xv)
        chipPolyParams = copy(chipPolyParams0[i, :, :])
        inparams = ChipPolyParams2Params(chipPolyParams)
        function nonlinear_loss_fit_partial(inparams)
            nonlinear_loss_fit!(chipPolyParams, inparams, xv[msk], yv[msk], chipIntv[msk];
                wporder = wporder, cporder = cporder, returnL2only = true)
        end
        res = optimize(
            nonlinear_loss_fit_partial, inparams, LBFGS(), Optim.Options(show_trace = false))
        nlParamsOpt = Optim.minimizer(res)
        linResid,
        linParamsOpt = nonlinear_loss_fit!(
            chipPolyParams, nlParamsOpt, xv[msk], yv[msk], chipIntv[msk];
            wporder = wporder, cporder = cporder, returnL2only = false)
        nlParams[i, :] = nlParamsOpt
        linParams[i, :] = linParamsOpt
        resid_vec[i, msk] = linResid
    end
    return linParams, nlParams, cavity_size, m_offset, resid_vec
end

function fit_poly_without_outliers(xvals, yvals, deg; nsigma = 3, max_repeat = 5, keep = nothing)
    summary = nanzeropercentile(yvals, percent_vec = [16.0, 50.0, 84.0])
    #change to median,-sigma,+sigma
    summary = [summary[2], summary[2] - summary[1], summary[3] - summary[2]]

    if isnothing(keep)
        keep = ones(Bool, size(xvals, 1))
    end
    keep = keep .&
           (yvals .>= summary[1] - nsigma * summary[2]) .&
           (yvals .<= summary[1] + nsigma * summary[3])

    for repeatInd in 1:max_repeat
        func = fit(xvals[keep], yvals[keep], deg)
        curr_resids = yvals .- func.(xvals)

        resid_summary = nanzeropercentile(curr_resids[keep], percent_vec = [16.0, 50.0, 84.0])
        resid_summary = [resid_summary[2],
            resid_summary[2] - resid_summary[1],
            resid_summary[3] - resid_summary[2]]

        new_good = (curr_resids .>= resid_summary[1] - nsigma * resid_summary[2]) .&
                   (curr_resids .<= resid_summary[1] + nsigma * resid_summary[3])

        if all(new_good .== keep)
            break
        end
        keep .= new_good
    end

    func = fit(xvals[keep], yvals[keep], deg)

    return func, keep
end

function interpolate_wave_params(fibInds, nlParams, linParams)
    mtp_inds = floor.((fibInds .- 1) ./ 30) .+ 1
    unique_mtp = 1:10

    interp_linParams = zeros(Float64, size(linParams))
    interp_nlParams = zeros(Float64, size(nlParams))

    #constant term in wave soln needs to fit each mtp separately
    #need a second order poly for zeroth order term
    for mtp_ind in unique_mtp
        keep_mtp = (mtp_inds .== mtp_ind)
        func, keep = fit_poly_without_outliers(fibInds[keep_mtp], linParams[keep_mtp, 1], 2)
        func_eval = func.(fibInds[keep_mtp])
        interp_linParams[keep_mtp, 1] .= func_eval .+
                                         nanmedian((linParams[keep_mtp, 1] .- func_eval)[keep])
    end

    #need a second order poly for first order term
    j = 2
    func, keep = fit_poly_without_outliers(fibInds, linParams[:, j], 2)
    func_eval = func.(fibInds)
    interp_linParams[:, j] .= func_eval .+ nanmedian((linParams[:, j] .- func_eval)[keep])

    #need a first order poly for higher order terms
    for j in 3:size(linParams, 2)
        func, keep = fit_poly_without_outliers(fibInds, linParams[:, j], 1)
        func_eval = func.(fibInds)
        interp_linParams[:, j] .= func_eval .+ nanmedian((linParams[:, j] .- func_eval)[keep])
    end

    #use a first order poly for the chip offset,scale parameters
    for j in 1:size(nlParams, 2)
        func, keep = fit_poly_without_outliers(fibInds, nlParams[:, j], 1)
        func_eval = func.(fibInds)
        interp_nlParams[:, j] .= func_eval .+ nanmedian((nlParams[:, j] .- func_eval)[keep])
    end

    return interp_nlParams, interp_linParams
end

function ingest_fpiLines_file(fileName)
    # Read in fpi line peaks data
    f = h5open(fileName, "r+")
    fpi_line_mat, fpi_line_cov_mat = try
        read(f["fpi_line_mat"]), read(f["fpi_line_cov_mat"])
    catch
        println(fileName)
        read(f["fpi_line_mat"]), read(f["fpi_line_cov_mat"])
    end
    close(f)

    fpi_line_xlst = fpi_line_mat[:, 2, :]
    fpi_line_xlst_errs = zeros(Float64, size(fpi_line_xlst))
    fill!(fpi_line_xlst_errs, NaN)
    good_ivars = (fpi_line_cov_mat[:, 2, 2, :] .> 0) .& (fpi_line_cov_mat[:, 2, 2, :] .< 100)
    fpi_line_xlst_errs[good_ivars] .= sqrt.(fpi_line_cov_mat[:, 2, 2, :][good_ivars])
    fpi_line_xlst[.!good_ivars] .= NaN
    return fpi_line_xlst, fpi_line_xlst_errs
end

function ingest_skyLines_file(fileName)
    # Read in sky line peaks data
    f = h5open(fileName, "r+")
    sky_line_mat_clean = try
        read(f["sky_line_mat_clean"])
    catch
        println(fileName)
        read(f["sky_line_mat_clean"])
    end
    # wavecal class per line (see WAVECAL_CLASS_CODES): only class B (code 2)
    # enters the wavelength solution; class-A guards are fit for association
    # protection and QA but NaN-masked here, which removes them from every
    # downstream fit (per-exposure, nightAve, dither) via the existing NaN
    # handling. Files written before the class column existed carry no
    # line_class dataset and are ingested unmasked (all lines were class B).
    line_class = if haskey(f, "line_class")
        read(f["line_class"])
    else
        fill(2, size(sky_line_mat_clean, 1))
    end
    close(f)
    sky_line_xlst = (sky_line_mat_clean[:, 1, :] .- (N_XPIX ÷ 2)) ./ N_XPIX
    sky_line_wlst = sky_line_mat_clean[:, 2, :]
    not_for_solution = line_class .!= 2
    sky_line_xlst[not_for_solution, :] .= NaN
    sky_line_wlst[not_for_solution, :] .= NaN
    return sky_line_xlst, sky_line_wlst
end

# this takes in a filename and replaces the chip index (make "a" default approx)
function ingest_skyLines_exp(fname; chip_lst = CHIP_LIST)
    # println("Ingesting sky lines for $fname")
    sky_line_uxlst = Matrix{Float64}[]
    sky_line_fwlst = Matrix{Float64}[]
    sky_line_chipInt = Matrix{Int}[]
    for chip in CHIP_LIST
        fnameloc = replace(fname, "_$(FIRST_CHIP)_" => "_$(chip)_")
        if isfile(fnameloc)
            sky_line_xlst, sky_line_wlst = ingest_skyLines_file(fnameloc)
            chipIndx = getChipIndx(chip)
            push!(sky_line_uxlst, sky_line_xlst)
            push!(sky_line_fwlst, sky_line_wlst)
            push!(sky_line_chipInt, chipIndx * ones(Int, size(sky_line_wlst)))
        else
            println("Sky line file $fnameloc does not exist")
            push!(sky_line_uxlst, Matrix{Float64}(undef, 0, 0))
            push!(sky_line_fwlst, Matrix{Float64}(undef, 0, 0))
            push!(sky_line_chipInt, Matrix{Int}(undef, 0, 0))
        end
    end
    sky_line_uxlst = vcat(sky_line_uxlst...)
    sky_line_fwlst = vcat(sky_line_fwlst...)
    sky_line_chipInt = vcat(sky_line_chipInt...)

    if size(sky_line_uxlst, 1) > 1
        msk_large_scatter = dropdims(nanzeroiqr(sky_line_uxlst, 2) .> 0.002, dims = 2)
        sky_line_uxlst[msk_large_scatter, :] .= NaN
        sky_line_fwlst[msk_large_scatter, :] .= NaN
    end
    # dims are num_sky_lines x num_fibers
    return sky_line_uxlst, sky_line_fwlst, sky_line_chipInt
end

# this takes in a filename and replaces the chip index (make "a" default approx)
function ingest_fpiLines_exp(fname)
    fpi_line_uxlst = Matrix{Float64}[]
    fpi_line_uxlst_errs = Matrix{Float64}[]
    fpi_line_chipInt = Matrix{Int}[]
    for chip in CHIP_LIST
        fnameloc = replace(fname, "_$(FIRST_CHIP)_" => "_$(chip)_")
        if isfile(fnameloc)
            fpi_line_xlst, fpi_line_xlst_errs = ingest_fpiLines_file(fnameloc)
            chipIndx = getChipIndx(chip)
            push!(fpi_line_uxlst, fpi_line_xlst)
            push!(fpi_line_uxlst_errs, fpi_line_xlst_errs)
            push!(fpi_line_chipInt, chipIndx .* ones(Int, size(fpi_line_xlst)))
        else
            println("$(fnameloc) is not a file")
            push!(fpi_line_uxlst, [])
            push!(fpi_line_uxlst_errs, [])
            push!(fpi_line_chipInt, [])
        end
    end
    fpi_line_uxlst = vcat(fpi_line_uxlst...)
    fpi_line_uxlst_errs = vcat(fpi_line_uxlst_errs...)
    fpi_line_chipInt = vcat(fpi_line_chipInt...)

    # dims are num_fpi_lines x num_fibers
    return fpi_line_uxlst, fpi_line_uxlst_errs, fpi_line_chipInt
end

function ingest_fpiLines(fname_list)
    fpi_line_uxlst = Matrix{Float64}[]
    fpi_line_uxlst_errs = Matrix{Float64}[]
    fpi_line_chipInt = Matrix{Int}[]
    fpi_line_expInt = Matrix{Int}[]
    fpi_line_peakInt = Matrix{Int}[]
    for fname_ind in 1:size(fname_list, 1)
        fname = fname_list[fname_ind]
        curr_out = ingest_fpiLines_exp(fname)
        push!(fpi_line_uxlst, curr_out[1])
        push!(fpi_line_uxlst_errs, curr_out[2])
        push!(fpi_line_chipInt, curr_out[3])
        push!(fpi_line_peakInt,
            (collect(1:size(curr_out[1], 1)) .- 1) .* ones(Int, size(curr_out[1])))
        push!(fpi_line_expInt, fname_ind .* ones(Int, size(curr_out[1])))
    end
    fpi_line_uxlst = vcat(fpi_line_uxlst...)
    fpi_line_uxlst_errs = vcat(fpi_line_uxlst_errs...)
    fpi_line_chipInt = vcat(fpi_line_chipInt...)
    fpi_line_expInt = vcat(fpi_line_expInt...)
    fpi_line_peakInt = vcat(fpi_line_peakInt...)

    # dims are num_fpi_lines x num_fibers
    return fpi_line_uxlst, fpi_line_uxlst_errs, fpi_line_chipInt, fpi_line_expInt, fpi_line_peakInt
end

## TODO local parallelization mode (which would also need to pamp inside of sky_wave_plots)
function skyline_medwavecal_skyline_dither(tele, mjd, mjd_list_wavecal, all1DObjectWavecal, all1DObjectSkyPeaks, checkpoint_mode; outdir = "../outdir")

    outname = joinpath(
        outdir, "wavecal", "wavecalNightAve_$(tele)_$(mjd).h5")
    mkpath(dirname(outname))
    
    if check_file(outname, mode = checkpoint_mode)
        #then skip making new file
        return outname
    end

    mskMJD = (mjd_list_wavecal .== mjd)
    if size(all1DObjectWavecal[mskMJD], 1) > 0
        all1DObjectWavecal_mjd = all1DObjectWavecal[mskMJD]
        # using all skyline wavelength solutions to determine median solution for SJD
        night_linParams, night_nlParams, night_wave_soln = get_ave_night_wave_soln(
            all1DObjectWavecal_mjd, fit_dither = true)

        # using skylines to measure dither offsets from nightly skyline wavelength solution for SJD
        get_and_save_sky_dither_per_fiber_partial(fname) = get_and_save_sky_dither_per_fiber(
                fname, night_linParams, night_nlParams; dporder = 1, wavetype = "sky", max_offset = 1.0)
        
        map(get_and_save_sky_dither_per_fiber_partial, all1DObjectSkyPeaks[mskMJD])

        # plotting skyline wavelength solution diagnostic figures
        # make sure plot directory exists
        mkpath(joinpath(outdir, "plots", string(mjd)))
        sky_wave_plots(
            all1DObjectWavecal[mskMJD], night_linParams, night_nlParams, night_wave_soln,
            dirNamePlots = joinpath(outdir, "plots/"),
            plot_fibers = (1, 50, 100, 150, 200, 250, 300),
            plot_pixels = (1, 512, 1024, 1536, 2048))

        curr_wave_type = "sky"

        h5open(outname, "w") do f
            #should include metadata at some point...
            f["best_wave_type"] = curr_wave_type
            f["$(curr_wave_type)/used_fnames"] = all1DObjectWavecal_mjd
            f["$(curr_wave_type)/nightAve_wave_soln"] = night_wave_soln
            f["$(curr_wave_type)/nightAve_nlParams"] = night_nlParams
            f["$(curr_wave_type)/nightAve_linParams"] = night_linParams
        end	

        return outname
    else
        return nothing
    end    
end

"""
    fpi_medwavecal_skyline_dither(mjd, mjd_list_fpi, mjd_list_wavecal,
        all1DfpiPeaks_a, all1DObjectSkyPeaks, wavecalNightAve_fname;
        verbose = true, checkpoint_mode = "commit_same", almanac_file = "")

Promote the night's wavelength solution from the sky-line solution to the
FPI solution -- but only if the FPI solution passes the acceptance gate in
`fpi_gate.jl`. See that file for the two guards, the thresholds and the corpus
they were measured against.

This promotion used to be unconditional and silent, which is how five LCO
nights with a dark FPI were delivered with a wavelength solution wrong by
~1.2 px RMS. When either guard fails the night is NOT failed: the sky solution
already written by `skyline_medwavecal_skyline_dither` is kept, no
`waveCalNightfpiDither` files are written for the night's object exposures (so
`reinterp_spectra` uses the sky solution), a warning naming the night, the
guard and its value is emitted, and the verdict is recorded under `fpi_qa/` in
the night's `wavecalNightAve` file.

`almanac_file` is needed to derive the night's FPI guide fibers. When it is
empty the primary guard falls back to the all-fiber statistic and says so.
"""
function fpi_medwavecal_skyline_dither(mjd, mjd_list_fpi, mjd_list_wavecal,
					     all1DfpiPeaks_a, all1DObjectSkyPeaks,
					     wavecalNightAve_fname;
                         verbose = true,
					     checkpoint_mode = "commit_same",
					     almanac_file = "")
    mskMJD_fpi = (mjd_list_fpi .== mjd)
    mskMJD_obj = (mjd_list_wavecal .== mjd)

    if isnothing(wavecalNightAve_fname)
        @warn "Could not find nightly average wave soln at $(wavecalNightAve_fname)"
    elseif (size(all1DfpiPeaks_a[mskMJD_fpi], 1) > 0)
        fpi_fnames = all1DfpiPeaks_a[mskMJD_fpi]
        tele = split(basename(fpi_fnames[1]), "_")[2]

        # ---- Guard 1 (primary): is there FPI light in the FPI fibers?
        # Runs on the fpiPeaks products, before any of the wavecal work, so a
        # dark FPI costs nothing and cannot reach the promotion below.
        fpi_fiberIDs = if isempty(almanac_file)
            Int[]
        else
            get_fpi_fiberIDs_from_almanac(almanac_file, tele, mjd)
        end
        fpi_fiberIndxs = try
            filter(i -> 1 <= i <= N_FIBERS, fiberID2fiberIndx.(fpi_fiberIDs))
        catch
            Int[]
        end
        light = fpi_light_check(fpi_fnames; fpi_fiberIndxs = fpi_fiberIndxs)

        qa = Dict{String, Any}(
            "tele" => String(tele),
            "mjd" => mjd,
            "n_fpi_exposures" => length(fpi_fnames),
            "n_fpiPeaks_files_read" => light.n_files,
            "fpi_fiberIDs" => collect(Int, fpi_fiberIDs),
            "fpi_fiberIDs_from_almanac" => !isempty(fpi_fiberIDs),
            "lit_fraction" => light.lit,
            "lit_fraction_all_fibers" => light.lit_all,
            "lit_fraction_fpi_fibers" => light.lit_fpifib,
            "lit_fraction_threshold" => FPI_LIT_FRACTION_MIN,
            "resid_rms_threshold" => FPI_RESID_RMS_MAX,
            "resid_rms" => NaN,
            "resid_mad" => NaN,
            "n_peaks_used" => -1,
            "frac_peaks_used" => NaN,
            "gate_pass" => false,
            "failed_guard" => "none",
            "best_wave_type" => "sky")

        if !light.pass
            @warn "FPI WAVECAL GATE FAILED for $(tele) MJD $(mjd): guard fpi_light. " *
                  "Lit fraction of fitted FPI peaks = $(round(light.lit, sigdigits = 3)) " *
                  "(all fibers $(round(light.lit_all, sigdigits = 3)), " *
                  "FPI guide fibers $(round(light.lit_fpifib, sigdigits = 3)), " *
                  "IDs $(isempty(fpi_fiberIDs) ? "unknown - almanac named none" : join(fpi_fiberIDs, "/"))), " *
                  "threshold $(FPI_LIT_FRACTION_MIN). The FPI produced no usable light; " *
                  "FALLING BACK to the nightly sky-line wavelength solution for this night. " *
                  "The night is NOT rejected."
            qa["failed_guard"] = "fpi_light"
            save_fpi_qa!(wavecalNightAve_fname, qa)
            return
        end

        f = h5open(wavecalNightAve_fname, "r")
        curr_best_wave_type = read(f["best_wave_type"])
	curr_linParams = read(f["$(curr_best_wave_type)/nightAve_linParams"])
	curr_nlParams = read(f["$(curr_best_wave_type)/nightAve_nlParams"])
	close(f)

        # using FPI exposures to measure high-precision nightly wavelength solution
        fpi_out = comb_exp_get_and_save_fpi_wavecal(
            fpi_fnames, curr_linParams, curr_nlParams, cporder = 1, wporder = 4, dporder = 2, n_sigma = 4, max_ang_sigma = 0.2, max_iter = 2, verbose = verbose, checkpoint_mode = checkpoint_mode)

        if isnothing(fpi_out)
            @warn "FPI WAVECAL GATE FAILED for $(tele) MJD $(mjd): guard fpi_solution. " *
                  "comb_exp_get_and_save_fpi_wavecal returned nothing. " *
                  "FALLING BACK to the nightly sky-line wavelength solution. The night is NOT rejected."
            qa["failed_guard"] = "fpi_solution"
            save_fpi_qa!(wavecalNightAve_fname, qa)
            return
        end
        outfname, night_linParams, night_nlParams, night_wave_soln = fpi_out

        # ---- Guard 2 (secondary): FPI fit residual RMS.
        resid = fpi_resid_stats(outfname)
        qa["resid_rms"] = resid.rms
        qa["resid_mad"] = resid.mad
        qa["n_peaks_used"] = resid.n_used
        qa["frac_peaks_used"] = resid.frac_used
        if resid.n_used == 0
            @warn "FPI WAVECAL GATE FAILED for $(tele) MJD $(mjd): guard fpi_no_peaks_used. " *
                  "The FPI wavelength fit used ZERO peaks, so $(basename(outfname)) holds " *
                  "the seed solution rather than a fit. " *
                  "FALLING BACK to the nightly sky-line wavelength solution for this night. " *
                  "The night is NOT rejected."
            qa["failed_guard"] = "fpi_no_peaks_used"
            save_fpi_qa!(wavecalNightAve_fname, qa)
            return
        elseif !resid.pass
            @warn "FPI WAVECAL GATE FAILED for $(tele) MJD $(mjd): guard fpi_resid_rms. " *
                  "FPI fit residual RMS = $(round(resid.rms, sigdigits = 3)) Angstroms " *
                  "(threshold $(FPI_RESID_RMS_MAX); $(round(100 * resid.frac_used, sigdigits = 3))% of peaks used). " *
                  "FALLING BACK to the nightly sky-line wavelength solution for this night. " *
                  "The night is NOT rejected."
            qa["failed_guard"] = "fpi_resid_rms"
            save_fpi_qa!(wavecalNightAve_fname, qa)
            return
        end

        qa["gate_pass"] = true
        qa["best_wave_type"] = "fpi"

        f = h5open(wavecalNightAve_fname, "r+")
	if haskey(f, "best_wave_type")
            delete_object(f, "best_wave_type")
        end

	curr_best_wave_type = "fpi"
	f["best_wave_type"] = curr_best_wave_type
	f["$(curr_best_wave_type)/used_fnames"] = fpi_fnames
	f["$(curr_best_wave_type)/nightAve_wave_soln"] = night_wave_soln
	f["$(curr_best_wave_type)/nightAve_nlParams"] = night_nlParams
	f["$(curr_best_wave_type)/nightAve_linParams"] = night_linParams
	close(f)

        save_fpi_qa!(wavecalNightAve_fname, qa)

        # using skylines to measure dither offsets from FPI-defined wavelength solution
        get_and_save_sky_dither_per_fiber_partial(fname) = get_and_save_sky_dither_per_fiber(
            fname, night_linParams, night_nlParams; dporder = 2,
            wavetype = "fpi", max_offset = 1.0)

        map(get_and_save_sky_dither_per_fiber_partial, all1DObjectSkyPeaks[mskMJD_obj])
    end
end

function sky_wave_plots(
        fname_list, night_linParams, night_nlParams, night_wave_soln;
        dirNamePlots = "../outdir/plots/", plot_fibers = (1, 50, 100, 150, 200, 250, 300),
        plot_pixels = (1, 512, 1024, 1536, 2048), wavetype = "sky")
    #guard against repetition
    for chip in CHIP_LIST
        fname_list .= replace.(fname_list, "_$(chip)_" => "_")
    end
    unique_fname_list = unique(fname_list)

    sname = split(split(unique_fname_list[1], "/")[end], "_")
    tele, mjd = sname[(end - 3):(end - 2)]

    n_fnames = size(unique_fname_list, 1)
    n_lin_coeffs = size(night_linParams, 2)
    n_nl_coeffs = size(night_nlParams, 2)
    fiber_inds = 1:N_FIBERS

    clims = (1, max(2,n_fnames))

    fig = Figure(size = (1200, 400 * n_lin_coeffs), fontsize = 22)
    axis_dict = Dict()

    for pind in 1:n_lin_coeffs
        if pind == 1
            axis_dict["$(pind)"] = Axis(fig[pind, 1],
                ylabel = "Linear Parameter $(pind)",
                title = "Tele: $(tele), MJD: $(mjd), Wave. Soln. Linear Params")
        elseif pind == n_lin_coeffs
            axis_dict["$(pind)"] = Axis(fig[pind, 1],
                xlabel = "FIBERID",
                ylabel = "Linear Parameter $(pind)")
        else
            axis_dict["$(pind)"] = Axis(fig[pind, 1],
                ylabel = "Linear Parameter $(pind)")
        end
    end

    for fname_ind in 1:n_fnames
        f = h5open(unique_fname_list[fname_ind], "r+")

        linParams = try
            read(f["linParams"])
        catch
            println(unique_fname_list[fname_ind])
            read(f["linParams"])
        end

        close(f)

        for pind in 1:n_lin_coeffs
            scatter!(axis_dict["$(pind)"], 301 .- fiber_inds, linParams[:, pind],
                color = fname_ind * ones(size(fiber_inds, 1)), colorrange = clims)
        end
    end

    sky_wave_linParams_Path = dirNamePlots *
                              "$(mjd)/$(wavetype)wave_linParams_$(tele)_$(mjd).png"
    save(sky_wave_linParams_Path, fig)

    fig = Figure(size = (1200, 400 * n_nl_coeffs), fontsize = 22)
    axis_dict = Dict()

    for pind in 1:n_nl_coeffs
        if pind == 1
            axis_dict["$(pind)"] = Axis(fig[pind, 1],
                ylabel = "Non-Linear Parameter $(pind)",
                title = "Tele: $(tele), MJD: $(mjd), Wave. Soln. Non-Linear Params")
        elseif pind == n_nl_coeffs
            axis_dict["$(pind)"] = Axis(fig[pind, 1],
                xlabel = "FIBERID",
                ylabel = "Non-Linear Parameter $(pind)")
        else
            axis_dict["$(pind)"] = Axis(fig[pind, 1],
                ylabel = "Non-Linear Parameter $(pind)")
        end
    end

    for fname_ind in 1:n_fnames
        f = h5open(unique_fname_list[fname_ind], "r+")

        nlParams = try
            read(f["nlParams"])
        catch
            println(unique_fname_list[fname_ind])
            read(f["nlParams"])
        end

        close(f)

        for pind in 1:n_nl_coeffs
            scatter!(axis_dict["$(pind)"], 301 .- fiber_inds, nlParams[:, pind],
                color = fname_ind * ones(size(fiber_inds, 1)), colorrange = clims)
        end
    end

    sky_wave_nlParams_Path = dirNamePlots *
                             "$(mjd)/$(wavetype)wave_nlParams_$(tele)_$(mjd).png"
    save(sky_wave_nlParams_Path, fig)

    n_plot_fibers = size(plot_fibers, 1)

    fig = Figure(size = (800 * 3, 400 * n_plot_fibers), fontsize = 22)
    axis_dict = Dict()

    for pind in 1:n_plot_fibers
        fibid = 301 - plot_fibers[pind]
        chip_a_diff = nanmedian(diff(night_wave_soln[:, plot_fibers[pind], 1]))
        chip_b_diff = nanmedian(diff(night_wave_soln[:, plot_fibers[pind], 2]))
        chip_c_diff = nanmedian(diff(night_wave_soln[:, plot_fibers[pind], 3]))
        if pind == 1
            axis_dict["b $(pind)"] = Axis(fig[pind, 2],
                title = "Tele: $(tele), MJD: $(mjd)\nChip $(CHIP_LIST[2]), Δλ = $(round(chip_b_diff,digits=3)) Å/pixel")
            axis_dict["a $(pind)"] = Axis(fig[pind, 1],
                ylabel = "FID$(fibid) Wavelength Offset (Å)",
                title = "\nChip $(CHIP_LIST[1]), Δλ = $(round(chip_a_diff,digits=3)) Å/pixel")
            axis_dict["c $(pind)"] = Axis(fig[pind, 3],
                title = "\nChip $(CHIP_LIST[3]), Δλ = $(round(chip_c_diff,digits=3)) Å/pixel")
        elseif pind == n_plot_fibers
            axis_dict["b $(pind)"] = Axis(fig[pind, 2],
                xlabel = "X Pixel")
            axis_dict["a $(pind)"] = Axis(fig[pind, 1],
                xlabel = "X Pixel",
                ylabel = "FID$(fibid) Wavelength Offset (Å)")
            axis_dict["c $(pind)"] = Axis(fig[pind, 3],
                xlabel = "X Pixel")
        else
            axis_dict["b $(pind)"] = Axis(fig[pind, 2])
            axis_dict["a $(pind)"] = Axis(fig[pind, 1],
                ylabel = "FID$(fibid) Wavelength Offset (Å)")
            axis_dict["c $(pind)"] = Axis(fig[pind, 3])
        end
    end

    dlam_dx = zeros(Float64, size(night_wave_soln))
    dlam_dx[(begin + 1):end, :, :] .= diff(night_wave_soln, dims = 1)
    dlam_dx[begin, :, :] .= dlam_dx[begin + 1, :, :]
    for fname_ind in 1:n_fnames
        f = h5open(unique_fname_list[fname_ind], "r+")

        chipWaveSoln = try
            read(f["chipWaveSoln"])
        catch
            println(unique_fname_list[fname_ind])
            read(f["chipWaveSoln"])
        end

        close(f)

        dx = nanmedian(nanmedian((night_wave_soln .- chipWaveSoln) ./ dlam_dx, 1), 3)[1, :, 1]

        for pind in 1:n_plot_fibers
            scatter!(axis_dict["a $(pind)"],
                1:N_XPIX,
                chipWaveSoln[:, plot_fibers[pind], 1] .- night_wave_soln[:, plot_fibers[pind], 1] .+
                dx[plot_fibers[pind]] .* dlam_dx[:, plot_fibers[pind], 1],
                color = fname_ind * ones(N_XPIX),
                colorrange = clims)
            scatter!(axis_dict["b $(pind)"],
                1:N_XPIX,
                chipWaveSoln[:, plot_fibers[pind], 2] .- night_wave_soln[:, plot_fibers[pind], 2] .+
                dx[plot_fibers[pind]] .* dlam_dx[:, plot_fibers[pind], 2],
                color = fname_ind * ones(N_XPIX),
                colorrange = clims)
            scatter!(axis_dict["c $(pind)"],
                1:N_XPIX,
                chipWaveSoln[:, plot_fibers[pind], 3] .- night_wave_soln[:, plot_fibers[pind], 3] .+
                dx[plot_fibers[pind]] .* dlam_dx[:, plot_fibers[pind], 3],
                color = fname_ind * ones(N_XPIX),
                colorrange = clims)
        end
    end

    sky_wave_fiber_vs_pixel_Path = dirNamePlots *
                                   "$(mjd)/$(wavetype)wave_per_fiber_vs_pixel_$(tele)_$(mjd).png"
    save(sky_wave_fiber_vs_pixel_Path, fig)

    n_plot_pixels = size(plot_pixels, 1)

    fig = Figure(size = (800 * 3, 400 * n_plot_pixels), fontsize = 22)
    axis_dict = Dict()

    chip_a_diff = nanmedian(diff(night_wave_soln[:, 150, 1]))
    chip_b_diff = nanmedian(diff(night_wave_soln[:, 150, 2]))
    chip_c_diff = nanmedian(diff(night_wave_soln[:, 150, 3]))
    for pind in 1:n_plot_pixels
        if pind == 1
            axis_dict["b $(pind)"] = Axis(fig[pind, 2],
                title = "Tele: $(tele), MJD: $(mjd)\nChip $(CHIP_LIST[2]), Δλ = $(round(chip_b_diff,digits=3)) Å/pixel")
            axis_dict["a $(pind)"] = Axis(fig[pind, 1],
                ylabel = "XPix$(plot_pixels[pind]) Wavelength Offset (Å)",
                title = "\nChip $(CHIP_LIST[1]), Δλ = $(round(chip_a_diff,digits=3)) Å/pixel")
            axis_dict["c $(pind)"] = Axis(fig[pind, 3],
                title = "\nChip $(CHIP_LIST[3]), Δλ = $(round(chip_c_diff,digits=3)) Å/pixel")
        elseif pind == n_plot_fibers
            axis_dict["b $(pind)"] = Axis(fig[pind, 2],
                xlabel = "FIBER ID")
            axis_dict["a $(pind)"] = Axis(fig[pind, 1],
                xlabel = "FIBER ID",
                ylabel = "XPix$(plot_pixels[pind]) Wavelength Offset (Å)")
            axis_dict["c $(pind)"] = Axis(fig[pind, 3],
                xlabel = "FIBER ID")
        else
            axis_dict["b $(pind)"] = Axis(fig[pind, 2])
            axis_dict["a $(pind)"] = Axis(fig[pind, 1],
                ylabel = "XPix$(plot_pixels[pind]) Wavelength Offset (Å)")
            axis_dict["c $(pind)"] = Axis(fig[pind, 3])
        end
    end

    for fname_ind in 1:n_fnames
        f = h5open(unique_fname_list[fname_ind], "r+")

        chipWaveSoln = try
            read(f["chipWaveSoln"])
        catch
            println(unique_fname_list[fname_ind])
            read(f["chipWaveSoln"])
        end

        close(f)

        dx = nanmedian(nanmedian((night_wave_soln .- chipWaveSoln) ./ dlam_dx, 1), 3)[1, :, 1]

        for pind in 1:n_plot_pixels
            scatter!(axis_dict["a $(pind)"],
                301 .- fiber_inds,
                chipWaveSoln[plot_pixels[pind], :, 1] .- night_wave_soln[plot_pixels[pind], :, 1] .+
                dx .* dlam_dx[plot_pixels[pind], :, 1],
                color = fname_ind * ones(N_FIBERS),
                colorrange = clims)
            scatter!(axis_dict["b $(pind)"],
                301 .- fiber_inds,
                chipWaveSoln[plot_pixels[pind], :, 2] .- night_wave_soln[plot_pixels[pind], :, 2] .+
                dx .* dlam_dx[plot_pixels[pind], :, 2],
                color = fname_ind * ones(N_FIBERS),
                colorrange = clims)
            scatter!(axis_dict["c $(pind)"],
                301 .- fiber_inds,
                chipWaveSoln[plot_pixels[pind], :, 3] .- night_wave_soln[plot_pixels[pind], :, 3] .+
                dx .* dlam_dx[plot_pixels[pind], :, 3],
                color = fname_ind * ones(N_FIBERS),
                colorrange = clims)
        end
    end

    sky_wave_pixel_vs_fiber_Path = dirNamePlots *
                                   "$(mjd)/$(wavetype)wave_per_pixel_vs_fiber_$(tele)_$(mjd).png"
    save(sky_wave_pixel_vs_fiber_Path, fig)

    return sky_wave_linParams_Path, sky_wave_nlParams_Path, sky_wave_fiber_vs_pixel_Path,
    sky_wave_pixel_vs_fiber_Path
end
