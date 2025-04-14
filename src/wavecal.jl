using HDF5
using Polynomials: Polynomial, fit

function linear_loss_fit(x, y; wporder = 2, returnL2only = false)
    A = positional_poly_mat(x, porder = wporder)
    linparam = A \ y
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
        return [ones(length(x)) x x .^ 2 x .^ 3 x .^ 4 x .^ 5 x .^ 6 x .^ 7 x .^ 8 x .^ 9 x .^ 10 x .^ 11]
    elseif porder == 12
        return [ones(length(x)) x x .^ 2 x .^ 3 x .^ 4 x .^ 5 x .^ 6 x .^ 7 x .^ 8 x .^ 9 x .^ 10 x .^ 11 x .^ 12]
    elseif porder == 13
        return [ones(length(x)) x x .^ 2 x .^ 3 x .^ 4 x .^ 5 x .^ 6 x .^ 7 x .^ 8 x .^ 9 x .^ 10 x .^ 11 x .^ 12 x .^ 13]
    elseif porder == 14
        return [ones(length(x)) x x .^ 2 x .^ 3 x .^ 4 x .^ 5 x .^ 6 x .^ 7 x .^ 8 x .^ 9 x .^ 10 x .^ 11 x .^ 12 x .^ 13 x .^ 14]
    elseif porder == 15
        return [ones(length(x)) x x .^ 2 x .^ 3 x .^ 4 x .^ 5 x .^ 6 x .^ 7 x .^ 8 x .^ 9 x .^ 10 x .^ 11 x .^ 12 x .^ 13 x .^ 14 x .^ 15]
    elseif porder == 16
        return [ones(length(x)) x x .^ 2 x .^ 3 x .^ 4 x .^ 5 x .^ 6 x .^ 7 x .^ 8 x .^ 9 x .^ 10 x .^ 11 x .^ 12 x .^ 13 x .^ 14 x .^ 15 x .^ 16]
    elseif porder == 17
        return [ones(length(x)) x x .^ 2 x .^ 3 x .^ 4 x .^ 5 x .^ 6 x .^ 7 x .^ 8 x .^ 9 x .^ 10 x .^ 11 x .^ 12 x .^ 13 x .^ 14 x .^ 15 x .^ 16 x .^ 17]
    elseif porder == 18
        return [ones(length(x)) x x .^ 2 x .^ 3 x .^ 4 x .^ 5 x .^ 6 x .^ 7 x .^ 8 x .^ 9 x .^ 10 x .^ 11 x .^ 12 x .^ 13 x .^ 14 x .^ 15 x .^ 16 x .^ 17 x .^ 18]
    elseif porder == 19
        return [ones(length(x)) x x .^ 2 x .^ 3 x .^ 4 x .^ 5 x .^ 6 x .^ 7 x .^ 8 x .^ 9 x .^ 10 x .^ 11 x .^ 12 x .^ 13 x .^ 14 x .^ 15 x .^ 16 x .^ 17 x .^ 18 x .^ 19]
    elseif porder == 20
        return [ones(length(x)) x x .^ 2 x .^ 3 x .^ 4 x .^ 5 x .^ 6 x .^ 7 x .^ 8 x .^ 9 x .^ 10 x .^ 11 x .^ 12 x .^ 13 x .^ 14 x .^ 15 x .^ 16 x .^ 17 x .^ 18 x .^ 19 x .^ 20]
    end
end

# make sure to pass this a copy so you can modify in place
function nonlinear_loss_fit!(
        chipPolyParams, inparams, x, y, chipInt; wporder = 2, cporder = 1, returnL2only = false)
    params2ChipPolyParams!(chipPolyParams, inparams, cporder)
    xt = zeros(Float64, length(x))
    for i in 1:3
        msk = chipInt .== i
        xt[msk] .= transform_x_chips(x[msk], chipPolyParams[i, :])
    end
    return linear_loss_fit(xt, y, wporder = wporder, returnL2only = returnL2only)
end

function transform_x_chips(x, chipPolyParams)
    porder = length(chipPolyParams) - 1
    if porder == 1
        return chipPolyParams[1] .+ chipPolyParams[2] .* x
    elseif porder == 2
        return chipPolyParams[1] .+ chipPolyParams[2] .* x .+ chipPolyParams[3] .* x .^ 2
    elseif porder == 3
        return chipPolyParams[1] .+ chipPolyParams[2] .* x .+ chipPolyParams[3] .* x .^ 2 .+
               chipPolyParams[4] .* x .^ 3
    end
end

function params2ChipPolyParams!(chipPolyParams, inparams, cporder)
    chipPolyParams[1, :] .= inparams[1:(cporder + 1), :]
    chipPolyParams[3, :] .= inparams[(cporder + 2):(2 * (cporder + 1)), :]
    return nothing
end

function ChipPolyParams2Params(chipPolyParams)
    return vcat(chipPolyParams[1, :], chipPolyParams[3, :])
end

# Sky line wavecal
function get_and_save_sky_wavecal(fname; cporder = 1, wporder = 2)
    # initial guess for the (low-order)chip polynomial parameters
    if occursin("_apo_",fname)
        # chipPolyParams0 = [-1.0716 1.00111
        #                    0 1
        #                    1.07009 0.98803]
        offset_func_chip1 = Polynomial([-1.07221500e+00,  4.52450367e-06])
        scale_func_chip1 = Polynomial([1.00093875e+00, -4.41886670e-07])
        offset_func_chip3 = Polynomial([1.06972294e+00, 2.77444782e-06])
        scale_func_chip3 = Polynomial([9.87857338e-01, 1.09350510e-06])
    elseif occursin("_lco_",fname)
        # chipPolyParams0 = [-1.0748 1.00168
        #                    0 1
        #                    1.07089 0.98763]
        offset_func_chip1 = Polynomial([-1.07456222e+00, -1.52100076e-07])
        scale_func_chip1 = Polynomial([1.00123795e+00, 5.30751281e-07])
        offset_func_chip3 = Polynomial([1.07199520e+00, -7.11920517e-06])
        scale_func_chip3 = Polynomial([9.87968936e-01, -2.76150881e-07])
    else
        # chipPolyParams0 = [-1.070 1
        #                    0 1
        #                    1.076 1]
        offset_func_chip1 = Polynomial([-1.070, 0.0])
        scale_func_chip1 = Polynomial([1.0, 0.0])
        offset_func_chip3 = Polynomial([1.076, 0.0])
        scale_func_chip3 = Polynomial([1.0, 0.0])
    end

    fibInds = 1:300
    chipPolyParams0 = zeros(Float64, (size(fibInds,1), 3, cporder+1))
    chipPolyParams0[:, 1, 1] .= offset_func_chip1.(fibInds)
    chipPolyParams0[:, 1, 2] .= scale_func_chip1.(fibInds)
    chipPolyParams0[:, 2, 1] .= 0.0
    chipPolyParams0[:, 2, 2] .= 1.0
    chipPolyParams0[:, 3, 1] .= offset_func_chip3.(fibInds)
    chipPolyParams0[:, 3, 2] .= scale_func_chip3.(fibInds)

    outname = replace(replace(fname, "skyLine_peaks" => "wavecal_skyline"), "_a_" => "_")
    sky_line_uxlst, sky_line_fwlst, sky_line_chipInt = ingest_skyLines_exp(fname)
    linParams, nlParams,
    resid_vec = get_sky_wavecal(
        sky_line_uxlst, sky_line_fwlst, sky_line_chipInt,
        chipPolyParams0; cporder = cporder, wporder = wporder)

    # iterpolate between fibers to
    # regularize the nlParams and linParams
    interp_nlParams, interp_linParams = interpolate_wave_params(fibInds,nlParams,linParams;linParam_deg=2,nlParam_deg=1)

    # do final pass with interp params to
    # get a new constant-term for wave soln
    # (which we do not think should be interpolated)

    interp_resid_vec = zeros(Float64, (size(fibInds,1), size(sky_line_uxlst, 1)))
    fill!(interp_resid_vec, NaN)
    for fibIndx in fibInds
        xv = sky_line_uxlst[:, fibIndx]
        yv = sky_line_fwlst[:, fibIndx]
        chipIntv = sky_line_chipInt[:, fibIndx]
        sky_msk = .!isnan.(xv)
        chipPolyParams = copy(chipPolyParams0[fibIndx, :, :])
        inparams = interp_nlParams[fibIndx, :]

        params2ChipPolyParams!(chipPolyParams, inparams, cporder)
        xt = zeros(Float64, length(xv[sky_msk]))
        for i in 1:3
            msk = chipIntv[sky_msk] .== i
            xt[msk] .= transform_x_chips(xv[sky_msk][msk], chipPolyParams[i, :])
        end

        A = positional_poly_mat(xt, porder = wporder)
	curr_resids = yv[sky_msk] .- A * interp_linParams[fibIndx, :]
	const_offset = mean(curr_resids) #could also be median
	interp_linParams[fibIndx, 1] += const_offset
	interp_resid_vec[fibIndx, sky_msk] .= curr_resids .- const_offset
    end

    chipWaveSoln = zeros(Float64, 2048, 300, 3)
    interp_chipWaveSoln = zeros(Float64, 2048, 300, 3)
    x = 1:2048
    ximport = (x .- 1024) ./ 2048
    for chip in ["a", "b", "c"]
        chipIndx = getChipIndx(chip)
        for fibIndx in 1:300
            params2ChipPolyParams!(chipPolyParams0[fibIndx, :, :], nlParams[fibIndx, :], cporder)
            xt = transform_x_chips(ximport, chipPolyParams0[fibIndx, chipIndx, :])
            Ax = positional_poly_mat(xt, porder = 2)
            yt = Ax * linParams[fibIndx, :]
            chipWaveSoln[:, fibIndx, chipIndx] .= yt

            params2ChipPolyParams!(chipPolyParams0[fibIndx, :, :], interp_nlParams[fibIndx, :], cporder)
            xt .= transform_x_chips(ximport, chipPolyParams0[fibIndx, chipIndx, :])
            Ax .= positional_poly_mat(xt, porder = 2)
            yt .= Ax * interp_linParams[fibIndx, :]
            interp_chipWaveSoln[:, fibIndx, chipIndx] .= yt
        end
    end

    # the best parameters are likely the interpolated ones
    # but return the raw measurements as well
    safe_jldsave(outname; linParams = interp_linParams, nlParams = interp_nlParams,
        resid_vec = interp_resid_vec, chipWaveSoln = interp_chipWaveSoln, 
	raw_linParams = linParams, raw_nlParams = nlParams,
        raw_resid_vec = resid_vec, raw_chipWaveSoln = chipWaveSoln)
end

function get_sky_wavecal(
        sky_line_uxlst, sky_line_fwlst, sky_line_chipInt, chipPolyParams0; cporder = 1, wporder = 2)
    linParams = zeros(Float64, 300, wporder + 1)
    nlParams = zeros(Float64, 300, 2 * (cporder + 1))
    resid_vec = zeros(Float64, 300, size(sky_line_uxlst, 1))
    fill!(resid_vec, NaN)
    for i in 1:300
        xv = sky_line_uxlst[:, i]
        yv = sky_line_fwlst[:, i]
        chipIntv = sky_line_chipInt[:, i]
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
    return linParams, nlParams, resid_vec
end

function fit_poly_without_outliers(xvals,yvals,deg;nsigma=5,max_repeat=3)
    summary = nanzeropercentile(yvals, [16,50,84])
    #change to median,-sigma,+sigma
    summary = [summary[2],summary[2]-summary[1],summary[3]-summary[2]]

    keep = (yvals .>= summary[1]-nsigma*summary[2]) .&
           (yvals .<= summary[1]+nsigma*summary[3])

    for repeatInd in 1:max_repeat
        func = fit(xvals[keep],yvals[keep],deg)
	curr_resids = yvals .- func.(xvals)

        resid_summary = nanzeropercentile(curr_resids[keep], [16,50,84])
        resid_summary = [resid_summary[2],
			  resid_summary[2]-resid_summary[1],
			  resid_summary[3]-resid_summary[2]]

        new_good = (curr_resids .>= resid_summary[1]-nsigma*resid_summary[2]) .&
                   (curr_resids .<= resid_summary[1]+nsigma*resid_summary[3])

	if all(new_good .== keep)
	    break
	end
	keep .= new_good
    end
    return func, keep
end

function interpolate_wave_params(fibInds,nlParams,linParams;linParam_deg=2,nlParam_deg=1)
    mtp_inds = floor.((fibInds .- 1) ./ 30) .+ 1
    unique_mtp = 1:10

    interp_linParams = zeros(Float64,size(linParams))
    interp_nlParams = zeros(Float64,size(nlParams))

    #constant term in wave soln needs to fit each mtp separately
    for mtp_ind in unique_mtp
        keep_mtp = (mtp_inds .== mtp_ind)
        func,keep = fit_poly_without_outliers(fibInds,linParams[1,keep_mtp],deg=linParam_deg)
	func_eval = func.(fibInds[keep_mtp])
	interp_linParams[1,keep_mtp] .= func_eval .+ nanmedian((linParams[1,keep_mtp].-func_eval)[keep])
    end

    for j in 2:size(linParams,1)
        func,keep = fit_poly_without_outliers(fibInds,linParams[j,:],deg=linParam_deg)
	func_eval = func.(fibInds)
	interp_linParams[j,keep_mtp] .= func_eval .+ nanmedian((linParams[j,:].-func_eval)[keep])
    end

    for j in 1:size(nlParams,1)
        func,keep = fit_poly_without_outliers(fibInds,nlParams[j,:],deg=nlParam_deg)
	func_eval = func.(fibInds)
	interp_nlParams[j,keep_mtp] .= func_eval .+ nanmedian((nlParams[j,:].-func_eval)[keep])
    end

    return interp_nlParams, interp_linParams
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
    close(f)
    sky_line_xlst = (sky_line_mat_clean[:, 1, :] .- 1024) ./ 2048
    sky_line_wlst = sky_line_mat_clean[:, 2, :]
    return sky_line_xlst, sky_line_wlst
end

# this takes in a filename and replaces the chip index (make "a" default approx)
function ingest_skyLines_exp(fname)
    sky_line_uxlst = Matrix{Float64}[]
    sky_line_fwlst = Matrix{Float64}[]
    sky_line_chipInt = Matrix{Int}[]
    for chip in ["a", "b", "c"]
        fnameloc = replace(fname, "_a_" => "_$(chip)_")
        if isfile(fnameloc)
            sky_line_xlst, sky_line_wlst = ingest_skyLines_file(fnameloc)
            chipIndx = getChipIndx(chip)
            push!(sky_line_uxlst, sky_line_xlst)
            push!(sky_line_fwlst, sky_line_wlst)
            push!(sky_line_chipInt, chipIndx * ones(Int, size(sky_line_wlst)))
        else
            push!(sky_line_uxlst, [])
            push!(sky_line_fwlst, [])
            push!(sky_line_chipInt, [])
        end
    end
    sky_line_uxlst = vcat(sky_line_uxlst...)
    sky_line_fwlst = vcat(sky_line_fwlst...)
    sky_line_chipInt = vcat(sky_line_chipInt...)

    msk_large_scatter = dropdims(nanzeroiqr(sky_line_uxlst, 2) .> 0.002, dims = 2)
    sky_line_uxlst[msk_large_scatter, :] .= NaN
    sky_line_fwlst[msk_large_scatter, :] .= NaN
    # dims are num_sky_lines x num_fibers
    return sky_line_uxlst, sky_line_fwlst, sky_line_chipInt
end
