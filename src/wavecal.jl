using HDF5
using Polynomials: Polynomial, fit

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

# make sure to pass this a copy so you can modify in place
function nonlinear_loss_fit!(
        chipPolyParams, inparams, x, y, chipInt; wporder = 2, cporder = 1,
        returnL2only = false, linparam = nothing)
    params2ChipPolyParams!(chipPolyParams, inparams, cporder)
    xt = zeros(Float64, length(x))
    for i in 1:3
        msk = chipInt .== i
        xt[msk] .= transform_x_chips(x[msk], chipPolyParams[i, :])
    end
    return linear_loss_fit(
        xt, y, wporder = wporder, returnL2only = returnL2only, linparam = linparam)
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
function get_and_save_sky_wavecal(fname; cporder = 1, wporder = 2)
    # initial guess for the (low-order)chip polynomial parameters
    if occursin("_apo_", fname)
        # chipPolyParams0 = [-1.0716 1.00111
        #                    0 1
        #                    1.07009 0.98803]
        offset_func_chip1 = Polynomial([-1.07221500e+00, 4.52450367e-06])
        scale_func_chip1 = Polynomial([1.00093875e+00, -4.41886670e-07])
        offset_func_chip3 = Polynomial([1.06972294e+00, 2.77444782e-06])
        scale_func_chip3 = Polynomial([9.87857338e-01, 1.09350510e-06])
    elseif occursin("_lco_", fname)
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

    #don't use the polynomials above for now...
    offset_func_chip1 = Polynomial([-1.070, 0.0])
    scale_func_chip1 = Polynomial([1.0, 0.0])
    offset_func_chip3 = Polynomial([1.076, 0.0])
    scale_func_chip3 = Polynomial([1.0, 0.0])

    fibInds = 1:N_FIBERS
    chipPolyParams0 = zeros(Float64, (size(fibInds, 1), 3, cporder + 1))
    chipPolyParams0[:, 1, 1] .= offset_func_chip1.(fibInds)
    chipPolyParams0[:, 1, 2] .= scale_func_chip1.(fibInds)
    chipPolyParams0[:, 2, 1] .= 0.0
    chipPolyParams0[:, 2, 2] .= 1.0
    chipPolyParams0[:, 3, 1] .= offset_func_chip3.(fibInds)
    chipPolyParams0[:, 3, 2] .= scale_func_chip3.(fibInds)

    outname = replace(replace(fname, "skyLinePeaks" => "waveCalSkyLine"), "_a_" => "_")
    sky_line_uxlst, sky_line_fwlst, sky_line_chipInt = ingest_skyLines_exp(fname)
    linParams, nlParams,
    resid_vec = get_sky_wavecal(
        sky_line_uxlst, sky_line_fwlst, sky_line_chipInt,
        chipPolyParams0; cporder = cporder, wporder = wporder)

    # iterpolate between fibers to
    # regularize the nlParams and linParams
    interp_nlParams, interp_linParams = interpolate_wave_params(
        fibInds, nlParams, linParams)

    # do final pass with interp params to
    # get a new constant-term for wave soln
    # (which we do not think should be interpolated)

    interp_resid_vec = zeros(Float64, (size(fibInds, 1), size(sky_line_uxlst, 1)))
    fill!(interp_resid_vec, NaN)
    for fibIndx in fibInds
        xv = sky_line_uxlst[:, fibIndx]
        yv = sky_line_fwlst[:, fibIndx]
        chipIntv = sky_line_chipInt[:, fibIndx]
        sky_msk = .!isnan.(xv)
        if sum(sky_msk) == 0
            continue
        end

        chipPolyParams = copy(chipPolyParams0[fibIndx, :, :])
        inparams = interp_nlParams[fibIndx, :]

        params2ChipPolyParams!(chipPolyParams, inparams, cporder)
        xt = zeros(Float64, length(xv[sky_msk]))
        for i in 1:N_CHIPS
            msk = chipIntv[sky_msk] .== i
            xt[msk] .= transform_x_chips(xv[sky_msk][msk], chipPolyParams[i, :])
        end

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

    chipWaveSoln = zeros(Float64, N_XPIX, N_FIBERS, 3)
    interp_chipWaveSoln = zeros(Float64, N_XPIX, N_FIBERS, 3)
    x = 1:N_XPIX
    ximport = (x .- (N_XPIX ÷ 2)) ./ N_XPIX
    for chip in ["a", "b", "c"]
        chipIndx = getChipIndx(chip)
        for fibIndx in 1:N_FIBERS
            params2ChipPolyParams!(chipPolyParams0, nlParams, cporder, fibIndx = fibIndx)
            xt = transform_x_chips(ximport, chipPolyParams0[fibIndx, chipIndx, :])
            Ax = positional_poly_mat(xt, porder = 2)
            yt = Ax * linParams[fibIndx, :]
            chipWaveSoln[:, fibIndx, chipIndx] .= yt

            params2ChipPolyParams!(chipPolyParams0, interp_nlParams, cporder, fibIndx = fibIndx)
            xt = transform_x_chips(ximport, chipPolyParams0[fibIndx, chipIndx, :])
            Ax = positional_poly_mat(xt, porder = 2)
            yt = Ax * interp_linParams[fibIndx, :]
            interp_chipWaveSoln[:, fibIndx, chipIndx] .= yt
        end
    end

    # the best parameters are likely the interpolated ones
    # but return the raw measurements as well
    safe_jldsave(outname; linParams = interp_linParams, nlParams = interp_nlParams,
        resid_vec = interp_resid_vec, chipWaveSoln = interp_chipWaveSoln,
        raw_linParams = linParams, raw_nlParams = nlParams,
        raw_resid_vec = resid_vec, raw_chipWaveSoln = chipWaveSoln)
    return outname
end

#read in a list of fnames, then take the 
#average/median of the wavelength solution
#parameters (removing dither differences)
#to have a stable average for the night
function get_ave_night_wave_soln(fname_list)
    #open first to get order of linear and non-linear parameters

    for chip in ["a", "b", "c"]
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
    n_fibers = size(linParams, 1)

    all_linParams = zeros(Float64, (n_fibers, n_lin_coeffs, n_fnames))
    all_nlParams = zeros(Float64, (n_fibers, n_nl_coeffs, n_fnames))
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

    med_linParams = nanmedian((all_linParams .- linParam_offsets), 3)[:, :, 1]
    med_nlParams = nanmedian((all_nlParams .- nlParam_offsets), 3)[:, :, 1]

    wporder = n_lin_coeffs - 1
    cporder = Int(n_nl_coeffs / 2) - 1
    chipPolyParams0 = zeros(Float64, (n_fibers, 3, cporder + 1))
    chipPolyParams0[:, 2, 2] .= 1.0 #chip b scale
    chipWaveSoln = zeros(Float64, 2048, n_fibers, 3)
    x = 1:2048
    ximport = (x .- 1024) ./ 2048
    for chip in ["a", "b", "c"]
        chipIndx = getChipIndx(chip)
        for fibIndx in 1:n_fibers
            params2ChipPolyParams!(chipPolyParams0, med_nlParams, cporder, fibIndx = fibIndx)
            xt = transform_x_chips(ximport, chipPolyParams0[fibIndx, chipIndx, :])
            Ax = positional_poly_mat(xt, porder = wporder)
            yt = Ax * med_linParams[fibIndx, :]
            chipWaveSoln[:, fibIndx, chipIndx] .= yt
        end
    end

    return med_linParams, med_nlParams, chipWaveSoln
end

function get_sky_wavecal(
        sky_line_uxlst, sky_line_fwlst, sky_line_chipInt, chipPolyParams0; cporder = 1, wporder = 2)
    linParams = zeros(Float64, N_FIBERS, wporder + 1)
    nlParams = zeros(Float64, N_FIBERS, 2 * (cporder + 1))
    resid_vec = zeros(Float64, N_FIBERS, size(sky_line_uxlst, 1))
    fill!(resid_vec, NaN)
    for i in 1:N_FIBERS
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

# FPI line wavecal
function get_and_save_fpi_wavecal(
        fname, initial_linParams, initial_nlParams; cporder = 1, wporder = 4)

    #first guess FPI cavity size, Angstroms
    if occursin("_lco_", fname)
        cavity_size = 3.73610e7
    elseif occursin("_apo_", fname)
        cavity_size = 3.736125e7
    else
        println("ERROR: No initial cavity size determined for filename $(fname)")
        return nothing
    end

    initial_n_lin_coeffs = size(initial_linParams, 2)
    initial_n_nl_coeffs = size(initial_nlParams, 2)
    n_fibers = size(initial_linParams, 1)

    initial_wporder = initial_n_lin_coeffs - 1
    initial_cporder = Int(initial_n_nl_coeffs / 2) - 1

    outname = replace(replace(fname, "fpi_peaks" => "wavecal_fpi"), "_a_" => "_")

    #read in the peak positions
    fpi_line_uxlst, fpi_line_uxlst_errs, fpi_line_chipInt = ingest_fpiLines_exp(fname)
    n_peaks = size(fpi_line_uxlst, 1)

    #use the initial wavelength solution and the 
    #chip b peaks to define the m integers 
    #(m = 2*cavity_size/wavelength) of the 
    #peaks on all the chips and for all fibers
    #(they might be offset from each other if
    # some peaks were missed during fitting)

    #transform peak positions to wavelengths
    fiber_m0s = zeros(Float64, n_fibers)
    peak_waves = zeros(Float64, n_peaks, n_fibers)
    peak_ints = zeros(Float64, n_peaks, n_fibers) #integer m values
    chipWaveSoln = zeros(Float64, 2048, n_fibers, 3)
    fpi_ximport = (fpi_line_uxlst .- 1024) ./ 2048
    chipPolyParams0 = zeros(Float64, (n_fibers, 3, initial_cporder + 1))
    chipPolyParams0[:, 2, 2] .= 1.0 #chip b scaling
    for chip in ["a", "b", "c"]
        chipIndx = getChipIndx(chip)
        for fibIndx in 1:n_fibers
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
    med_m0 = nanmedian(peak_ints[1, :])

    #best fit non-integer offset to m0
    m_offset = nanmedian(2 * cavity_size ./ peak_waves .- peak_ints)
    #update cavity size guess
    cavity_size = nanmedian(peak_waves .* (peak_ints .+ m_offset) ./ 2.0)

    #need to fit additional parameters
    #to the chip offset/scaling and wave soln,
    #specifically, the FPI cavity size and
    #the non-integer offset to m0 integer
    #(should be smaller than 1 in size)

    linParams, nlParams, cavity_size, m_offset,
    resid_vec = get_fpi_wavecal(
        fpi_line_uxlst, fpi_line_uxlst_errs, peak_ints, fpi_line_chipInt,
        cavity_size, m_offset, initial_linParams, initial_wporder,
        initial_nlParams, initial_cporder, chipPolyParams0;
        cporder = cporder, wporder = wporder)
    slksjhsjkl

    chipWaveSoln = zeros(Float64, 2048, 300, 3)
    x = 1:2048
    ximport = (x .- 1024) ./ 2048
    for chip in ["a", "b", "c"]
        chipIndx = getChipIndx(chip)
        for fibIndx in 1:300
            params2ChipPolyParams!(chipPolyParams0, nlParams, cporder, fibIndx = fibIndx)
            xt = transform_x_chips(ximport, chipPolyParams0[fibIndx, chipIndx, :])
            Ax = positional_poly_mat(xt, porder = 2)
            yt = Ax * linParams[fibIndx, :]
            chipWaveSoln[:, fibIndx, chipIndx] .= yt
        end
    end

    # the best parameters are likely the interpolated ones
    # but return the raw measurements as well
    safe_jldsave(outname; linParams = linParams, nlParams = nlParams,
        resid_vec = resid_vec, chipWaveSoln = chipWaveSoln)
end

function get_fpi_wavecal(
        fpi_line_uxlst, fpi_line_uxlst_errs, peak_ints,
        fpi_line_chipInt, cavity_size, m_offset,
        initial_linParams, initial_wporder,
        initial_nlParams, initial_cporder,
        chipPolyParams0; cporder = 1, wporder = 4)
    linParams = zeros(Float64, 300, wporder + 1)
    nlParams = zeros(Float64, 300, 2 * (cporder + 1) + 2) #extra 2 for integer offset and cavity size
    resid_vec = zeros(Float64, 300, size(fpi_line_uxlst, 1))
    fill!(resid_vec, NaN)

    linParams[:, begin:(begin + size(initial_linParams, 2))] .= initial_linParams
    nlParams[:, begin:(begin + size(initial_nlParams, 2))] .= initial_nlParams

    #best fit non-integer offset to m0
    m_offset = nanmedian(2 * cavity_size ./ peak_waves .- peak_ints)
    #update cavity size guess
    cavity_size = nanmedian(peak_waves .* (peak_ints .+ m_offset) ./ 2.0)

    for i in 1:300
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
    close(f)
    sky_line_xlst = (sky_line_mat_clean[:, 1, :] .- (N_XPIX ÷ 2)) ./ N_XPIX
    sky_line_wlst = sky_line_mat_clean[:, 2, :]
    return sky_line_xlst, sky_line_wlst
end

# this takes in a filename and replaces the chip index (make "a" default approx)
function ingest_skyLines_exp(fname; chip_lst = ["a", "b", "c"])
    # println("Ingesting sky lines for $fname")
    sky_line_uxlst = Matrix{Float64}[]
    sky_line_fwlst = Matrix{Float64}[]
    sky_line_chipInt = Matrix{Int}[]
    for chip in chip_lst
        fnameloc = replace(fname, "_a_" => "_$(chip)_")
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

    msk_large_scatter = dropdims(nanzeroiqr(sky_line_uxlst, 2) .> 0.002, dims = 2)
    sky_line_uxlst[msk_large_scatter, :] .= NaN
    sky_line_fwlst[msk_large_scatter, :] .= NaN
    # dims are num_sky_lines x num_fibers
    return sky_line_uxlst, sky_line_fwlst, sky_line_chipInt
end

# this takes in a filename and replaces the chip index (make "a" default approx)
function ingest_fpiLines_exp(fname)
    fpi_line_uxlst = Matrix{Float64}[]
    fpi_line_uxlst_errs = Matrix{Float64}[]
    fpi_line_chipInt = Matrix{Int}[]
    for chip in ["a", "b", "c"]
        fnameloc = replace(fname, "_a_" => "_$(chip)_")
        if isfile(fnameloc)
            fpi_line_xlst, fpi_line_xlst_errs = ingest_fpiLines_file(fnameloc)
            chipIndx = getChipIndx(chip)
            push!(fpi_line_uxlst, fpi_line_xlst)
            push!(fpi_line_uxlst_errs, fpi_line_xlst_errs)
            push!(fpi_line_chipInt, chipIndx * ones(Int, size(fpi_line_xlst)))
        else
            push!(fpi_line_uxlst, [])
            push!(fpi_line_fwlst, [])
            push!(fpi_line_chipInt, [])
        end
    end
    fpi_line_uxlst = vcat(fpi_line_uxlst...)
    fpi_line_uxlst_errs = vcat(fpi_line_uxlst_errs...)
    fpi_line_chipInt = vcat(fpi_line_chipInt...)

    # dims are num_fpi_lines x num_fibers
    return fpi_line_uxlst, fpi_line_uxlst_errs, fpi_line_chipInt
end

function sky_wave_plots(
        fname_list, night_linParams, night_nlParams, night_wave_soln;
        dirNamePlots = "../outdir/plots/", plot_fibers = (1, 50, 100, 150, 200, 250, 300),
        plot_pixels = (1, 512, 1024, 1536, 2048), chip_names = ("a", "b", "c"))
    for chip in chip_names
        fname_list .= replace.(fname_list, "_$(chip)_" => "_")
    end
    unique_fname_list = unique(fname_list)

    sname = split(split(unique_fname_list[1], "/")[end], "_")
    tele, mjd = sname[(end - 3):(end - 2)]

    n_fnames = size(unique_fname_list, 1)
    n_lin_coeffs = size(night_linParams, 2)
    n_nl_coeffs = size(night_nlParams, 2)
    n_fibers = size(night_linParams, 1)
    fiber_inds = 1:n_fibers

    clims = (1, n_fnames)

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
                              "skywave_linParams_$(tele)_$(mjd).png"
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
                             "skywave_nlParams_$(tele)_$(mjd).png"
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
                title = "Tele: $(tele), MJD: $(mjd)\nChip $(chip_names[2]), Δλ = $(round(chip_b_diff,digits=3)) Å/pixel")
            axis_dict["a $(pind)"] = Axis(fig[pind, 1],
                ylabel = "FID$(fibid) Wavelength Offset (Å)",
                title = "\nChip $(chip_names[1]), Δλ = $(round(chip_a_diff,digits=3)) Å/pixel")
            axis_dict["c $(pind)"] = Axis(fig[pind, 3],
                title = "\nChip $(chip_names[3]), Δλ = $(round(chip_c_diff,digits=3)) Å/pixel")
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
                1:2048,
                chipWaveSoln[:, plot_fibers[pind], 1] .- night_wave_soln[:, plot_fibers[pind], 1] .+
                dx[plot_fibers[pind]] .* dlam_dx[:, plot_fibers[pind], 1],
                color = fname_ind * ones(2048),
                colorrange = clims)
            scatter!(axis_dict["b $(pind)"],
                1:2048,
                chipWaveSoln[:, plot_fibers[pind], 2] .- night_wave_soln[:, plot_fibers[pind], 2] .+
                dx[plot_fibers[pind]] .* dlam_dx[:, plot_fibers[pind], 2],
                color = fname_ind * ones(2048),
                colorrange = clims)
            scatter!(axis_dict["c $(pind)"],
                1:2048,
                chipWaveSoln[:, plot_fibers[pind], 3] .- night_wave_soln[:, plot_fibers[pind], 3] .+
                dx[plot_fibers[pind]] .* dlam_dx[:, plot_fibers[pind], 3],
                color = fname_ind * ones(2048),
                colorrange = clims)
        end
    end

    sky_wave_fiber_vs_pixel_Path = dirNamePlots *
                                   "skywave_per_fiber_vs_pixel_$(tele)_$(mjd).png"
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
                title = "Tele: $(tele), MJD: $(mjd)\nChip $(chip_names[2]), Δλ = $(round(chip_b_diff,digits=3)) Å/pixel")
            axis_dict["a $(pind)"] = Axis(fig[pind, 1],
                ylabel = "XPix$(plot_pixels[pind]) Wavelength Offset (Å)",
                title = "\nChip $(chip_names[1]), Δλ = $(round(chip_a_diff,digits=3)) Å/pixel")
            axis_dict["c $(pind)"] = Axis(fig[pind, 3],
                title = "\nChip $(chip_names[3]), Δλ = $(round(chip_c_diff,digits=3)) Å/pixel")
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
                color = fname_ind * ones(n_fibers),
                colorrange = clims)
            scatter!(axis_dict["b $(pind)"],
                301 .- fiber_inds,
                chipWaveSoln[plot_pixels[pind], :, 2] .- night_wave_soln[plot_pixels[pind], :, 2] .+
                dx .* dlam_dx[plot_pixels[pind], :, 2],
                color = fname_ind * ones(n_fibers),
                colorrange = clims)
            scatter!(axis_dict["c $(pind)"],
                301 .- fiber_inds,
                chipWaveSoln[plot_pixels[pind], :, 3] .- night_wave_soln[plot_pixels[pind], :, 3] .+
                dx .* dlam_dx[plot_pixels[pind], :, 3],
                color = fname_ind * ones(n_fibers),
                colorrange = clims)
        end
    end

    sky_wave_pixel_vs_fiber_Path = dirNamePlots *
                                   "skywave_per_pixel_vs_fiber_$(tele)_$(mjd).png"
    save(sky_wave_pixel_vs_fiber_Path, fig)

    return sky_wave_linParams_Path, sky_wave_nlParams_Path, sky_wave_fiber_vs_pixel_Path,
    sky_wave_pixel_vs_fiber_Path
end
