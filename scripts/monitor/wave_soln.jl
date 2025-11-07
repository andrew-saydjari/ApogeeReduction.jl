using Pkg;
Pkg.instantiate();
using JLD2, ProgressMeter, ArgParse, SlackThreads, Glob, StatsBase, Optim, HDF5
using Polynomials: Polynomial

using ApogeeReduction
using ApogeeReduction: nanmedian, nanzeropercentile, safe_jldsave, params2ChipPolyParams!
src_dir = "../../src/"
include(src_dir * "/makie_plotutils.jl")

## Parse command line arguments
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--outdir"
        required = false
        help = "outdir of large data run (many nights)"
        arg_type = String
        default = "/mnt/home/sdssv/ceph/work/asaydjari/2025_10_03/outdir/"
        "--plotdir"
        required = false
        help = "plotdir where output figures will be saved"
        arg_type = String
        default = "/mnt/home/kmckinnon/ceph/scratch/20250923/monitor/"
        "--max_plot"
        required = false
        help = "max number of MJD to use when plotting outputs for FPI"
        arg_type = Int
        default = 10000000
        "--min_mjd"
        required = false
        help = "min MJD to use when finding wave solutions"
        arg_type = Int
        default = 0
    end
    return parse_args(s)
end

parg = parse_commandline()

max_plot = parg["max_plot"]
min_mjd = parg["min_mjd"]


function dither_from_params_loss(dither_params,orig_poly,target_params;
				  return_new=false)
    new_coeffs = orig_poly(Polynomial(dither_params)).coeffs
    if return_new
        return new_coeffs
    else
        l2 = sum((new_coeffs[begin:size(target_params,1)] .- target_params) .^ 2)
        if size(new_coeffs,1) > size(target_params,1)
            l2 += sum((new_coeffs[size(target_params,1)+1:end]) .^ 2)
        end
        return l2
    end
end

function params_from_dither_loss(target_params,dither_params,orig_params;
				  return_new=false)
    new_coeffs = Polynomial(target_params)(Polynomial(dither_params)).coeffs
    if return_new
        return new_coeffs
    else
        l2 = sum((new_coeffs[begin:size(target_params,1)] .- orig_params) .^ 2)
        if size(new_coeffs,1) > size(orig_params,1)
            l2 += sum((new_coeffs[size(orig_params,1)+1:end]) .^ 2)
        end
        return l2
    end
end


function save_and_score_all_wave_solns(tele,mjd_list;
			      dporder=2,plotdir=parg["plotdir"],outdir=parg["outdir"])

    waveChi2_outname = joinpath(plotdir, "bulkWaveSoln_chi2_summary.h5")
    if isfile(waveChi2_outname)
        io_mode = "r+"
    else
        io_mode = "w"
    end
    f = h5open(waveChi2_outname, io_mode)


    for mjd in mjd_list
        curr_best_wave_type,ditherCorr_linParams,night_nlParams,ditherCorrParams,linParam_chi2,nlParam_chi2,ditherParam_chi2 = score_wave_soln(tele,mjd;dporder=2,plotdir=parg["plotdir"],outdir=parg["outdir"])

	savedict = Dict(
			"best_wave_type" => curr_best_wave_type,
			"ditherCorr_linParams" => ditherCorr_linParams,
			"night_nlParams" => night_nlParams,
			"ditherCorrParams" => ditherCorrParams,
			"linParam_chi2s" => linParam_chi2,
			"nlParam_chi2s" => nlParam_chi2,
			"ditherParam_chi2s" => ditherParam_chi2,
			)

	for (key,val) in savedict
            savekey = "$(tele)/$(mjd)/$(key)"
            if haskey(f, "$(tele)")
	        if haskey(f, "$(tele)/$(mjd)")
                    if haskey(f, savekey)
                        delete_object(f, savekey)
                    end
		end
	    end
            f[savekey] = val
	end
    end
    close(f)

    return waveChi2_outname

end


function score_wave_soln(tele_loc,mjd_loc;
			      dporder=2,plotdir=parg["plotdir"],outdir=parg["outdir"])

    outname = joinpath(plotdir, "fpiWaveSoln_summary.h5")
    med_linParams = load(outname, "$(tele_loc)/med_linParams")
    med_nlParams = load(outname, "$(tele_loc)/med_nlParams")
    med_ditherParams = load(outname, "$(tele_loc)/med_ditherParams")
    med_fpi_m0 = load(outname, "$(tele_loc)/med_fpi_m0")
    med_fpi_m0_offset = load(outname, "$(tele_loc)/med_fpi_m0_offset")
    med_fpi_cavity_size = load(outname, "$(tele_loc)/med_fpi_cavity_size")
    med_linParam_errs = load(outname, "$(tele_loc)/med_linParam_errs")
    med_nlParam_errs = load(outname, "$(tele_loc)/med_nlParam_errs")
    med_ditherParam_errs = load(outname, "$(tele_loc)/med_ditherParam_errs")
    med_fpi_m0_err = load(outname, "$(tele_loc)/med_fpi_m0_err")
    med_fpi_m0_offset_err = load(outname, "$(tele_loc)/med_fpi_m0_offset_err")
    med_fpi_cavity_size_err = load(outname, "$(tele_loc)/med_fpi_cavity_size_err")
    med_n_linParams = load(outname,"$(tele_loc)/n_linParams")
    med_n_nlParams = load(outname,"$(tele_loc)/n_nlParams")
    med_n_ditherParams = load(outname,"$(tele_loc)/n_ditherParams")

    fname_loc = joinpath(outdir, "wavecal/wavecalNightAve_$(tele_loc)_$(mjd_loc).h5")
    curr_best_wave_type = load(fname_loc, "best_wave_type")
    night_linParams = load(fname_loc, "$(curr_best_wave_type)/nightAve_linParams")
    night_nlParams = load(fname_loc, "$(curr_best_wave_type)/nightAve_nlParams")

    if length(size(night_linParams)) == 3
       night_linParams = dropdims(night_linParams,dims=3)
    end
    if length(size(night_nlParams)) == 3
       night_nlParams = dropdims(night_nlParams,dims=3)
    end

    n_linParams = size(night_linParams,2)
    n_nlParams = size(night_nlParams,2)
    if n_nlParams < med_n_nlParams
        new_night_nlParams = zeros(size(med_nlParams))

        med_n_nlParams_per_chip = med_n_nlParams ÷ 2
        n_nlParams_per_chip = n_nlParams ÷ 2
        med_cporder = med_n_nlParams_per_chip - 1
        cporder = n_nlParams_per_chip - 1

	new_night_nlParams[:,2] .= 1.0
	new_night_nlParams[:,med_n_nlParams_per_chip+2] .= 1.0
	new_night_nlParams[:,begin:n_nlParams_per_chip] .= night_nlParams[:,begin:n_nlParams_per_chip]
	new_night_nlParams[:,med_n_nlParams_per_chip+1:med_n_nlParams_per_chip+1+n_nlParams_per_chip] .= night_nlParams[:,n_nlParams_per_chip+1:end]
        
	night_nlParams = new_night_nlParams
        n_nlParams = size(night_nlParams,2)
    end
    ditherCorr_linParams = zeros(N_FIBERS,n_linParams)
    ditherCorrParams = zeros(N_FIBERS,dporder+1)
    fill!(ditherCorr_linParams,NaN)
    if dporder > 0
        ditherCorrParams[:,2] .= 1.0
    end
    #measure the dither parameters that get best agreement between
    #the observed linParams
    for fiber_ind in 1:N_FIBERS
        orig_poly = Polynomial(night_linParams[fiber_ind,:])

        function dither_nonlinear_loss_fit_partial(inparams)
            dither_from_params_loss(inparams,orig_poly,med_linParams[fiber_ind,:])
        end
 
        res = optimize(
            dither_nonlinear_loss_fit_partial, 
            ditherCorrParams[fiber_ind,:], 
     	    LBFGS(), Optim.Options(show_trace = false))
    
        ditherParamsOpt = Optim.minimizer(res)
        ditherCorrParams[fiber_ind,:] .= ditherParamsOpt
        ditherCorr_linParams[fiber_ind,:] .= dither_from_params_loss(
                                                     ditherParamsOpt,orig_poly,
						     med_linParams[fiber_ind,:],
						     return_new=true)[begin:n_linParams]
    end

    n_lin_param_comp = min(med_n_linParams,n_linParams)
    n_nl_param_comp = min(med_n_nlParams,n_nlParams)
    n_dither_param_comp = min(med_n_ditherParams,dporder+1)

    linParam_chi2 = ((ditherCorr_linParams[:,begin:n_lin_param_comp] .- med_linParams[:,begin:n_lin_param_comp]) ./ med_linParam_errs[:,begin:n_lin_param_comp]) .^ 2
    if med_n_linParams > n_lin_param_comp
        linParam_chi2 .+= (med_linParams[:,n_lin_param_comp:end] ./ med_linParam_errs[:,n_lin_param_comp:end]) .^ 2
    end
    nlParam_chi2 = ((night_nlParams[:,begin:n_nl_param_comp] .- med_nlParams[:,begin:n_nl_param_comp]) ./ med_nlParam_errs[:,begin:n_nl_param_comp]) .^ 2
    if med_n_nlParams > n_nl_param_comp 
        nlParam_chi2 .+= (med_nlParams[:,n_nl_param_comp:end] ./ med_nlParam_errs[:,n_nl_param_comp:end]) .^ 2
    end
    ditherParam_chi2 = ((ditherCorrParams[:,begin:n_dither_param_comp] .- med_ditherParams[:,begin:n_dither_param_comp]) ./ med_ditherParam_errs[:,begin:n_dither_param_comp]) .^ 2
    if med_n_ditherParams > n_dither_param_comp 
        ditherParam_chi2 .+= (med_ditherParams[:,n_dither_param_comp:end] ./ med_ditherParam_errs[:,n_dither_param_comp:end]) .^ 2
    end

    return curr_best_wave_type,ditherCorr_linParams,night_nlParams,ditherCorrParams,linParam_chi2,nlParam_chi2,ditherParam_chi2
end



function summarize_fpi_wave_solns(tele_loc,mjd_list;
			      dporder=2,plotdir=parg["plotdir"],outdir=parg["outdir"])

    comp_mjd_ind = 1
    comp_mjd_ind = size(mjd_list,1)
    fname_loc = joinpath(outdir, "apred/$(mjd_list[comp_mjd_ind])/waveCalFPI_$(tele_loc)_$(mjd_list[comp_mjd_ind])_arclamp.h5") 
    linParams = load(fname_loc, "linParams")
    nlParams = load(fname_loc, "nlParams")
    ditherParams = load(fname_loc, "ditherParams")
    fpi_m0 = load(fname_loc, "fpi_m0")
    fpi_m0_offset = load(fname_loc, "fpi_m0_offset")
    fpi_cavity_size = load(fname_loc, "fpi_cavity_size")

    n_mjd = size(mjd_list,1)

    n_fibers = size(linParams,1)
    fiber_inds = collect(1:n_fibers)
    n_linParams = size(linParams,2)
    n_nlParams = size(nlParams,2)
    n_ditherParams = size(ditherParams,2)
    all_linParams = zeros(n_mjd,n_fibers,n_linParams)
    all_nlParams = zeros(n_mjd,n_fibers,n_nlParams)
    all_fpiCavityParams = zeros(n_mjd, 3)

    fill!(all_linParams,NaN)
    fill!(all_nlParams,NaN)
    fill!(all_fpiCavityParams,NaN)

    med_linParams = zeros(n_fibers,n_linParams) 
    med_nlParams = zeros(n_fibers,n_nlParams)
    med_ditherParams = zeros(n_fibers,dporder+1)

    all_ditherCorr_linParams = zeros(n_mjd,n_fibers,n_linParams)
    all_ditherCorrParams = zeros(n_mjd,n_fibers,dporder+1)
    fill!(all_ditherCorr_linParams,NaN)
    if dporder > 0
        all_ditherCorrParams[:,:,2] .= 1.0
    end
    #save a copy of the linear parameters to measure
    #dithers from for the other dates
    comp_linParams = copy(linParams)
    true_med_linParams = zeros(n_fibers,n_linParams)
    all_ditherCorr_linParams[comp_mjd_ind,:,:] .= comp_linParams

    n_repeat = 3
    for r_ind in 1:n_repeat
        for (mjd_ind,mjd_loc) in enumerate(mjd_list)
            fname_loc = joinpath(outdir, "apred/$(mjd_loc)/waveCalFPI_$(tele_loc)_$(mjd_loc)_arclamp.h5") 
            curr_m0 = load(fname_loc, "fpi_m0")
            if (curr_m0 < 3500) | (curr_m0 > 4500)
                continue
            end

            if r_ind == 1
                all_linParams[mjd_ind,:,:] .= load(fname_loc, "linParams")
                all_nlParams[mjd_ind,:,:] .= load(fname_loc, "nlParams")
                all_fpiCavityParams[mjd_ind,1] = load(fname_loc, "fpi_m0")
                all_fpiCavityParams[mjd_ind,2] = load(fname_loc, "fpi_m0_offset")
                all_fpiCavityParams[mjd_ind,3] = load(fname_loc, "fpi_cavity_size")
            end

            if (mjd_ind == comp_mjd_ind) & (r_ind == 1)
                continue
            end

            #measure the dither parameters that get best agreement between
            #the observed linParams
            for fiber_ind in 1:n_fibers
                orig_poly = Polynomial(all_linParams[mjd_ind,fiber_ind,:])

                function dither_nonlinear_loss_fit_partial(inparams)
                    dither_from_params_loss(inparams,orig_poly,comp_linParams[fiber_ind,:])
                end
 
                res = optimize(
                    dither_nonlinear_loss_fit_partial, 
        	        all_ditherCorrParams[mjd_ind,fiber_ind,:], 
     	        LBFGS(), Optim.Options(show_trace = false))
    
                ditherParamsOpt = Optim.minimizer(res)
                all_ditherCorrParams[mjd_ind,fiber_ind,:] .= ditherParamsOpt
                all_ditherCorr_linParams[mjd_ind,fiber_ind,:] .= dither_from_params_loss(
         						ditherParamsOpt,orig_poly,
    						comp_linParams[fiber_ind,:],return_new=true)[begin:n_linParams]
            end
        end

        med_linParams .= dropdims(nanmedian(all_ditherCorr_linParams,1),dims=1)
        med_nlParams .= dropdims(nanmedian(all_nlParams,1),dims=1)
        med_ditherParams .= dropdims(nanmedian(all_ditherCorrParams,1),dims=1)

        if r_ind == n_repeat
            continue
        end

        #use the med dither and med linParams to define the original med_linParams
        for fiber_ind in 1:n_fibers
            function orig_param_nonlinear_loss_fit_partial(inparams)
               params_from_dither_loss(inparams,med_ditherParams[fiber_ind,:],med_linParams[fiber_ind,:])
            end
 
            res = optimize(
                orig_param_nonlinear_loss_fit_partial,
    	        med_linParams[fiber_ind,:],
                LBFGS(), Optim.Options(show_trace = false))
     
            true_med_linParams[fiber_ind,:] .= Optim.minimizer(res)
        end
        comp_linParams .= true_med_linParams
    end

    good_mjds = isfinite.(all_fpiCavityParams[:,1])
    good_mjd_list = mjd_list[good_mjds]
    good_mjd_range = minimum(good_mjd_list),maximum(good_mjd_list)
    final_fpiCavityParams = dropdims(nanmedian(all_fpiCavityParams,1),dims=1)
    med_fpi_m0 = final_fpiCavityParams[1]
    med_fpi_m0_offset = final_fpiCavityParams[2]
    med_fpi_cavity_size = final_fpiCavityParams[3]

    med_linParam_errs = dropdims(0.5 .* diff(nanzeropercentile(all_ditherCorr_linParams,dims=1,percent_vec=[16,84]),dims=1),dims=1)
    med_nlParam_errs = dropdims(0.5 .* diff(nanzeropercentile(all_nlParams,dims=1,percent_vec=[16,84]),dims=1),dims=1)
    med_ditherParam_errs = dropdims(0.5 .* diff(nanzeropercentile(all_ditherCorrParams,dims=1,percent_vec=[16,84]),dims=1),dims=1)
    final_fpiCavityParam_errs = dropdims(0.5 .* diff(nanzeropercentile(all_fpiCavityParams,dims=1,percent_vec=[16,84]),dims=1),dims=1)

    med_fpi_m0_err = final_fpiCavityParam_errs[1]
    med_fpi_m0_offset_err = final_fpiCavityParam_errs[2]
    med_fpi_cavity_size_err = final_fpiCavityParam_errs[3]

#    linParam_chi2 = ((all_ditherCorr_linParams .- med_linParams) ./ med_linParam_errs) .^ 2
#    nlParam_chi2 = ((all_nlParams .- med_nlParams) ./ med_nlParam_errs) .^ 2
#    ditherParam_chi2 = ((all_ditherCorrParams .- med_ditherParams) ./ med_ditherParam_errs) .^ 2

    n_nlParams_per_chip = n_nlParams ÷ 2
    cporder = n_nlParams_per_chip - 1
    wporder = n_linParams - 1
    final_porder = wporder * cporder
    
    linParams_per_chip = zeros(Float64, (N_FIBERS, N_CHIPS, final_porder + 1))

    chipPolyParams = zeros(Float64, (N_CHIPS, cporder + 1))
    if cporder > 0
        chipPolyParams[:,2] .= 1.0
    end

    for fiber_ind in 1:n_fibers
        curr_nlParams = nlParams[fiber_ind, :]
        params2ChipPolyParams!(chipPolyParams, curr_nlParams, cporder)
        for chip_ind in 1:N_CHIPS
            new_coeffs = Polynomial(med_linParams[fiber_ind,:])(Polynomial(chipPolyParams[chip_ind,:])).coeffs
            n_new_coeffs = size(new_coeffs,1)
            linParams_per_chip[fiber_ind,chip_ind,begin:n_new_coeffs] .= new_coeffs
        end
    end

    savedict = Dict(
		    "med_linParams" => med_linParams,
		    "med_nlParams" => med_nlParams,
		    "med_ditherParams" => med_ditherParams,
		    "med_fpi_m0" => med_fpi_m0,
		    "med_fpi_m0_offset" => med_fpi_m0_offset,
		    "med_fpi_cavity_size" => med_fpi_cavity_size,
		    "med_linParam_errs" => med_linParam_errs,
		    "med_nlParam_errs" => med_nlParam_errs,
		    "med_ditherParam_errs" => med_ditherParam_errs,
		    "med_fpi_m0_err" => med_fpi_m0_err,
		    "med_fpi_m0_offset_err" => med_fpi_m0_offset_err,
		    "med_fpi_cavity_size_err" => med_fpi_cavity_size_err,
		    "n_used_mjds" => sum(good_mjds),
		    "used_mjd_min" => good_mjd_range[1],
		    "used_mjd_max" => good_mjd_range[2],
		    "used_mjds" => good_mjd_list,
		    "data_outdir" => outdir,
                    "n_linParams" => n_linParams,
                    "n_nlParams" => n_nlParams,
                    "n_ditherParams" => dporder+1,
		    )

    for chip_ind in 1:N_CHIPS
        savedict["chip$(CHIP_LIST[chip_ind])_waveParams"] = linParams_per_chip[:,chip_ind,:]
    end

    outname = joinpath(plotdir, "fpiWaveSoln_summary.h5")
    if isfile(outname)
        io_mode = "r+"
    else
        io_mode = "w"
    end
    f = h5open(outname, io_mode)

    for (key,val) in savedict
        savekey = "$(tele_loc)/$(key)"
        if haskey(f, savekey)
            delete_object(f, savekey)
        end
        f[savekey] = val
    end

    close(f)

    clims = (minimum(mjd_list), maximum(mjd_list))

    fig = Figure(size = (1200, 400 * size(all_fpiCavityParams,2)), fontsize = 22)

    for pind in 1:size(all_fpiCavityParams,2)
        yvals = @view all_fpiCavityParams[:, pind]

        if pind == 1
            med_m0 = nanmedian(yvals)
            if all(yvals .== med_m0)
                limits = (nothing,(med_m0-1,med_m0+1))
            else
                limits = (nothing,nothing)
            end
            ax = Axis(fig[pind, 1],
                ylabel = "FPI Peak m0",
                title = "Tele: $(tele_loc), FPI Cavity Parameters",
		limits = limits)
        elseif pind == size(all_fpiCavityParams,2)
            ax = Axis(fig[pind, 1],
                xlabel = "MJD",
		ylabel = "Cavity Size (Å)")
        else
            ax = Axis(fig[pind, 1],
                ylabel = "Peak Non-integer Offset")
        end

        scatter!(ax, mjd_list, yvals, color=mjd_list, colorrange=clims)
        med_val = nanmedian(yvals)
        hlines!(ax, med_val, linestyle = :dash)
    end

    fpiParams_framePath = joinpath(plotdir, "fpiCavityParams_summary_$(tele_loc).png")
    fpiParams_string = "FPI cavity parameters for $(tele_loc) from MJD $(mjd_list[begin]) to $(mjd_list[end])"
    save(fpiParams_framePath, fig)

    fig = Figure(size = (1200, 400 * (dporder+1)), fontsize = 22)

    for pind in 1:dporder+1
        yvals = nanmedian(all_ditherCorrParams[:,:,pind],2)[:,1]

        mult = 1
        if pind == 1
            ax = Axis(fig[pind, 1],
                ylabel = "Dither Offset (pixels)",
                title = "Tele: $(tele_loc), Dither Correction Parameters")
            mult = 2048
        elseif pind == size(all_ditherCorrParams,2)
            ax = Axis(fig[pind, 1],
                xlabel = "MJD",
                ylabel = "Dither Parameter $(pind)")
        else
            ax = Axis(fig[pind, 1],
                ylabel = "Dither Parameter $(pind)")
        end

        scatter!(ax, mjd_list, yvals .* mult, color=mjd_list, colorrange=clims)
        med_val = nanmedian(yvals)
        hlines!(ax, med_val, linestyle = :dash)
    end

    medDitherParams_framePath = joinpath(plotdir, "medDitherParams_summary_$(tele_loc).png")
    medDitherParams_string = "Nightly median dither parameter corrections for $(tele_loc) from MJD $(mjd_list[begin]) to $(mjd_list[end])"
    save(medDitherParams_framePath, fig)


    fig = Figure(size = (1200, 400 * n_linParams), fontsize = 22)

    for pind in 1:n_linParams
        yvals = @view all_linParams[:, :, pind]
#        yvals = @view all_ditherCorr_linParams[:, :, pind]
        summary = nanzeropercentile(yvals, percent_vec = [16, 50, 84], dims=(1,2))
        n_sigma = 8
        ylim = (summary[2] - n_sigma * (summary[2] - summary[1]),
            summary[2] + n_sigma * (summary[3] - summary[2]))
        outside_limits = (yvals .< ylim[1]) .| (yvals .> ylim[2])
        if sum(outside_limits) == 0
            ylim = nothing
        end
        limits = (nothing, ylim)

        if pind == 1
            ax = Axis(fig[pind, 1],
                limits = limits,
                ylabel = "Linear Parameter $(pind)",
                title = "Tele: $(tele_loc), Wave. Soln. Linear Params")
        elseif pind == n_linParams
            ax = Axis(fig[pind, 1],
                limits = limits,
                xlabel = "FIBERID",
                ylabel = "Linear Parameter $(pind)")
        else
            ax = Axis(fig[pind, 1],
                limits = limits,
                ylabel = "Linear Parameter $(pind)")
        end

        for (mjd_ind,mjd_loc) in enumerate(mjd_list)
            scatter!(ax, 301 .- fiber_inds, all_linParams[mjd_ind, :, pind],
	     	     color=mjd_loc,colorrange=clims)
        end
    end

    linParams_framePath = joinpath(plotdir, "linParam_waveSoln_summary_$(tele_loc).png")
    save(linParams_framePath, fig)
    linParams_string = "Linear wavelength solution parameters for $(tele_loc) from MJD $(mjd_list[begin]) to $(mjd_list[end])"

    fig = Figure(size = (1200, 400 * n_linParams), fontsize = 22)

    for pind in 1:n_linParams
        yvals = @view all_ditherCorr_linParams[:, :, pind]
        summary = nanzeropercentile(yvals, percent_vec = [16, 50, 84], dims=(1,2))
        n_sigma = 8
        ylim = (summary[2] - n_sigma * (summary[2] - summary[1]),
            summary[2] + n_sigma * (summary[3] - summary[2]))
        outside_limits = (yvals .< ylim[1]) .| (yvals .> ylim[2])
        if sum(outside_limits) == 0
            ylim = nothing
        end
        limits = (nothing, ylim)

        if pind == 1
            ax = Axis(fig[pind, 1],
                limits = limits,
                ylabel = "Linear Parameter $(pind)",
                title = "Tele: $(tele_loc), Wave. Soln. Linear Params")
        elseif pind == n_linParams
            ax = Axis(fig[pind, 1],
                limits = limits,
                xlabel = "FIBERID",
                ylabel = "Linear Parameter $(pind)")
        else
            ax = Axis(fig[pind, 1],
                limits = limits,
                ylabel = "Linear Parameter $(pind)")
        end

        for (mjd_ind,mjd_loc) in enumerate(mjd_list)
            scatter!(ax, 301 .- fiber_inds, all_ditherCorr_linParams[mjd_ind, :, pind],
	     	     color=mjd_loc,colorrange=clims)
        end
    end

    linParamsCorr_framePath = joinpath(plotdir, "linParamCorr_waveSoln_summary_$(tele_loc).png")
    save(linParamsCorr_framePath, fig)
    linParamsCorr_string = "Dither-corrected linear wavelength solution parameters for $(tele_loc) from MJD $(mjd_list[begin]) to $(mjd_list[end])"

    fig = Figure(size = (1200, 400 * n_nlParams), fontsize = 22)

    for pind in 1:n_nlParams
        yvals = @view all_nlParams[:, :, pind]
        summary = nanzeropercentile(yvals, percent_vec = [16, 50, 84], dims=(1,2))
        n_sigma = 8
        ylim = (summary[2] - n_sigma * (summary[2] - summary[1]),
            summary[2] + n_sigma * (summary[3] - summary[2]))
        outside_limits = (yvals .< ylim[1]) .| (yvals .> ylim[2])
        if sum(outside_limits) == 0
            ylim = nothing
        end
        limits = (nothing, ylim)

        if pind == 1
            ax = Axis(fig[pind, 1],
                limits = limits,
                ylabel = "Non-Linear Parameter $(pind)",
                title = "Tele: $(tele_loc), Non-Linear Chip Params")
        elseif pind == n_nlParams
            ax = Axis(fig[pind, 1],
                limits = limits,
                xlabel = "FIBERID",
                ylabel = "Non-Linear Parameter $(pind)")
        else
            ax = Axis(fig[pind, 1],
                limits = limits,
                ylabel = "Non-Linear Parameter $(pind)")
        end

        for (mjd_ind,mjd_loc) in enumerate(mjd_list)
            scatter!(ax, 301 .- fiber_inds, all_nlParams[mjd_ind, :, pind],
	     	     color=mjd_loc,colorrange=clims)
        end
    end

    nlParams_framePath = joinpath(plotdir, "nlParam_waveSoln_summary_$(tele_loc).png")
    save(nlParams_framePath, fig)
    nlParams_string = "Non-linear chip parameters for $(tele_loc) from MJD $(mjd_list[begin]) to $(mjd_list[end])"

    fig = Figure(size = (1200, 400 * (dporder+1)), fontsize = 22)

    for pind in 1:(dporder+1)
        yvals = @view all_ditherCorrParams[:, :, pind]
        summary = nanzeropercentile(yvals, percent_vec = [16, 50, 84], dims=(1,2))
        n_sigma = 5
        ylim = (summary[2] - n_sigma * (summary[2] - summary[1]),
            summary[2] + n_sigma * (summary[3] - summary[2]))
        outside_limits = (yvals .< ylim[1]) .| (yvals .> ylim[2])
        if sum(outside_limits) == 0
            ylim = nothing
        end
        limits = (nothing, ylim)

        mult = 1
        if pind == 1
            mult = 2048
            if !isnothing(ylim)
                limits = (nothing, ylim .* mult)
            end
            ax = Axis(fig[pind, 1],
                limits = limits,
                ylabel = "Dither Offset (pixels)",
                title = "Tele: $(tele_loc), Dither Params")
        elseif pind == (dporder+1)
            ax = Axis(fig[pind, 1],
                limits = limits,
                xlabel = "FIBERID",
                ylabel = "Dither Parameter $(pind)")
        else
            ax = Axis(fig[pind, 1],
                limits = limits,
                ylabel = "Dither Parameter $(pind)")
        end

        for (mjd_ind,mjd_loc) in enumerate(mjd_list)
            scatter!(ax, 301 .- fiber_inds, all_ditherCorrParams[mjd_ind, :, pind] .* mult,
	     	     color=mjd_loc,colorrange=clims)
        end
    end

    ditherParams_framePath = joinpath(plotdir, "ditherParam_waveSoln_summary_$(tele_loc).png")
    save(ditherParams_framePath, fig)
    ditherParams_string = "Dither parameters that put all wavelength solutions of similar reference for $(tele_loc) from MJD $(mjd_list[begin]) to $(mjd_list[end])"

    return ((fpiParams_framePath,fpiParams_string),
	    (medDitherParams_framePath,medDitherParams_string),
	    (linParams_framePath,linParams_string),
	    (linParamsCorr_framePath,linParamsCorr_string),
	    (nlParams_framePath,nlParams_string),
	    (ditherParams_framePath,ditherParams_string),
	    )
end

poss_mjds = sort(readdir(joinpath(parg["outdir"], "apred/")))
good_mjds = zeros(Bool, size(poss_mjds,1))
good_mjds_fpi = zeros(Bool, size(poss_mjds,1))


#-loop over the FPI exposures in APO and LCO to summarize high quality wavelength solutions
#-combine median wavelength solution parameters into one file
#-loop over all the night average wave solutions (sky and fpi) and score them based on distance 
#     from expectation to flag the outliers
for tele in ["apo", "lco"]
    good_mjds[:] .= false
    good_mjds_fpi[:] .= false
    for (mjd_ind,mjd) in enumerate(poss_mjds)
        if isfile(joinpath(parg["outdir"], "apred/$(mjd)/waveCalFPI_$(tele)_$(mjd)_arclamp.h5"))
            good_mjds_fpi[mjd_ind] = true
	end
        if isfile(joinpath(parg["outdir"], "wavecal/wavecalNightAve_$(tele)_$(mjd).h5"))
            good_mjds[mjd_ind] = true
        end
    end
    good_mjd_ints = sort(map(x -> parse(Int, x), poss_mjds[good_mjds]))
    good_mjd_ints = good_mjd_ints[good_mjd_ints .>= min_mjd]

    good_mjd_ints_fpi = sort(map(x -> parse(Int, x), poss_mjds[good_mjds_fpi]))
    good_mjd_ints_fpi = good_mjd_ints_fpi[good_mjd_ints_fpi .>= min_mjd]
    

    println("Tele $(tele) has $(size(good_mjd_ints,1)) night average wavelength solutions along outdir=$(parg["outdir"])")
    println("Tele $(tele) has $(size(good_mjd_ints_fpi,1)) FPI-based wavelength solutions along outdir=$(parg["outdir"])")
    if size(good_mjd_ints_fpi,1) == 0
        continue
    end
    if size(good_mjd_ints_fpi,1) > max_plot
        step_size = ceil(Int,round(size(good_mjd_ints_fpi,1)/max_plot))
        good_mjd_ints_fpi = good_mjd_ints_fpi[begin:step_size:end]
        println("Only plotting $(size(good_mjd_ints_fpi,1)) MJDs with FPI data to be careful with memory.")
    end
    if size(good_mjd_ints,1) > max_plot
        step_size = ceil(Int,round(size(good_mjd_ints,1)/max_plot))
        good_mjd_ints = good_mjd_ints[begin:step_size:end]
        println("Only plotting $(size(good_mjd_ints,1)) MJDs with night wavelength solutions to be careful with memory.")
    end

    out_plot_fnames = summarize_fpi_wave_solns(tele,good_mjd_ints_fpi;plotdir=parg["plotdir"],outdir=parg["outdir"])
    thread = SlackThread()
    thread("Wavelength solution time stability summary for $(tele) using data at $(parg["outdir"])")
    thread("Used $(size(good_mjd_ints_fpi,1)) unique dates (in the MJD range of $(good_mjd_ints_fpi[begin]) to $(good_mjd_ints_fpi[end])) with useful FPI-measured wavelength solutions")

    for (framePath,frameString) in out_plot_fnames
        thread(frameString, framePath)
    end

    waveChi2_outname = save_and_score_all_wave_solns(tele,good_mjd_ints;
		                               dporder=2,plotdir=parg["plotdir"],outdir=parg["outdir"])

end



