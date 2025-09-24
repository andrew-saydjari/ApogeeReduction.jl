using Pkg;
Pkg.instantiate();
using JLD2, ProgressMeter, ArgParse, SlackThreads, Glob, StatsBase, Optim
using Polynomials: Polynomial

using ApogeeReduction
using ApogeeReduction: nanmedian, nanzeropercentile
src_dir = "../../src/"
include(src_dir * "/makie_plotutils.jl")

## Parse command line arguments
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--outdir"
        required = true
        help = "outdir of large data run (many nights)"
        arg_type = String
        default = ""
    end
    return parse_args(s)
end

parg = parse_commandline()

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


function summarize_wave_solns(tele_loc,mjd_list,plotdir;
			      dporder=2)

    fname_loc = joinpath(parg["outdir"], "apred/$(mjd_list[1])/waveCalFPI_$(tele_loc)_$(mjd_list[1])_arclamp.h5") 
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

    all_ditherCorr_linParams = zeros(n_mjd,n_fibers,n_linParams)
    all_ditherCorrParams = zeros(n_mjd,n_fibers,dporder+1)
    fill!(all_ditherCorr_linParams,NaN)
    if dporder > 0
        all_ditherCorrParams[:,:,2] .= 1.0
    end
    #save a copy of the linear parameters to measure
    #dithers from for the other dates
    comp_linParams = copy(linParams)
    all_ditherCorr_linParams[1,:,:] .= comp_linParams

    for (mjd_ind,mjd_loc) in enumerate(mjd_list)
        fname_loc = joinpath(parg["outdir"], "apred/$(mjd_loc)/waveCalFPI_$(tele_loc)_$(mjd_loc)_arclamp.h5") 

        all_linParams[mjd_ind,:,:] .= load(fname_loc, "linParams")
        all_nlParams[mjd_ind,:,:] .= load(fname_loc, "nlParams")
        all_fpiCavityParams[mjd_ind,1] = load(fname_loc, "fpi_m0")
        all_fpiCavityParams[mjd_ind,2] = load(fname_loc, "fpi_m0_offset")
        all_fpiCavityParams[mjd_ind,3] = load(fname_loc, "fpi_cavity_size")

        if mjd_ind == 1
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

    clims = (minimum(mjd_list), maximum(mjd_list))

    fig = Figure(size = (1200, 400 * size(all_fpiCavityParams,2)), fontsize = 22)

    for pind in 1:size(all_fpiCavityParams,2)
        yvals = @view all_fpiCavityParams[:, pind]

        if pind == 1
            ax = Axis(fig[pind, 1],
                ylabel = "FPI Peak m0",
                title = "Tele: $(tele_loc), FPI Cavity Parameters")
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

for tele in ["apo", "lco"]
    good_mjds[:] .= false
    for (mjd_ind,mjd) in enumerate(poss_mjds)
        if isfile(joinpath(parg["outdir"], "apred/$(mjd)/waveCalFPI_$(tele)_$(mjd)_arclamp.h5"))
            good_mjds[mjd_ind] = true
        end
    end
    good_mjd_ints = map(x -> parse(Int, x), poss_mjds[good_mjds])

    println("Tele $(tele) has $(size(good_mjd_ints,1)) FPI-based wavelength solutions along outdir=$(parg["outdir"])")
    if size(good_mjd_ints,1) == 0
        continue
    end

#    plotdir = joinpath(parg["outdir"], "monitor") 
    plotdir = "../outdir/monitor" 
    out_plot_fnames = summarize_wave_solns(tele,good_mjd_ints,plotdir)
    thread = SlackThread()
    thread("Wavelength solution time stability summary for $(tele) using data at $(parg["outdir"])")
    thread("Found $(size(good_mjd_ints,1)) unique dates (in the MJD range of $(good_mjd_ints[begin]) to $(good_mjd_ints[end])) with useful FPI-measured wavelength solutions")

    for (framePath,frameString) in out_plot_fnames
        thread(frameString, framePath)
    end
end




