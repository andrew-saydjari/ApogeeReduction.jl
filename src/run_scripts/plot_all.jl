# This is a script that makes the plots for a nightly processing and posts them to the slack channel.
# called by run_all.sh
## TODO add sky flux relFlux consistency check plots
using JLD2, ProgressMeter, ArgParse, Glob, StatsBase, Random, HDF5, DataFrames

src_dir = "../"
include(src_dir * "/fileNameHandling.jl")
include(src_dir * "/utils.jl")
include(src_dir * "/makie_plotutils.jl")
include(src_dir * "/ar1D.jl")
include(src_dir * "/ApogeeReduction.jl")

## Parse command line arguments
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--chips"
        required = false
        help = "chip name (a, b, c)"
        arg_type = String
        "--tele"
        required = true
        help = "telescope name (apo or lco)"
        arg_type = String
        default = ""
        "--mjd"
        required = false
        help = "mjd"
        arg_type = Int
        default = 0
        "--expid"
        required = false
        help = "exposure number to be run"
        arg_type = Int
        default = 1
        "--runlist"
        required = false
        help = "path name to hdf5 file with keys specifying list of exposures to run"
        arg_type = String
        default = ""
        "--outdir"
        required = true
        help = "directory where data is stored"
        arg_type = String
        default = "../outdir/"
        "--runname"
        required = true
        help = "name of the run (specifically almanac file)"
        arg_type = String
        default = "test"
    end
    return parse_args(s)
end

parg = parse_commandline()

println("Making plots for telescope: $(parg["tele"])")
println("Output directory: $(parg["outdir"])")

CHIP_LIST = ApogeeReduction.CHIP_LIST
FIRST_CHIP = ApogeeReduction.FIRST_CHIP
N_CHIPS = ApogeeReduction.N_CHIPS
N_FIBERS = ApogeeReduction.N_FIBERS

unique_mjds = if parg["runlist"] != ""
    subDic = load(parg["runlist"])
    unique(subDic["mjd"])
else
    [parg["mjd"]]
end

# Get the save directory for a plot (subfoldered by mjd)
function get_save_dir(mjd; outdir = parg["outdir"])
    dir = joinpath(outdir, "plots", string(mjd))
    mkpath(dir)
    return dir * "/"
end

if length(unique_mjds) == 0
    println(stderr, "ERROR: No MJDs found for plotting.")
    exit(1)
end

println("\nFound $(length(unique_mjds)) unique MJDs to process")

expid_list = if parg["runlist"] != ""
    subDic = load(parg["runlist"])
    subDic["expid"]
else
    [parg["expid"]]
end

list1DexpObject = []
list1DexpFPI = []
list1DexpArclamp = []
for mjd in unique_mjds
    df = read_almanac_exp_df(
        joinpath(parg["outdir"], "almanac/$(parg["runname"]).h5"), parg["tele"], mjd)
    function get_1d_name_partial(expid)
        if df.imagetyp[expid] == "Object"
            return parg["outdir"] * "/apred/$(mjd)/" * get_1d_name(expid, df, cal = true) * ".h5"
        else
            return nothing
        end
    end
    function get_1d_name_ARCLAMP_partial(expid)
        if (df.imagetyp[expid] == "ArcLamp") &
           ((df.lampthar[expid] == "T") | (df.lampune[expid] == "T"))
            return parg["outdir"] * "/apred/$(mjd)/" * get_1d_name(expid, df, cal = true) * ".h5"
        else
            return nothing
        end
    end
    function get_1d_name_FPI_partial(expid)
        if (df.imagetyp[expid] == "ArcLamp") & (df.lampthar[expid] == "F") &
           (df.lampune[expid] == "F")
            return parg["outdir"] * "/apred/$(mjd)/" * get_1d_name(expid, df, cal = true) * ".h5"
        else
            return nothing
        end
    end
    local1D = get_1d_name_partial.(expid_list)
    push!(list1DexpObject, filter(!isnothing, local1D))
    local1D_fpi = get_1d_name_FPI_partial.(expid_list)
    push!(list1DexpFPI, filter(!isnothing, local1D_fpi))
    local1D_arclamp = get_1d_name_ARCLAMP_partial.(expid_list)
    push!(list1DexpArclamp, filter(!isnothing, local1D_arclamp))
end
all1DObjecta = vcat(list1DexpObject...)
all1DFPIa = vcat(list1DexpFPI...)
all1DArclampa = vcat(list1DexpArclamp...)

println("\nFound:")
println("- $(length(all1DObjecta)) object exposures")
println("- $(length(all1DFPIa)) FPI exposures")
println("- $(length(all1DArclampa)) arc lamp exposures")

all1DObjectSkyPeaks = replace.(
    replace.(all1DObjecta, "ar1Dcal" => "skyLinePeaks"), "ar1D" => "skyLinePeaks")
all1DObjectSkyDither = replace.(
    replace.(all1DObjectSkyPeaks, "skyLinePeaks" => "waveCalNightskyDither"), "_$(FIRST_CHIP)_" => "_")
all1DObjectFPIDither = replace.(
    replace.(all1DObjectSkyPeaks, "skyLinePeaks" => "waveCalNightfpiDither"), "_$(FIRST_CHIP)_" => "_")

all1DObject_expid_strings = map(x -> x[end - 1], split.(all1DObjectSkyDither, "_"))
all1DObject_expids = map(x -> parse(Int, x), all1DObject_expid_strings)

function dither_plotter(fname_list, fname_expid_strings, mjd, tele)
    n_fnames = size(fname_list, 1)
    fname_ind = 1
    fname = fname_list[fname_ind]

    f = h5open(fname, "r+")

    linParams, nlParams, ditherParams, resid_vec, resid_xt = try
        read(f["linParams"]), read(f["nlParams"]), read(f["ditherParams"]), read(f["resid_vec"]),
        read(f["resid_xt"])
    catch
        println(fname)
        read(f["linParams"]), read(f["nlParams"]), read(f["ditherParams"]), read(f["resid_vec"]),
        read(f["resid_xt"])
    end

    close(f)

    n_dither_coeffs = size(ditherParams, 2)
    all_ditherParams = zeros(Float64, (n_fnames, N_FIBERS, n_dither_coeffs))
    fill!(all_ditherParams, NaN)
    all_ditherParams[fname_ind, :, :] .= ditherParams

    resid_plot_fnames = copy(fname_list)

    fiber_inds = collect(1:N_FIBERS)
    y_vals = 301 .- fiber_inds
    for (fname_ind, fname) in enumerate(fname_list)
        f = h5open(fname, "r+")

        linParams, nlParams, ditherParams, resid_vec, resid_xt = try
            read(f["linParams"]),
            read(f["nlParams"]), read(f["ditherParams"]), read(f["resid_vec"]), read(f["resid_xt"])
        catch
            println(fname)
            read(f["linParams"]),
            read(f["nlParams"]), read(f["ditherParams"]), read(f["resid_vec"]), read(f["resid_xt"])
        end

        close(f)

        all_ditherParams[fname_ind, :, :] .= ditherParams

        fig = Figure(size = (1200, 1200), fontsize = 22)
        x, y, z = Vector{Float64}(), Vector{Float64}(), Vector{Float64}()
        for fibIndx in 1:N_FIBERS
            x = vcat(x, resid_xt[fibIndx, :])
            y = vcat(y, y_vals[fibIndx] .* ones(size(resid_xt[fibIndx, :], 1)))
            z = vcat(z, resid_vec[fibIndx, :] .* 1000)
        end

        resid_summary = nanzeropercentile(z, percent_vec = [16, 84])
        resid_sigma = 0.5 * (resid_summary[2] - resid_summary[1])
        ax = Axis(fig[1, 1],
            ylabel = "FIBERID",
            xlabel = "Transformed X Peak Position",
            limits = ((-1.8, 1.8), (-10 + 1, 300 + 10)),
            title = "Sky Peak Residuals, σ = $(round(resid_sigma; digits = 3)) mÅ\nTele: $(tele), MJD: $(mjd), ExpID: $(fname_expid_strings[fname_ind])")

        clim = (-10.0, 10.0)
        z .= clamp.(z, clim[1], clim[2])
        cmap = :diverging_bkr_55_10_c35_n256
        scatter!(ax, x, y, color = z,
            markersize = 3, colorrange = clim,
            colormap = cmap)
        Colorbar(fig[1, 2], colorrange = clim, colormap = cmap,
            width = 20, height = Relative(1.0), label = "data - model (mÅ)")
        fpiPeakResiduals_Path = get_save_dir(mjd) *
                                "skyPeakResiduals_$(tele)_$(mjd)_$(fname_expid_strings[fname_ind]).png"
        save(fpiPeakResiduals_Path, fig)
        resid_plot_fnames[fname_ind] = fpiPeakResiduals_Path
    end
    all_ditherParams[:, :, 1] .*= 2048

    fig = Figure(size = (1200, 400 * n_dither_coeffs), fontsize = 22)

    for pind in 1:n_dither_coeffs
        yvals = @view all_ditherParams[:, :, pind]
        summary = nanzeropercentile(yvals, percent_vec = [16, 50, 84], dims = (1, 2))
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
                ylabel = "Dither Parameter $(pind)",
                title = "Tele: $(tele), MJD: $(mjd), Wave. Soln. Dither Params",
                limits = limits)
        elseif pind == n_dither_coeffs
            ax = Axis(fig[pind, 1],
                xlabel = "FIBERID",
                limits = limits,
                ylabel = "Dither Parameter $(pind)")
        else
            ax = Axis(fig[pind, 1],
                limits = limits,
                ylabel = "Dither Parameter $(pind)")
        end

        for fname_ind in 1:n_fnames
            scatter!(ax, 301 .- fiber_inds, all_ditherParams[fname_ind, :, pind],
                color = fname_ind .* ones(N_FIBERS), colorrange = (1, n_fnames))
        end
        if pind == 1
        end
    end
    ditherParams_Path = get_save_dir(mjd) *
                        "nightSkywave_ditherParams_$(tele)_$(mjd).png"
    save(ditherParams_Path, fig)

    return ditherParams_Path, resid_plot_fnames
end

if length(unique_mjds) > 1
    min_mjd, max_mjd = extrema(unique_mjds)
    println("\nGenerating wavelength solution stability plots for MJDs $min_mjd to $max_mjd")
else
    println("\nGenerating wavelength solution stability plots for MJD $(unique_mjds[1])")
end

for mjd_ind in 1:size(unique_mjds, 1)
    tele = parg["tele"]
    mjd = unique_mjds[mjd_ind]
    waveSoln_labels = ["Sky Line"]
    waveSoln_types = ["sky"]

    for j in 1:size(waveSoln_labels, 1)
        savePath = get_save_dir(unique_mjds[mjd_ind]) *
                   "$(waveSoln_types[j])wave_linParams_$(parg["tele"])_$(unique_mjds[mjd_ind]).png"
        if !isfile(savePath)
            println("Could not find any $(waveSoln_labels[j]) wavelength solution summary figures.")
            continue
        end
    end

    ditherParams_Path, resid_plot_fnames = dither_plotter(
        all1DObjectSkyDither, all1DObject_expid_strings, mjd, tele)

    outname = parg["outdir"] *
              "/apred/$(unique_mjds[mjd_ind])/waveCalFPI_$(parg["tele"])_$(unique_mjds[mjd_ind])_ARCLAMP.h5"
    if !isfile(outname)
        println("Could not find nightly FPI wavelength solution at $(outname)")
        continue
    end

    ditherParams_Path, resid_plot_fnames = dither_plotter(
        all1DObjectFPIDither, all1DObject_expid_strings, mjd, tele)

    f = h5open(outname, "r+")

    linParams, nlParams, ditherParams, exposure_names = try
        read(f["linParams"]), read(f["nlParams"]), read(f["ditherParams"]),
        read(f["exposure_names"])
    catch
        println(outname)
        read(f["linParams"]), read(f["nlParams"]), read(f["ditherParams"]),
        read(f["exposure_names"])
    end

    resid_vec, resid_chipInts, resid_exp_ints, resid_peak_ints,
    resid_used_in_fit, resid_xt = try
        read(f["resid_vec"]), read(f["resid_chipInts"]), read(f["resid_exp_ints"]),
        read(f["resid_peak_ints"]), read(f["resid_used_in_fit"]), read(f["resid_xt"])
    catch
        println(outname)
        read(f["resid_vec"]), read(f["resid_chipInts"]), read(f["resid_exp_ints"]),
        read(f["resid_peak_ints"]), read(f["resid_used_in_fit"]), read(f["resid_xt"])
    end

    close(f)

    n_fnames = size(exposure_names, 1)
    n_lin_coeffs = size(linParams, 2)
    n_nl_coeffs = size(nlParams, 2)
    n_dither_coeffs = size(ditherParams, 2)
    fiber_inds = 1:N_FIBERS

    fname_expid_strings = map(x -> x[end - 2], split.(exposure_names, "_"))
    fname_expids = map(x -> parse(Int, x), fname_expid_strings)

    #    clims = (minimum(fname_expids), maximum(fname_expids))
    clims = (1, n_fnames)

    fig = Figure(size = (1200, 400 * n_lin_coeffs), fontsize = 22)

    for pind in 1:n_lin_coeffs
        yvals = @view linParams[:, pind]
        summary = nanzeropercentile(yvals, percent_vec = [16, 50, 84])
        n_sigma = 5
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
                title = "Tele: $(parg["tele"]), MJD: $(unique_mjds[mjd_ind]), Wave. Soln. Linear Params")
        elseif pind == n_lin_coeffs
            ax = Axis(fig[pind, 1],
                limits = limits,
                xlabel = "FIBERID",
                ylabel = "Linear Parameter $(pind)")
        else
            ax = Axis(fig[pind, 1],
                limits = limits,
                ylabel = "Linear Parameter $(pind)")
        end

        scatter!(ax, 301 .- fiber_inds, linParams[:, pind])
    end

    sky_wave_linParams_Path = get_save_dir(unique_mjds[mjd_ind]) *
                              "nightFPIwave_linParams_$(parg["tele"])_$(unique_mjds[mjd_ind]).png"
    save(sky_wave_linParams_Path, fig)

    fig = Figure(size = (1200, 400 * n_nl_coeffs), fontsize = 22)

    for pind in 1:n_nl_coeffs
        yvals = @view nlParams[:, pind]
        summary = nanzeropercentile(yvals, percent_vec = [16, 50, 84])
        n_sigma = 5
        ylim = (summary[2] - n_sigma * (summary[2] - summary[1]),
            summary[2] + n_sigma * (summary[3] - summary[2]))
        outside_limits = (yvals .< ylim[1]) .| (yvals .> ylim[2])
        if sum(outside_limits) == 0
            ylim = nothing
        end
        limits = (nothing, ylim)

        if pind == 1
            ax = Axis(fig[pind, 1],
                ylabel = "Non-Linear Parameter $(pind)",
                limits = limits,
                title = "Tele: $(parg["tele"]), MJD: $(unique_mjds[mjd_ind]), Wave. Soln. Non-Linear Params")
        elseif pind == n_nl_coeffs
            ax = Axis(fig[pind, 1],
                xlabel = "FIBERID",
                limits = limits,
                ylabel = "Non-Linear Parameter $(pind)")
        else
            ax = Axis(fig[pind, 1],
                limits = limits,
                ylabel = "Non-Linear Parameter $(pind)")
        end

        scatter!(ax, 301 .- fiber_inds, nlParams[:, pind])
    end
    sky_wave_nlParams_Path = get_save_dir(unique_mjds[mjd_ind]) *
                             "nightFPIwave_nlParams_$(parg["tele"])_$(unique_mjds[mjd_ind]).png"
    save(sky_wave_nlParams_Path, fig)

    fig = Figure(size = (1200, 400 * n_dither_coeffs), fontsize = 22)

    for pind in 1:n_dither_coeffs
        yvals = @view ditherParams[:, pind]
        summary = nanzeropercentile(yvals, percent_vec = [16, 50, 84])
        n_sigma = 5
        ylim = (summary[2] - n_sigma * (summary[2] - summary[1]),
            summary[2] + n_sigma * (summary[3] - summary[2]))
        outside_limits = (yvals .< ylim[1]) .| (yvals .> ylim[2])
        if sum(outside_limits) == 0
            ylim = nothing
        end
        limits = (nothing, ylim)

        if pind == 1
            ax = Axis(fig[pind, 1],
                ylabel = "Dither Parameter $(pind)",
                limits = limits,
                title = "Tele: $(parg["tele"]), MJD: $(unique_mjds[mjd_ind]), Wave. Soln. Dither Params")
        elseif pind == n_dither_coeffs
            ax = Axis(fig[pind, 1],
                limits = limits,
                xlabel = "FIBERID",
                ylabel = "Dither Parameter $(pind)")
        else
            ax = Axis(fig[pind, 1],
                limits = limits,
                ylabel = "Dither Parameter $(pind)")
        end

        scatter!(ax, 301 .- fiber_inds, ditherParams[:, pind])
    end
    sky_wave_ditherParams_Path = get_save_dir(unique_mjds[mjd_ind]) *
                                 "nightFPIwave_ditherParams_$(parg["tele"])_$(unique_mjds[mjd_ind]).png"
    save(sky_wave_ditherParams_Path, fig)

    n_fnames = maximum(resid_exp_ints)
    fpi_resid_per_exp_per_fiber = zeros(Float64, (n_fnames, N_FIBERS))
    y_vals = 301 .- fiber_inds
    for fname_ind in 1:n_fnames
        fig = Figure(size = (1200, 1200), fontsize = 22)
        x, y, z = Vector{Float64}(), Vector{Float64}(), Vector{Float64}()
        for fibIndx in 1:N_FIBERS
            in_exp = findall(resid_exp_ints[:, fibIndx] .== fname_ind)
            msk = resid_used_in_fit[in_exp, fibIndx]
            #            x = vcat(x,resid_peak_ints[in_exp[msk],fibIndx])
            x = vcat(x, resid_xt[fibIndx, in_exp[msk]])
            y = vcat(y, y_vals[fibIndx] .* ones(sum(msk)))
            z = vcat(z, resid_vec[fibIndx, in_exp[msk]] .* 1000)
            curr_resid_summary = nanzeropercentile(
                resid_vec[fibIndx, in_exp[msk]] .* 1000, percent_vec = [16, 84])
            fpi_resid_per_exp_per_fiber[fname_ind, fibIndx] = 0.5 * (curr_resid_summary[2] -
                                                               curr_resid_summary[1])
        end

        resid_summary = nanzeropercentile(z, percent_vec = [16, 84])
        resid_sigma = 0.5 * (resid_summary[2] - resid_summary[1])
        ax = Axis(fig[1, 1],
            ylabel = "FIBERID",
            limits = ((-1.8, 1.8), (-10 + 1, 300 + 10)),
            xlabel = "Transformed X Peak Position",
            title = "FPI Peak Residuals, σ = $(round(resid_sigma; digits = 3)) mÅ\nTele: $(parg["tele"]), MJD: $(unique_mjds[mjd_ind]), ExpID: $(fname_expid_strings[fname_ind])")

        clim = (-10.0, 10.0)
        z .= clamp.(z, clim[1], clim[2])
        cmap = :diverging_bkr_55_10_c35_n256
        scatter!(ax, x, y, color = z,
            markersize = 3, colorrange = clim,
            colormap = cmap)
        Colorbar(fig[1, 2], colorrange = clim, colormap = cmap,
            width = 20, height = Relative(1.0), label = "data - model (mÅ)")
        fpiPeakResiduals_Path = get_save_dir(unique_mjds[mjd_ind]) *
                                "fpiPeakResiduals_$(parg["tele"])_$(unique_mjds[mjd_ind])_$(fname_expid_strings[fname_ind]).png"
        save(fpiPeakResiduals_Path, fig)
    end

    fig = Figure(size = (1200, 400), fontsize = 22)
    clims = (1, n_fnames)
    n_sigma = 8
    summary = nanzeropercentile(
        fpi_resid_per_exp_per_fiber, percent_vec = [16, 50, 84], dims = (1, 2))
    med_val = summary[2]
    ylim = (summary[2] - n_sigma * (summary[2] - summary[1]),
        summary[2] + n_sigma * (summary[3] - summary[2]))
    outside_limits = (fpi_resid_per_exp_per_fiber .< ylim[1]) .|
                     (fpi_resid_per_exp_per_fiber .> ylim[2])
    if sum(outside_limits) == 0
        ylim = nothing
    end
    limits = (nothing, ylim)
    ax = Axis(fig[1, 1],
        xlabel = "FIBERID",
        ylabel = "Residual Scatter (mÅ)",
        limits = limits,
        title = "FPI Peak Residual Scatter Per Fiber\nTele: $(parg["tele"]), MJD: $(unique_mjds[mjd_ind])")

    for fname_ind in 1:n_fnames
        scatter!(ax, 301 .- fiber_inds, fpi_resid_per_exp_per_fiber[fname_ind, :],
            color = fname_ind * ones(N_FIBERS), colorrange = clims)
    end
    hlines!(ax, med_val, linestyle = :dash)

    fpiResidual_scatter_Path = get_save_dir(unique_mjds[mjd_ind]) *
                               "fpiPeakResidualScatterPerFiber_$(parg["tele"])_$(unique_mjds[mjd_ind]).png"
    save(fpiResidual_scatter_Path, fig)

    vel_mult = 3e8 / 1000 / 16000
    fig = Figure(size = (1200, 400), fontsize = 22)
    n_sigma = 8
    summary = nanzeropercentile(
        fpi_resid_per_exp_per_fiber .* vel_mult, percent_vec = [16, 50, 84], dims = (1, 2))
    med_val = summary[2]
    ylim = (summary[2] - n_sigma * (summary[2] - summary[1]),
        summary[2] + n_sigma * (summary[3] - summary[2]))
    outside_limits = (fpi_resid_per_exp_per_fiber .* vel_mult .< ylim[1]) .|
                     (fpi_resid_per_exp_per_fiber .* vel_mult .> ylim[2])
    if sum(outside_limits) == 0
        ylim = nothing
    end
    limits = (nothing, ylim)
    ax = Axis(fig[1, 1],
        xlabel = "FIBERID",
        ylabel = "Approx. Velocity Uncertainty (m/s)",
        limits = limits,
        title = "FPI Peak Residual Scatter Per Fiber\nTele: $(parg["tele"]), MJD: $(unique_mjds[mjd_ind])")

    for fname_ind in 1:n_fnames
        scatter!(ax, 301 .- fiber_inds, fpi_resid_per_exp_per_fiber[fname_ind, :] .* vel_mult,
            color = fname_ind * ones(N_FIBERS), colorrange = clims)
    end
    hlines!(ax, med_val, linestyle = :dash)

    fpiResidual_scatter_Path = get_save_dir(unique_mjds[mjd_ind]) *
                               "fpiPeakResidualVelErrPerFiber_$(parg["tele"])_$(unique_mjds[mjd_ind]).png"
    save(fpiResidual_scatter_Path, fig)
end

list2Dexp = []
for mjd in unique_mjds
    df = read_almanac_exp_df(parg["outdir"] * "almanac/$(parg["runname"]).h5", parg["tele"], mjd)
    df.cartidInt = parseCartID.(df.cartid)
    df.exposure_int = if typeof(df.exposure) <: Array{Int}
        df.exposure
    else
        parse.(Int, df.exposure)
    end
    df.exposure_str = if typeof(df.exposure) <: Array{String}
        df.exposure
    else
        lpad.(string.(df.exposure), 8, "0")
    end

    function get_2d_name_partial(expid)
        parg["outdir"] * "/apred/$(mjd)/" *
        replace(get_1d_name(expid, df), "ar1D" => "ar2Dresidualscal") * ".h5"
    end
    local2D = get_2d_name_partial.(expid_list)
    push!(list2Dexp, local2D)
end
all2Da = vcat(list2Dexp...)

# per chip example 2D flux residuals from 1D extraction
for chip in CHIP_LIST
    println("\nGenerating 2D residual plots for chip $chip")
    all2D = replace.(all2Da, "_$(FIRST_CHIP)_" => "_$(chip)_")

    rng = MersenneTwister(351 + unique_mjds[1])

    # TODO parallelize plotting
    # we should customize this to the exposures we want to see and types of stars we want
    nsamp = minimum([length(all2D), 10])
    sample_exposures = sample(rng, all2D, nsamp, replace = false)
    f = h5open(parg["outdir"] * "almanac/$(parg["runname"]).h5")
    for exp_fname in sample_exposures
        sname = split(split(split(exp_fname, "/")[end], ".h5")[1], "_")
        fnameType, tele, mjd, expnum, chiploc, exptype = sname[(end - 5):end]

        flux_2d = load(exp_fname, "resid_flux")
        ivar_2d = load(exp_fname, "resid_ivar")
        resid_zscore = flux_2d .* (ivar_2d .^ 0.5)

        fig = Figure(size = (840, 800), fontsize = 24)
        ax = Axis(fig[1, 1])
        hm = heatmap!(ax, resid_zscore,
            colormap = :diverging_bkr_55_10_c35_n256,
            colorrange = (-5, 5),
            interpolate = false
        )

        text!(ax,
            0.5, 1.05,
            text = "Zscore Residuals\nTele: $(tele), MJD: $(mjd), ExpNum: $(expnum), Chip: $(chiploc)",
            align = (:center, :bottom),
            space = :relative
        )

        Colorbar(fig[1, 2], hm, width = 20, height = Relative(1.0))
        colgap!(fig.layout, 1, 20)  # Set spacing between image 1 and colorbar 1
        data_aspect = diff(hm[1][])[1] / (diff(hm[2][])[1])
        colsize!(fig.layout, 1, Aspect(1, data_aspect))
        resize_to_layout!(fig)

        savePath = get_save_dir(mjd) *
                   "ar2Dresidualscal_$(tele)_$(mjd)_$(expnum)_$(chiploc)_$(exptype).png"
        save(savePath, fig, px_per_unit = 3)
    end
end

all1Da = String[] # all 1D files for chip a
for mjd in unique_mjds
    df = read_almanac_exp_df(parg["outdir"] * "almanac/$(parg["runname"]).h5", parg["tele"], mjd)
    function get_1d_name_partial(expid)
        parg["outdir"] * "apred/$(mjd)/" * get_1d_name(expid, df) * ".h5"
    end

    file_list = get_1d_name_partial.(expid_list)
    append!(all1Da, file_list)
end

allExptype = convert.(String, map(x -> split(split(split(x, "/")[end], ".")[1], "_")[end], all1Da))

# Define custom sorting order for exposure types
function get_exptype_priority(exptype::AbstractString)
    priority_map = Dict(
        "OBJECT" => 1,
        "DOMEFLAT" => 2,
        "INTERNALFLAT" => 3,
        "QUARTZFLAT" => 4,
        "DARK" => 5,
        "ARCLAMP" => 6
    )
    return get(priority_map, exptype, 999)  # Unknown types get high number
end

# Sort allExptype using the custom priority function
sorted_exptypes = sort(unique(allExptype), by = get_exptype_priority)

# per chip example spectra
for chip in string.(collect(parg["chips"]))
    println("\nGenerating example spectra plots for chip $chip")
    all1D = replace.(all1Da, "_$(FIRST_CHIP)_" => "_$(chip)_")

    rng = MersenneTwister(351 + unique_mjds[1])

    # TODO parallelize plotting
    # we should customize this to the types of objects we want
    # for example, we want to be looking at the tellurics and pure sky
    for exptype2plot in sorted_exptypes
        msk_exptype = allExptype .== exptype2plot
        if any(msk_exptype)
            nsamp = minimum([count(msk_exptype), 3])
            sample_exposures = sample(rng, all1D[msk_exptype], nsamp, replace = false)

            f = h5open(parg["outdir"] * "almanac/$(parg["runname"]).h5")
            for exp_fname in sample_exposures
                sname = split(split(split(exp_fname, "/")[end], ".h5")[1], "_")
                fnameType, tele, mjd, expnum, chiploc, exptype = sname[(end - 5):end]

                flux_1d = load(exp_fname, "flux_1d")
                mask_1d = load(exp_fname, "mask_1d")
                msk_loc = (mask_1d .& bad_pix_bits .== 0)
                # msk_loc = (mask_1d .&
                #         (bad_pix_bits + bad_1d_failed_extract + bad_1d_no_good_pix + bad_1d_neff) .== 0)

                fibtargDict = get_fibTargDict(f, tele, mjd, expnum)
                sample_fibers = sample(rng, 1:300, 3, replace = false)
                for fib in sample_fibers
                    fibID = fiberIndx2fiberID(fib)
                    fibType = fibtargDict[fibID]
                    dat = flux_1d[:, fib]
                    mskt = msk_loc[:, fib]
                    dat = nanify(flux_1d[mskt, fib], mskt)
                    datbad = nanify(flux_1d[.!mskt, fib], .!mskt)

                    fig = Figure(size = (800, 800))
                    ax1 = Axis(fig[1, 1])
                    lines!(ax1, 1:2048, dat, color = :dodgerblue)
                    scatter!(ax1, 1:2048, datbad, color = :red, markersize = 2)
                    xlims!(ax1, 0, 2049)
                    ax1.ylabel = "ADU"

                    ax2 = Axis(fig[2, 1])
                    lines!(ax2, 1:2048, dat, color = :dodgerblue)
                    scatter!(ax2, 1:2048, datbad, color = :red, markersize = 2)
                    ymaxt = 2 * nanzeromedian(dat)
                    ymaxv = if .!isnan(ymaxt)
                        ymaxt
                    else
                        1
                    end
                    ylims!(ax2, 0, ymaxv)
                    xlims!(ax2, 0, 2049)
                    ax2.xlabel = "Pixel Index"
                    ax2.ylabel = "ADU"

                    savePath = get_save_dir(mjd) *
                               "ar1D_$(tele)_$(mjd)_$(expnum)_$(chiploc)_$(fib)_$(fibType)_$(exptype).png"
                    save(savePath, fig)
                end
            end
        end
    end
end

# 1D reinterpolated spectra examples

println("\nGenerating reinterpolated spectra examples")

function plot_1d_uni(
        fib, fibtargDict, outflux, outmsk, bname, tele, mjd, expnum, exptype, expuni_fname)
    fibID = fiberIndx2fiberID(fib)
    fibType = fibtargDict[fibID]

    fig = Figure(size = (600, 1200), fontsize = 12)

    # Top panel showing full wavelength range
    ax1 = Axis(fig[1, 1:2], title = "Full Spectrum", ylabel = "Flux")
    # hidexdecorations!(ax1, grid=false)
    nan_flux = copy(outflux[:, fib])
    nan_flux[isnanorzero.(nan_flux)] .= NaN
    # TODO: want to add red coloring for the pixels we are dropping because of masking
    nan_flux[outmsk[:, fib] .== 0] .= NaN

    # Get full wavelength limits and add padding
    full_xlims = if all(isnan, nan_flux)
        (minimum(logUniWaveAPOGEE), maximum(logUniWaveAPOGEE))
    else
        (minimum(logUniWaveAPOGEE[.!isnan.(nan_flux)]) - 2,
            maximum(logUniWaveAPOGEE[.!isnan.(nan_flux)]) + 2)
    end
    full_mask = full_xlims[1] .<= logUniWaveAPOGEE .<= full_xlims[2]
    filtered_dat = filter(!isnan, nan_flux[full_mask])
    full_ylims = if isempty(filtered_dat)
        (0, 1)
    else
        extrema(filtered_dat)
    end
    full_yrange = full_ylims[2] - full_ylims[1]
    full_ylims = (full_ylims[1] - 0.05 * full_yrange, full_ylims[2] + 0.05 * full_yrange)

    limits!(ax1, full_xlims..., full_ylims...)
    lines!(ax1, logUniWaveAPOGEE, nan_flux, color = :dodgerblue)

    # Middle panel showing zoomed region
    ax2 = Axis(fig[2, 1:2], title = "Strong DIB Region", ylabel = "Flux")
    # hidexdecorations!(ax2, grid=false)
    xlims = (15220, 15340)
    mask = xlims[1] .<= logUniWaveAPOGEE .<= xlims[2]
    filtered_dat = filter(!isnan, nan_flux[mask])
    ylims = if isempty(filtered_dat)
        (0, 1)
    else
        extrema(filtered_dat)
    end
    yrange = ylims[2] - ylims[1]
    ylims = (ylims[1] - 0.05 * yrange, ylims[2] + 0.05 * yrange)

    limits!(ax2, xlims..., ylims...)
    lines!(ax2, logUniWaveAPOGEE, nan_flux, color = :dodgerblue)

    # Bottom left panel showing left edge
    ax3 = Axis(fig[3, 1], title = "Blue Edge", ylabel = "Flux")
    # hidexdecorations!(ax3, grid=false)
    left_xlims = (minimum(full_xlims), minimum(full_xlims) + 100)
    left_mask = left_xlims[1] .<= logUniWaveAPOGEE .<= left_xlims[2]
    filtered_dat = filter(!isnan, nan_flux[left_mask])
    left_ylims = if isempty(filtered_dat)
        (0, 1)
    else
        extrema(filtered_dat)
    end
    left_yrange = left_ylims[2] - left_ylims[1]
    left_ylims = (left_ylims[1] - 0.05 * left_yrange, left_ylims[2] + 0.05 * left_yrange)

    limits!(ax3, left_xlims..., left_ylims...)
    lines!(ax3, logUniWaveAPOGEE, nan_flux, color = :dodgerblue)

    # Bottom right panel showing right edge
    ax4 = Axis(fig[3, 2], title = "Red Edge")
    hideydecorations!(ax4, grid = false)
    # hidexdecorations!(ax4, grid=false)
    right_xlims = (maximum(full_xlims) - 100, maximum(full_xlims))
    right_mask = right_xlims[1] .<= logUniWaveAPOGEE .<= right_xlims[2]
    filtered_dat = filter(!isnan, nan_flux[right_mask])
    right_ylims = if isempty(filtered_dat)
        (0, 1)
    else
        extrema(filtered_dat)
    end
    right_yrange = right_ylims[2] - right_ylims[1]
    right_ylims = (right_ylims[1] - 0.05 * right_yrange, right_ylims[2] + 0.05 * right_yrange)

    limits!(ax4, right_xlims..., right_ylims...)
    lines!(ax4, logUniWaveAPOGEE, nan_flux, color = :dodgerblue)

    # Bottom row showing chip boundaries
    ax5 = Axis(
        fig[4, 1], title = "Chip Blue-Green Edge", xlabel = "Wavelength (Å)", ylabel = "Flux")
    chip_b_xlims = (15834 - 70, 15834 + 70)
    chip_b_mask = chip_b_xlims[1] .<= logUniWaveAPOGEE .<= chip_b_xlims[2]
    filtered_dat = filter(!isnan, nan_flux[chip_b_mask])
    chip_b_ylims = if isempty(filtered_dat)
        (0, 1)
    else
        extrema(filtered_dat)
    end
    chip_b_yrange = chip_b_ylims[2] - chip_b_ylims[1]
    chip_b_ylims = (
        chip_b_ylims[1] - 0.05 * chip_b_yrange, chip_b_ylims[2] + 0.05 * chip_b_yrange)

    limits!(ax5, chip_b_xlims..., chip_b_ylims...)
    lines!(ax5, logUniWaveAPOGEE, nan_flux, color = :dodgerblue)

    ax6 = Axis(fig[4, 2], title = "Chip Green-Red Edge", xlabel = "Wavelength (Å)")
    hideydecorations!(ax6, grid = false)
    chip_a_xlims = (16454 - 70, 16454 + 70)
    chip_a_mask = chip_a_xlims[1] .<= logUniWaveAPOGEE .<= chip_a_xlims[2]
    filtered_dat = filter(!isnan, nan_flux[chip_a_mask])
    chip_a_ylims = if isempty(filtered_dat)
        (0, 1)
    else
        extrema(filtered_dat)
    end
    chip_a_yrange = chip_a_ylims[2] - chip_a_ylims[1]
    chip_a_ylims = (
        chip_a_ylims[1] - 0.05 * chip_a_yrange, chip_a_ylims[2] + 0.05 * chip_a_yrange)

    limits!(ax6, chip_a_xlims..., chip_a_ylims...)
    lines!(ax6, logUniWaveAPOGEE, nan_flux, color = :dodgerblue)

    resize_to_layout!(fig)

    savePath = get_save_dir(mjd) *
               "$(bname)_$(tele)_$(mjd)_$(expnum)_$(fib)_$(fibType)_$(exptype).png"
    save(savePath, fig)
    return
end

rng = MersenneTwister(536 + unique_mjds[1])

# TODO: parallelize plotting
for exptype2plot in sorted_exptypes
    msk_exptype = allExptype .== exptype2plot
    if any(msk_exptype)
        nsamp = minimum([count(msk_exptype), 3])
        sample_exposures = sample(rng, all1Da[msk_exptype], nsamp, replace = false)
        f = h5open(parg["outdir"] * "almanac/$(parg["runname"]).h5")
        for exp_fname in sample_exposures
            sname = split(split(split(exp_fname, "/")[end], ".h5")[1], "_")
            fnameType, tele, mjd, expnum, chiploc, exptype = sname[(end - 5):end]
            expuni_fname = replace(
                replace(exp_fname, "ar1D" => "ar1Duni"), "_$(FIRST_CHIP)_" => "_")
            outflux = load(expuni_fname, "flux_1d")
            outmsk = load(expuni_fname, "mask_1d")
            expunical_fname = replace(
                replace(exp_fname, "ar1D" => "ar1Dunical"), "_$(FIRST_CHIP)_" => "_")
            outfluxcal = load(expunical_fname, "flux_1d")
            outmskcal = load(expunical_fname, "mask_1d")
            # need to switch this back when the masking is updated
            # msk_loc = (outmsk .& bad_pix_bits .== 0)

            fibtargDict = get_fibTargDict(f, tele, mjd, expnum)
            sample_fibers = sample(rng, 1:300, 3, replace = false)
            for fib in sample_fibers
                plot_1d_uni(fib, fibtargDict, outflux, outmsk, "ar1Duni",
                    tele, mjd, expnum, exptype, expuni_fname)
                plot_1d_uni(fib, fibtargDict, outfluxcal, outmskcal, "ar1Dunical",
                    tele, mjd, expnum, exptype, expunical_fname)
            end
        end
    end
end
