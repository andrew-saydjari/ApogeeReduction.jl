# This is a script that makes the plots for a nightly processing and posts them to the slack channel.
# called by run_all.sh
using JLD2, ProgressMeter, ArgParse, SlackThreads, Glob, StatsBase, Random, HDF5, DataFrames

src_dir = "../"
include(src_dir * "/fileNameHandling.jl")
include(src_dir * "/utils.jl")
include(src_dir * "/makie_plotutils.jl")
include(src_dir * "/ar1D.jl")

## TODO add sky flux relFlux consistency check plots

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

dirNamePlots = parg["outdir"] * "plots/"
mkpath(dirNamePlots) # will work even if it already exists

unique_mjds = if parg["runlist"] != ""
    subDic = load(parg["runlist"])
    unique(subDic["mjd"])
else
    [parg["mjd"]]
end

expid_list = if parg["runlist"] != ""
    subDic = load(parg["runlist"])
    subDic["expid"]
else
    [parg["expid"]]
end

sky_thread = SlackThread()
if length(unique_mjds) > 1
    min_mjd, max_mjd = extrema(unique_mjds)
    sky_thread("Here are the wavelength solution stability plots from $(parg["tele"]) for SJD $(min_mjd) to $(max_mjd)")
else
    sky_thread("Here are the wavelength solution stability plots from $(parg["tele"]) for SJD $(unique_mjds[1])")
end

for mjd_ind in 1:size(unique_mjds, 1)
    savePath = dirNamePlots * "skywave_linParams_$(parg["tele"])_$(unique_mjds[mjd_ind]).png"
    sky_thread(
        "Sky Line Wavelength Solution Linear Parameters: MJD $(unique_mjds[mjd_ind])", savePath)
    savePath = dirNamePlots * "skywave_nlParams_$(parg["tele"])_$(unique_mjds[mjd_ind]).png"
    sky_thread(
        "Sky Line Wavelength Solution Non-linear Parameters: MJD $(unique_mjds[mjd_ind])", savePath)
    savePath = dirNamePlots *
               "skywave_per_fiber_vs_pixel_$(parg["tele"])_$(unique_mjds[mjd_ind]).png"
    sky_thread("Sky Line Wavelength Solution per fiber: MJD $(unique_mjds[mjd_ind])", savePath)
    savePath = dirNamePlots *
               "skywave_per_pixel_vs_fiber_$(parg["tele"])_$(unique_mjds[mjd_ind]).png"
    sky_thread("Sky Line Wavelength Solution per pixel: MJD $(unique_mjds[mjd_ind])", savePath)
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
for chip in ["a", "b", "c"]
    all2D = replace.(all2Da, "_a_" => "_$(chip)_")

    thread = SlackThread()
    if length(unique_mjds) > 1
        min_mjd, max_mjd = extrema(unique_mjds)
        thread("Here are some example 2D residuals after 1D extraction on chip $(chip) from $(parg["tele"]) for SJD $(min_mjd) to $(max_mjd)")
    else
        thread("Here are some example 2D residuals after 1D extraction on chip $(chip) from $(parg["tele"]) for SJD $(unique_mjds[1])")
    end
    rng = MersenneTwister(351 + unique_mjds[1])

    # TODO parallelize plotting
    # we should customize this to the exposures we want to see and types of stars we want
    nsamp = minimum([length(all2D), 10])
    sample_exposures = sample(rng, all2D, nsamp, replace = false)
    f = h5open(parg["outdir"] * "almanac/$(parg["runname"]).h5")
    for exp_fname in sample_exposures
        sname = split(split(split(exp_fname, "/")[end],".h5")[1], "_")
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

        savePath = dirNamePlots *
                   "ar2Dresidualscal_$(tele)_$(mjd)_$(expnum)_$(chiploc)_$(exptype).png"
        save(savePath, fig, px_per_unit = 3)
        thread("$(exp_fname)", savePath)
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
    all1D = replace.(all1Da, "_a_" => "_$(chip)_")

    thread = SlackThread()
    if length(unique_mjds) > 1
        min_mjd, max_mjd = extrema(unique_mjds)
        thread("Here are some example spectra on chip $(chip) from $(parg["tele"]) for SJD $(min_mjd) to $(max_mjd)")
    else
        thread("Here are some example spectra on chip $(chip) from $(parg["tele"]) for SJD $(unique_mjds[1])")
    end
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
                sname = split(split(split(exp_fname, "/")[end],".h5")[1], "_")
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

                    savePath = dirNamePlots *
                               "ar1D_$(tele)_$(mjd)_$(expnum)_$(chiploc)_$(fib)_$(fibType)_$(exptype).png"
                    save(savePath, fig)

                    thread("Fiberindex: $(fib) $(fibType), $(exp_fname)", savePath)
                end
            end
        end
    end
end

# 1D reinterpolated spectra examples

thread = SlackThread();
if length(unique_mjds) > 1
    min_mjd, max_mjd = extrema(unique_mjds)
    thread("Here are some example reinterpolated spectra from $(parg["tele"]) for SJD $(min_mjd) to $(max_mjd)")
else
    thread("Here are some example reinterpolated spectra from $(parg["tele"]) for SJD $(unique_mjds[1])")
end

function plot_1d_uni(fib, fibtargDict, outflux, outmsk, thread, bname,
        tele, mjd, expnum, chiploc, exptype, expuni_fname)
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

    savePath = dirNamePlots *
               "$(bname)_$(tele)_$(mjd)_$(expnum)_$(chiploc)_$(fib)_$(fibType)_$(exptype).png"
    save(savePath, fig)

    thread("Fiberindex: $(fib) $(fibType), $(expuni_fname)", savePath)
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
            sname = split(split(split(exp_fname, "/")[end],".h5")[1], "_")
            fnameType, tele, mjd, expnum, chiploc, exptype = sname[(end - 5):end]
            expuni_fname = replace(replace(exp_fname, "ar1D" => "ar1Duni"), "_a_" => "_")
            outflux = load(expuni_fname, "flux_1d")
            outmsk = load(expuni_fname, "mask_1d")
            expunical_fname = replace(replace(exp_fname, "ar1D" => "ar1Dunical"), "_a_" => "_")
            outfluxcal = load(expunical_fname, "flux_1d")
            outmskcal = load(expunical_fname, "mask_1d")
            # need to switch this back when the masking is updated
            # msk_loc = (outmsk .& bad_pix_bits .== 0)

            fibtargDict = get_fibTargDict(f, tele, mjd, expnum)
            sample_fibers = sample(rng, 1:300, 3, replace = false)
            for fib in sample_fibers
                plot_1d_uni(fib, fibtargDict, outflux, outmsk, thread, "ar1Duni",
                    tele, mjd, expnum, chiploc, exptype, expuni_fname)
                plot_1d_uni(fib, fibtargDict, outfluxcal, outmskcal, thread, "ar1Dunical",
                    tele, mjd, expnum, chiploc, exptype, expunical_fname)
            end
        end
    end
end
