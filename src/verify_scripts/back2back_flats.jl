using JLD2, ProgressMeter, ArgParse, SlackThreads, Glob, StatsBase, BinnedStatistics

src_dir = "../"
include(src_dir * "/fileNameHandling.jl")
include(src_dir * "/utils.jl")
include(src_dir * "/makie_plotutils.jl")

## TODO Make sure this function is extended to check BOTH ar2D and ar2Dcal

## Parse command line arguments
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--tele"
        required = true
        help = "telescope name (apo or lco)"
        arg_type = String
        default = ""
        "--mjd"
        required = true
        help = "start mjd"
        arg_type = Int
        default = 0
        "--expid-start"
        required = true
        help = "start expid"
        arg_type = Int
        default = 0
        "--expid-end"
        required = true
        help = "end expid"
        arg_type = Int
        default = 0
        "--data_dir"
        required = true
        help = "base directory where 2D extractions of flats are stored"
        arg_type = String
        default = ""
        "--verify_dir"
        required = false
        help = "base directory where to store verification plots (defaults to data_dir if not specified)"
        arg_type = String
        default = ""
    end
    return parse_args(s)
end

parg = parse_commandline()

# Set verify_dir to data_dir if not provided
if parg["verify_dir"] == ""
    parg["verify_dir"] = parg["data_dir"]
end

thread = SlackThread();
thread("Back2Back flats ivar test for Tele: $(parg["tele"]), MJD: $(parg["mjd"])")

dirNamePlots = joinpath(parg["verify_dir"], "plots/")
if !ispath(dirNamePlots)
    mkpath(dirNamePlots)
end

function get_IQR(x, y; cnts_cut = 100, nbin_med = 100)
    edges, centers,
    binStat = binnedStatistic(
        y,
        x,
        nbins = nbin_med,
        statistic = :f,
        f = nanzeroiqr
    )

    edges, centers, binStat1 = binnedStatistic(
        y,
        x,
        nbins = nbin_med,
        statistic = :median
    )

    edges, centers, binCnts = binnedStatistic(
        y,
        ones(size(y)),
        nbins = nbin_med,
        statistic = :sum
    )

    med_line = binStat1
    med_line[binCnts .< cnts_cut] .= NaN
    iqr_line = binStat
    iqr_line[binCnts .< cnts_cut] .= NaN
    return centers, iqr_line, med_line
end

# why not also ar2Dcal?
# sorting might be reverse compared to normal because of abc to RGB (come back to this next b2b run)
flist = glob("ar2D_*.h5", joinpath(parg["data_dir"], "apred/$(parg["mjd"])/"));
expnum_v = map(x -> parse(Int, split(basename(x), "_")[end - 2]), flist)
chip_v = map(x -> split(basename(x), "_")[end - 1], flist);
mskCal = (parg["expid-start"] % 10000) .<= expnum_v .<= (parg["expid-end"] % 10000);

p = sortperm(expnum_v[mskCal])

@showprogress for indoff in 1:(length(flist[mskCal][p]) ÷ 6 - 1)
    fig = Figure(size = (2000, 800), fontsize = 32)

    expid1 = expnum_v[mskCal][p][1 + 6 * indoff]
    for (ind, chipn) in enumerate(CHIP_LIST)
        ax = Axis(fig[1, ind], title = "Tele: $(parg["tele"]), Chip: $chipn",
            xlabel = "Z-Score", ylabel = "Log Average Flux")

        dimage1 = load(flist[mskCal][p][ind + 6 * indoff], "dimage")
        ivarimage1 = load(flist[mskCal][p][ind + 6 * indoff], "ivarimage")
        dimage2 = load(flist[mskCal][p][ind + 3 + 6 * indoff], "dimage")
        ivarimage2 = load(flist[mskCal][p][ind + 3 + 6 * indoff], "ivarimage")

        datdiff = (dimage1 .- dimage2)
        avgFlux = (dimage1 .+ dimage2) ./ 2
        logavgFlux = log10n.(avgFlux)
        differr = sqrt.(1 ./ ivarimage1 .+ 1 ./ ivarimage2)
        zmat = datdiff ./ differr

        msk = (abs.(zmat) .< 15) .& (-3 .<= logavgFlux .<= 4)

        h = hexbin!(ax, zmat[msk], logavgFlux[msk],
            bins = 500,
            colormap = :linear_bmy_10_95_c71_n256,
            colorscale = log10
        )

        if count(msk) > 0
            centers, iqr_line, med_line = get_IQR(zmat[msk], logavgFlux[msk], cnts_cut = 300)

            lines!(ax, med_line, centers, color = :white, linewidth = 2)
            lines!(ax, iqr_line, centers, color = :limegreen, linewidth = 2)
        end

        vlines!(ax, 0, color = :white, linewidth = 2, linestyle = :dash)
        vlines!(ax, -1, color = :limegreen, linewidth = 2, linestyle = :dash)
        vlines!(ax, 1, color = :limegreen, linewidth = 2, linestyle = :dash)

        xlims!(ax, -10, 10)
        ylims!(ax, -3, 4)
    end

    resize_to_layout!(fig)
    # need a better naming scheme that keeps the exposure ids
    plot_name = "back2backFlat_$(parg["tele"])_$(parg["mjd"])_$(expid1).png"
    save(joinpath(dirNamePlots, plot_name), fig, px_per_unit = 3)
    thread("back2backFlat_$(parg["tele"])_$(parg["mjd"])_$(expid1)", joinpath(dirNamePlots, plot_name))
end
