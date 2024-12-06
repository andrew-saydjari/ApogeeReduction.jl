# This is a script that makes the plots for a nightly processing and posts them to the slack channel.
using JLD2, ProgressMeter, ArgParse, SlackThreads, Glob, StatsBase, Random, HDF5, DataFrames

src_dir = "../"
include(src_dir * "/fileNameHandling.jl")
include(src_dir * "/utils.jl")
include(src_dir * "/makie_plotutils.jl")
include(src_dir * "/ap1D.jl")

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

all1Da = String[] # all 1D files for chip a
for mjd in unique_mjds
    df = h5open(parg["outdir"] * "almanac/$(parg["runname"]).h5") do f
        df = DataFrame(read(f["$(parg["tele"])/$(mjd)/exposures"]))
    end
    function get_1d_name_partial(expid)
        parg["outdir"] * "apred/$(mjd)/" * get_1d_name(expid, df) * ".jld2"
    end

    file_list = get_1d_name_partial.(expid_list)
    append!(all1Da, file_list)
end

all1Dperchip = []
for chip in string.(collect(parg["chips"]))
    all1Dchip = replace.(all1Da, "_a_" => "_$(chip)_")
    push!(all1Dperchip, all1Dchip)
end
all1D = vcat(all1Dperchip...)

thread = SlackThread();
if length(unique_mjds) > 1
    min_mjd, max_mjd = extrema(unique_mjds)
    thread("Here are some example spectra on chip(s) $(parg["chips"]) from $(parg["tele"]) for SJD $(min_mjd) to $(max_mjd)")
else
    thread("Here are some example spectra on chip(s) $(parg["chips"]) from $(parg["tele"]) for SJD $(unique_mjds[1])")
end
rng = MersenneTwister(351 + unique_mjds[1])

# TODO parallelize plotting
# we should customize this to the exposures we want to see and types of stars we want
# for example, we want to be looking at the tellurics and pure sky
nsamp = minimum([length(all1D), 5])
sample_exposures = sample(rng, all1D, nsamp, replace = false)
f = h5open(parg["outdir"] * "almanac/$(parg["runname"]).h5")
for exp_fname in sample_exposures
    sname = split(exp_fname, "_")
    tele, mjd, chiploc, expid = sname[(end - 4):(end - 1)]
    expid_num = parse(Int, last(expid, 4))
    flux_1d = load(exp_fname, "flux_1d")
    mask_1d = load(exp_fname, "mask_1d")
    msk_loc = (mask_1d .& bad_pix_bits .== 0)

    fibtargDict = get_fibTargDict(f, tele, parse(Int, mjd), expid_num)
    sample_fibers = sample(rng, 1:300, 5, replace = false)
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
        ylims!(ax2, 0, 2 * nanzeromedian(dat))
        xlims!(ax2, 0, 2049)
        ax2.xlabel = "Pixel Index"
        ax2.ylabel = "ADU"

        savePath = dirNamePlots *
                   "ap1D_$(tele)_$(mjd)_$(chiploc)_$(expid)_$(fib)_$(fibType).png"
        save(savePath, fig)

        thread("Fiberindex: $(fib) $(fibType), $(exp_fname)", savePath)
    end
end
