# This is a script that makes the plots for a nightly processing and posts them to the slack channel.
using JLD2, ProgressMeter, ArgParse, SlackThreads, Glob, StatsBase, Random, HDF5, DataFrames

src_dir = "../"
include(src_dir * "/fileNameHandling.jl")
include(src_dir * "/utils.jl")
include(src_dir * "/plotutils.jl")

## Parse command line arguments
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--chip"
        required = false
        help = "chip name (a, b, c)"
        arg_type = String
        "--tele"
        required = true
        help = "telescope name (apo or lco)"
        arg_type = String
        default = ""
        "--mjd"
        required = true
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
chip = parg["chip"]

thread = SlackThread();
thread("Here are some example spectra from $(parg["tele"]) $(chip) for $(parg["mjd"])")

dirNamePlots = parg["outdir"] * "plots/"
if !ispath(dirNamePlots)
    mkpath(dirNamePlots)
end

rng = MersenneTwister(351)

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

list1Dfiles = []
for mjd in unique_mjds
    df = h5open(parg["outdir"] * "almanac/$(parg["runname"]).h5") do f
        df = DataFrame(read(f["$(parg["tele"])/$(mjd)/exposures"]))
    end
    function get_1d_name_partial(expid)
        parg["outdir"] * "apred/$(mjd)/" * get_1d_name(expid, df) * ".jld2"
    end

    file_list = get_1d_name_partial.(expid_list)
    push!(list1Dfiles, file_list)
end
all1D = vcat(list1Dfiles...)

# we should customize this to the exposures we want to see and types of stars we want
# for example, we want to be looking at the tellurics and pure sky
nsamp = minimum([length(all1D), 5])
sample_exposures = sample(rng, all1D, nsamp, replace = false)
for exp_fname in sample_exposures
    sname = split(exp_fname, "_")
    tele, mjd, chiploc, expid = sname[(end - 4):(end - 1)]
    flux_1d = load(exp_fname, "flux_1d")

    sample_fibers = sample(rng, 1:300, 5, replace = false)
    for fib in sample_fibers
        dat = flux_1d[:, fib]
        fig = PythonPlot.figure(figsize = (8, 8), dpi = 150)
        ax = fig.add_subplot(2, 1, 1)
        ax.plot(1:2048, dat, color = "dodgerblue")
        ax.set_xlim(0, 2049)
        ax.set_ylabel("ADU")

        ax = fig.add_subplot(2, 1, 2)
        ax.plot(1:2048, dat, color = "dodgerblue")
        ax.set_ylim(0, 2 * median(dat))
        ax.set_xlim(0, 2049)
        ax.set_xlabel("Pixel Index")
        ax.set_ylabel("ADU")
        savePath = dirNamePlots * "ap1D_$(tele)_$(mjd)_$(chiploc)_$(expid)_$(fib).png"
        fig.savefig(savePath, bbox_inches = "tight", pad_inches = 0.1)
        PythonPlot.plotclose(fig)
        thread("Fiberindex: $(fib), $(exp_fname)", savePath)
    end
end
