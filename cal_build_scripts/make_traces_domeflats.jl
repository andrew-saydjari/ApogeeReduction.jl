using JLD2, ProgressMeter, ArgParse, SlackThreads, Glob, StatsBase, Einsum

src_dir = "./src"
include(src_dir * "/utils.jl")
include(src_dir * "/fileNameHandling.jl")
include(src_dir * "/../plotutils.jl")
include(src_dir * "/cal_build/traceExtract.jl")

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
        "--mjd-start"
        required = true
        help = "start mjd"
        arg_type = Int
        default = 0
        "--mjd-end"
        required = true
        help = "end mjd"
        arg_type = Int
        default = 0
        "--trace_dir"
        required = true
        help = "directory where 2D extractions of traces are stored"
        arg_type = String
        default = ""
        "--runlist"
        required = true
        help = "path name to hdf5 file with keys specifying list of exposures to run"
        arg_type = String
        default = ""
    end
    return parse_args(s)
end

# I am really in a pickle now. These are really really nightly things, not something
# for a range of MJD. I will need to change the way I am doing this.

parg = parse_commandline()
chip = parg["chip"]

# make summary plots and gifs and send to slack in outer loop
thread = SlackThread();
thread("DOMEFLAT TRACES for $(parg["tele"]) $(chip) from $(parg["mjd-start"]) to $(parg["mjd-end"])")

dirNamePlots = parg["trace_dir"] * "plots/"
if !ispath(dirNamePlots)
    mkpath(dirNamePlots)
end

# Load in the exact set of exposures
mjd = load(parg["runlist"], "mjd")
expid = load(parg["runlist"], "expid");
flist = get_cal_file.(
    Ref(parg["trace_dir"]), Ref(parg["tele"]), mjd, expid, Ref(chip), Ref("DOMEFLAT"))

fpifib1, fpifib2 = get_fpi_guide_fiberID(parg["tele"])

@showprogress for (indx, fname) in enumerate(flist)
    sname = split(fname, "_")
    teleloc, mjdloc, chiploc, expidloc = sname[(end - 4):(end - 1)]
    mjdfps2plate = get_fps_plate_divide(teleloc)
    f = jldopen(fname)
    image_data = f["dimage"][1:2048, 1:2048]
    ivar_image = f["ivarimage"][1:2048, 1:2048]
    close(f)

    trace_params = trace_extract(
        image_data, ivar_image, teleloc, mjdloc, chiploc, expidloc; image_mask = nothing)
    jldsave(
        parg["trace_dir"] *
        "dome_flats/domeTrace_$(teleloc)_$(mjdloc)_$(expidloc)_$(chiploc).jld2";
        trace_params = trace_params)

    cut = 750
    fig = PythonPlot.figure(figsize = (8, 8), dpi = 150)
    ax = fig.add_subplot(1, 1, 1)
    y = dropdims(nanzeromedian(trace_params[:, :, 1], 1), dims = 1)
    ax.scatter(301 .- (1:300), y)
    if parse(Int, mjdloc) > mjdfps2plate
        ax.scatter(fpifib1, y[301 - fpifib1], color = "red")
        ax.scatter(fpifib2, y[301 - fpifib2], color = "red")
    end
    ax.axhline(cut, linestyle = "--")

    ax.set_xlabel("FIBERID")
    ax.set_ylabel("Fit Height")

    plt.text(0.5,
        1.01,
        "Dome Flat Fit, Tele: $(parg["tele"]), MJD: $(mjdloc), Chip: $(chiploc), Expid: $(expidloc)",
        ha = "center",
        va = "bottom",
        transform = ax.transAxes)

    tracePlotPath = dirNamePlots *
                    "domeTrace_$(parg["tele"])_$(teleloc)_$(mjdloc)_$(expidloc)_$(chiploc).png"
    fig.savefig(tracePlotPath, bbox_inches = "tight", pad_inches = 0.1)
    thread("Trace extraction for $(parg["tele"]) $(chiploc) $(mjdloc) $(expidloc) done.")
    PythonPlot.plotclose(fig)
    thread("Here is the median flux per fiber", tracePlotPath)
end
thread("DomeFlat traces done.")
