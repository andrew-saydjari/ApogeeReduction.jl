using ArgParse, Distributed, SlurmClusterManager, SlackThreads

## Parse command line arguments
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--chips"
        required = false
        help = "chip names, i.e. abc"
        default = "abc"
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

# memory usage is under 3GB per worker in testing
if parg["runlist"] != "" # only multiprocess if we have a list of exposures
    if "SLURM_NTASKS" in keys(ENV)
        using SlurmClusterManager
        addprocs(SlurmManager(), exeflags = ["--project=./"])
    else
        addprocs(16)
    end
end

@everywhere begin
    using JLD2, ProgressMeter, ArgParse, Glob, StatsBase, ParallelDataTransfer
    src_dir = "../"
    include(src_dir * "/makie_plotutils.jl")
    include(src_dir * "/utils.jl")
    include(src_dir * "/fileNameHandling.jl")
    include(src_dir * "/cal_build/traceExtract_GH.jl")
end

@passobj 1 workers() parg # make it available to all workers

@everywhere begin
    chips = collect(String, split(parg["chips"], ""))

    # make summary plots and gifs and send to slack in outer loop

    dirNamePlots = parg["trace_dir"] * "plots/"
    if !ispath(dirNamePlots)
        mkpath(dirNamePlots)
    end

    # Load in the exact set of exposures
    mjd = load(parg["runlist"], "mjd")
    expid = load(parg["runlist"], "expid")
    flist = [get_cal_file(parg["trace_dir"], parg["tele"], mjd[i],
                 expid[i], chip, "DOMEFLAT", use_cal = true)
             for chip in chips, i in eachindex(mjd)]

    fpifib1, fpifib2 = get_fpi_guide_fiberID(parg["tele"])
end

desc = "trace extract for $(parg["tele"]) $(chips)"
plot_paths = @showprogress desc=desc pmap(flist) do fname
    sname = split(fname, "_")
    teleloc, mjdloc, chiploc, expidloc = sname[(end - 4):(end - 1)]
    mjdfps2plate = get_fps_plate_divide(teleloc)
    f = jldopen(fname)
    image_data = f["dimage"][1:2048, 1:2048]
    ivar_image = f["ivarimage"][1:2048, 1:2048]
    pix_bitmask_image = f["pix_bitmask"][1:2048, 1:2048]
    close(f)

    med_center_to_fiber_func, x_prof_min, x_prof_max_ind, n_sub, min_prof_fib, max_prof_fib,
    all_y_prof, all_y_prof_deriv = gh_profiles(
        teleloc, mjdloc, chiploc, expidloc; n_sub = 100)

    #    trace_params, trace_param_covs = trace_extract(
    #        image_data, ivar_image, teleloc, mjdloc, chiploc, expidloc; good_pixels = nothing)
    good_pixels = ((pix_bitmask_image .& bad_pix_bits) .== 0)
    trace_params,
    trace_param_covs = trace_extract(
        image_data, ivar_image, teleloc, mjdloc, chiploc, expidloc,
        med_center_to_fiber_func, x_prof_min, x_prof_max_ind, n_sub, min_prof_fib, max_prof_fib, all_y_prof, all_y_prof_deriv
        ; good_pixels = good_pixels)

    safe_jldsave(
        parg["trace_dir"] *
        "dome_flats/domeTrace_$(teleloc)_$(mjdloc)_$(expidloc)_$(chiploc).jld2";
        trace_params = trace_params, trace_param_covs = trace_param_covs)

    return trace_plots("dome", trace_params, teleloc, mjdloc, chiploc, expidloc, mjdfps2plate, fpifib1, fpifib2)
end

thread = SlackThread()
thread("DOMEFLAT TRACES for $(parg["tele"]) $(chips) from $(parg["mjd-start"]) to $(parg["mjd-end"])")
for (filename, heights_widths_path) in zip(flist, plot_paths)
    thread("Here is the median flux and width per fiber for $(filename)", heights_widths_path)
end
thread("DomeFlat traces done.")
