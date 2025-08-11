using ArgParse, Distributed, SlurmClusterManager, SlackThreads

## Parse command line arguments
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--chips"
        required = false
        help = "chip names, i.e. RGB"
        default = "RGB"
        arg_type = String
        "--tele"
        required = true
        help = "telescope name (apo or lco)"
        arg_type = String
        default = ""
        "--trace_dir"
        required = true
        help = "directory where 2D extractions of traces are stored"
        arg_type = String
        default = ""
        "--runlist"
        required = true
        help = "path name to hdf5 file with keys specifying list of exposures to run"
        "--flat_type"
        required = true
        help = "flat type, i.e. dome or quartz"
        arg_type = String
        default = "quartz"
    end
    return parse_args(s)
end

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
    using ApogeeReduction: get_cal_file, get_fpi_guide_fiberID, get_fps_plate_divide, gh_profiles, trace_extract, safe_jldsave, trace_plots, bad_pix_bits
    include("../../src/makie_plotutils.jl")
end

@passobj 1 workers() parg # make it available to all workers

# maybe the flist should be passed to everywhere?
@everywhere begin
    chips = collect(String, split(parg["chips"], ""))

    # make summary plots and gifs and send to slack in outer loop

    dirNamePlots = parg["trace_dir"] * "plots/"
    if !ispath(dirNamePlots)
        mkpath(dirNamePlots)
    end

    # Load in the exact set of exposures
    tele_list = load(parg["runlist"], "tele")
    unique_teles = unique(tele_list)
    mjd_list = load(parg["runlist"], "mjd")
    unique_mjds = unique(mjd_list)
    expid_list = load(parg["runlist"], "expid")
    cal_flist = String[] # all cal files for FIRST CHIP
    for tele in unique_teles
        expid2do = findall(tele_list .== tele)
        for chip in chips
            for i in expid2do
                push!(cal_flist, get_cal_file(parg["trace_dir"], tele_list[i], mjd_list[i], expid_list[i], chip, "$(uppercase(parg["flat_type"]))FLAT", use_cal = true))
            end
        end
    end
    flist = vcat(cal_flist...)

    fpifib1, fpifib2 = get_fpi_guide_fiberID(parg["tele"])
end

desc = "trace extract for $(parg["tele"]) $(chips)"
plot_paths = @showprogress desc=desc pmap(flist) do fname
    sname = split(split(split(fname, "/")[end],".h5")[1], "_")
    fnameType, teleloc, mjdloc, expnumloc, chiploc, exptype = sname[(end - 5):end]
    
    mjdfps2plate = get_fps_plate_divide(teleloc)
    f = jldopen(fname)
    image_data = f["dimage"][1:2048, 1:2048]
    ivar_image = f["ivarimage"][1:2048, 1:2048]
    pix_bitmask_image = f["pix_bitmask"][1:2048, 1:2048]
    close(f)

    med_center_to_fiber_func, x_prof_min, x_prof_max_ind, n_sub, min_prof_fib, max_prof_fib,
    all_y_prof, all_y_prof_deriv = gh_profiles(
        teleloc, mjdloc, expnumloc, chiploc; n_sub = 100)

    #    trace_params, trace_param_covs = trace_extract(
    #        image_data, ivar_image, teleloc, mjdloc, chiploc, expidloc; good_pixels = nothing)
    good_pixels = ((pix_bitmask_image .& bad_pix_bits) .== 0)
    trace_params,
    trace_param_covs = trace_extract(
        image_data, ivar_image, teleloc, mjdloc, expnumloc, chiploc,
        med_center_to_fiber_func, x_prof_min, x_prof_max_ind, n_sub, min_prof_fib, max_prof_fib, all_y_prof, all_y_prof_deriv
        ; good_pixels = good_pixels)

    safe_jldsave(
        parg["trace_dir"] *
        "$(parg["flat_type"])_flats/$(parg["flat_type"])Trace_$(teleloc)_$(mjdloc)_$(expnumloc)_$(chiploc).h5";
        trace_params = trace_params, trace_param_covs = trace_param_covs)

    return trace_plots(dirNamePlots, parg["flat_type"], trace_params, teleloc, mjdloc, expnumloc, chiploc, mjdfps2plate, fpifib1, fpifib2)
end

thread = SlackThread()
if length(unique_mjds) > 1
    thread("$(parg["flat_type"])Flat Traces for $(parg["tele"]) $(chips) from SJD $(minimum(unique_mjds)) to $(maximum(unique_mjds))")
    for (filename, heights_widths_path) in zip(flist, plot_paths)
        thread("Here is the median flux and width per fiber for $(filename)", heights_widths_path)
    end
else
    thread("$(parg["flat_type"])Flat Traces for $(parg["tele"]) $(chips) had no data")
end
thread("$(parg["flat_type"])Flat traces done.")
