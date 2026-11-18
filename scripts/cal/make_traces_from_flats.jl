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
        "--slack_quiet"
        required = false
        help = "if true, don't send slack notifications"
        arg_type = Bool
        default = false
        "--checkpoint_mode"
        required = false
        help = "checkpoint mode (clobber, commit_exists, commit_same)"
        arg_type = String
        default = "commit_same"
        "--trace_identify"
        required = false
        help = "fiber identification on the flat: \"robust\" (default: template-free whole-slit identification, refused with a named reason when ambiguous; unlit fibers placed from lit neighbours on both sides; see src/trace_qa.jl) or \"production\" (the pre-existing path). The per-column refinement of lit fibers is the same in both."
        arg_type = String
        default = "robust"
        "--almanac_file"
        required = false
        help = "almanac hdf5 file; used to identify the FPI guide fibers from the night's configurations for the dome-flat trace plots. If empty, FPI fibers are not annotated."
        arg_type = String
        default = ""
        "--trace_template_cal"
        required = false
        help = "per-epoch trace-template calibration file (HDF5, scripts/cal/trace_templates/; src/trace_template_cal.jl). OFF BY DEFAULT (none): traces exactly as without templates. A path turns it on: each flat's QA record carries the night's template identity (trace_template_id) and a checkpointed trace whose record has another identity, or none, is refitted. A path that does not exist or cannot be read is a hard error. Use the same file as pipeline_2d_1d.jl --trace_template_cal."
        arg_type = String
        default = "none"
    end
    return parse_args(s)
end

parg = parse_commandline()
# trace templates: off by default; a configured path that is missing is a hard error
const TRACE_TEMPLATE_OFF = strip(parg["trace_template_cal"]) in ("", "none")
TRACE_TEMPLATE_OFF || isfile(parg["trace_template_cal"]) ||
    error("trace-template calibration file $(parg["trace_template_cal"]) (--trace_template_cal) does not exist; " *
          "install it, or run with --trace_template_cal none (scripts/cal/trace_templates/README.md, Bootstrap).")

# memory usage is under 3GB per worker in testing
proj_path = dirname(Base.active_project()) * "/"
if parg["runlist"] != "" # only multiprocess if we have a list of exposures
    if "SLURM_NTASKS" in keys(ENV)
        using SlurmClusterManager
        addprocs(SlurmManager(), exeflags = ["--project=$proj_path"])
    else
        addprocs(parse(Int, get(ENV, "AR_LOCAL_WORKERS", "16")),
            exeflags = ["--project=$proj_path"])
    end
end

@everywhere begin
    using JLD2, ProgressMeter, ArgParse, Glob, StatsBase, ParallelDataTransfer
    using ApogeeReduction
    using ApogeeReduction: get_cal_file, get_fpi_fiberIDs_from_almanac, get_fps_plate_divide, trace_extract,
                           safe_jldsave, trace_plots, bad_pix_bits, check_file, nanzeropercentile, regularize_trace,
                           trace_lamp_metrics, flat_trace_qa, write_trace_qa, read_trace_qa, trace_qa_path,
                           trace_reference_positions, load_trace_template_cal, trace_template_metadata,
                           template_needs_recompute
end
@everywhere TRACE_TEMPLATE_CAL = $(TRACE_TEMPLATE_OFF) ? nothing :
                                 load_trace_template_cal($(parg["trace_template_cal"]))
println(TRACE_TEMPLATE_OFF ? "Trace templates: OFF (--trace_template_cal none)" :
        "Trace templates from $(TRACE_TEMPLATE_CAL) [$(ApogeeReduction.trace_template_cal_id(TRACE_TEMPLATE_CAL))]")

@passobj 1 workers() parg # make it available to all workers
@passobj 1 workers() proj_path

# maybe the flist should be passed to everywhere?
@everywhere begin
    include(joinpath(proj_path, "src/makie_plotutils.jl"))

    chips = collect(String, split(parg["chips"], ""))
    parg["trace_identify"] in ("production", "robust") ||
        error("--trace_identify must be production or robust, got $(parg["trace_identify"])")

    # make summary plots and gifs and send to slack in outer loop

    dirNamePlots = joinpath(parg["trace_dir"], "plots/")
    if !ispath(dirNamePlots)
        mkpath(dirNamePlots)
    end

    # Load in the exact set of exposures
    tele_list = load(parg["runlist"], "tele")
    unique_teles = if parg["tele"] == ""
        unique(tele_list)
    else
        [parg["tele"]]
    end
    mjd_list = load(parg["runlist"], "mjd")
    unique_mjds = unique(mjd_list)
    dfindx_list = load(parg["runlist"], "dfindx")
    cal_flist = String[] # all cal files for FIRST CHIP
    for tele in unique_teles
        expid2do = findall(tele_list .== tele)
        for chip in chips
            for i in expid2do
                push!(cal_flist, get_cal_file(parg["trace_dir"], tele_list[i], mjd_list[i], dfindx_list[i], chip, "$(parg["flat_type"])flat", use_cal = true))
            end
        end
    end
    flist = vcat(cal_flist...)

    # FPI guide fibers are derived per night from the configurations in the
    # almanac ("bonus" stubs), not from a hardcoded per-telescope pair: the
    # constant is wrong on the LCO FPS commissioning nights (MJD 59820/59826/
    # 59827, where the FPI is on fiber_id 142/153, not 82/213). Flats carry no
    # configuration of their own, so the night's object configurations are used;
    # the FPI feed does not move within a night. Cached per MJD -- one small
    # HDF5 read per night rather than per file.
    _fpi_cache = Dict{String, Vector{Int}}()
    function fpi_fiberIDs_for(teleloc, mjdloc)
        parg["almanac_file"] == "" && return Int[]
        get!(_fpi_cache, "$(teleloc)_$(mjdloc)") do
            get_fpi_fiberIDs_from_almanac(parg["almanac_file"], teleloc, mjdloc)
        end
    end

    # Every flat gets a QA record (`<flat_type>TraceQA_<tele>_<mjd>_<exp>_<chip>.h5`, next to
    # the trace file): the outcome ("ok", or the NAMED reason it was refused: lamp_off,
    # trace_count, no_flux, id_*, fit_error) and the gate metrics (src/trace_qa.jl). A refusal is
    # still a refusal -- the same flats are dropped as before in the default mode -- but it is no
    # longer silent: create_traceMain and anyone auditing a night read the reason from the record.
    function write_flat_qa(qaname, image_data, ivar_image, pix_bitmask_image, teleloc, chiploc,
            trace_params, outcome, extra)
        qa = flat_trace_qa(image_data, ivar_image, pix_bitmask_image, teleloc, chiploc;
            trace_params = trace_params, outcome = outcome,
            refpos = trace_reference_positions(teleloc, chiploc; data_path = joinpath(proj_path, "data")),
            extra = extra)
        write_trace_qa(qaname, qa)
        return qa
    end

    function make_traces(fname, flat_type; checkpoint_mode = "commit_same", identify = "robust",
            template_cal = nothing)
        sname = split(split(split(fname, "/")[end], ".h5")[1], "_")
        fnameType, teleloc, mjdloc, expnumloc, chiploc, image_type = sname[(end - 5):end]
        mjdfps2plate = get_fps_plate_divide(teleloc)
        # Lamp-off gate (unchanged; now `trace_lamp_metrics` in src/trace_qa.jl): flux_68p must
        # reach ~20% of the typical value (of the smallest-flux chip, and smallest flux from dome
        # vs quartz) from days when lamps were on: 40 (APO), 10 (LCO).
        # n.b. the median test (16 APO / 10 LCO) is DELIBERATELY NOT USED: turning it on is held
        # pending the post-pass-1 classifier work (AKS 2026-09-08, "let's see how far we can get
        # with the classifier after pass1"), because on the testbed corpus it fixes a confirmed
        # false acceptance but also false-rejects at least one genuinely good faint LCO
        # quartzflat. Do not enable it as a drive-by cleanup; do not delete it either. It is
        # recorded in the QA record as `flux_med`.

        savename = joinpath(parg["trace_dir"], "$(flat_type)_flats", "$(mjdloc)", "$(flat_type)Trace_$(teleloc)_$(mjdloc)_$(expnumloc)_$(chiploc).h5")
        qaname = trace_qa_path(parg["trace_dir"], flat_type, teleloc, mjdloc, expnumloc, chiploc)
        extra = Dict{String, Any}("flat_type" => String(flat_type), "mjd" => parse(Int, mjdloc),
            "expnum" => parse(Int, expnumloc), "image_file" => String(fname),
            "trace_file" => "", "identify" => String(identify))
        # the night's template identity ("none" when off) goes into the QA record; a trace whose
        # record carries another identity (or none, or has no record while templates are on) is
        # refitted
        merge!(extra, trace_template_metadata(template_cal, teleloc, chiploc, parse(Int, mjdloc)))
        trace_done = check_file(savename, mode = checkpoint_mode)
        if trace_done
            qa_done = check_file(qaname, mode = checkpoint_mode)
            rec = qa_done ? read_trace_qa(qaname) : Dict{String, Any}()
            if template_needs_recompute(rec, template_cal, teleloc, chiploc, parse(Int, mjdloc))
                trace_done = false
            elseif qa_done
                return trace_plots(dirNamePlots, flat_type, savename, teleloc, mjdloc, expnumloc, chiploc, mjdfps2plate, fpi_fiberIDs_for(teleloc, mjdloc); checkpoint_mode = checkpoint_mode)
            end
        end

        f = jldopen(fname, "r")
        image_data = f["dimage"][1:N_XPIX, 1:N_XPIX]
        ivar_image = f["ivarimage"][1:N_XPIX, 1:N_XPIX]
        pix_bitmask_image = f["pix_bitmask"][1:N_XPIX, 1:N_XPIX]
        close(f)

        if trace_done
            # trace kept by the checkpoint (e.g. made before QA records existed): add its QA record
            extra["trace_file"] = savename
            extra["identify"] = "unknown(checkpointed)"
            write_flat_qa(qaname, image_data, ivar_image, pix_bitmask_image, teleloc, chiploc,
                load(savename, "trace_params"), "ok", extra)
            return trace_plots(dirNamePlots, flat_type, savename, teleloc, mjdloc, expnumloc, chiploc, mjdfps2plate, fpi_fiberIDs_for(teleloc, mjdloc); checkpoint_mode = checkpoint_mode)
        end

        good_pixels = ((pix_bitmask_image .& bad_pix_bits) .== 0)
        lamp = trace_lamp_metrics(image_data, pix_bitmask_image, teleloc)
        flux_med, flux_68p = lamp.flux_med, lamp.flux_68p
        if !lamp.lamp_ok
            @warn "File $(fname) (med_flux = $(round(flux_med,digits=2)), 68% interval = $(round(flux_68p,digits=2))) was skipped for appearing to have the lamp turned off. Used $(round(sum(good_pixels)/length(good_pixels)*100,digits=2))% of the pixels for this calculation."
            write_flat_qa(qaname, image_data, ivar_image, pix_bitmask_image, teleloc, chiploc,
                nothing, "lamp_off", extra)
            return nothing
        end

        (med_center_to_fiber_func, x_prof_min, x_prof_max_ind, n_sub, min_prof_fib, max_prof_fib,
        all_y_prof, all_y_prof_deriv) = ApogeeReduction.get_default_trace_hyperparams(teleloc, chiploc, profile_path = joinpath(proj_path, "data"), plot_path = joinpath(parg["trace_dir"], "plots/"))

        tqa = Dict{String, Any}()
        trace_params, trace_param_covs = try
            trace_extract(
                image_data, ivar_image, teleloc, mjdloc, expnumloc, chiploc,
                med_center_to_fiber_func, x_prof_min, x_prof_max_ind, n_sub, 
		min_prof_fib, max_prof_fib, all_y_prof, all_y_prof_deriv, fname; 
                good_pixels = good_pixels, median_trace_pos_path = joinpath(proj_path, "data"),
                identify = identify, qa = tqa)
        catch e
            println("Trace fitting failed for $(fname)")
            println("Flux info is (med_flux = $(round(flux_med,digits=2)), 68% interval = $(round(flux_68p,digits=2)), percent_good_pixels = $(round(sum(good_pixels)/length(good_pixels)*100,digits=2))%)")
            extra["fit_error"] = first(sprint(showerror, e), 500)
            write_flat_qa(qaname, image_data, ivar_image, pix_bitmask_image, teleloc, chiploc,
                nothing, "fit_error", extra)
            # retried once as before (a second failure propagates, as before)
            trace_extract(
                image_data, ivar_image, teleloc, mjdloc, expnumloc, chiploc,
                med_center_to_fiber_func, x_prof_min, x_prof_max_ind, n_sub, 
		min_prof_fib, max_prof_fib, all_y_prof, all_y_prof_deriv, fname;
                good_pixels = good_pixels, median_trace_pos_path = joinpath(proj_path, "data"),
                identify = identify, qa = tqa)
        end
        for k in ("n_detected", "n_good_throughput", "n_unlit")
            haskey(tqa, k) && (extra[k] = tqa[k])
        end

        if size(trace_params,1) == 0
	    # KEVIN should look at what is going on with these weird LCO early MJD (e.g. 57801) cases more closely in future
            @warn "File $(fname) (med_flux = $(round(flux_med,digits=2)), 68% interval = $(round(flux_68p,digits=2)), percent_good_pixels = $(round(sum(good_pixels)/length(good_pixels)*100,digits=2))%) was skipped: $(get(tqa, "outcome", "not the correct amount of traces"))."
            write_flat_qa(qaname, image_data, ivar_image, pix_bitmask_image, teleloc, chiploc,
                nothing, get(tqa, "outcome", "trace_count"), extra)
            return nothing
	end

        regularized_trace_params = regularize_trace(trace_params)

        mkpath(dirname(savename))
        safe_jldsave(savename; trace_params = trace_params, trace_param_covs = trace_param_covs, regularized_trace_params = regularized_trace_params, no_metadata = true)
        extra["trace_file"] = savename
        write_flat_qa(qaname, image_data, ivar_image, pix_bitmask_image, teleloc, chiploc,
            trace_params, "ok", extra)
    
        return trace_plots(dirNamePlots, flat_type, savename, teleloc, mjdloc, expnumloc, chiploc, mjdfps2plate, fpi_fiberIDs_for(teleloc, mjdloc); checkpoint_mode = checkpoint_mode)
    end
end

desc = "trace extract for $(parg["tele"]) $(chips)"
@everywhere make_traces_partial(fname) = make_traces(fname, parg["flat_type"], checkpoint_mode = parg["checkpoint_mode"],
    identify = parg["trace_identify"], template_cal = TRACE_TEMPLATE_CAL)
plot_paths = @showprogress desc=desc pmap(make_traces_partial, flist)

if !parg["slack_quiet"]
    thread = SlackThread()
    if length(unique_mjds) > 1
        thread("$(parg["flat_type"])Flat Traces for $(parg["tele"]) $(chips) from SJD $(minimum(unique_mjds)) to $(maximum(unique_mjds))")
        for (filename, heights_widths_path) in zip(flist, plot_paths)
            if isnothing(heights_widths_path)
                sname = split(split(basename(filename), ".h5")[1], "_")
                qa = read_trace_qa(trace_qa_path(parg["trace_dir"], parg["flat_type"], sname[end - 4], sname[end - 3], sname[end - 2], sname[end - 1]))
                why = isnothing(qa) ? "unknown reason" : get(qa, "outcome", "unknown reason")
                thread("WARNING: File $(filename) was skipped: $(why).")
                continue
            end
            thread("Here is the median flux and width per fiber for $(filename)", heights_widths_path)
        end
    else
        thread("$(parg["flat_type"])Flat Traces for $(parg["tele"]) $(chips) had no data")
    end
    thread("$(parg["flat_type"])Flat traces done.")
end
