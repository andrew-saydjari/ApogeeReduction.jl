## This is a reduction pipeline for APOGEE
using InteractiveUtils
versioninfo();
@time "Package activation" begin
    import Pkg
    Pkg.instantiate()
    Pkg.precompile() # no need for Pkg.activate("./") because of invocation w/ environment
    using Distributed, ArgParse, TimerOutputs
end

## Parse command line arguments
# we have forced this one to run over all chips for the sake of the wavelength solution
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--tele"
        required = true
        help = "telescope name (apo or lco)"
        arg_type = String
        default = ""
        "--mjd"
        required = false
        help = "mjd of the exposure(s) to be run"
        arg_type = Int
        default = 1
        # probably want to add in defaults that loops over them all
        "--dfindx"
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
        required = false
        help = "output directory"
        arg_type = String
        default = "../outdir/"
        "--runname"
        required = true
        help = "name of the run (specifically almanac file)"
        arg_type = String
        default = "test"
        "--extraction"
        required = false
        help = "extraction method (boxcar or optimal)"
        arg_type = String
        default = "optimal"
        "--workers_per_node"
        required = false
        help = "number of workers per node"
        arg_type = Int
        default = -1 # -1 means use all the cores on the node
        "--checkpoint_mode"
        required = false
        help = "checkpoint mode (clobber, commit_exists, commit_same)"
        arg_type = String
        default = "commit_same"
        "--relFlux"
        required = false
        help = "use relFluxing (true or false)"
        arg_type = Bool
        default = true
        "--waveSoln"
        required = false
        help = "do the wavelength solution (true or false)"
        arg_type = Bool
        default = true
        "--doUncals"
        required = false
        help = "move uncalibrated (no flats/darks etc.) 2d data forward"
        arg_type = Bool
        default = false
        "--arclamp_wavecal"
        required = false
        help = "write the stage-1 arc-lamp night product arclampWaveCalNight_<tele>_<mjd>.h5 (read by nothing yet; default false)"
        arg_type = Bool
        default = false
        "--arclamp_linelist"
        required = false
        help = "arc line list for --arclamp_wavecal (default: data/arclamp_linelist_v1.csv in the repository)"
        arg_type = String
        default = ""
        "--arclamp_linelist_selection"
        required = false
        help = "selection column of the arc line list (use_default, use_with_neon, use_with_blended, prelim_classB)"
        arg_type = String
        default = "use_default"
        "--chip_transform_cal"
        required = false
        help = "comb-measured chip-transform calibration file (HDF5, built by scripts/cal/chip_transform/build_chip_transform.jl). A calibration input configured by path like --exp_class_model; ON BY DEFAULT with the pinned file below (airflow/dags/ar_common.py sets it for production). The sky-line wavelength fit then fixes the R/B chip scales to the file's values instead of 1. Pass the word none (or an empty string) to turn it off deliberately: chip scale 1, the pre-calibration behaviour, and the bootstrap mode used to rebuild the file from scratch (scripts/cal/chip_transform/README.md, Bootstrap). A path that does not exist or cannot be read is a hard error, never a silent fallback to scale 1."
        arg_type = String
        default = "/mnt/ceph/users/sdssv/work/asaydjari/cal_ref/chip_transform/chiptransform_pass1dr21_v2_allblock.h5"
        "--wavecal_order"
        required = false
        help = "order of the wavelength chain: legacy (sky -> FPI comb -> night average -> per-exposure sky/FPI dither; the default, unchanged) or arc_first (arcs -> arc-seeded comb -> night assembly -> per-exposure FPI-fiber/sky terms, each stage once and in order; src/wavecal_arcfirst.jl). Nights without a usable arc fit run the legacy chain and are flagged WAVE_NOT_ARCFIRST."
        arg_type = String
        default = "legacy"
        "--arcfirst_shape_template"
        required = false
        help = "arc_first only: comb shape template (HDF5, scripts/refit/build_comb_shape_template.jl) for FPS-era exposures on the arc route (no night comb, or no usable FPI fibers), as in the archive refit. none (or empty) = arc shape (WAVE_SHAPE_ARC). A path that does not exist is a hard error."
        arg_type = String
        default = "/mnt/ceph/users/sdssv/work/asaydjari/2026_11_17/wavecal_refit/cal/template_shape_v1.h5"
        "--arcfirst_fps_shape"
        required = false
        help = "arc_first only: shape of FPS-era arc-route exposures, both telescopes: template (the refit default) or arc. The plate era always uses the arc shape (PI 2026-11-18)."
        arg_type = String
        default = "template"
        "--arcfirst_fiber_fill"
        required = false
        help = "arc_first only: fill_only (fibers without an arc fit get the across-slit interpolate_wave_params value; fitted fibers unchanged; WAVE_FIBER_FILLED) or none (such fibers have no solution and the exposure falls back at resampling)"
        arg_type = String
        default = "fill_only"
        "--arcfirst_max_reanchor"
        required = false
        help = "arc_first only: cap on conditional re-anchors per night when the end-of-night check (median sky-line offset of the FPI-fiber exposures against the final solution, per chip, and comb mode slips) fails; each re-anchor takes the anchor rule's next candidate and reruns comb -> exposures for that night (WAVE_REANCHORED; cap reached: WAVE_ANCHOR_CHECK_FAILED, best candidate kept)"
        arg_type = Int
        default = 2
        "--arcfirst_sky_check_ms"
        required = false
        help = "arc_first only: end-of-night check threshold on |per-chip median sky offset| (m/s); a wrong anchor block is a dither step (2-3 km/s), the arc-vs-sky ensemble offset 40-140 m/s"
        arg_type = Float64
        default = 1000.0
        "--arcfirst_noarc_fallback"
        required = false
        help = "arc_first only: what a night without a usable arc fit runs. Only legacy (the legacy sky/FPI chain for that night, WAVE_NOT_ARCFIRST) is implemented."
        arg_type = String
        default = "legacy"
    end
    return parse_args(s)
end
parg = parse_commandline()
dirNamePlots = parg["outdir"] * "plots/"
mkpath(dirNamePlots) # will work even if it already exists

# if we are not refluxing, we should not be doing the wavelength solution
if !parg["relFlux"] && parg["waveSoln"]
    error("Should not perform wavelength solution when relFlux is false. Set --waveSoln=false or --relFlux=true")
end

# The chip-transform calibration is ON by default. "none" (or "") turns it off: the sky fit keeps
# chip scale 1. That is also the BOOTSTRAP mode for (re)building the file from scratch, which needs
# no chip-transform file at all (scripts/cal/chip_transform/README.md, "Bootstrap"). A configured
# file that is missing is a hard error here (an unreadable one when it is loaded), never a silent
# chip scale 1.
const CHIP_TRANSFORM_OFF = strip(parg["chip_transform_cal"]) in ("", "none")
if !CHIP_TRANSFORM_OFF && !isfile(parg["chip_transform_cal"])
    error("chip-transform calibration file $(parg["chip_transform_cal"]) (--chip_transform_cal) does not exist. Either\n" *
          "  (1) install it at that path (or point --chip_transform_cal / CHIP_TRANSFORM_CAL in airflow/dags/ar_common.py /\n" *
          "      AR_CHIP_TRANSFORM_CAL for the SLURM scripts at an existing file), or\n" *
          "  (2) run with --chip_transform_cal none (chip scale 1; bootstrap mode) -- this is how the file is rebuilt\n" *
          "      from scratch: see the \"Bootstrap\" section of scripts/cal/chip_transform/README.md.")
end

# wavelength-chain order (legacy = unchanged default; arc_first = src/wavecal_arcfirst.jl)
parg["wavecal_order"] in ("legacy", "arc_first") ||
    error("--wavecal_order must be legacy or arc_first, got $(parg["wavecal_order"])")
if parg["wavecal_order"] == "arc_first"
    strip(parg["arcfirst_shape_template"]) in ("", "none") && (parg["arcfirst_shape_template"] = "")
    (isempty(parg["arcfirst_shape_template"]) || isfile(parg["arcfirst_shape_template"])) ||
        error("--arcfirst_shape_template $(parg["arcfirst_shape_template"]) does not exist (pass none for the arc shape)")
    parg["arcfirst_fps_shape"] in ("template", "arc") || error("--arcfirst_fps_shape must be template or arc")
    parg["arcfirst_fiber_fill"] in ("fill_only", "none") || error("--arcfirst_fiber_fill must be fill_only or none")
    parg["arcfirst_noarc_fallback"] == "legacy" || error("--arcfirst_noarc_fallback: only legacy is implemented")
end

workers_per_node = parg["workers_per_node"]
proj_path = dirname(Base.active_project()) * "/"

@time "Allocating workers" begin
    if parg["runlist"] != "" # only multiprocess if we have a list of exposures
        if "SLURM_NTASKS" in keys(ENV)
            using SlurmClusterManager
            addprocs(SlurmManager(), exeflags = ["--project=$proj_path"])
            if workers_per_node != -1
                ntasks = parse(Int, ENV["SLURM_NTASKS"])
                nnodes = parse(Int, ENV["SLURM_NNODES"])
                cpus_per_node = parse(Int, ENV["SLURM_CPUS_ON_NODE"])
                total_workers = nnodes * workers_per_node
                workers_to_keep = []
                for node in 0:(nnodes - 1)
                    node_start = 1 + node * cpus_per_node
                    spacing = cpus_per_node ÷ workers_per_node
                    append!(workers_to_keep,
                        [node_start + spacing * i for i in 0:(workers_per_node - 1)])
                end
                rmprocs(setdiff(1:ntasks, workers_to_keep))
            end
        else
            addprocs(workers_per_node, exeflags = ["--project=$proj_path"])
        end
    end
end
println("Running Main on ", gethostname());
flush(stdout);

@time "Setting up worker packages" @everywhere begin
    using LinearAlgebra
    BLAS.set_num_threads(1)
    using FITSIO, HDF5, FileIO, JLD2, Glob, CSV
    using DataFrames, EllipsisNotation, StatsBase
    using AstroTime # can remove after Adam merges the PR to recast as Float
    using ParallelDataTransfer, ProgressMeter
    using ApogeeReduction
    using ApogeeReduction: read_almanac_exp_df, get_1d_name, get_and_save_sky_wavecal,
                           get_and_save_sky_dither_per_fiber, get_and_save_sky_peaks,
                           get_ave_night_wave_soln, sky_wave_plots, reinterp_spectra,
                           get_and_save_arclamp_peaks, get_and_save_fpi_peaks,
                           comb_exp_get_and_save_fpi_wavecal, skyline_medwavecal_skyline_dither,
                           fpi_medwavecal_skyline_dither,
                           get_and_save_arclamp_night_wavecal, load_arclamp_linelist,
                           safe_jldsave, process_1D, check_file, modify_saved_dict,
                           load_chip_transform_cal

    ###decide order to look for traces from these cal types (i.e. dome or quartz flats)
    trace_type_order = ["quartz","dome"]
end

println(BLAS.get_config());
flush(stdout);

@time "Passing objects to workers" begin
    # Is this really causing 3 min of overhead?
    @passobj 1 workers() parg
    @everywhere proj_path = dirname(Base.active_project()) * "/"
end

##### 1D stage
@time "include makie_plotutils" @everywhere include(joinpath(proj_path, "src/makie_plotutils.jl"))

@time "Generating file lists" begin
    # Load everywhere so we can make list in parallel
    @everywhere begin
        tele_list = if parg["runlist"] != ""
            load(parg["runlist"], "tele")
        else
            [parg["tele"]]
        end
        # wow this is really not robust, you HAVE to pass a runlist and a tele on the command line
        unique_teles = unique(tele_list)
        mskTele = tele_list .== parg["tele"]

        mjd_list = if parg["runlist"] != ""
            load(parg["runlist"], "mjd")
        else
            [parg["mjd"]]
        end
        unique_mjds = unique(mjd_list[mskTele])

        dfindx_list = if parg["runlist"] != ""
            load(parg["runlist"], "dfindx")
        else
            [parg["dfindx"]]
        end
    end

    @everywhere function get_2d_names_for_mjd(mjd)
        df = read_almanac_exp_df(
            joinpath(parg["outdir"], "almanac/$(parg["runname"]).h5"), parg["tele"], mjd)
        function get_2d_name_partial(expid)
            parg["outdir"] * "/apred/$(mjd)/" *
            replace(get_1d_name(expid, df), "ar1D" => "ar2D") * ".h5"
        end
        mskMJD = (mjd_list .== mjd) .& mskTele
        return get_2d_name_partial.(dfindx_list[mskMJD])
    end
    
    list2Dexp = pmap(get_2d_names_for_mjd, unique_mjds)
    all2Da = vcat(list2Dexp...)

    all2Dperchip = []
    for chip in CHIP_LIST
        all2Dchip = replace.(all2Da, "_$(FIRST_CHIP)_" => "_$(chip)_")
        push!(all2Dperchip, all2Dchip)
    end
    all2D = vcat(all2Dperchip...)
end

### This should move to ar1D.jl

# make trace param file copies for each date
# searching outwards from the given mjd, 
# looking for cal types in the order of trace_type_order
@everywhere function create_traceMain(mjd,chip,trace_type_order;
				       tele = parg["tele"],
				       outdir = parg["outdir"],
				       checkpoint_mode = "commit_same")

    traceMain_fname = outdir * "apred/$(mjd)/traceMain_$(tele)_$(mjd)_$(chip).h5"
    if check_file(traceMain_fname, mode = checkpoint_mode)
        return
    end
    orig_trace_param_fname = nothing
    trace_type = nothing
    found_match = false

    for check_trace_type in trace_type_order
        # check for different cal types in preferential
        # order given by trace_type_order

        traceList = sort(glob("$(check_trace_type)Trace_$(tele)_$(mjd)_*_$(chip).h5",
            outdir * "$(check_trace_type)_flats/$(mjd)/"))
        if length(traceList) >= 1
            orig_trace_param_fname = traceList[end]
            trace_type = check_trace_type
            found_match = true
            break
        end
    end

    if !found_match
        orig_trace_param_fname = joinpath(proj_path, "data", "fallback_traces", "traceMain_$(tele)_60933_$(chip).h5")
        trace_type = "quartz_fallback"
        if !isfile(orig_trace_param_fname)
            error("Fallback trace file $(orig_trace_param_fname) does not exist")
        end
    end

    # TODO (probably for KM): 
    #     -use all the night's flats (quartz and dome) to define
    #          the time-evolution of trace params, then save
    #          that output (likely using this function)
    # KM says: this should be okay because ultimately we will
	# be using information about how the trace parameters change as a function
    # of time over a night, whose outputs we will want to save
    # (ie final version will produce a new trace param file per tele/mjd/chip combo) 
    
    # come back to why this symlink does not work
    # is this causing memory bloat?
    cp(orig_trace_param_fname, traceMain_fname, force = true)
    trace_metadata = Dict{String, Any}()
    trace_metadata["trace_type"] = trace_type
    trace_metadata["trace_found_match"] = found_match
    trace_metadata["trace_orig_param_fname"] = orig_trace_param_fname

    modify_saved_dict(traceMain_fname, "metadata", trace_metadata)
    return
end

# we should do somthing smart to assemble the traces from a night into a single file
# that gives us the trace of a fiber as a funciton of time or something
# for now, for each MJD, take the first one (or do that in the trace cal step)
# I think dome flats needs to swtich to dome_flats/mjd/
@time "Finding traceMain files" begin
    # Create all (mjd, chip) combinations for parallel processing
    mjd_chip_combinations = [(mjd, chip) for mjd in unique_mjds for chip in CHIP_LIST]
    
    @everywhere create_traceMain_wrapper(mjd_chip) = create_traceMain(
        mjd_chip[1], mjd_chip[2], trace_type_order;
        tele = parg["tele"],
        outdir = parg["outdir"],
        checkpoint_mode = parg["checkpoint_mode"])
    
    pmap(create_traceMain_wrapper, mjd_chip_combinations)
end
flush(stdout);

# extract the 2D to 1D
@everywhere process_1D_wrapper(fname) = process_1D(
    fname,
    outdir = parg["outdir"],
    runname = parg["runname"],
    extraction = parg["extraction"],
    relFlux = parg["relFlux"],
    chip_list = CHIP_LIST,
    profile_path = joinpath(proj_path, "data"),
    plot_path = joinpath(parg["outdir"], "plots/"),
    checkpoint_mode = parg["checkpoint_mode"]
)
if parg["doUncals"]
    desc = "Extracting 2D to 1D (uncals):"
    @showprogress desc=desc pmap(process_1D_wrapper, all2D)
end
all2Dcal = replace.(all2D, "ar2D" => "ar2Dcal")
desc = "Extracting 2Dcal to 1Dcal:"
@showprogress desc=desc pmap(process_1D_wrapper, all2Dcal)

### Only do the wavelength solution if we are relFluxing
if parg["relFlux"]
    ## get all OBJECT files (happy to add any other types that see sky?)
    ## also get FPI and arclamp files
    list1DexpObject = []
    list1DexpFPI = []
    list1DexpArclamp = []
    # this needs to be parallelized
    for mjd in unique_mjds
        df = read_almanac_exp_df(
            joinpath(parg["outdir"], "almanac/$(parg["runname"]).h5"), parg["tele"], mjd)
        function get_1d_name_partial(expid)
            if df.image_type[expid] == "object"
                return parg["outdir"] * "/apred/$(mjd)/" * get_1d_name(expid, df, cal = true) *
                       ".h5"
            else
                return nothing
            end
        end
        function get_1d_name_ARCLAMP_partial(expid)
            if (df.image_type[expid] == "arclamp") &
               ((df.lamp_thar[expid] == 1) | (df.lamp_une[expid] == 1))
                return parg["outdir"] * "/apred/$(mjd)/" * get_1d_name(expid, df, cal = true) *
                       ".h5"
            else
                return nothing
            end
        end
        function get_1d_name_FPI_partial(expid)
            if (df.image_type[expid] == "arclamp") & (df.lamp_thar[expid] == 0) &
               (df.lamp_une[expid] == 0)
                return parg["outdir"] * "/apred/$(mjd)/" * get_1d_name(expid, df, cal = true) *
                       ".h5"
            else
                return nothing
            end
        end
        mskMJD = (mjd_list .== mjd) .& mskTele
        local1D = get_1d_name_partial.(dfindx_list[mskMJD])
        push!(list1DexpObject, filter(!isnothing, local1D))
        local1D_fpi = get_1d_name_FPI_partial.(dfindx_list[mskMJD])
        push!(list1DexpFPI, filter(!isnothing, local1D_fpi))
        local1D_arclamp = get_1d_name_ARCLAMP_partial.(dfindx_list[mskMJD])
        push!(list1DexpArclamp, filter(!isnothing, local1D_arclamp))
    end
    all1DObjecta = vcat(list1DexpObject...)
    all1DFPIa = vcat(list1DexpFPI...)
    all1DArclampa = vcat(list1DexpArclamp...)
    all1DfpiPeaks_a = replace.(replace.(all1DFPIa, "ar1Dcal" => "fpiPeaks"), "ar1D" => "fpiPeaks")

    all1DObjectperchip = []
    all1DArclampperchip = []
    all1DFPIperchip = []
    for chip in CHIP_LIST
        all1DObjectchip = replace.(all1DObjecta, "_$(FIRST_CHIP)_" => "_$(chip)_")
        push!(all1DObjectperchip, all1DObjectchip)
        all1DArclampchip = replace.(all1DArclampa, "_$(FIRST_CHIP)_" => "_$(chip)_")
        push!(all1DArclampperchip, all1DArclampchip)
        all1DFPIchip = replace.(all1DFPIa, "_$(FIRST_CHIP)_" => "_$(chip)_")
        push!(all1DFPIperchip, all1DFPIchip)
    end
    all1DObject = vcat(all1DObjectperchip...)
    all1DArclamp = vcat(all1DArclampperchip...)
    all1DFPI = vcat(all1DFPIperchip...)

    ## load rough wave dict and sky lines list
    @everywhere begin
        roughwave_dict = load(joinpath(proj_path, "data", "roughwave_dict.jld2"), "roughwave_dict")
        df_sky_lines = CSV.read(joinpath(proj_path, "data", "APOGEE_lines.csv"), DataFrame)
        df_sky_lines.linindx = 1:size(df_sky_lines, 1)
    end

    ## get sky line peaks
    @everywhere get_and_save_sky_peaks_partial(fname) = get_and_save_sky_peaks(
        fname, roughwave_dict, df_sky_lines, checkpoint_mode = parg["checkpoint_mode"])
    desc = "Fitting sky line peaks: "
    @showprogress desc=desc pmap(get_and_save_sky_peaks_partial, all1DObject)

    # get arclamp peaks
    if size(all1DArclamp, 1) > 0
        ## get (non-fpi) arclamp peaks
        desc = "Fitting arclamp peaks: "
        @everywhere get_and_save_arclamp_peaks_partial(fname) = get_and_save_arclamp_peaks(
            fname, checkpoint_mode = parg["checkpoint_mode"])
        @showprogress desc=desc pmap(get_and_save_arclamp_peaks_partial, all1DArclamp)
    end

    all1DfpiPeaks_a = replace.(
        replace.(all1DFPIa, "ar1Dcal" => "fpiPeaks"), "ar1D" => "fpiPeaks")
    all1DfpiPeaks = if size(all1DFPI, 1) > 0
        ## get FPI peaks
        desc = "Fitting FPI peaks: "
        @everywhere get_and_save_fpi_peaks_partial(fname) = get_and_save_fpi_peaks(
            fname, data_path = joinpath(proj_path, "data"),
            checkpoint_mode = parg["checkpoint_mode"])
        @showprogress desc=desc pmap(get_and_save_fpi_peaks_partial, all1DFPI)
    else
        []
    end
    all1DfpiPeaks_out = reshape(all1DfpiPeaks, length(all1DfpiPeaks_a), length(CHIP_LIST))
    mskFPInothing = .!any.(isnothing.(eachrow(all1DfpiPeaks_out)))
    println("FPI peaks found for $(sum(mskFPInothing)) of $(length(mskFPInothing)) exposures")

    ## comb-measured chip transform (calibration input; nothing = chip scale fixed to 1)
    @everywhere CHIP_TRANSFORM_CAL = $(CHIP_TRANSFORM_OFF) ? nothing :
                                     load_chip_transform_cal($(parg["chip_transform_cal"]))
    println(CHIP_TRANSFORM_OFF ? "Sky wavecal chip scales: OFF (chip scale 1, --chip_transform_cal none)" :
            "Sky wavecal chip scales from $(CHIP_TRANSFORM_CAL) [$(ApogeeReduction.chip_transform_id(CHIP_TRANSFORM_CAL))]")

    ## ---- arc-first order (--wavecal_order arc_first; skipped entirely with legacy) ----------
    ## Each stage once, in order, reading only earlier stages (src/wavecal_arcfirst.jl):
    ## arcs -> arc-seeded comb -> FPI fibers -> night assembly + per-exposure terms -> waveCalExp.
    ## Nights without a usable arc fit continue below through the legacy chain, restricted to
    ## those nights; their exposures are flagged WAVE_NOT_ARCFIRST at resampling.
    if parg["wavecal_order"] == "arc_first"
        @everywhere ARCFIRST_CFG = ApogeeReduction.RefitConfig(
            chip_transform_file = $(CHIP_TRANSFORM_OFF) ? "" : $(parg["chip_transform_cal"]),
            shape_template_file = $(parg["arcfirst_shape_template"]),
            fps_shape_apo = $(parg["arcfirst_fps_shape"]), fps_shape_lco = $(parg["arcfirst_fps_shape"]),
            anchor_samepos = "none")
        almanac_file_af = joinpath(parg["outdir"], "almanac/$(parg["runname"]).h5")
        sendto(workers(), almanac_file_af = almanac_file_af)
        af_mjd(f) = parse(Int, split(basename(f), "_")[3])
        all1DcalFirst = replace.(all2Dperchip[1], "ar2D" => "ar1Dcal")
        af_jobs = [(mjd, ApogeeReduction.arcfirst_exps_from_names(all1DcalFirst, parg["tele"], mjd),
                    filter(f -> af_mjd(f) == mjd, all1DfpiPeaks_a)) for mjd in unique_mjds]
        @everywhere arcfirst_partial(job) = try
            ApogeeReduction.arcfirst_night(parg["tele"], job[1]; outdir = parg["outdir"],
                almanac_file = almanac_file_af, exps = job[2], fpi_fnames = job[3], cfg = ARCFIRST_CFG,
                checkpoint_mode = parg["checkpoint_mode"], fiber_fill = parg["arcfirst_fiber_fill"],
                noarc_fallback = parg["arcfirst_noarc_fallback"], max_reanchor = parg["arcfirst_max_reanchor"],
                sky_check_ms = parg["arcfirst_sky_check_ms"])
        catch e
            @warn "ARC-FIRST CHAIN FAILED for $(parg["tele"]) $(job[1]): the night falls back to the legacy chain (WAVE_NOT_ARCFIRST)" exception = (e, catch_backtrace())
            (; route = :legacy, arc = "error", comb = "", fpifib = "", assemble = "", n_exp_written = 0)
        end
        desc = "Arc-first wavelength chain (arcs -> comb -> exposures): "
        af_res = @showprogress desc=desc pmap(arcfirst_partial, af_jobs)
        for (job, r) in zip(af_jobs, af_res)
            println("arc_first $(parg["tele"]) $(job[1]): route=$(r.route) arc=$(r.arc) comb=$(r.comb) ",
                "fpifib=$(r.fpifib) assemble=$(r.assemble) waveCalExp=$(r.n_exp_written)/$(length(job[2])) ",
                "anchor_attempts=$(get(r, :n_attempts, 0))")
        end
        legacy_mjds = [job[1] for (job, r) in zip(af_jobs, af_res) if r.route == :legacy]
        println("arc_first: $(length(unique_mjds) - length(legacy_mjds)) of $(length(unique_mjds)) nights in order; ",
            "legacy chain for $(length(legacy_mjds)) nights without a usable arc fit: $(legacy_mjds)")
        af_keep(f) = af_mjd(f) in legacy_mjds
        all1DObjecta = filter(af_keep, all1DObjecta)
        all1DFPIa = filter(af_keep, all1DFPIa)
        all1DfpiPeaks_a = filter(af_keep, all1DfpiPeaks_a)
        all1DFPI = filter(af_keep, all1DFPI)
        all1DArclamp = filter(af_keep, all1DArclamp)
        unique_mjds = legacy_mjds
    end

    ## get wavecal from sky line peaks
    #only need to give one chip's list because internal
    #logic handles finding other chips when ingesting data
    #Andrew says that that is a bit worrisome and would should revisit that logic
    all1DObjectSkyPeaks = replace.(
        replace.(all1DObjecta, "ar1Dcal" => "skyLinePeaks"), "ar1D" => "skyLinePeaks")
    desc = "Skyline wavelength solutions:"
    @everywhere get_and_save_sky_wavecal_partial(fname) = get_and_save_sky_wavecal(
        fname, checkpoint_mode = parg["checkpoint_mode"], chip_transform_cal = CHIP_TRANSFORM_CAL)
    all1DObjectWavecal = @showprogress desc=desc pmap(
        get_and_save_sky_wavecal_partial, all1DObjectSkyPeaks)
    good_fnames = .!isnothing.(all1DObjectWavecal)
    all1DObjectWavecal = convert(Vector{String}, all1DObjectWavecal[good_fnames])
    all1DObjectSkyPeaks = all1DObjectSkyPeaks[good_fnames]

    # putting this parallelized within each mjd is really not good in the bulk run context
    mjd_list_wavecal = map(x -> parse(Int, split(basename(x), "_")[3]), all1DObjectWavecal)

    sendto(workers(), mjd_list_wavecal = mjd_list_wavecal)
    sendto(workers(), all1DObjectWavecal = all1DObjectWavecal)
    sendto(workers(), all1DObjectSkyPeaks = all1DObjectSkyPeaks)

    # the almanac provides the commanded dither positions (DITHPIX labels) for
    # the deterministic night-frame anchor, and the FPI wavecal gate's guide
    # fibers (they are not a constant: LCO MJD 59810-59850 used 142/153)
    almanac_file = joinpath(parg["outdir"], "almanac/$(parg["runname"]).h5")
    sendto(workers(), almanac_file = almanac_file)

    desc = "Skyline medwave/skyline dither: "
    @everywhere skyline_medwavecal_skyline_dither_partial(mjd) = skyline_medwavecal_skyline_dither(
        parg["tele"], mjd, mjd_list_wavecal, all1DObjectWavecal,
	all1DObjectSkyPeaks, parg["checkpoint_mode"]; outdir = parg["outdir"],
        almanac_file = almanac_file)
    wavecalNightAve_fnames = @showprogress desc=desc pmap(skyline_medwavecal_skyline_dither_partial, unique_mjds)
    
    sendto(workers(), wavecalNightAve_fnames = wavecalNightAve_fnames)
    sendto(workers(), unique_mjds = unique_mjds)
    unique_mjd_inds = collect(1:size(unique_mjds,1))

    # the FPI/arclamp version of wavecal is still a TODO from Kevin McKinnon
    if size(all1DFPI, 1) > 0
        mjd_list_fpi = map(x -> parse(Int, split(basename(x), "_")[3]), all1DfpiPeaks_a)
        desc = "FPI medwave/skyline dither: "
        sendto(workers(), mjd_list_fpi = mjd_list_fpi)
        sendto(workers(), all1DfpiPeaks_a = all1DfpiPeaks_a)
        sendto(workers(), all1DObjectSkyPeaks = all1DObjectSkyPeaks)
        # almanac_file (hoisted above the skyline night average) is needed by
        # the FPI wavecal gate to derive the night's FPI guide fibers
        @everywhere fpi_medwavecal_skyline_dither_partial(mjd_ind) = fpi_medwavecal_skyline_dither(
            unique_mjds[mjd_ind], mjd_list_fpi, mjd_list_wavecal, all1DfpiPeaks_a, all1DObjectSkyPeaks,
            wavecalNightAve_fnames[mjd_ind], verbose = false, checkpoint_mode = parg["checkpoint_mode"],
            almanac_file = almanac_file)
        @showprogress desc=desc pmap(fpi_medwavecal_skyline_dither_partial, unique_mjd_inds)
    end

    ## stage-1 arc-lamp night product (off by default; nothing downstream reads it yet).
    ## It needs the arc ar1Dcal files and, for the sequence gate, the per-exposure sky/FPI
    ## wavelength solutions written above.
    if parg["arclamp_wavecal"] && size(all1DArclamp, 1) > 0
        almanac_file_arc = joinpath(parg["outdir"], "almanac/$(parg["runname"]).h5")
        sendto(workers(), almanac_file_arc = almanac_file_arc)
        @everywhere arclamp_linelist = load_arclamp_linelist(
            parg["arclamp_linelist"] == "" ? ApogeeReduction.ARCLAMP_LINELIST_DEFAULT :
            parg["arclamp_linelist"]; selection = parg["arclamp_linelist_selection"])
        @everywhere arclamp_night_partial(mjd) = try
            get_and_save_arclamp_night_wavecal(parg["tele"], mjd;
                apred_dir = joinpath(parg["outdir"], "apred"), almanac_file = almanac_file_arc,
                linelist = arclamp_linelist, checkpoint_mode = parg["checkpoint_mode"])
        catch e
            @warn "arc-lamp night product failed for $(parg["tele"]) $(mjd)" exception = e
            nothing
        end
        desc = "Arc-lamp night products (stage 1): "
        @showprogress desc=desc pmap(arclamp_night_partial, unique_mjds)
    end

    # ## TODO when are we going to split into individual fiber files? Then we should be writing fiber type to the file name

    ## combine chips for single exposure onto loguniform wavelength grid
    ## pushing off the question of dither combinations for now (to apMADGICS stage)
    @everywhere reinterp_spectra_partial(fname) = reinterp_spectra(
        fname, roughwave_dict, 
        checkpoint_mode = parg["checkpoint_mode"], outdir = parg["outdir"],
        wavecal_order = parg["wavecal_order"])
    if parg["doUncals"]
        all1Da = replace.(all2Dperchip[1], "ar2D" => "ar1D")
        desc = "Reinterp exposure spectra (uncals):"
        @showprogress desc=desc pmap(reinterp_spectra_partial, all1Da)
    end

    # very cyclic CPU usage? MEM log-jam? Some slow step I didn't notice?
    all1Da = replace.(all2Dperchip[1], "ar2D" => "ar1Dcal")
    desc = "Reinterp exposure spectra:"
    @showprogress desc=desc pmap(reinterp_spectra_partial, all1Da)
end
