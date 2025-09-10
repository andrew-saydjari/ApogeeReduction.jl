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
                           safe_jldsave, process_1D, check_file

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

    list2Dexp = []
    for mjd in unique_mjds
        df = read_almanac_exp_df(
            joinpath(parg["outdir"], "almanac/$(parg["runname"]).h5"), parg["tele"], mjd)
        function get_2d_name_partial(expid)
            parg["outdir"] * "/apred/$(mjd)/" *
            replace(get_1d_name(expid, df), "ar1D" => "ar2D") * ".h5"
        end
        mskMJD = (mjd_list .== mjd) .& mskTele
        local2D = get_2d_name_partial.(dfindx_list[mskMJD])
        push!(list2Dexp, local2D)
    end
    all2Da = vcat(list2Dexp...)

    all2Dperchip = []
    for chip in CHIP_LIST
        all2Dchip = replace.(all2Da, "_$(FIRST_CHIP)_" => "_$(chip)_")
        push!(all2Dperchip, all2Dchip)
    end
    all2D = vcat(all2Dperchip...)
end

# make trace param file copies for each date
# searching outwards from the given mjd, 
# looking for cal types in the order of trace_type_order
@everywhere function create_traceMain(mjd,chip,trace_type_order;
				       max_mjd_offset = 7,
				       tele = parg["tele"],
				       outdir = parg["outdir"],
				       clobber = true)

    new_trace_param_fname = outdir * "apred/$(mjd)/traceMain_$(tele)_(mjd)_$(chip).h5"
    if isfile(new_trace_param_fname) & (!clobber)
        curr_trace_type,curr_mjd_offset,found_first_choice,found_on_mjd,found_match,orig_trace_param_fname = h5open(new_trace_param_fname, "r") do f
            (f["trace_type"],f["trace_mjd_offset"],f["trace_found_best_choice"],f["trace_match_mjd"],f["trace_found_match"],f["trace_orig_param_fname"])
        end

        return curr_trace_type,curr_mjd_offset,found_first_choice,found_on_mjd,found_match,orig_trace_param_fname,new_trace_param_fname
    end
    orig_trace_param_fname = nothing

    curr_trace_type = trace_type_order[1]
    curr_mjd_offset = 0 # number of mjd that trace parameters are away from given mjd
    found_first_choice = true # ie same mjd, using first entry of trace_type_order
    found_on_mjd = true # ie same mjd
    found_match = false

    # best case scenario will have:
    #    curr_trace_type = trace_type_order[1]
    #    curr_mjd_offset = 0
    #    found_first_choice = true
    #    found_on_mjd = true
    #    found_match = true

    #work outwards from current MJD
    check_mjd_offsets = zeros(Int,max_mjd_offset*2+1) 
    check_mjd_offsets[begin+1:2:end-1] .= collect(-1:-1:-max_mjd_offset)
    check_mjd_offsets[begin+2:2:end] .= collect(1:1:max_mjd_offset)

    for mjd_ind in 1:length(check_mjd_offsets)
        curr_mjd_offset = check_mjd_offsets[mjd_ind]
        check_mjd = mjd + curr_mjd_offset

        for check_trace_type in trace_type_order
            # check for different cal types in preferential
            # order given by trace_type_order
            curr_trace_type = check_trace_type

            traceList = sort(glob("$(check_trace_type)Trace_$(tele)_$(check_mjd)_*_$(chip).h5",
                outdir * "$(check_trace_type)_flats/$(check_mjd)/"))
            if length(traceList) >= 1
                found_match = true
                orig_trace_param_fname = traceList[1]
            end
        end

        if found_match
            break
        end
    end

    if !found_match
#        @error "No trace files found for $(tele) $(mjd) $(chip). Looked in $(outdir) using +/- $(max_mjd_offset) days around given MJD looking for trace types of $(trace_type_order)."
        @warn "No trace files found for $(tele) $(mjd) $(chip). Looked in $(outdir) using +/- $(max_mjd_offset) days around given MJD looking for trace types of $(trace_type_order)."
        #POSSIBLE TODO: have some fallback trace parameter file 
        return nothing
    end

    # TODO (probably for KM): 
    #     -use all the night's flats (quartz and dome) to define
    #          the time-evolution of trace params, then save
    #          that output (likely using this function)

    found_on_mjd = (curr_mjd_offset == 0)
    found_first_choice = found_on_mjd & (curr_trace_type == trace_type_order[1])
    
    if (!isfile(new_trace_param_fname)) | clobber
        # come back to why this symlink does not work
        # is this causing memory bloat?
        cp(orig_trace_param_fname, new_trace_param_fname)

	#KM says: this should be okay because ultimately we will
	#be using information about how the trace parameters change as a function
        #of time over a night, whose outputs we will want to save
        #(ie final version will produce a new trace param file per tele/mjd/chip combo) 

        #add in some additional information to final trace file
        h5open(new_trace_param_fname, "r+") do f
            write(f, "trace_type", curr_trace_type)
            write(f, "trace_mjd_offset", curr_mjd_offset)
            write(f, "trace_found_best_choice", found_first_choice)
            write(f, "trace_match_mjd", found_on_mjd)
            write(f, "trace_found_match", found_match)
            write(f, "trace_orig_param_fname", orig_trace_param_fname)
        end
    end

    return curr_trace_type,curr_mjd_offset,found_first_choice,found_on_mjd,found_match,orig_trace_param_fname,new_trace_param_fname
end

# we should do somthing smart to assemble the traces from a night into a single file
# that gives us the trace of a fiber as a funciton of time or something
# for now, for each MJD, take the first one (or do that in run_trace_cal.sh)
# I think dome flats needs to swtich to dome_flats/mjd/
@time "Processing trace files" for mjd in unique_mjds
    for chip in CHIP_LIST
        create_traceMain(mjd,chip,trace_type_order;
#                         max_mjd_offset = 7,
                         max_mjd_offset = 100000, #setting to a large number for now, to ensure traceMain files exist
                         tele = parg["tele"],
                         outdir = parg["outdir"],
                         clobber = false)
    end
end

# extract the 2D to 1D
@everywhere process_1D_wrapper(fname) = process_1D(
    fname,
    outdir = parg["outdir"],
    runname = parg["runname"],
    extraction = parg["extraction"],
    relFlux = parg["relFlux"],
    chip_list = CHIP_LIST,
    profile_path = joinpath(proj_path, "data"),
    plot_path = joinpath(parg["outdir"], "plots/")
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
    for mjd in unique_mjds
        df = read_almanac_exp_df(
            joinpath(parg["outdir"], "almanac/$(parg["runname"]).h5"), parg["tele"], mjd)
        function get_1d_name_partial(expid)
            if df.imagetyp[expid] == "Object"
                return parg["outdir"] * "/apred/$(mjd)/" * get_1d_name(expid, df, cal = true) *
                       ".h5"
            else
                return nothing
            end
        end
        function get_1d_name_ARCLAMP_partial(expid)
            if (df.imagetyp[expid] == "ArcLamp") &
               ((df.lampthar[expid] == "T") | (df.lampune[expid] == "T"))
                return parg["outdir"] * "/apred/$(mjd)/" * get_1d_name(expid, df, cal = true) *
                       ".h5"
            else
                return nothing
            end
        end
        function get_1d_name_FPI_partial(expid)
            if (df.imagetyp[expid] == "ArcLamp") & (df.lampthar[expid] == "F") &
               (df.lampune[expid] == "F")
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

    ## get wavecal from sky line peaks
    #only need to give one chip's list because internal
    #logic handles finding other chips when ingesting data
    #Andrew says that that is a bit worrisome and would should revisit that logic
    all1DObjectSkyPeaks = replace.(
        replace.(all1DObjecta, "ar1Dcal" => "skyLinePeaks"), "ar1D" => "skyLinePeaks")
    desc = "Skyline wavelength solutions:"
    @everywhere get_and_save_sky_wavecal_partial(fname) = get_and_save_sky_wavecal(
        fname, checkpoint_mode = parg["checkpoint_mode"])
    all1DObjectWavecal = @showprogress desc=desc pmap(
        get_and_save_sky_wavecal_partial, all1DObjectSkyPeaks)
    all1DObjectWavecal = filter(!isnothing, all1DObjectWavecal)

    # putting this parallelized within each mjd is really not good in the bulk run context
    mjd_list_wavecal = map(x -> parse(Int, split(basename(x), "_")[3]), all1DObjectWavecal)

    sendto(workers(), mjd_list_wavecal = mjd_list_wavecal)
    sendto(workers(), all1DObjectWavecal = all1DObjectWavecal)
    sendto(workers(), all1DObjectSkyPeaks = all1DObjectSkyPeaks)

    @everywhere outname = joinpath(
        parg["outdir"], "wavecal", "skyline_wavecal_$(parg["tele"])_$(parg["runname"])_dict.jld2")
    if !check_file(outname, mode = parg["checkpoint_mode"])
        desc = "Skyline medwave/skyline dither: "
        @everywhere skyline_medwavecal_skyline_dither_partial(mjd) = skyline_medwavecal_skyline_dither(
            mjd, mjd_list_wavecal, all1DObjectWavecal, all1DObjectSkyPeaks; outdir = parg["outdir"])
        pout = @showprogress desc=desc pmap(skyline_medwavecal_skyline_dither_partial, unique_mjds)

        night_wave_soln_dict = Dict(unique_mjds .=> map(x -> x[1], pout))
        night_linParams_dict = Dict(unique_mjds .=> map(x -> x[2], pout))
        night_nlParams_dict = Dict(unique_mjds .=> map(x -> x[3], pout))

        mkpath(dirname(outname))
        safe_jldsave(outname, night_wave_soln_dict = night_wave_soln_dict,
            night_linParams_dict = night_linParams_dict,
            night_nlParams_dict = night_nlParams_dict, no_metadata = true)
    end
    @everywhere begin
        night_wave_soln_dict, night_linParams_dict, night_nlParams_dict = load(
            outname, "night_wave_soln_dict", "night_linParams_dict", "night_nlParams_dict")
    end

    # the FPI/arclamp version of wavecal is still a TODO from Kevin McKinnon
    if size(all1DFPI, 1) > 0
        mjd_list_fpi = map(x -> parse(Int, split(basename(x), "_")[3]), all1DfpiPeaks_a)
        desc = "FPI medwave/skyline dither: "
        sendto(workers(), mjd_list_fpi = mjd_list_fpi)
        sendto(workers(), all1DfpiPeaks_a = all1DfpiPeaks_a)
        sendto(workers(), all1DObjectSkyPeaks = all1DObjectSkyPeaks)
        @everywhere fpi_medwavecal_skyline_dither_partial(mjd) = fpi_medwavecal_skyline_dither(
            mjd, mjd_list_fpi, mjd_list_wavecal, all1DfpiPeaks_a, all1DObjectSkyPeaks,
            night_linParams_dict, night_nlParams_dict, checkpoint_mode = parg["checkpoint_mode"])
        @showprogress desc=desc pmap(fpi_medwavecal_skyline_dither_partial, unique_mjds)
    end

    ## TODO when are we going to split into individual fiber files? Then we should be writing fiber type to the file name

    ## combine chips for single exposure onto loguniform wavelength grid
    ## pushing off the question of dither combinations for now (to apMADGICS stage)
    @everywhere reinterp_spectra_partial(fname) = reinterp_spectra(
        fname, roughwave_dict, backupWaveSoln = night_wave_soln_dict,
        checkpoint_mode = parg["checkpoint_mode"])
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
