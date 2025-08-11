## This is a reduction pipeline for APOGEE
import Pkg;
using Dates;
t0 = now();
t_then = t0;
using InteractiveUtils;
versioninfo();
Pkg.instantiate();
Pkg.precompile(); # no need for Pkg.activate("./") because of invocation w/ environment

using Distributed, ArgParse, TimerOutputs
t_now = now();
dt = Dates.canonicalize(Dates.CompoundPeriod(t_now - t_then));
println("Package activation took $dt");
t_then = t_now;
flush(stdout);

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
        default = 32
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
if parg["runlist"] != "" # only multiprocess if we have a list of exposures
    if "SLURM_NTASKS" in keys(ENV)
        using SlurmClusterManager
        addprocs(SlurmManager(), exeflags = ["--project=$proj_path"])
        ntasks = parse(Int, ENV["SLURM_NTASKS"])
        nnodes = ntasks ÷ 64  # Each node has 64 cores
        total_workers = nnodes * workers_per_node
        workers_to_keep = []
        for node in 0:(nnodes - 1)
            node_start = 1 + node * 64
            spacing = 64 ÷ workers_per_node
            append!(workers_to_keep, [node_start + spacing * i for i in 0:(workers_per_node - 1)])
        end
        rmprocs(setdiff(1:ntasks, workers_to_keep))
    else
        addprocs(workers_per_node, exeflags = ["--project=$proj_path"])
    end
end
t_now = now();
dt = Dates.canonicalize(Dates.CompoundPeriod(t_now - t_then));
println("Allocating $(nworkers()) workers took $dt");
t_then = t_now;
flush(stdout);
println("Running Main on ", gethostname());
flush(stdout);

@everywhere begin
    using LinearAlgebra
    BLAS.set_num_threads(1)
    using FITSIO, HDF5, FileIO, JLD2, Glob, CSV
    using DataFrames, EllipsisNotation, StatsBase
    using AstroTime # can remove after Adam merges the PR to recast as Float
    using ParallelDataTransfer, ProgressMeter
    using ApogeeReduction
    using ApogeeReduction: read_almanac_exp_df, gh_profiles, read_metadata, regularize_trace,
                           extract_boxcar, extract_optimal_iter, safe_jldsave, get_fluxing_file,
                           get_fibTargDict, get_1d_name, get_and_save_sky_wavecal,
                           get_and_save_sky_dither_per_fiber, get_and_save_sky_peaks,
                           get_ave_night_wave_soln, sky_wave_plots, reinterp_spectra,
                           get_and_save_arclamp_peaks, get_and_save_fpi_peaks,
                           comb_exp_get_and_save_fpi_wavecal

    include("src/makie_plotutils.jl")

    ###decide which type of cal to use for traces (i.e. dome or quartz flats)
    # trace_type = "dome"
    trace_type = "quartz" #probably should use this!
end

println(BLAS.get_config());
flush(stdout);

##### 1D stage
@everywhere begin
    function process_1D(fname)
        sname = split(split(split(fname, "/")[end], ".h5")[1], "_")
        fnameType, tele, mjd, expnum, chip, exptype = sname[(end - 5):end]

        # how worried should I be about loading this every time?
        falm = h5open(joinpath(parg["outdir"], "almanac/$(parg["runname"]).h5"))
        dfalmanac = read_almanac_exp_df(falm, parg["tele"], mjd)

        med_center_to_fiber_func, x_prof_min, x_prof_max_ind, n_sub, min_prof_fib, max_prof_fib,
        all_y_prof, all_y_prof_deriv = gh_profiles(tele, mjd, expnum, chip; n_sub = 100)

        fnamecal = if (fnameType == "ar2D")
            replace(fname, "ar2D" => "ar2Dcal")
        else
            fname
        end

        dimage = load(fname, "dimage")
        ivarimage = load(fname, "ivarimage")
        pix_bitmask = load(fnamecal, "pix_bitmask")
        metadata = read_metadata(fname)

        # this seems annoying to load so often if we know we are doing a daily... need to ponder
        traceList = sort(glob("$(trace_type)TraceMain_$(tele)_$(mjd)_*_$(chip).h5",
            parg["outdir"] * "apred/$(mjd)/"))
        trace_params = load(traceList[1], "trace_params")

        # adam: should this be saved somewhere?  It's fairly simple to reproduce, but that's true of
        # everything to some degree
        regularized_trace_params = regularize_trace(trace_params)

        flux_1d, ivar_1d,
        mask_1d,
        dropped_pixels_mask_1d,
        resid_flux,
        resid_ivar = if parg["extraction"] == "boxcar"
            extract_boxcar(
                dimage, ivarimage, pix_bitmask, regularized_trace_params, return_resids = true)
        elseif parg["extraction"] == "optimal"
            #            extract_optimal(dimage, ivarimage, pix_bitmask, regularized_trace_params)
            extract_optimal_iter(dimage, ivarimage, pix_bitmask, regularized_trace_params,
                med_center_to_fiber_func, x_prof_min, x_prof_max_ind, n_sub,
                min_prof_fib, max_prof_fib, all_y_prof, all_y_prof_deriv, return_resids = true)
        else
            error("Extraction method $(parg["extraction"]) not recognized")
        end

        outfname = replace(fname, "ar2D" => "ar1D")
        resid_outfname = replace(fname, "ar2D" => "ar2Dresiduals")
        safe_jldsave(resid_outfname, metadata; resid_flux, resid_ivar)
        if parg["relFlux"]
            # relative fluxing (using "B" only for now)
            # this is the path to the underlying fluxing file.
            # it is symlinked below to an exposure-specific file (linkPath).
            calPath = get_fluxing_file(
                dfalmanac, parg["outdir"], tele, mjd, expnum, fluxing_chip = CHIP_LIST[end])
            expid_num = parse(Int, last(expnum, 4)) #this is silly because we translate right back
            fibtargDict = get_fibTargDict(falm, tele, mjd, expid_num)
            fiberTypeList = map(x -> fibtargDict[x], 1:300)

            if isnothing(calPath)
                @warn "No fluxing file available for $(tele) $(mjd) $(expnum) $(chip)"
                relthrpt = ones(size(flux_1d, 2))
                bitmsk_relthrpt = 2^2 * ones(Int, size(flux_1d, 2))
            elseif !isfile(calPath)
                error("Fluxing file $(calPath) for $(tele) $(mjd) $(expnum) $(chip) does not exist")
            else
                linkPath = abspath(joinpath(
                    dirname(fname), "relFlux_$(tele)_$(mjd)_$(expnum)_$(chip).h5"))
                if !islink(linkPath)
                    symlink(abspath(calPath), linkPath)
                end
                relthrpt = load(linkPath, "relthrpt")
                relthrptr = reshape(relthrpt, (1, length(relthrpt)))
                bitmsk_relthrpt = load(linkPath, "bitmsk_relthrpt")
            end

            # don't flux broken fibers (don't use warn fibers for sky)
            msk_goodwarn = (bitmsk_relthrpt .== 0) .| (bitmsk_relthrpt .& 2^0) .== 2^0
            if any(msk_goodwarn)
                flux_1d[:, msk_goodwarn] ./= relthrptr[:, msk_goodwarn]
                ivar_1d[:, msk_goodwarn] .*= relthrptr[:, msk_goodwarn] .^ 2
            end

            # we probably want to append info from the fiber dictionary from alamanac into the file name
            safe_jldsave(outfname, metadata; flux_1d, ivar_1d, mask_1d, dropped_pixels_mask_1d,
                extract_trace_centers = regularized_trace_params[:, :, 2],
                relthrpt, bitmsk_relthrpt, fiberTypeList)
        else
            outfname_norelflux = replace(outfname, "apred" => "apredrelflux")
            dirName = dirname(outfname_norelflux)
            if !ispath(dirName)
                mkpath(dirName)
            end
            safe_jldsave(
                outfname_norelflux, metadata; flux_1d, ivar_1d, mask_1d, dropped_pixels_mask_1d,
                extract_trace_centers = regularized_trace_params[:, :, 2])
        end
        close(falm)
    end
end
t_now = now();
dt = Dates.canonicalize(Dates.CompoundPeriod(t_now - t_then));
println("Function definitions took $dt");
t_then = t_now;
flush(stdout);

@passobj 1 workers() parg

# Find the 2D calibration files for the relevant MJDs
tele_list = if parg["runlist"] != ""
    load(parg["runlist"], "tele")
else
    [parg["tele"]]
end
unique_teles = unique(tele_list)
mskTele = tele_list .== parg["tele"]

mjd_list = if parg["runlist"] != ""
    load(parg["runlist"], "mjd")
else
    [parg["mjd"]]
end
unique_mjds = unique(mjd_list[mskTele])

# make file name list
expid_list = if parg["runlist"] != ""
    load(parg["runlist"], "expid")
else
    [parg["expid"]]
end

# need to be building a msk on the expid for the different MJDs
# or a list of lists

list2Dexp = []
for mjd in unique_mjds
    df = read_almanac_exp_df(
        joinpath(parg["outdir"], "almanac/$(parg["runname"]).h5"), parg["tele"], mjd)
    function get_2d_name_partial(expid)
        parg["outdir"] * "/apred/$(mjd)/" *
        replace(get_1d_name(expid, df), "ar1D" => "ar2D") * ".h5"
    end
    mskMJD = (mjd_list .== mjd) .& mskTele
    local2D = get_2d_name_partial.(expid_list[mskMJD])
    push!(list2Dexp, local2D)
end
all2Da = vcat(list2Dexp...)

all2Dperchip = []
for chip in CHIP_LIST
    all2Dchip = replace.(all2Da, "_$(FIRST_CHIP)_" => "_$(chip)_")
    push!(all2Dperchip, all2Dchip)
end
all2D = vcat(all2Dperchip...)

# we should do somthing smart to assemble the traces from a night into a single file
# that gives us the trace of a fiber as a funciton of time or something
# for now, for each MJD, take the first one (or do that in run_trace_cal.sh)
# I think dome flats needs to swtich to dome_flats/mjd/
for mjd in unique_mjds
    for chip in CHIP_LIST
        traceList = sort(glob("$(trace_type)Trace_$(parg["tele"])_$(mjd)_*_$(chip).h5",
            parg["outdir"] * "$(trace_type)_flats/"))
        if length(traceList) > 1
            @warn "Multiple $(trace_type) trace files found for $(parg["tele"]) $(mjd) $(chip): $(traceList)"
        elseif length(traceList) == 0
            @error "No $(trace_type) trace files found for $(parg["tele"]) $(mjd) $(chip). Looked in $(parg["outdir"] * "$(trace_type)_flats/")."
        end
        calPath = traceList[1]
        linkPath = parg["outdir"] * "apred/$(mjd)/" *
                   replace(basename(calPath), "$(trace_type)Trace" => "$(trace_type)TraceMain")
        if !isfile(linkPath)
            # come back to why this symlink does not work
            cp(calPath, linkPath)
        end
    end
end

# extract the 2D to 1D, ideally the calibrated files
# need to think hard about batching daily versus all data for cal load in
# someday we might stop doing the uncal extractions, but very useful for testing

###### HELLO A HACK HERE WE NEED TO REMOVE
# desc = "Extracting 2D to 1D"
# @showprogress desc=desc pmap(process_1D, all2D)

# desc = "Extracting 2Dcal to 1Dcal"
# all2Dcal = replace.(all2D, "ar2D" => "ar2Dcal")
# @showprogress desc=desc pmap(process_1D, all2Dcal)

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
        local1D = get_1d_name_partial.(expid_list[mskMJD])
        push!(list1DexpObject, filter(!isnothing, local1D))
        local1D_fpi = get_1d_name_FPI_partial.(expid_list[mskMJD])
        push!(list1DexpFPI, filter(!isnothing, local1D_fpi))
        local1D_arclamp = get_1d_name_ARCLAMP_partial.(expid_list[mskMJD])
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
        roughwave_dict = load("data/roughwave_dict.jld2", "roughwave_dict")
        df_sky_lines = CSV.read("data/APOGEE_lines.csv", DataFrame)
        df_sky_lines.linindx = 1:size(df_sky_lines, 1)
    end

    ##### HARD HACK HERE, MUST REMOVE
    # ## get sky line peaks
    # @everywhere get_and_save_sky_peaks_partial(fname) = get_and_save_sky_peaks(
    #     fname, roughwave_dict, df_sky_lines)
    # desc = "Fitting sky line peaks: "
    # @showprogress desc=desc pmap(get_and_save_sky_peaks_partial, all1DObject)

    # # get arclamp peaks
    # if size(all1DArclamp, 1) > 0
    #     ## get (non-fpi) arclamp peaks
    #     desc = "Fitting arclamp peaks: "
    #     @showprogress desc=desc pmap(get_and_save_arclamp_peaks, all1DArclamp)
    # end

    all1DfpiPeaks_a = replace.(
        replace.(all1DFPIa, "ar1Dcal" => "fpiPeaks"), "ar1D" => "fpiPeaks")
    all1DfpiPeaks = if size(all1DFPI, 1) > 0
        ## get FPI peaks
        desc = "Fitting FPI peaks: "
        @showprogress desc=desc pmap(get_and_save_fpi_peaks, all1DFPI)
    else
        []
    end
    all1DfpiPeaks_out = reshape(all1DfpiPeaks,length(all1DfpiPeaks_a),length(CHIP_LIST))
    mskFPInothing = .!any.(isnothing.(eachrow(all1DfpiPeaks_out)))
    println("FPI peaks found for $(sum(mskFPInothing)) of $(length(mskFPInothing)) exposures")

    ## get wavecal from sky line peaks
    println("Solving skyline wavelength solution:")
    flush(stdout)
    #only need to give one chip's list because internal
    #logic handles finding other chips when ingesting data
    #Andrew says that that is a bit worrisome and would should revisit that logic
    all1DObjectSkyPeaks = replace.(
        replace.(all1DObjecta, "ar1Dcal" => "skyLinePeaks"), "ar1D" => "skyLinePeaks")
    all1DObjectWavecal = @showprogress pmap(get_and_save_sky_wavecal, all1DObjectSkyPeaks)
    all1DObjectWavecal = filter(x -> !isnothing(x), all1DObjectWavecal)

    # putting this parallelized within each mjd is really not good in the bulk run context
    night_wave_soln_dict = Dict{Int, Any}()
    night_linParams_dict = Dict{Int, Any}()
    night_nlParams_dict = Dict{Int, Any}()
    mjd_list_wavecal = map(x -> parse(Int, split(basename(x), "_")[3]), all1DObjectWavecal)
    for mjd in unique_mjds
        mskMJD = (mjd_list_wavecal .== mjd)
        if size(all1DObjectWavecal[mskMJD], 1) > 0
            println("Using all skyline wavelength solutions to determine median solution for MJD $mjd.")
            flush(stdout)
            all1DObjectWavecal_mjd = all1DObjectWavecal[mskMJD]
            night_linParams, night_nlParams, night_wave_soln = get_ave_night_wave_soln(
                all1DObjectWavecal_mjd, fit_dither = true)
            night_wave_soln_dict[mjd] = night_wave_soln
            night_linParams_dict[mjd] = night_linParams
            night_nlParams_dict[mjd] = night_nlParams
            sendto(workers(), night_wave_soln = night_wave_soln)
            sendto(workers(), night_linParams = night_linParams)
            sendto(workers(), night_nlParams = night_nlParams)

            println("Using skylines to measure dither offsets from nightly skyline average wavelength solution for MJD $mjd.")
            @everywhere get_and_save_sky_dither_per_fiber_partial(fname) = get_and_save_sky_dither_per_fiber(
                fname, night_linParams, night_nlParams; dporder = 1, wavetype = "sky", max_offset = 1.0)

            @showprogress pmap(
                get_and_save_sky_dither_per_fiber_partial, all1DObjectSkyPeaks[mskMJD])

            println("Plotting skyline wavelength solution diagnostic figures.")

            # I am very worried about the lack of parallelization here.
            sky_wave_plots(
                all1DObjectWavecal[mskMJD], night_linParams, night_nlParams, night_wave_soln,
                dirNamePlots = dirNamePlots,
                plot_fibers = (1, 50, 100, 150, 200, 250, 300),
                plot_pixels = (1, 512, 1024, 1536, 2048))
        else
            night_wave_soln_dict[mjd] = nothing
            night_linParams_dict[mjd] = nothing
            night_nlParams_dict[mjd] = nothing
        end
    end
    sendto(workers(), night_wave_soln_dict = night_wave_soln_dict)
    sendto(workers(), night_linParams_dict = night_linParams_dict)
    sendto(workers(), night_nlParams_dict = night_nlParams_dict)

    if size(all1DFPI, 1) > 0
        mjd_list_fpi = map(x -> parse(Int, split(basename(x), "_")[3]), all1DfpiPeaks_a)
        #change the condition once there are wavelength 
        #solutions from the ARCLAMPs as well
        for mjd in unique_mjds
            mskMJD_fpi = (mjd_list_fpi .== mjd)
            mskMJD_obj = (mjd_list_wavecal .== mjd)
            if (!isnothing(night_linParams_dict[mjd])) & (size(all1DfpiPeaks_a[mskMJD_fpi], 1) > 0)
                println("Using $(size(all1DfpiPeaks_a[mskMJD_fpi],1)) FPI exposures to measure high-precision nightly wavelength solution")
                outfname, night_linParams, night_nlParams, night_wave_soln = comb_exp_get_and_save_fpi_wavecal(
                    all1DfpiPeaks_a[mskMJD_fpi], night_linParams_dict[mjd], night_nlParams_dict[mjd], cporder = 1, wporder = 4,
                    dporder = 2, n_sigma = 4, max_ang_sigma = 0.2, max_iter = 2)
                sendto(workers(), night_wave_soln = night_wave_soln)
                sendto(workers(), night_linParams = night_linParams)
                sendto(workers(), night_nlParams = night_nlParams)

                println("Using skylines to measure dither offsets from FPI-defined wavelength solution")
                @everywhere get_and_save_sky_dither_per_fiber_partial(fname) = get_and_save_sky_dither_per_fiber(
                    fname, night_linParams, night_nlParams; dporder = 2,
                    wavetype = "fpi", max_offset = 1.0)

                @showprogress pmap(
                    get_and_save_sky_dither_per_fiber_partial, all1DObjectSkyPeaks[mskMJD_obj])
            end
        end
    end

    ## TODO when are we going to split into individual fiber files? Then we should be writing fiber type to the file name
    ## combine chips for single exposure onto loguniform wavelength grid
    ## pushing off the question of dither combinations for now (to apMADGICS stage)
    all1Da = replace.(all2Dperchip[1], "ar2D" => "ar1D")
    println("Reinterpolating exposure spectra:")
    flush(stdout)
    @everywhere reinterp_spectra_partial(fname) = reinterp_spectra(
        fname, roughwave_dict, backupWaveSoln = night_wave_soln)
    @showprogress pmap(reinterp_spectra_partial, all1Da)

    all1Da = replace.(all2Dperchip[1], "ar2D" => "ar1Dcal")
    println("Reinterpolating calibrated exposure spectra:")
    flush(stdout)
    @showprogress pmap(reinterp_spectra_partial, all1Da)
end