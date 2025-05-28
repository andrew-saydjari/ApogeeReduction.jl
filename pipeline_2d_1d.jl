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
    end
    return parse_args(s)
end

parg = parse_commandline()
workers_per_node = parg["workers_per_node"]
if parg["runlist"] != "" # only multiprocess if we have a list of exposures
    if "SLURM_NTASKS" in keys(ENV)
        using SlurmClusterManager
        addprocs(SlurmManager(), exeflags = ["--project=./"])
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
        addprocs(workers_per_node, exeflags = ["--project=./"])
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

    src_dir = "./"
    include(src_dir * "src/ar1D.jl")
    include(src_dir * "src/spectraInterpolation.jl")
    include(src_dir * "src/fileNameHandling.jl")
    include(src_dir * "src/utils.jl")
    include(src_dir * "src/skyline_peaks.jl")
    include(src_dir * "src/wavecal.jl")
    include(src_dir * "src/cal_build/traceExtract_GH.jl")

    ###decide which type of cal to use for traces (i.e. dome or quartz flats)
    # trace_type = "dome"
    trace_type = "quartz" #probably should use this!
end

println(BLAS.get_config());
flush(stdout);

##### 1D stage
@everywhere begin
    function process_1D(fname)
        sname = split(split(split(fname, "/")[end],".h5")[1], "_")
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
            # relative fluxing (using "c" only for now)
            # this is the path to the underlying fluxing file.
            # it is symlinked below to an exposure-specific file (linkPath).
            calPath = get_fluxing_file(
                dfalmanac, parg["outdir"], tele, mjd, expnum, fluxing_chip = "c")
            expid_num = parse(Int, last(expnum, 4)) #this is silly because we translate right back
            fibtargDict = get_fibTargDict(falm, tele, mjd, expid_num)
            fiberTypeList = map(x -> fibtargDict[x], 1:300)

            if isnothing(calPath)
                # TODO uncomment this
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
            safe_jldsave(outfname, metadata; flux_1d, ivar_1d, mask_1d,
                relthrpt, bitmsk_relthrpt, fiberTypeList)
        else
            safe_jldsave(outfname, metadata; flux_1d, ivar_1d, mask_1d)
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
unique_mjds = if parg["runlist"] != ""
    subDic = load(parg["runlist"])
    unique(subDic["mjd"])
else
    [parg["mjd"]]
end

# make file name list
expid_list = if parg["runlist"] != ""
    subDic = load(parg["runlist"])
    subDic["expid"]
else
    [parg["expid"]]
end

list2Dexp = []
for mjd in unique_mjds
    df = read_almanac_exp_df(joinpath(parg["outdir"], "almanac/$(parg["runname"]).h5"), parg["tele"], mjd)
    function get_2d_name_partial(expid)
        parg["outdir"] * "/apred/$(mjd)/" *
        replace(get_1d_name(expid, df), "ar1D" => "ar2D") * ".h5"
    end
    local2D = get_2d_name_partial.(expid_list)
    push!(list2Dexp, local2D)
end
all2Da = vcat(list2Dexp...)

all2Dperchip = []
for chip in CHIP_LST
    all2Dchip = replace.(all2Da, "_$(CHIP_LST[1])_" => "_$(chip)_")
    push!(all2Dperchip, all2Dchip)
end
all2D = vcat(all2Dperchip...)

# we should do somthing smart to assemble the traces from a night into a single file
# that gives us the trace of a fiber as a funciton of time or something
# for now, for each MJD, take the first one (or do that in run_trace_cal.sh)
# I think dome flats needs to swtich to dome_flats/mjd/
for mjd in unique_mjds
    for chip in CHIP_LST
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
println("Extracting 2D to 1D:");
flush(stdout);
@showprogress pmap(process_1D, all2D)
println("Extracting 2Dcal to 1Dcal:");
flush(stdout);
all2Dcal = replace.(all2D, "ar2D" => "ar2Dcal")
@showprogress pmap(process_1D, all2Dcal)

## get all OBJECT files (happy to add any other types that see sky?)
list1DexpObject = []
for mjd in unique_mjds
    df = read_almanac_exp_df(joinpath(parg["outdir"], "almanac/$(parg["runname"]).h5"), parg["tele"], mjd)
    function get_1d_name_partial(expid)
        if df.imagetyp[expid] == "Object"
            return parg["outdir"] * "/apred/$(mjd)/" * get_1d_name(expid, df, cal = true) * ".h5"
        else
            return nothing
        end
    end
    local1D = get_1d_name_partial.(expid_list)
    push!(list1DexpObject, filter(!isnothing, local1D))
end
all1DObjecta = vcat(list1DexpObject...)

all1DObjectperchip = []
for chip in CHIP_LST
    all1DObjectchip = replace.(all1DObjecta, "_$(CHIP_LST[1])_" => "_$(chip)_")
    push!(all1DObjectperchip, all1DObjectchip)
end
all1DObject = vcat(all1DObjectperchip...)

## load rough wave dict and sky lines list
@everywhere begin
    roughwave_dict = load(src_dir * "data/roughwave_dict.jld2", "roughwave_dict")
    df_sky_lines = CSV.read(src_dir * "data/APOGEE_lines.csv", DataFrame)
    df_sky_lines.linindx = 1:size(df_sky_lines, 1)
end

## get sky line peaks
println("Fitting sky line peaks:");
flush(stdout);
@everywhere get_and_save_sky_peaks_partial(fname) = get_and_save_sky_peaks(
    fname, roughwave_dict, df_sky_lines)
@showprogress pmap(get_and_save_sky_peaks_partial, all1DObject)

## get wavecal from sky line peaks
println("Solving skyline wavelength solution:");
flush(stdout);
all1DObjectSkyPeaks = replace.(
replace.(all1DObject, "ar1Dcal" => "skyLinePeaks"), "ar1D" => "skyLinePeaks")
@showprogress pmap(get_and_save_sky_wavecal, all1DObjectSkyPeaks)

## TODO when are we going to split into individual fiber files? Then we should be writing fiber type to the file name
## combine chips for single exposure onto loguniform wavelength grid
## pushing off the question of dither combinations for now (to apMADGICS stage)
all1Da = replace.(all2Dperchip[1], "ar2D" => "ar1D")
println("Reinterpolating exposure spectra:");
flush(stdout);
@showprogress pmap(reinterp_spectra, all1Da)

all1Da = replace.(all2Dperchip[1], "ar2D" => "ar1Dcal")
println("Reinterpolating calibrated exposure spectra:");
flush(stdout);
@showprogress pmap(reinterp_spectra, all1Da)

## I should probably add some slack plots from the wavecal skylines
