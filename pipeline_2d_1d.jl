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
    end
    return parse_args(s)
end

parg = parse_commandline()
if parg["runlist"] != "" # only multiprocess if we have a list of exposures
    if "SLURM_NTASKS" in keys(ENV)
        using SlurmClusterManager
        addprocs(SlurmManager(), exeflags = ["--project=./"])
    else
        addprocs(32)
    end
end
t_now = now();
dt = Dates.canonicalize(Dates.CompoundPeriod(t_now - t_then));
println("Worker allocation took $dt");
t_then = t_now;
flush(stdout);
println("Running Main on ", gethostname());
flush(stdout);

@everywhere begin
    using LinearAlgebra
    BLAS.set_num_threads(1)
    using FITSIO, HDF5, FileIO, JLD2, Glob
    using DataFrames, EllipsisNotation, StatsBase
    using ParallelDataTransfer, SIRS, ProgressMeter

    src_dir = "./"
    include(src_dir * "src/ap1D.jl")
    include(src_dir * "src/fileNameHandling.jl")
    include(src_dir * "src/utils.jl")
end

println(BLAS.get_config());
flush(stdout);
using LibGit2;
git_branch, git_commit = initalize_git(src_dir);
@passobj 1 workers() git_branch;
@passobj 1 workers() git_commit;
## some hard coded parameters

##### 1D stage
@everywhere begin
    function process_1D(fname)
        sname = split(fname, "_")
        fnameType, tele, mjd, chip, expid = sname[(end - 5):(end - 1)]

        fnamecal = if (fnameType == "ap2D")
            replace(fname, "ap2D" => "ap2Dcal")
        else
            fname
        end

        dimage = load(fname, "dimage")
        ivarimage = load(fname, "ivarimage")
        pix_bitmask = load(fnamecal, "pix_bitmask") #strip out the replace once we are happy with ap2Dcal

        # this seems annoying to load so often if we know we are doing a daily... need to ponder
        traceList = sort(glob("domeTraceMain_$(tele)_$(mjd)_*_$(chip).jld2",
            parg["outdir"] * "apred/$(mjd)/"))
        trace_params = load(traceList[1], "trace_params")

        # adam: should this be saved somewhere?  It's fairly simple to reproduce, but that's true of 
        # everything to some degree
        regularized_trace_params = regularize_trace(trace_params)

        flux_1d, ivar_1d, mask_1d = if parg["extraction"] == "boxcar"
            extract_boxcar(dimage, ivarimage, pix_bitmask, regularized_trace_params)
        elseif parg["extraction"] == "optimal"
            extract_optimal(dimage, ivarimage, pix_bitmask, regularized_trace_params)
        else
            error("Extraction method $(parg["extraction"]) not recognized")
        end

        # we probably want to append info from the fiber dictionary from alamanac into the file name
        outfname = replace(fname, "ap2D" => "ap1D")
        jldsave(outfname; flux_1d, ivar_1d, mask_1d, git_branch, git_commit)
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
    f = h5open(parg["outdir"] * "almanac/$(parg["runname"]).h5")
    df = DataFrame(read(f["$(parg["tele"])/$(mjd)/exposures"]))
    close(f)
    function get_2d_name_partial(expid)
        parg["outdir"] * "/apred/$(mjd)/" *
        replace(get_1d_name(expid, df), "ap1D" => "ap2D") * ".jld2"
    end
    local2D = get_2d_name_partial.(expid_list)
    push!(list2Dexp, local2D)
end
all2Da = vcat(list2Dexp...)

all2Dperchip = []
for chip in ["a","b","c"]
    all2Dchip = replace.(all2Da, "_a_" => "_$(chip)_")
    push!(all2Dperchip, all2Dchip)
end
all2D = vcat(all2Dperchip...)

# we should do somthing smart to assemble the traces from a night into a single file 
# that gives us the trace of a fiber as a funciton of time or something
# for now, for each MJD, take the first one (or do that in run_trace_cal.sh)
# I think dome flats needs to swtich to dome_flats/mjd/
for mjd in unique_mjds
    for chip in ["a","b","c"]
        traceList = sort(glob("domeTrace_$(parg["tele"])_$(mjd)_*_$(chip).jld2",
            parg["outdir"] * "dome_flats/"))
        if length(traceList) > 1
            @warn "Multiple dome trace files found for $(parg["tele"]) $(mjd) $(chip): $(traceList)"
        elseif length(traceList) == 0
            @error "No dome trace files found for $(parg["tele"]) $(mjd) $(chip). Looked in $(parg["outdir"] * "dome_flats/")."
        end
        calPath = traceList[1]
        linkPath = parg["outdir"] * "apred/$(mjd)/" *
                   replace(basename(calPath), "domeTrace" => "domeTraceMain")
        if !isfile(linkPath)
            # come back to why this symlink does not work
            cp(calPath, linkPath)
        end
    end
end

# extract the 2D to 1D, ideally the calibrated files
# need to think hard about batching daily versus all data for cal load in
println("Extracting 2D to 1D:")
@showprogress pmap(process_1D, all2D)
println("Extracting 2Dcal to 1Dcal:")
all2Dcal = replace.(all2D, "ap2D" => "ap2Dcal")
@showprogress pmap(process_1D, all2Dcal)

## get all OBJECT files (happy to add any other types that see sky?)
list1DexpObject = []
for mjd in unique_mjds
    f = h5open(parg["outdir"] * "almanac/$(parg["runname"]).h5")
    df = DataFrame(read(f["$(parg["tele"])/$(mjd)/exposures"]))
    close(f)
    function get_1d_name_partial(expid)
        if df.imagetyp[expid] == "Object"
            return parg["outdir"] * "/apred/$(mjd)/" * get_1d_name(expid, df) * ".jld2"
        else
            return []
        end
    end
    local1D = get_1d_name_partial.(expid_list)
    push!(list1DexpObject, local1D)
end
all1DObjecta = vcat(list1DexpObject...)

all1DObjectperchip = []
for chip in ["a","b","c"]
    all1DObjectchip = replace.(all1DObjecta, "_a_" => "_$(chip)_")
    push!(all1DObjectperchip, all1DObjectchip)
end
all1DObject = vcat(all1DObjectperchip...)

## load rough wave dict and sky lines list
roughwave_dict = load(src_dir * "data/roughwave_dict.jld2","roughwave_dict")
df_sky_lines = CSV.read(src_dir * "data/df_sky_lines/APOGEE_lines.csv", DataFrame);
df_sky_lines.linindx = 1:size(df_sky_lines,1);

## get sky line peaks
println("Fitting sky line peaks:")
get_and_save_sky_peaks_partial(fname) = get_and_save_sky_peaks(fname,roughwave_dict,df_sky_lines)
@showprogress pmap(get_and_save_sky_peaks_partial, all1DObject)

## get wavecal from sky line peaks
println("Solving skyline wavelength solution:")
@showprogress pmap(get_and_save_sky_wavecal, all1DObjecta) is 

## add a plot to plot all to just show the chips together
## I should probably add some slack plots from the wavecal skylines
## then do a dither combination (on Fri)