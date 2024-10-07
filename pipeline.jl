## This is a reduction pipeline for APOGEE

import Pkg;
using Dates;
t0 = now();
t_then = t0;
using InteractiveUtils;
versioninfo();

# Pkg.activate("./") # just call with is the environment already activated
# Pkg seems to default to the nonexistant "master" branch for these packages and throw an error
# if I don't specify main

# When this (https://github.com/JuliaLang/Pkg.jl/pull/3783) is in stable Julia Pkg,
# these won't be required.
#Pkg.add(url = "https://github.com/nasa/SIRS.gite#main")
#Pkg.add(url = "https://github.com/andrew-saydjari/SlackThreads.jl.git#main")
Pkg.instantiate()
Pkg.precompile()

@timeit "imports" using Distributed, ArgParse
t_now = now();
dt = Dates.canonicalize(Dates.CompoundPeriod(t_now - t_then));
println("Package activation took $dt");
t_then = t_now;
flush(stdout);

## Parse command line arguments
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
        "--chip"
        required = true
        help = "chip to run, usually a, b, or c"
        arg_type = String
        default = "a"
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
    end
    return parse_args(s)
end

parg = parse_commandline()
if parg["runlist"] != "" # only multiprocess if we have a list of exposures
    if "SLURM_JOB_ID" in keys(ENV)
        using SlurmClusterManager
        addprocs(SlurmManager(), exeflags = ["--project=./"])
    else
        addprocs(16)
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
    using FITSIO, HDF5, FileIO, JLD2
    using DataFrames, EllipsisNotation, StatsBase
    using ParallelDataTransfer, SIRS, ProgressMeter

    src_dir = "./"
    include(src_dir * "src/ap3D.jl")
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

##### 3D stage
@everywhere begin
    # firstind overriden for APO dome flats
    function process_3D(outdir, caldir, runname, mjd, expid, chip; firstind = 3,
            cor1fnoise = true, extractMethod = "sutr_tb")
        dirName = outdir * "/apred/$(mjd)/"
        if !ispath(dirName)
            mkpath(dirName)
        end

        f = h5open(outdir * "almanac/$(runname).h5")
        df = DataFrame(read(f["$(parg["tele"])/$(mjd)/exposures"]))
        close(f)

        # check if chip is in the llist of chips in df.something[expid] (waiting on Andy Casey to update alamanc)
        rawpath = build_raw_path(
            df.observatory[expid], df.mjd[expid], chip, df.exposure[expid])
        # decompress and convert apz data format to a standard 3D cube of reads
        cubedat, hdr = apz2cube(rawpath)

        # ADD? reset anomaly fix (currently just drop first ind as our "fix")
        # REMOVES FIRST READ (as a view)
        # might need to adjust for the few read cases (2,3,4,5)
        firstind_loc = if ((df.exptype[expid] == "DOMEFLAT") &
                           (df.observatory[expid] == "apo")) # NREAD 5, and lamp gets shutoff too soon (needs to be DCS)
            2
        else
            firstind
        end

        tdat = @view cubedat[:, :, firstind_loc:end]

        ## remove 1/f correlated noise (using SIRS.jl) [some preallocs would be helpful]
        if cor1fnoise
            in_data = Float64.(tdat[1:2048, :, :])
            outdat = zeros(Float64, size(tdat, 1), size(tdat, 2), size(in_data, 3)) #obviously prealloc...
            sirssub!(sirs4amps, in_data, f_max = 95.0)
            outdat[1:2048, :, :] .= in_data

            # fixing the 1/f in the reference array is a bit of a hack right now (IRRC might fix)
            in_data = Float64.(tdat[1:2048, :, :])
            in_data[513:1024, :, :] .= tdat[2049:end, :, :]
            sirssub!(sirsrefas2, in_data, f_max = 95.0)
            outdat[2049:end, :, :] .= in_data[513:1024, :, :]
        else
            outdat = Float64.(tdat)
        end

        ## should this be done on the difference cube, or is this enough?
        # vert_ref_edge_corr!(outdat)
        refarray_zpt!(outdat)
        vert_ref_edge_corr_amp!(outdat)

        # ADD? reference array-based masking/correction

        # ADD? nonlinearity correction

        # extraction 3D -> 2D
        dimage, ivarimage, chisqimage = if extractMethod == "dcs"
            dcs(outdat, gainMat, readVarMat)
        elseif extractMethod == "sutr_tb"
            sutr_tb(outdat, gainMat, readVarMat)
        else
            error("Extraction method not recognized")
        end

        # need to clean up exptype to account for FPI versus ARCLAMP
        outfname = join(
            ["ap2D", df.observatory[expid], df.mjd[expid],
                chip, df.exposure[expid], df.exptype[expid]],
            "_")
        # probably change to FITS to make astronomers happy (this JLD2, which is HDF5, is just for debugging)
        jldsave(
            outdir * "/apred/$(mjd)/" * outfname * ".jld2"; dimage, ivarimage, chisqimage)
    end
end
t_now = now();
dt = Dates.canonicalize(Dates.CompoundPeriod(t_now - t_then));
println("Function definitions took $dt");
t_then = t_now;
flush(stdout);

@passobj 1 workers() parg
@everywhere caldir = "/uufs/chpc.utah.edu/common/home/u6039752/scratch1/working/2024_08_14/outdir/cal/" # hard coded for now

## load these based on the chip keyword to the pipeline parg
# load gain and readnoise calibrations
# currently globals, should pass and wrap in the partial
@everywhere gainMat = 1.9 * ones(Float32, 2560, 2048)
@everywhere readVarMat = 25 * ones(Float32, 2560, 2048)
# ADD load the dark currrent map
# load SIRS.jl models
@everywhere sirs4amps = SIRS.restore(caldir * "sirs_test_d12_r60_n15.jld"); # these really are too big... we need to work on reducing their size
@everywhere sirsrefas2 = SIRS.restore(caldir * "sirs_test_ref2_d12_r60_n15.jld");
# write out sym links in the level of folder that MUST be uniform in their cals? or a billion symlinks with expid

try
    if parg["runlist"] != ""
        subDic = load(parg["runlist"])
        subiter = Iterators.zip(subDic["mjd"], subDic["expid"])
        @everywhere process_3D_partial((mjd, expid)) = process_3D(
            parg["outdir"], caldir, parg["runname"], mjd, expid, parg["chip"]) # does Julia LRU cache this?
        @showprogress pmap(process_3D_partial, subiter)
    else
        process_3D(parg["outdir"], caldir, parg["runname"],
            parg["mjd"], parg["expid"], parg["chip"])
    end
finally
    rmprocs(workers())
end
