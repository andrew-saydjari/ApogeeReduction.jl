## This is a reduction pipeline for APOGEE
using InteractiveUtils;
versioninfo();
@time "Package instantiation and precompilation" begin
    import Pkg
    Pkg.instantiate()
    Pkg.precompile() # no need for Pkg.activate("./") because of invocation w/ environment
end

@time "Package activation" using Distributed, ArgParse, TimerOutputs
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
        "--chips"
        required = false
        help = "chip(s) to run, usually R, G, or B"
        arg_type = String
        default = "RGB"
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
        "--caldir_darks"
        required = false
        help = "outdir where to look for the dark cals"
        arg_type = String
        default = "../outdir/"
        "--caldir_flats"
        required = false
        help = "outdir where to look for the flat cals"
        arg_type = String
        default = "../outdir/"
        "--doCal2d"
        required = false
        help = "run the 2D calibration step"
        arg_type = Bool
        default = true
        "--workers_per_node"
        required = false
        help = "number of workers per node"
        arg_type = Int
        default = -1 # -1 means use all the cores on the node
        "--cluster"
        required = false
        help = "cluster name (sdss or cca or path)"
        arg_type = String
        default = "sdss"
        "--checkpoint_mode"
        required = false
        help = "checkpoint mode (clobber, commit_exists, commit_same)"
        arg_type = String
        default = "commit_same"
        "--suppress_cluster_path_warning"
        required = false
        help = "suppress cluster path warnings"
        arg_type = Bool
        default = false
        "--gain_read_cal_dir"
        required = false
        help = "path to the gain and read noise calibration directory"
        arg_type = String
        default = "/uufs/chpc.utah.edu/common/home/u6039752/scratch1/working/2025_06_03/pass_clean/"
    end
    return parse_args(s)
end
parg = parse_commandline()

@time "Allocating workers" begin
    workers_per_node = parg["workers_per_node"]
    proj_path = dirname(Base.active_project()) * "/"
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

@time "Loading worker packages" @everywhere begin
    using LinearAlgebra
    BLAS.set_num_threads(1)
    using FITSIO, HDF5, FileIO, JLD2, Glob
    using DataFrames, EllipsisNotation, StatsBase
    using ParallelDataTransfer, ProgressMeter
    using AstroTime: TAIEpoch, modified_julian, days, value
    using ApogeeReduction: load_read_var_maps, load_gain_maps, load_saturation_maps, process_3D,
                           process_2Dcal, cal2df, get_cal_path, TAIEpoch
end
@passobj 1 workers() parg
@passobj 1 workers() proj_path
println(BLAS.get_config());
flush(stdout);

##### 3D stage

## load these based on the chip keyword to the pipeline parg
# load gain and readnoise calibrations
# currently globals, should pass and wrap in the partial
# read noise is DN/read
@everywhere begin
    readVarMatDict = load_read_var_maps(parg["gain_read_cal_dir"], parg["tele"], parg["chips"])
    # gain is e-/DN
    gainMatDict = load_gain_maps(parg["gain_read_cal_dir"], parg["tele"], parg["chips"])
    saturationMatDict = load_saturation_maps(
        parg["tele"], parg["chips"], datadir = proj_path * "data/saturation_maps")
end

# write out sym links in the level of folder that MUST be uniform in their cals? or a billion symlinks with expid

# setup the (sjd, expid, chip) tuples to iterate over
# if we have a runlist, we iterate over the mjd and expid in the runlist
# otherwise we iterate over the mjd, expid, and chips specified on the command line
subiter = if parg["runlist"] != ""
    subDic = load(parg["runlist"])
    msk = subDic["tele"] .== parg["tele"]
    # add a filter on tele
    Iterators.product(
        Iterators.zip(subDic["mjd"][msk], subDic["expid"][msk]),
        string.(collect(parg["chips"]))
    )
else
    Iterators.product(
        [(parg["mjd"], parg["expid"])],
        string.(collect(parg["chips"]))
    )
end

t_now = now();
dt = Dates.canonicalize(Dates.CompoundPeriod(t_now - t_then));
println("Setting up the iterator took $dt");
t_then = t_now;
flush(stdout);
# partially apply the process_3D function to everything except the (sjd, expid, chip) values
@everywhere process_3D_partial(((mjd, expid), chip)) = process_3D(
    parg["outdir"], parg["runname"], parg["tele"], mjd, expid, chip,
    gainMatDict, readVarMatDict, saturationMatDict, cluster = parg["cluster"],
    suppress_warning = parg["suppress_cluster_path_warning"],
    checkpoint_mode = parg["checkpoint_mode"])
desc = "3D->2D for $(parg["tele"]) $(parg["chips"])"
ap2dnamelist = @showprogress desc=desc pmap(process_3D_partial, subiter)

# Find the 2D calibration files for the relevant MJDs
unique_mjds = if parg["runlist"] != ""
    subDic = load(parg["runlist"])
    unique(subDic["mjd"])
else
    [parg["mjd"]]
end

# probably need to capture that calFlag somehow, write a meta cal file?
all2D = vcat(ap2dnamelist...)
if parg["doCal2d"]
    darkFlist = sort(glob("darkRate*.h5", parg["caldir_darks"] * "darks/"))
    df_dark = cal2df(darkFlist)

    flatFlist = sort(glob("flatFraction*.h5", parg["caldir_flats"] * "flats/"))
    df_flat = cal2df(flatFlist)

    for mjd in unique_mjds
        for chip in parg["chips"]
            calPath, calFlag = get_cal_path(df_dark, parg["tele"], mjd, string(chip))
            linkPath = parg["outdir"] * "/apred/$(mjd)/" * basename(calPath)
            if !isfile(linkPath)
                mkpath(dirname(linkPath))
                symlink(calPath, linkPath)
            end

            calPath, calFlag = get_cal_path(df_flat, parg["tele"], mjd, string(chip))
            linkPath = parg["outdir"] * "/apred/$(mjd)/" * basename(calPath)
            if !isfile(linkPath)
                mkpath(dirname(linkPath))
                symlink(calPath, linkPath)
            end
        end
    end

    # process the 2D calibration for all exposures
    @everywhere process_2Dcal_partial(fname) = process_2Dcal(
        fname, checkpoint_mode = parg["checkpoint_mode"])
    @showprogress desc="2D Calibration" pmap(process_2Dcal_partial, all2D)
end
