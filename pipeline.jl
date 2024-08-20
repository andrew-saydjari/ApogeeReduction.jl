## This is a reduction pipeline for APOGEE

import Pkg; using Dates; t0 = now(); t_then = t0;
using InteractiveUtils; versioninfo()
Pkg.activate("./")
Pkg.add("ProgressMeter")
Pkg.rm("SIRS")
Pkg.add(url="https://github.com/nasa/SIRS.git")
Pkg.resolve(); Pkg.instantiate(); Pkg.precompile()

using Distributed, ArgParse
t_now = now(); dt = Dates.canonicalize(Dates.CompoundPeriod(t_now-t_then)); println("Package activation took $dt"); t_then = t_now; flush(stdout)
if "SLURM_JOB_ID" in keys(ENV)
    using SlurmClusterManager
    addprocs(SlurmManager(),exeflags=["--project=./"])
else
    addprocs(64)
end
t_now = now(); dt = Dates.canonicalize(Dates.CompoundPeriod(t_now-t_then)); println("Worker allocation took $dt"); t_then = t_now; flush(stdout)
println("Running Main on ", gethostname()); flush(stdout)

@everywhere begin
    using LinearAlgebra
    BLAS.set_num_threads(1)
    using FITSIO, HDF5, FileIO, JLD2
    using DataFrames, EllipsisNotation, StatsBase
    using ParallelDataTransfer, SIRS, ProgressMeter

    src_dir = "./"
    include(src_dir*"src/ap3D.jl")
    include(src_dir*"src/fileNameHandling.jl")
    include(src_dir*"src/utils.jl")
end

println(BLAS.get_config()); flush(stdout)
using LibGit2; git_branch, git_commit = initalize_git(src_dir); @passobj 1 workers() git_branch; @passobj 1 workers() git_commit
## some hard coded parameters

##### 3D stage
@everywhere begin
    function process_3D(release_dir,outdir,runname,mjd,expid;firstind=2,cor1fnoise=true)
	caldir = "/uufs/chpc.utah.edu/common/home/u6039752/scratch1/working/2024_08_14/outdir/cal/"
        dirName = outdir*"/ap2D/"
        if !ispath(dirName)
            mkpath(dirName)
        end

        f = h5open(outdir*"almanac/$(runname).h5")
        df = DataFrame(read(f["apo/$(mjd)/exposures"]))
        close(f)

        # load gain and readnoise calibrations (is there a good way to move these outside?)
        gainMat = 1.9*ones(Float32,2560,2048)
        readVarMat = 25*ones(Float32,2560,2048)
        # ADD load the dark currrent map
        # load SIRS.jl models
        sirs4amps = SIRS.restore(caldir*"sirs_test_d12_r60_n15.jld");
        sirsrefas2 = SIRS.restore(caldir*"sirs_test_ref2_d12_r60_n15.jld");
        # write out sym links in the level of folder that MUST be uniform in their cals? or a billion symlinks with expid

        # this is only doing chip A for now (because of almanac)
        rawpath = build_raw_path(release_dir,df.observatory[expid],df.mjd[expid],df.chip[expid],df.exposure[expid])
        # decompress and convert apz data format to a standard 3D cube of reads
        cubedat, hdr = apz2cube(rawpath);

        # ADD? reset anomaly fix (currently just drop first ind as our "fix")
        tdat = @view cubedat[:,:,firstind:end]

        ## remove 1/f correlated noise (using SIRS.jl) [some preallocs would be helpful]
        if cor1fnoise
            in_data = Float64.(tdat[1:2048,:,:]);
            outdat = zeros(Float64,size(tdat,1),size(tdat,2),size(in_data,3))
            sirssub!(sirs4amps, in_data, f_max=95.);
            outdat[1:2048,:,:] .= in_data
            
            # fixing the 1/f in the reference array is a bit of a hack right now (IRRC might fix)
            in_data = Float64.(tdat[1:2048,:,:]);
            in_data[513:1024,:,:].=tdat[2049:end,:,:]
            sirssub!(sirsrefas2, in_data, f_max=95.);
            outdat[2049:end,:,:] .= in_data[513:1024,:,:]
        else
            outdat = Float64.(tdat)
        end

        ## wondering if we need to add back in the vert_ref_edge_corr?

        # ADD? reference array-based masking/correction

        # ADD? nonlinearity correction

        # extraction 3D -> 2D
        #dimage, ivarimage = dcs(outdat,gainMat,readVarMat,firstind=1);
        dimage, ivarimage = sutr_tb(outdat,gainMat,readVarMat,firstind=1);

	print(mean(dimage),"max",maximum(dimage))
        # dark current subtraction

        # need to clean up exptype to account for FPI versus ARCLAMP
        outfname = join(["ap2D",df.observatory[expid],df.mjd[expid],df.chip[expid],df.exposure[expid],df.exptype[expid]],"_")
        # probably change to FITS to make astronomers happy (this JLD2, which is HDF5, is just for debugging)
        jldsave(outdir*"/ap2D/"*outfname*".jld2"; dimage, ivarimage)
    end
end

## Parse command line arguments
function parse_commandline()
    s=ArgParseSettings()
    @add_arg_table s begin
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
        # 3D raws never move so we could leave release_dir out
        "--release_dir"
            required = false
            help = "path to the release directory where the raw data is stored"
            arg_type = String
            default = "sdsswork"
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
t_now = now(); dt = Dates.canonicalize(Dates.CompoundPeriod(t_now-t_then)); println("Function definitions took $dt"); t_then = t_now; flush(stdout)

parg = parse_commandline()
@passobj 1 workers() parg
if parg["runlist"] != ""
    subDic = load(parg["runlist"])
    subiter = Iterators.zip(subDic["mjd"],subDic["expid"])
    @everywhere process_3D_partial((mjd,expid)) = process_3D(parg["release_dir"],parg["outdir"],parg["runname"],mjd,expid)
    @showprogress pmap(process_3D_partial,subiter)
else
    process_3D(parg["release_dir"],parg["outdir"],parg["runname"],parg["mjd"],parg["expid"])
end
rmprocs(workers())
