## This is a reduction pipeline for APOGEE

import Pkg; using Dates; t0 = now(); t_then = t0;
using InteractiveUtils; versioninfo()
Pkg.activate("./"); Pkg.instantiate(); Pkg.precompile()
using Distributed, ArgParse
t_now = now(); dt = Dates.canonicalize(Dates.CompoundPeriod(t_now-t_then)); println("Package activation took $dt"); t_then = t_now; flush(stdout)

@everywhere begin
    using LinearAlgebra
    BLAS.set_num_threads(1)
    using FITSIO, HDF5, FileIO, JLD2
    using DataFrames, EllipsisNotation, StatsBase
    using ParallelDataTransfer, SIRS

    src_dir = "./"
    include(src_dir*"src/ap3D.jl")
    include(src_dir*"src/fileNameHandling.jl")
    include(src_dir*"src/utils.jl")
end

println(BLAS.get_config()); flush(stdout)
using LibGit2; git_branch, git_commit = initalize_git(src_dir); @passobj 1 workers() git_branch; @passobj 1 workers() git_commit
## some hard coded parameters

##### 3D stage
function process_3D(release_dir,mjd,expid,outdir;firstind=2)
    dirName = outdir*"/ap2D/"
    if !ispath(dirName)
        mkpath(dirName)
    end

    f = h5open(outdir*"almanac/$(mjd).h5")
    df = DataFrame(read(f["apo/$(mjd)/exposures"]))
    close(f)

    # load gain and readnoise calibrations
    gainMat = 1.9*ones(Float32,2560,2048)
    readVarMat = 25*ones(Float32,2560,2048)
    # ADD load the dark currrent map
    # load SIRS.jl models
    sirs4amps = SIRS.restore(outdir*"cal/sirs_test_d12_r60_n15.jld");
    sirsrefas2 = SIRS.restore(outdir*"cal/sirs_test_ref2_d12_r60_n15.jld");
    # write out sym links in the level of folder that MUST be uniform in their cals? or a billion symlinks with expid

    # this is only doing chip A for now (because of almanac)
    rawpath = build_raw_path(release_dir,df.observatory[expid],df.mjd[expid],df.chip[expid],df.exposure[expid])
    # decompress and convert apz data format to a standard 3D cube of reads
    cubedat, hdr = apz2cube(rawpath);

    # ADD? reset anomaly fix (currently just drop first ind as our "fix")
    tdat = @view cubedat[:,:,firstind:end]

    ## remove 1/f correlated noise (using SIRS.jl) [some preallocs would be helpful]
    in_data = Float64.(tdat[1:2048,:,:]);
    outdat = zeros(Float64,size(tdat,1),size(tdat,2),size(in_data,3))
    sirssub!(sirs4amps, in_data, f_max=95.);
    outdat[1:2048,:,:] .= in_data
    
    # fixing the 1/f in the reference array is a bit of a hack right now (IRRC might fix)
    in_data = Float64.(tdat[1:2048,:,:]);
    in_data[513:1024,:,:].=tdat[2049:end,:,:]
    sirssub!(sirsrefas2, in_data, f_max=95.);
    outdat[2049:end,:,:] .= in_data[513:1024,:,:]

    # ADD? reference array-based masking/correction

    # ADD? nonlinearity correction

    # extraction 3D -> 2D
    dimage, ivarimage = dcs(outdat,gainMat,readVarMat,firstind=1);

    # dark current subtraction

    # need to clean up exptype to account for FPI versus ARCLAMP
    outfname = join(["ap2D",df.observatory[expid],df.mjd[expid],df.chip[expid],df.exposure[expid],df.exptype[expid]],"_")
    # probably change to FITS to make astronomers happy (this JLD2, which is HDF5, is just for debugging)
    jldsave(outdir*"/ap2D/"*outfname*".jld2"; dimage, ivarimage)
end

## Parse command line arguments
function parse_commandline()
    s=ArgParseSettings()
    @add_arg_table s begin
        "--mjd"
            required = true
            help = "mjd of the exposure(s) to be run"
            arg_type = Int
            default = 1

        # probably want to add in defaults that loops over them all
        "--expid"
            required = true
            help = "exposure number to be run"
            arg_type = Int
            default = 1
        "--outdir"
            required = true
            help = "output directory"
            arg_type = String
            default = "../outdir/"
    end
    return parse_args(s)
end

parg = parse_commandline()
release_dir = "sdsswork" # hardcoding for now, 3D raws never move so we could leave this out
process_3D(release_dir,parg["mjd"],parg["expid"],parg["outdir"])

