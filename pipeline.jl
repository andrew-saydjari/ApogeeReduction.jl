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
    using ParallelDataTransfer

    src_dir = "./"
    include(src_dir*"src/ap3D.jl")
    include(src_dir*"src/fileNameHandling.jl")
    include(src_dir*"src/utils.jl")
end

println(BLAS.get_config()); flush(stdout)
using LibGit2; git_branch, git_commit = initalize_git(src_dir); @passobj 1 workers() git_branch; @passobj 1 workers() git_commit
## some hard coded parameters

##### 3D stage
function process_3D(release_dir,mjd,expid)
    f = h5open("$(mjd).h5")
    df = DataFrame(read(f["apo/$(mjd)/exposures"]))
    close(f)

    # this is only doing chip A for now (because of almanac)
    rawpath = build_raw_path(release_dir,df.observatory[expid],df.mjd[expid],df.chip[expid],df.exposure[expid])
    # decompress and convert apz data format to a standard 3D cube of reads
    cubedat, hdr = apz2cube(rawpath);
    # drop first read, then convert cube to differences
    dcubedat = diff(cubedat[:,:,2:end],dims=3);
    dcubedat_out = refcorr(dcubedat);
    dcubedat_out_v = vert_ref_edge_corr(dcubedat_out);
    # dcubedat_out_vh = horz_ref_edge_corr(dcubedat_out_v);

    # need to clean up exptype to account for FPI versus ARCLAMP
    outfname = join(["ap3D",df.observatory[expid],df.mjd[expid],df.chip[expid],df.exposure[expid],df.exptype[expid]],"_")
    # probably change to FITS to make astronomers happy (this JLD2, which is HDF5, is just for debugging)
    jldsave(outfname*".jld2"; dcubedat_out_v)
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
    end
    return parse_args(s)
end

parg = parse_commandline()
release_dir = "sdsswork" # hardcoding for now, 3D raws never move so we could leave this out
process_3D(release_dir,parg["mjd"],parg["expid"])

