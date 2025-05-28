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
        help = "chip(s) to run, usually a, b, or c"
        arg_type = String
        default = "abc"
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
        default = 58
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
println("Allocating $(nworkers()) workers took $dt")
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
    using AstroTime: TAIEpoch, modified_julian, days, value

    src_dir = "./"
    include(src_dir * "src/ar3D.jl")
    include(src_dir * "src/ar2Dcal.jl")
    include(src_dir * "src/fileNameHandling.jl")
    include(src_dir * "src/utils.jl")
end

println(BLAS.get_config());
flush(stdout);

##### 3D stage
@everywhere begin
    # firstind overriden for APO dome flats
    function process_3D(outdir, sirscaldir, runname, mjd, expid, chip; firstind = 3,
            cor1fnoise = true, extractMethod = "sutr_wood")
        dirName = joinpath(outdir, "apred/$(mjd)/")
        if !ispath(dirName)
            mkpath(dirName)
        end

        df = read_almanac_exp_df(joinpath(outdir, "almanac/$(runname).h5"), parg["tele"], mjd)
        #        println(expid,chip,size(df.observatory),size(df.mjd),size(df.exposure_int))
        # check if chip is in the llist of chips in df.something[expid] (waiting on Andy Casey to update alamanc)
        rawpath = build_raw_path(
            df.observatory[expid], chip, df.mjd[expid], lpad(df.exposure_int[expid], 8, "0"))
        cartid = df.cartidInt[expid]
        # decompress and convert apz data format to a standard 3D cube of reads
        cubedat, hdr_dict = apz2cube(rawpath)

        nread_total = size(cubedat, 3)

        # ADD? reset anomaly fix (currently just drop first ind or two as our "fix")
        # REMOVES FIRST READ (as a view)
        # might need to adjust for the few read cases (2,3,4,5)
        firstind_loc,
        extractMethod_loc = if ((df.exptype[expid] == "DOMEFLAT") &
                                (df.observatory[expid] == "apo")) # NREAD 5, and lamp gets shutoff too soon (needs to be DCS)
            2, "dcs"
        elseif ((df.exptype[expid] == "QUARTZFLAT") & (nread_total == 3))
            2, "dcs"
        elseif (nread_total == 3)
            #catch some weird cases (like nreads=3 with Darks)
            #but still reduce them to prevent errors later in pipeline_2d_1d
            #ULTIMATELY want to make it so these exposures are removed from runlist earlier
            2, "dcs"
        else
            firstind, extractMethod
        end

        # Might want lastind_loc as well to truncate when all reads are
        # saturated (important for calculating the exposure mid time)
        lastind_loc = size(cubedat, 3)

        tdat = @view cubedat[:, :, firstind_loc:lastind_loc]
        nread_used = size(tdat, 3) - 1 #this is actually nread_used-1 (ie it is truly ndiff_used)

        n_read_dropped = firstind_loc - 1
        first_image_start_time = TAIEpoch(hdr_dict[firstind_loc]["DATE-OBS"])
        last_image_start_time = TAIEpoch(hdr_dict[lastind_loc]["DATE-OBS"])
        dtime_read = (hdr_dict[firstind_loc]["INTOFF"] / 1000 / 3600 / 24)days #dt_read, orig in ms, convert to days
        dtime_delay = (hdr_dict[firstind_loc]["INTDELAY"] / 3600 / 24)days #orig in seconds, convert to days
        exptime_est = (hdr_dict[firstind_loc]["EXPTIME"] / 3600 / 24)days

        # NOT using dtime_delay because we start directly at first_image_start_time
        # (though we might need to think about this more: not sure what dtime_delay does)
        # REMEMBER to add half of a dtime_read to shift to center of exposure

        # Like DRP outputs (we think)
        mjd_mid_exposure_old = modified_julian(first_image_start_time
                                               +
                                               0.5 * exptime_est * size(tdat, 3) / nread_total)
        # Using dread_time*nread_USED
        mjd_mid_exposure_rough = modified_julian(first_image_start_time
                                                 +
                                                 dtime_read *
                                                 (0.5 + 0.5 * (size(tdat, 3) - 1)))
        # Using times directly from header
        mjd_mid_exposure_precise = modified_julian(first_image_start_time + 0.5 * dtime_read
                                                   +
                                                   0.5 *
                                                   (last_image_start_time - first_image_start_time))
        mjd_mid_exposure = mjd_mid_exposure_precise

        ## remove 1/f correlated noise (using SIRS.jl) [some preallocs would be helpful]
        if cor1fnoise
            in_data = Float64.(tdat[1:2048, :, :])
            # sirsub! modifies in_data, but make a copy so that it's faster the second time.
            # better not to need to do this at all
            copied_in_data = copy(in_data)
            outdat = zeros(
                Float64, size(tdat, 1), size(tdat, 2), size(in_data, 3)) #obviously prealloc...
            sirssub!(sirs4amps, in_data, f_max = 95.0)
            outdat[1:2048, :, :] .= in_data

            # fixing the 1/f in the reference array is a bit of a hack right now (IRRC might fix)
            in_data = copied_in_data
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
        dimage, ivarimage,
        chisqimage,
        CRimage,
        saturation_image = if extractMethod_loc == "dcs"
            # TODO some kind of outlier rejection, this just keeps all diffs
            images = dcs(outdat, gainMatDict[chip], readVarMatDict[chip])
            images..., zeros(Int, size(images[1])), zeros(Int, size(images[1]))
        elseif extractMethod_loc == "sutr_wood"
            # n.b. this will mutate outdat
            sutr_wood!(outdat, gainMatDict[chip], readVarMatDict[chip])
        else
            error("Extraction method not recognized")
        end

        # need to clean up exptype to account for FPI versus ARCLAMP
        outfname = join(
            ["ar2D", df.observatory[expid], df.mjd[expid],
                last(df.exposure_str[expid], 4), chip, df.exptype[expid]],
            "_")
        # probably change to FITS to make astronomers happy (this JLD2, which is HDF5, is just for debugging)

        metadata = Dict(
            "cartid" => cartid,
            "nread_used" => nread_used,
            "nread_total" => nread_total,
            "mjd_mid_exposure_old" => value(mjd_mid_exposure_old),
            "mjd_mid_exposure_rough" => value(mjd_mid_exposure_rough),
            "mjd_mid_exposure_precise" => value(mjd_mid_exposure_precise),
            "mjd_mid_exposure" => value(mjd_mid_exposure)
        )
        fname = joinpath(outdir, "apred/$(mjd)/" * outfname * ".h5")
        safe_jldsave(fname, metadata; dimage, ivarimage, chisqimage, CRimage, saturation_image)
        return fname
    end

    # come back to tuning the chi2perdofcut once more rigorously establish noise model
    function process_2Dcal(fname; chi2perdofcut = 100)
        sname = split(split(split(fname, "/")[end], ".h5")[1], "_")
        fnameType, tele, mjd, expnum, chip, exptype = sname[(end - 5):end]

        dimage = load(fname, "dimage")
        ivarimage = load(fname, "ivarimage")
        CRimage = load(fname, "CRimage")
        chisqimage = load(fname, "chisqimage")
        saturation_image = load(fname, "saturation_image")

        metadata = read_metadata(fname)
        nread_used = metadata["nread_used"]

        ### dark current subtraction
        darkRateflst = sort(glob("darkRate_$(tele)_$(chip)_*.h5", dirname(fname)))
        if length(darkRateflst) != 1
            error("I didn't just find one darkRate file for mjd $mjd, I found $(length(darkRateflst))")
        end
        darkRate = load(darkRateflst[1], "dark_rate")
        pix_bitmask = load(darkRateflst[1], "dark_pix_bitmask")

        dimage .-= darkRate
        # should I be modifying ivarimage? (uncertainty on dark rate in quad... but dark subtraction has bigger sys)

        ### flat fielding
        flatFractionflst = sort(glob("flatFraction_$(tele)_$(chip)_*.h5", dirname(fname)))
        if length(flatFractionflst) != 1
            error("I didn't just find one flatFraction file for mjd $mjd, I found $(length(flatFractionflst))")
        end
        flat_im = load(flatFractionflst[1], "flat_im")
        flat_pix_bitmask = load(flatFractionflst[1], "flat_pix_bitmask")
        dimage[5:2044, 5:2044] ./= flat_im
        ivarimage[5:2044, 5:2044] .*= flat_im .^ 2
        pix_bitmask[5:2044, 5:2044] .|= flat_pix_bitmask

        pix_bitmask .|= (CRimage .== 1) * 2^7
        pix_bitmask .|= (CRimage .> 1) * 2^8
        pix_bitmask .|= ((chisqimage ./ nread_used) .> chi2perdofcut) * 2^9
        pix_bitmask .|= saturation_image * 2^13

        outfname = replace(fname, "ar2D" => "ar2Dcal")
        safe_jldsave(outfname, metadata; dimage, ivarimage, pix_bitmask)
    end
end
t_now = now();
dt = Dates.canonicalize(Dates.CompoundPeriod(t_now - t_then));
println("Function definitions took $dt");
t_then = t_now;
flush(stdout);

@passobj 1 workers() parg
@everywhere sirscaldir = "/uufs/chpc.utah.edu/common/home/u6039752/scratch1/working/2024_08_14/outdir/cal/" # hard coded for now
@everywhere gainReadCalDir = "/uufs/chpc.utah.edu/common/home/u6039752/scratch1/working/2025_02_03/"

## load these based on the chip keyword to the pipeline parg
# load gain and readnoise calibrations
# currently globals, should pass and wrap in the partial
# read noise is DN/read
readVarMatDict = load_read_var_maps(gainReadCalDir, parg["tele"], parg["chips"])
@passobj 1 workers() readVarMatDict
# gain is e-/DN
gainMatDict = load_gain_maps(gainReadCalDir, parg["tele"], parg["chips"])
@passobj 1 workers() gainMatDict

# ADD load the dark currrent map
# load SIRS.jl models
@everywhere sirs4amps = SIRS.restore(sirscaldir * "sirs_test_d12_r60_n15.jld"); # these really are too big... we need to work on reducing their size
@everywhere sirsrefas2 = SIRS.restore(sirscaldir * "sirs_test_ref2_d12_r60_n15.jld");
# write out sym links in the level of folder that MUST be uniform in their cals? or a billion symlinks with expid

# clean up this statement to have less replication
desc = "3D->2D for $(parg["tele"]) $(parg["chips"])"
if parg["runlist"] != ""
    subDic = load(parg["runlist"])

    subiter = Iterators.product(
        Iterators.zip(subDic["mjd"], subDic["expid"]),
        string.(collect(parg["chips"])))
    @everywhere process_3D_partial(((mjd, expid), chip)) = process_3D(
        parg["outdir"], sirscaldir, parg["runname"], mjd, expid, chip) # does Julia LRU cache this?
    ap2dnamelist = @showprogress desc=desc pmap(process_3D_partial, subiter)
else
    subiter = string.(collect(parg["chips"]))
    @everywhere process_3D_partial(chip) = process_3D(
        parg["outdir"], sirscaldir, parg["runname"], parg["mjd"], parg["expid"], chip)
    ap2dnamelist = @showprogress desc=desc pmap(process_3D_partial, subiter)
end

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
                symlink(calPath, linkPath)
            end

            calPath, calFlag = get_cal_path(df_flat, parg["tele"], mjd, string(chip))
            linkPath = parg["outdir"] * "/apred/$(mjd)/" * basename(calPath)
            if !isfile(linkPath)
                symlink(calPath, linkPath)
            end
        end
    end

    # process the 2D calibration for all exposures
    @showprogress desc="2D Calibration" pmap(process_2Dcal, all2D)
end
