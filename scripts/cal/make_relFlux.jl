using ArgParse, Distributed, SlurmClusterManager, SlackThreads

## Parse command line arguments
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--chips"
        required = false
        help = "chip names, i.e. RGB"
        default = "RGB"
        arg_type = String
        "--tele"
        required = false
        help = "telescope name (apo or lco)"
        arg_type = String
        default = ""
        "--mjd"
        required = false
        help = "mjd of the exposure(s) to be run (overridden by runlist)"
        arg_type = Int
        default = 1
        "--expid"
        required = false
        help = "exposure number to be run (overridden by runlist)"
        arg_type = Int
        default = 1
        "--trace_dir"
        required = true
        help = "directory where extractions of traces are stored"
        arg_type = String
        default = ""
        "--runlist"
        required = true
        help = "path name to hdf5 file with keys specifying list of exposures to run"
        arg_type = String
        default = ""
        "--runname"
        required = true
        help = "name of the run (specifically almanac file)"
        arg_type = String
        default = "test"
    end
    return parse_args(s)
end

parg = parse_commandline()

proj_path = dirname(Base.active_project()) * "/"
if parg["runlist"] != "" # only multiprocess if we have a list of exposures
    if "SLURM_NTASKS" in keys(ENV)
        using SlurmClusterManager
        addprocs(SlurmManager(), exeflags = ["--project=$proj_path"])
    else
        addprocs(16)
    end
end

@everywhere begin
    using JLD2, ProgressMeter, ArgParse, Glob, StatsBase, ParallelDataTransfer
    using HDF5, DataFrames, SlackThreads
    using ApogeeReduction
    using ApogeeReduction: safe_jldsave, read_almanac_exp_df, get_1d_name, get_relFlux,
                           read_metadata
end

@passobj 1 workers() parg # make it available to all workers
@passobj 1 workers() proj_path

@everywhere begin
    include(joinpath(proj_path, "src/makie_plotutils.jl"))

    chips = collect(String, split(parg["chips"], ""))

    dirNamePlots = joinpath(parg["trace_dir"], "plots/")
    if !ispath(dirNamePlots)
        mkpath(dirNamePlots)
    end
end

## TODO Confirm multiple MJD handling is correct/ideal
# using runlist for dome/quartz flats only
tele_list = if parg["runlist"] != ""
    load(parg["runlist"], "tele")
else
    [parg["tele"]]
end
unique_teles = unique(tele_list)

mjd_list = if parg["runlist"] != ""
    load(parg["runlist"], "mjd")
else
    [parg["mjd"]]
end
unique_mjds = unique(mjd_list)

expid_list = if parg["runlist"] != ""
    load(parg["runlist"], "expid")
else
    [parg["expid"]]
end

all1Da = String[] # all 1D files for chip a
for tele in unique_teles
    for mjd in unique_mjds
        mskMJD = (mjd_list .== mjd) .& (tele_list .== tele)
        if count(mskMJD) > 0
            df = read_almanac_exp_df(
                joinpath(parg["trace_dir"], "almanac/$(parg["runname"]).h5"), tele, mjd)
            function get_1d_name_partial(expid)
                joinpath(parg["trace_dir"], "apredrelflux/$(mjd)",
                    get_1d_name(expid, df, cal = true) * ".h5")
            end
            file_list = get_1d_name_partial.(expid_list[mskMJD])
            append!(all1Da, file_list)
        end
    end
end

if length(all1Da) > 0
    ## need to get cal_type from runlist
    exp_type_lst = map(x -> split(split(x, "FLAT")[1], "_")[end], all1Da)
    unique_exp_lst = unique(exp_type_lst)
    if length(unique_exp_lst) > 1
        error("Multiple cal types found in runlist")
    end
    cal_type = lowercase(unique_exp_lst[1])
    @passobj 1 workers() cal_type # make cal_type available to all workers

    flist_chips = []
    for chip in chips
        push!(flist_chips, replace.(all1Da, "_$(FIRST_CHIP)_" => "_$(chip)_"))
    end
    all1D = vcat(flist_chips...)

    @everywhere begin
        function get_and_save_relFlux(fname)
            absthrpt, relthrpt, bitmsk_relthrpt, metadata = get_relFlux(fname)
            cartid = Int(metadata["cartid"])
            outfname = replace(fname,
                "apredrelflux" => "$(cal_type)_flats",
                "apred" => "$(cal_type)_flats",
                "ar1Dcal" => "$(cal_type)Flux",
                ".h5" => "_$(cartid).h5"
            )
            dname = dirname(outfname)
            if !ispath(dname)
                mkpath(dname)
            end
            safe_jldsave(outfname, metadata; absthrpt, relthrpt, bitmsk_relthrpt)
            return outfname
        end
    end

    desc = "get and save relFlux from Flats"
    outfnames = @showprogress desc=desc pmap(get_and_save_relFlux, all1D)
    outfnamesa = outfnames[1:length(all1Da)]

    ## I want to take the average of relFlux for adjacent exposures with the same CartID
    ## average over the chips except b (or should we just take one chip?)
    # just going to use chip c for now, can revisit this question if we think it adds stability

    # add plotting
    thread = SlackThread()
    if length(unique_mjds) == 1
        thread("$(cal_type) relFluxing for $(unique_teles[begin]) for SJD $(unique_mjds[begin])")

        @everywhere begin
            function plot_relFlux(fname)
                sname = split(split(split(fname, "/")[end], ".h5")[1], "_")
                fnameType, tele, mjd, expnum, chiploc, exptype, cartid = sname[(end - 6):end]

                xvec = if tele == "apo"
                    1:300
                elseif tele == "lco"
                    301:600
                else
                    error("Unknown telescope: $(tele)")
                end
                absthrpt = zeros(300, length(chips))
                relthrpt = zeros(300, length(chips))
                bitmsk_relthrpt = zeros(Int, 300, length(chips))
                cartid = Int(read_metadata(fname)["cartid"])
                for (cindx, chip) in enumerate(chips)
                    local_fname = replace(fname, "_$(chiploc)_" => "_$(chip)_")
                    f = jldopen(local_fname)
                    absthrpt[:, cindx] = f["absthrpt"]
                    relthrpt[:, cindx] = f["relthrpt"]
                    bitmsk_relthrpt[:, cindx] = f["bitmsk_relthrpt"]
                    close(f)
                end

                # Initialize strings to collect fiber status
                status_str = "Fiber Status Report for $(tele) MJD=$(mjd) CartID=$(cartid):\n"

                # plot the relFlux
                fig = Figure(size = (1200, 400))
                for (cindx, chip) in enumerate(chips)
                    ax = Axis(fig[1, cindx], title = "RelFlux Chip $(chip)")
                    msk = bitmsk_relthrpt[:, cindx] .== 0
                    broken_msk = (bitmsk_relthrpt[:, cindx] .& 2^1) .== 2^1
                    warn_msk = (bitmsk_relthrpt[:, cindx] .& 2^0) .== 2^0
                    warn_only_msk = warn_msk .& .!broken_msk

                    scatter!(ax, xvec[msk], relthrpt[msk, cindx], color = "limegreen")

                    # Broken fibers
                    scatter!(ax, xvec[broken_msk], relthrpt[broken_msk, cindx], color = "red")
                    broken_fibers = xvec[broken_msk]
                    if !isempty(broken_fibers) && (chip == LAST_CHIP) # hardcoded for now
                        status_str *= "\nChip $(chip) Broken Fibers (adjfiberindx):\n    $(join(broken_fibers, "\n    "))"
                    end

                    # Warn fibers
                    scatter!(
                        ax, xvec[warn_only_msk], relthrpt[warn_only_msk, cindx], color = "orange")
                    warn_fibers = xvec[warn_only_msk]
                    if !isempty(warn_fibers) && (chip == LAST_CHIP) # hardcoded for now
                        status_str *= "\nChip $(chip) Warning Fibers (adjfiberindx):\n    $(join(warn_fibers, "\n    "))"
                    end
                end

                # we should be more uniform about the naming convention
                savePath = dirNamePlots *
                           "$(mjd)/relFlux_$(cal_type)_$(tele)_$(mjd)_$(expnum)_$(cartid).png"
                save(savePath, fig)

                return savePath, status_str
            end
        end

        desc = "plotting relFlux"
        pout = @showprogress desc=desc pmap(plot_relFlux, outfnamesa)
        status_str_lst = map(x -> x[2], pout)
        unique_status_strs = unique(status_str_lst)
        savePath_lst = map(x -> x[1], pout)
        for status_str in unique_status_strs
            thread("Relfluxing: $(status_str)")
        end
        for savePath in savePath_lst
            thread("Relfluxing: $(savePath)", savePath)
        end
    elseif length(unique_mjds) > 1
        ## need to think harder about the plotting we want to do in this case
        thread("$(cal_type) relFluxing for $(unique_teles) with multiple SJDs from $(unique_mjds[begin]) to $(unique_mjds[end])")
    end
else
    thread = SlackThread()
    thread("No files of at least one cal type found for relFluxing")
end
