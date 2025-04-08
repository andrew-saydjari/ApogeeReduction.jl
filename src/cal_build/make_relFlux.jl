using ArgParse, Distributed, SlurmClusterManager, SlackThreads

## Parse command line arguments
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--chips"
        required = false
        help = "chip names, i.e. abc"
        default = "abc"
        arg_type = String
        "--tele"
        required = true
        help = "telescope name (apo or lco)"
        arg_type = String
        default = ""
        "--mjd-start"
        required = true
        help = "start mjd"
        arg_type = Int
        default = 0
        "--mjd-end"
        required = true
        help = "end mjd"
        arg_type = Int
        default = 0
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

if parg["runlist"] != "" # only multiprocess if we have a list of exposures
    if "SLURM_NTASKS" in keys(ENV)
        using SlurmClusterManager
        addprocs(SlurmManager(), exeflags = ["--project=./"])
    else
        addprocs(16)
    end
end

@everywhere begin
    using JLD2, ProgressMeter, ArgParse, Glob, StatsBase, ParallelDataTransfer
    using HDF5, DataFrames, SlackThreads
    src_dir = "../"
    include(src_dir * "/makie_plotutils.jl")
    include(src_dir * "/utils.jl")
    include(src_dir * "/fileNameHandling.jl")
    include(src_dir * "/ar1D.jl")
end

@passobj 1 workers() parg # make it available to all workers

@everywhere begin
    chips = collect(String, split(parg["chips"], ""))

    dirNamePlots = parg["trace_dir"] * "plots/"
    if !ispath(dirNamePlots)
        mkpath(dirNamePlots)
    end
end

## TODO Confirm multiple MJD handling is correct/ideal
# using runlist for dome/quartz flats only
unique_mjds = if parg["runlist"] != ""
    subDic = load(parg["runlist"])
    unique(subDic["mjd"])
else
    [parg["mjd"]]
end

expid_list = if parg["runlist"] != ""
    subDic = load(parg["runlist"])
    subDic["expid"]
else
    [parg["expid"]]
end

all1Da = String[] # all 1D files for chip a
for mjd in unique_mjds
    df = h5open(parg["trace_dir"] * "almanac/$(parg["runname"]).h5") do f
        DataFrame(read(f["$(parg["tele"])/$(mjd)/exposures"]))
    end
    function get_1d_name_partial(expid)
        parg["trace_dir"] * "apred/$(mjd)/" * get_1d_name(expid, df, cal=true) * ".jld2"
    end

    file_list = get_1d_name_partial.(expid_list)
    append!(all1Da, file_list)
end

## need to get cal_type from runlist
exp_type_lst = map(x->split(split(x,"FLAT")[1],"_")[end],all1Da)
unique_exp_lst = unique(exp_type_lst)
if length(unique_exp_lst) > 1
    error("Multiple cal types found in runlist")
end
cal_type = lowercase(unique_exp_lst[1])
@passobj 1 workers() cal_type # make cal_type available to all workers

flist_chips  = []
for chip in chips
    push!(flist_chips, replace.(all1Da, "_a_"=>"_$(chip)_"))
end
all1D = vcat(flist_chips...)

all1Daout = map(x->replace(replace(x, "apred"=>"$(cal_type)_flats"),"ar1Dcal" => "$(cal_type)Flux"),all1Da)
dname = dirname(all1Daout[1])
if !ispath(dname)
    mkpath(dname)
end

@everywhere begin
    function get_and_save_relFlux(fname)
        absthrpt, relthrpt, bitmsk_relthrpt = get_relFlux(fname)
        outfname = replace(replace(fname, "apred"=>"$(cal_type)_flats"),"ar1Dcal" => "$(cal_type)Flux")
        jldsave(outfname; absthrpt, relthrpt, bitmsk_relthrpt)
    end
end

desc = "get and save relFlux from DomeFlats"
@showprogress desc=desc pmap(get_and_save_relFlux,all1D)

# add plotting
thread = SlackThread()
thread("$(cal_type) relFluxing")

@everywhere begin
    function plot_relFlux(fname)
        sname = split(fname, "_")
        tele, mjd, chiploc, expid = sname[(end - 4):(end - 1)]
        xvec = if tele == "apo"
            1:300
        elseif tele == "lco"
            301:600
        else
            error("Unknown telescope: $(tele)")
        end
        absthrpt = zeros(300,length(chips))
        relthrpt = zeros(300,length(chips))
        bitmsk_relthrpt = zeros(Int,300,length(chips))
        for (cindx, chip) in enumerate(chips)
            local_fname = replace(fname, "_a_"=>"_$(chip)_")
            f = jldopen(local_fname)
            absthrpt[:,cindx] = f["absthrpt"]
            relthrpt[:,cindx] = f["relthrpt"]
            bitmsk_relthrpt[:,cindx] = f["bitmsk_relthrpt"]
            close(f)
        end
        # plot the relFlux
        fig = Figure(size=(1200, 400))
        for (cindx, chip) in enumerate(chips)
            ax = Axis(fig[1, cindx], title="RelFlux Chip $(chip)")
            msk = bitmsk_relthrpt[:,cindx].==0
            scatter!(ax, xvec[msk], relthrpt[msk,cindx], color="limegreen")
            msk = (bitmsk_relthrpt[:,cindx].& 2^0).==2^0
            scatter!(ax, xvec[msk], relthrpt[msk,cindx], color="orange")
            msk = (bitmsk_relthrpt[:,cindx].& 2^1).==2^1
            scatter!(ax, xvec[msk], relthrpt[msk,cindx], color="red")
        end
        # we should be more uniform about the naming convention
        savePath = dirNamePlots * "relFlux_$(cal_type)_$(tele)_$(mjd)_$(expid).png"
        save(savePath, fig)
        return savePath
    end
end

desc = "plotting relFlux"
savePaths = @showprogress desc=desc pmap(plot_relFlux,all1Daout)
for savePath in savePaths
    thread("Relfluxing", savePath)
end