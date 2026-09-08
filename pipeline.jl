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
        "--dfindx"
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
        default = "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_06/pass_clean/"
        "--exp_class_model"
        required = false
        help = "exposure-type classifier artifact (JLD2). A calibration input, configured by path like --caldir_darks / --gain_read_cal_dir; airflow/dags/ar_common.py sets it for production. The post-2D check RUNS BY DEFAULT. Pass an empty string to turn it off deliberately (1D products then record exp_class_status=\"notrun\"). A path that does not exist is a hard error."
        arg_type = String
        default = "/mnt/ceph/users/sdssv/work/asaydjari/cal_ref/exposure_classifier/exposure_classifier_rf_v6.jld2"
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

# The classifier artifact is a calibration input like the dark/flat/gain
# directories: a path supplied by the caller, defaulted in the arg table above,
# and set explicitly by airflow/dags/ar_common.py for production. Version
# pinning lives in that path (v6), NOT in a glob or a "newest wins" lookup.
if parg["exp_class_model"] == ""
    println("Exposure-type check: DISABLED by explicit --exp_class_model \"\". " *
            "1D products will record exp_class_status=\"notrun\".")
elseif !isfile(parg["exp_class_model"])
    # Fail fast and loudly. Silently skipping the check would hand a bulk run
    # thousands of nights of "notrun" metadata that look indistinguishable from
    # a deliberate opt-out, which is precisely the failure this wiring exists to
    # remove. If you meant to turn it off, say so with --exp_class_model "".
    error("--exp_class_model points at a file that does not exist:\n  " *
          parg["exp_class_model"] *
          "\nPass --exp_class_model \"\" to run without the exposure-type check.")
else
    println("Exposure-type check: ENABLED, model = ", parg["exp_class_model"])
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

# write out sym links in the level of folder that MUST be uniform in their cals? or a billion symlinks with dfindx
@time "Setting up the iterator" begin
    # setup the (sjd, dfindx, chip) tuples to iterate over
    # if we have a runlist, we iterate over the mjd and dfindx in the runlist
    # otherwise we iterate over the mjd, dfindx, and chips specified on the command line
    subiter = if parg["runlist"] != ""
        subDic = load(parg["runlist"])
        msk = subDic["tele"] .== parg["tele"]
        # add a filter on tele
        Iterators.product(
            Iterators.zip(subDic["mjd"][msk], subDic["dfindx"][msk]),
            string.(collect(parg["chips"]))
        )
    else
        Iterators.product(
            [(parg["mjd"], parg["dfindx"])],
            string.(collect(parg["chips"]))
        )
    end
end

# partially apply the process_3D function to everything except the (sjd, dfindx, chip) values
@everywhere process_3D_partial(((mjd, dfindx), chip)) = process_3D(
    parg["outdir"], parg["runname"], parg["tele"], mjd, dfindx, chip,
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

# ar2D chip-file path => exposure-type classifier feature vector, handed forward
# from the 2D calibration pass so the exposure-type check does not re-read the
# same images. Empty when the check is disabled, when doCal2d is false, or for
# any chip the checkpoint let the 2D pass skip; every one of those cases falls
# back to reading the file, so a miss costs time, never correctness.
exp_class_featcache = Dict{String, Vector{Float64}}()

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

    # process the 2D calibration for all exposures.
    # When the exposure-type check is enabled, this pass ALSO returns each
    # chip's classifier features, computed from the ar2D array it already holds
    # in memory. The check below then needs no second read of the same files.
    @everywhere process_2Dcal_partial(fname) = process_2Dcal(
        fname, checkpoint_mode = parg["checkpoint_mode"],
        exp_class_features = (parg["exp_class_model"] != ""))
    cal2d_feats = @showprogress desc="2D Calibration" pmap(process_2Dcal_partial, all2D)
    # chip-file path => feature vector, for the chips this pass actually
    # recomputed. Entries are absent for chips the checkpoint skipped, and the
    # exposure-type check reads those from disk instead.
    for (fname, f) in zip(all2D, cal2d_feats)
        isnothing(f) || (exp_class_featcache[fname] = f)
    end
end

##### Exposure-type check (POST-2D, PRE-1D) #####
# Classify each exposure from its ar2D images alone and compare to the
# commanded image_type + lamp flags. Mislabeled cals (ThAr/UNe swaps, FPI vs
# arclamp, lamp-on "darks") poison downstream calibrations; this writes a
# per-mjd audit table and warns on confident disagreements. Warning-only:
# labels are never changed automatically and NO exposure is dropped from the
# reduction here — engineering and known-bad frames are still reduced. The
# verdict is advisory metadata (carried into the 1D products, see
# `process_1D`) plus an input to the calibration runlists (see
# `make_runlist_fiber_flats.jl`), never a reduction filter.
#
# Placement matters: this block is the last thing pipeline.jl does, and
# pipeline.jl is the 3D->2D stage, so the verdict is written before any caller
# runs pipeline_2d_1d.jl. Every DAG path (bulk, daily, cal, regression) invokes
# pipeline.jl strictly before pipeline_2d_1d.jl, so "after 2D, before 1D" holds
# by construction rather than by convention.
if parg["exp_class_model"] != ""
    @everywhere begin
        using ApogeeReduction: exposure_class_features, load_exposure_classifier,
                               classify_exposure_type, exposure_class_label,
                               exposure_check_category, exposure_class_metadata,
                               read_almanac_exp_df,
                               EXPFLAG_PREDICTED_BAD, EXPFLAG_NOTRUN,
                               CHIP_LIST
        const EXP_CLF = Ref{Any}(nothing)
        function get_exp_clf()
            if EXP_CLF[] === nothing
                EXP_CLF[] = load_exposure_classifier(parg["exp_class_model"])
            end
            EXP_CLF[]
        end
        # `work` is (chip => ar2D path, chip => precomputed features). Features
        # present for a chip are used as-is; absent ones are read from the file,
        # so a checkpointed rerun that skipped the 2D pass still gets a verdict
        # rather than "notrun".
        function exp_type_check_one(work)
            fnames_by_chip, feats_by_chip = work
            # any chip filename parses to (tele, mjd, expnum, imtype)
            sname = split(split(basename(first(values(fnames_by_chip))), ".h5")[1], "_")
            _, tele, mjdstr, expnumstr, _, _ = sname[(end - 5):end]
            mjd, expnum = parse(Int, mjdstr), parse(Int, expnumstr)
            # the check must never break a reduction run
            try
                df = read_almanac_exp_df(
                    joinpath(parg["outdir"], "almanac/$(parg["runname"]).h5"), tele, mjd)
                erow = df[expnum, :]
                labeled = exposure_class_label(erow.image_type,
                    get(erow, :lamp_quartz, "?"), get(erow, :lamp_thar, "?"),
                    get(erow, :lamp_une, "?"))
                clf = get_exp_clf()
                # All three chips are required: the design row concatenates
                # R, G and B and the colour features are inter-chip ratios, so
                # a per-chip partial verdict would be meaningless. The grouping
                # below only submits exposures with a complete chip set.
                chip_features = Dict(
                    chip => get(feats_by_chip, chip, nothing) !== nothing ?
                            feats_by_chip[chip] :
                            exposure_class_features(load(fnames_by_chip[chip], "dimage"))
                for chip in keys(fnames_by_chip))
                res = classify_exposure_type(clf, chip_features, tele)
                flag = exposure_check_category(labeled, res, clf.flag_tau)
                # persistence prior: a dark taken right after a bright exposure
                # may carry afterglow (recorded in the table, not warned)
                if flag == "ok" && res.pred == "dark_q0t0u0" &&
                   labeled == "dark_q0t0u0" && expnum > 1
                    prevtype = lowercase(string(df.image_type[expnum - 1]))
                    if prevtype in
                       ApogeeReduction.CLASSIFIER_PERSIST_SOURCES
                        flag = "persistence_prior"
                    end
                end
                (tele = String(tele), mjd = mjd, expnum = expnum, labeled = labeled,
                    pred = res.pred, prob = res.prob, flag = flag)
            catch e
                @warn "Exposure-type check failed for $tele $mjd exp $expnum" exception=e
                (tele = String(tele), mjd = mjd, expnum = expnum, labeled = "checkfail",
                    pred = "checkfail", prob = NaN, flag = "checkfail")
            end
        end
    end

    @time "Exposure-type check" begin
        # group ar2D filenames by exposure (drop the chip distinction)
        groups = Dict{String, Dict{String, String}}()
        for fname in all2D
            isfile(fname) || continue
            sname = split(split(basename(fname), ".h5")[1], "_")
            chip = String(sname[end - 1])
            key = join(vcat(sname[(end - 5):(end - 2)], sname[end]), "_")
            get!(groups, key, Dict{String, String}())[chip] = fname
        end
        complete = [g for g in values(groups) if length(g) == length(CHIP_LIST)]
        # Pair each exposure's chip files with whatever the 2D pass already
        # computed. Feature vectors are 21 Float64s, so shipping them to the
        # workers is free next to re-reading a 32 MB image.
        work = map(complete) do g
            fx = Dict{String, Vector{Float64}}()
            for (chip, path) in g
                f = get(exp_class_featcache, path, nothing)
                isnothing(f) || (fx[chip] = f)
            end
            (g, fx)
        end
        n_chip_cached = sum(length(w[2]) for w in work; init = 0)
        n_chip_total = sum(length(w[1]) for w in work; init = 0)
        println("Exposure-type check: $(n_chip_cached)/$(n_chip_total) chip images reused " *
                "from the 2D stage; $(n_chip_total - n_chip_cached) re-read from disk")
        checks = @showprogress desc="Exposure-type check" pmap(exp_type_check_one, work)

        if !isempty(checks)
            dfc = DataFrame(checks)
            # Resolve the downstream masking policy HERE, once, and store it, so
            # that the 1D products, the almanac decoration and the calibration
            # runlists all read one answer rather than each re-deriving it from
            # (labeled, pred, flag) and risking three subtly different ones.
            #
            # One column: EXPFLAG_NOTRUN inside the mask carries "no verdict",
            # so a zero byte means judged-and-fine and needs no companion.
            # UInt8[...] on purpose: the Dict is Dict{String,Any}, so an
            # unannotated comprehension would build a Vector{Any} column and
            # write an untyped HDF5 dataset.
            dfc.exposure_flags = UInt8[exposure_class_metadata(
                                           r.labeled, r.pred, r.prob, r.flag)["exposure_flags"]
                                       for r in eachrow(dfc)]
            for sub in groupby(dfc, :mjd)
                mjd = sub.mjd[1]
                safe_jldsave(
                    joinpath(parg["outdir"],
                        "apred/$(mjd)/exposureTypeCheck_$(parg["tele"])_$(mjd).h5"),
                    Dict{String, Any}();
                    tele = String.(sub.tele), mjd = collect(sub.mjd),
                    expnum = collect(sub.expnum), labeled = String.(sub.labeled),
                    pred = String.(sub.pred), prob = collect(sub.prob),
                    flag = String.(sub.flag),
                    exposure_flags = collect(sub.exposure_flags))
            end
            # persistence_prior is informational (recorded, not warned)
            warnable = (dfc.flag .!= "ok") .& (dfc.flag .!= "persistence_prior")
            nbad = sum(warnable)
            if nbad > 0
                bad = dfc[warnable, :]
                for r in eachrow(bad)
                    @warn "Exposure-type check: $(r.tele) $(r.mjd) exp $(r.expnum) labeled '$(r.labeled)' classified '$(r.pred)' (p=$(round(r.prob, digits = 3)), $(r.flag))"
                end
                nflag_frac = nbad / nrow(dfc)
                nflag_frac > 0.10 &&
                    @warn "Exposure-type check flagged $(round(100 * nflag_frac, digits = 1))% of exposures (> 10% prior on mislabel rate) — inspect before trusting the flags."
            end
            npbad = sum((dfc.exposure_flags .& EXPFLAG_PREDICTED_BAD) .!= 0x00)
            nunk = sum((dfc.exposure_flags .& EXPFLAG_NOTRUN) .!= 0x00)
            println("Exposure-type check: $(nrow(dfc)) exposures checked, $nbad flagged, " *
                    "$npbad predicted_bad, $nunk unjudged (check failed)")
            println("Exposure-type check: NOTHING was dropped from the reduction here; " *
                    "the verdict is advisory 1D metadata plus a calibration-runlist input.")
        end
    end
end
