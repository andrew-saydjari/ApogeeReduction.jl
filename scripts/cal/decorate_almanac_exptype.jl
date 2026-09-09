using Pkg;
Pkg.instantiate();
using HDF5, ArgParse, DataFrames, JLD2
using ApogeeReduction: exposure_class_label, exposure_class_metadata, initalize_git,
                       exposure_class_verdict, EXP_CLASS_UNKNOWN_STR,
                       exposure_predicted_bad, exposure_flag_bits,
                       exposure_engineering_from_almanac,
                       ENGINEERING_CARTON_PREFIXES, ENGINEERING_CARTON_PURITY,
                       ENGINEERING_CHECK_IMAGE_TYPES,
                       ENGINEERING_FLAG_SCIENCELESS_POSITION_STARS,
                       apply_designless_clause, ENGINEERING_FLAG_DESIGNLESS,
                       ENGINEERING_DESIGNLESS_DESIGN_ID,
                       EXPFLAG_PREDICTED_BAD, EXPFLAG_ENGINEERING, EXPFLAG_NOTRUN

# recompute at runtime: the module-level git consts are frozen at precompile
# time and can go stale (see comment in src/utils.jl)
git_branch, git_commit, git_clean = initalize_git(dirname(dirname(@__DIR__)) * "/")

## Decorate an almanac file with the exposure-type classifier verdicts,
## in the same spirit as `almanac metadata` decorating with photometry:
## a sibling top-level group `exposure_class/<tele>/<mjd>` is (re)written
## with one dataset per column, aligned row-for-row with
## `raw/<tele>/<mjd>/exposures`. Nothing under `raw/` is ever touched.
##
## Two input modes:
##   --apred_dir <outdir>/apred   in-pipeline: gathers the per-MJD
##                                exposureTypeCheck_*.h5 tables pipeline.jl
##                                writes between the 2D and 1D stages. This is
##                                the mode the DAGs use.
##   --results_file <sweep.h5>    offline: a classifier sweep audit table.
##
## Columns per (tele, mjd):
##   exposure               exposure number (join key, matches exposures group)
##   exposure_flags         UInt8 bitmask, written by BOTH producers (PR #397 +
##                          PR #398): 2^0 predicted_bad, 2^1 engineering,
##                          2^2 notrun (no verdict formed). exposure_flags == 0
##                          therefore means JUDGED AND FINE; "never judged" is
##                          2^2, a distinct value. 2^0 and 2^2 are mutually
##                          exclusive. Every row that had no classifier verdict
##                          — including exposures with no reduced 2D data —
##                          gets 2^2. 2^1 is INDEPENDENT of the classifier and
##                          can therefore co-occur with 2^2: an engineering
##                          exposure the classifier never judged is 2^1|2^2.
##   exposure_class_pred    predicted content class ("unknown" if no verdict)
##   exposure_class_prob    max forest probability (NaN if no verdict)
##   exposure_class_status  ok / mislabel_candidate / lamp_off_candidate /
##                          persistence_risk / faint_twilight / unknown /
##                          rare_label / nofiles / unclassified
##   engineering            UInt8 0/1 — carton check. Clause 1: ALL of the
##                          configuration's science fibers carry an engineering
##                          carton (ENGINEERING_CARTON_PREFIXES, purity rule, see
##                          ENGINEERING_CARTON_PURITY). Clause 2: the config has
##                          ZERO science fibers and ANY fiber carries one (see
##                          ENGINEERING_FLAG_SCIENCELESS_POSITION_STARS)
##   engineering_frac       Float64 fraction of science fibers matching (NaN
##                          when there is no configuration: plate era, cals)
##   engineering_carton     the matched carton name ("" if none), for audit
##   engineering_basis      which fibers the fraction was computed over:
##                          "science" (clause 1: purity over the science fibers),
##                          "scienceless_position_stars" (clause 2: zero science
##                          fibers and any fiber carries the carton), or "none"
##                          (clause 3, designless, reports basis "designless")
## The `exposure_class` group carries git branch/commit/clean and the source
## path (results_file / apred_dir) as attributes, so the model + policy version
## is pinned to the pipeline git hash. The engineering policy constants are
## written as attributes too, so the meaning of the bit travels with the file.
##
## NOTE ON SCOPE: every column here is ADVISORY METADATA for downstream
## consumers (prior builds, cal runlists, catalog construction). Nothing in the
## 3D->2D->1D reduction reads it: engineering exposures are still reduced to 1D
## in full, exactly as before. `scripts/bulk/make_runlist_all.jl` deliberately
## does not consult these columns.

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--almanac_file"
        required = true
        help = "path to the almanac file to decorate (modified in place)"
        arg_type = String
        "--results_file"
        required = false
        help = "offline classifier results table (sweep audit): obs/mjd/expnum/pred/maxp/status"
        arg_type = String
        default = ""
        "--apred_dir"
        required = false
        help = "in-pipeline mode: reduction apred/ directory; gathers every apred/<mjd>/exposureTypeCheck_*.h5 written by pipeline.jl between the 2D and 1D stages"
        arg_type = String
        default = ""
    end
    return parse_args(s)
end

parg = parse_commandline()

if (parg["results_file"] == "") == (parg["apred_dir"] == "")
    error("give exactly one of --results_file (offline sweep) or --apred_dir (in-pipeline exposureTypeCheck tables)")
end

# classifier verdicts keyed by (obs, mjd, expnum) => (pred, maxp, status)
verdict = if parg["results_file"] != ""
    res = h5open(parg["results_file"], "r") do f
        DataFrame(obs = read(f["obs"]), mjd = read(f["mjd"]), expnum = read(f["expnum"]),
            pred = read(f["pred"]), maxp = read(f["maxp"]), status = read(f["status"]))
    end
    println("source: offline results file ", parg["results_file"])
    Dict(zip(zip(res.obs, res.mjd, res.expnum),
        zip(res.pred, res.maxp, res.status)))
else
    # In-pipeline mode: the per-MJD tables pipeline.jl writes after the 2D stage.
    # Column names differ from the offline sweep (tele/flag vs obs/status).
    d = Dict{Tuple{String, Int, Int}, Tuple{String, Float64, String}}()
    files = String[]
    for mjddir in readdir(parg["apred_dir"]; join = true)
        isdir(mjddir) || continue
        append!(files, filter(p -> occursin("exposureTypeCheck_", basename(p)) &&
                                  endswith(p, ".h5"),
            readdir(mjddir; join = true)))
    end
    println("source: $(length(files)) in-pipeline exposureTypeCheck table(s) under ",
        parg["apred_dir"])
    for p in files
        t = load(p)
        for i in eachindex(t["expnum"])
            d[(String(t["tele"][i]), Int(t["mjd"][i]), Int(t["expnum"][i]))] = (
                String(t["pred"][i]), Float64(t["prob"][i]), String(t["flag"][i]))
        end
    end
    d
end
println("classifier verdicts: ", length(verdict))
isempty(verdict) &&
    @warn "no classifier verdicts found — every exposure will be decorated as UNJUDGED (exposure_flags bit 2^2 notrun), which is correct but means nothing downstream will be filtered. Was --exp_class_model set on the pipeline.jl call?"

nbad = 0
nunknown = 0
neng = 0
nobj = 0
nobj_nocfg = 0
ntot = 0
eng_carton_counts = Dict{String, Int}()
h5open(parg["almanac_file"], "r+") do f
    rawgrp = haskey(f, "raw") ? "raw" : ""
    haskey(f, "exposure_class") && delete_object(f, "exposure_class")
    g = create_group(f, "exposure_class")
    attrs(g)["exposure_flags_bits"] = "2^0=predicted_bad,2^1=engineering,2^2=notrun"
    attrs(g)["git_branch"] = git_branch
    attrs(g)["git_commit"] = string(git_commit)
    attrs(g)["git_clean"] = string(git_clean)
    attrs(g)["results_file"] = parg["results_file"] == "" ? "" : abspath(parg["results_file"])
    attrs(g)["apred_dir"] = parg["apred_dir"] == "" ? "" : abspath(parg["apred_dir"])
    # engineering-carton policy travels with the file so the bit is self-describing
    attrs(g)["engineering_carton_prefixes"] = join(ENGINEERING_CARTON_PREFIXES, ",")
    attrs(g)["engineering_carton_purity"] = ENGINEERING_CARTON_PURITY
    attrs(g)["engineering_check_image_types"] = join(ENGINEERING_CHECK_IMAGE_TYPES, ",")
    attrs(g)["engineering_flag_scienceless_position_stars"] =
        string(ENGINEERING_FLAG_SCIENCELESS_POSITION_STARS)
    attrs(g)["engineering_flag_designless"] = string(ENGINEERING_FLAG_DESIGNLESS)
    attrs(g)["engineering_designless_design_id"] = ENGINEERING_DESIGNLESS_DESIGN_ID
    for tele in keys(f[rawgrp == "" ? "/" : rawgrp])
        tele in ("exposure_class", "meta") && continue
        tele_out = create_group(g, tele)
        for mjd in keys(f[joinpath(rawgrp, tele)])
            haskey(f[joinpath(rawgrp, tele, mjd)], "exposures") || continue
            exp_grp = f[joinpath(rawgrp, tele, mjd, "exposures")]
            expnum = read(exp_grp["exposure"])
            imtype = read(exp_grp["image_type"])
            lq, lt, lu = read(exp_grp["lamp_quartz"]), read(exp_grp["lamp_thar"]),
            read(exp_grp["lamp_une"])
            # plate-era files may not carry config_id at all; -1 == no configuration
            cfgid = haskey(exp_grp, "config_id") ? read(exp_grp["config_id"]) :
                    fill(-1, length(expnum))
            # plate-era files carry no design_id; `nothing` reads as "unknown",
            # which is_designless_design deliberately does NOT treat as -999
            designid = haskey(exp_grp, "design_id") ? read(exp_grp["design_id"]) :
                       fill(nothing, length(expnum))
            n = length(expnum)
            pred = fill(EXP_CLASS_UNKNOWN_STR, n)
            prob = fill(NaN, n)
            status = fill("unclassified", n)
            # Default every row to NOTRUN, not to zero: a row we never formed a
            # verdict for must not decorate as "judged and fine". Rows with a
            # verdict overwrite the classifier bits below; the engineering bit
            # is ORed in afterwards and is independent of the classifier.
            flags = fill(EXPFLAG_NOTRUN, n)
            eng = falses(n)
            engfrac = fill(NaN, n)
            engcarton = fill("", n)
            engbasis = fill("none", n)
            # exposures on the same night share configurations; read each fiber
            # table at most once (the DR21 corpus has ~48k configurations)
            engcache = Dict{Int, NamedTuple}()
            for i in 1:n
                # classifier verdict (only for exposures that were classified)
                v = get(verdict, (tele, parse(Int, mjd), expnum[i]), nothing)
                if !isnothing(v)
                    pred[i], prob[i], status[i] = v
                    labeled = exposure_class_label(imtype[i], lq[i], lt[i], lu[i])
                    # a checkfail verdict comes back with NOTRUN still set
                    flags[i] = UInt8(exposure_class_metadata(
                        labeled, pred[i], prob[i], status[i])["exposure_flags"])
                end
                # engineering carton check — INDEPENDENT of the image classifier:
                # it runs on every object exposure whether or not a 2D prediction
                # exists, because it is derived from targeting, not from pixels.
                is_obj = lowercase(strip(String(imtype[i]))) in ENGINEERING_CHECK_IMAGE_TYPES
                e = if !is_obj
                    (engineering = false, frac = NaN, carton = "", nsci = 0,
                        basis = "none")
                else
                    # design_id is per EXPOSURE; clause 3 is applied outside the
                    # per-config cache so one exposure's design cannot leak onto
                    # its config-mates
                    apply_designless_clause(
                        get!(engcache, Int(cfgid[i])) do
                            exposure_engineering_from_almanac(f, tele, mjd, cfgid[i],
                                imtype[i]; root = rawgrp)
                        end, designid[i])
                end
                engbasis[i] = e.basis
                eng[i] = e.engineering
                engfrac[i] = e.frac
                engcarton[i] = e.carton
                if lowercase(strip(String(imtype[i]))) == "object"
                    global nobj += 1
                    e.nsci == 0 && (global nobj_nocfg += 1)
                end
                if e.engineering
                    eng_carton_counts[e.carton] = get(eng_carton_counts, e.carton, 0) + 1
                end
                # The engineering bit is INDEPENDENT of the classifier, so it is
                # ORed in last rather than passed through exposure_class_metadata:
                # on the checkfail/notrun path that helper routes to
                # exposure_class_unknown_metadata, which hardcodes
                # engineering = false and would silently drop the bit. 2^1 may
                # therefore co-occur with 2^2, which is correct and intended.
                e.engineering && (flags[i] |= EXPFLAG_ENGINEERING)
            end
            out = create_group(tele_out, mjd)
            out["exposure"] = expnum
            out["exposure_flags"] = flags
            out["exposure_class_pred"] = pred
            out["exposure_class_prob"] = prob
            out["exposure_class_status"] = status
            out["engineering"] = UInt8.(eng)
            out["engineering_frac"] = engfrac
            out["engineering_carton"] = engcarton
            out["engineering_basis"] = engbasis
            global nbad += sum((flags .& EXPFLAG_PREDICTED_BAD) .!= 0x00)
            global nunknown += sum((flags .& EXPFLAG_NOTRUN) .!= 0x00)
            global neng += sum(eng)
            global ntot += n
        end
    end
end
println("decorated $(parg["almanac_file"]): $ntot exposures, $nbad predicted_bad, " *
        "$nunknown unknown (no classifier verdict), " *
        "$neng engineering (of $nobj object exposures; $nobj_nocfg had no readable " *
        "configuration/carton table and are therefore NOT engineering by construction)")
if !isempty(eng_carton_counts)
    println("engineering exposures by dominant carton:")
    for (k, v) in sort(collect(eng_carton_counts), by = last, rev = true)
        println("  $(k): $(v)")
    end
end
