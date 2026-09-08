using Pkg;
Pkg.instantiate();
using HDF5, ArgParse, DataFrames
using ApogeeReduction: exposure_class_label, exposure_predicted_bad, initalize_git,
                       exposure_engineering_from_almanac, exposure_flag_bits,
                       ENGINEERING_CARTON_PREFIXES, ENGINEERING_CARTON_PURITY,
                       ENGINEERING_CHECK_IMAGE_TYPES, EXPFLAG_PREDICTED_BAD,
                       EXPFLAG_ENGINEERING,
                       ENGINEERING_FLAG_SCIENCELESS_POSITION_STARS

# recompute at runtime: the module-level git consts are frozen at precompile
# time and can go stale (see comment in src/utils.jl)
git_branch, git_commit, git_clean = initalize_git(dirname(dirname(@__DIR__)) * "/")

## Decorate an almanac file with the exposure-type classifier verdicts,
## in the same spirit as `almanac metadata` decorating with photometry:
## a sibling top-level group `exposure_class/<tele>/<mjd>` is (re)written
## with one dataset per column, aligned row-for-row with
## `raw/<tele>/<mjd>/exposures`. Nothing under `raw/` is ever touched.
##
## Columns per (tele, mjd):
##   exposure               exposure number (join key, matches exposures group)
##   predicted_bad          UInt8 0/1 mask for downstream cal/wavecal runlists
##                          (policy: ApogeeReduction.exposure_predicted_bad)
##   exposure_class_pred    predicted content class ("" if unclassified)
##   exposure_class_prob    max forest probability (NaN if unclassified)
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
##   exposure_flags         UInt8 bitmask: 2^0 predicted_bad, 2^1 engineering
##                          (see README "Exposure-Level Flag Bits")
## The `exposure_class` group carries git branch/commit/clean and the results
## file path as attributes, so the model + policy version is pinned to the
## pipeline git hash. The engineering policy constants are written as
## attributes too, so the meaning of the bit travels with the file.
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
        required = true
        help = "classifier results table (sweep audit or merged exposureTypeCheck): obs/mjd/expnum/pred/maxp/status"
        arg_type = String
    end
    return parse_args(s)
end

parg = parse_commandline()

# classifier verdicts keyed by (obs, mjd, expnum)
res = h5open(parg["results_file"], "r") do f
    DataFrame(obs = read(f["obs"]), mjd = read(f["mjd"]), expnum = read(f["expnum"]),
        pred = read(f["pred"]), maxp = read(f["maxp"]), status = read(f["status"]))
end
verdict = Dict(zip(zip(res.obs, res.mjd, res.expnum),
    zip(res.pred, res.maxp, res.status)))
println("classifier verdicts: ", length(verdict))

nbad = 0
neng = 0
nobj = 0
nobj_nocfg = 0
ntot = 0
eng_carton_counts = Dict{String, Int}()
h5open(parg["almanac_file"], "r+") do f
    rawgrp = haskey(f, "raw") ? "raw" : ""
    haskey(f, "exposure_class") && delete_object(f, "exposure_class")
    g = create_group(f, "exposure_class")
    attrs(g)["git_branch"] = git_branch
    attrs(g)["git_commit"] = string(git_commit)
    attrs(g)["git_clean"] = string(git_clean)
    attrs(g)["results_file"] = abspath(parg["results_file"])
    # engineering-carton policy travels with the file so the bit is self-describing
    attrs(g)["engineering_carton_prefixes"] = join(ENGINEERING_CARTON_PREFIXES, ",")
    attrs(g)["engineering_carton_purity"] = ENGINEERING_CARTON_PURITY
    attrs(g)["engineering_check_image_types"] = join(ENGINEERING_CHECK_IMAGE_TYPES, ",")
    attrs(g)["engineering_flag_scienceless_position_stars"] =
        string(ENGINEERING_FLAG_SCIENCELESS_POSITION_STARS)
    attrs(g)["exposure_flags_bits"] = "2^0=predicted_bad,2^1=engineering"
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
            n = length(expnum)
            pred = fill("", n)
            prob = fill(NaN, n)
            status = fill("unclassified", n)
            bad = falses(n)
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
                    bad[i] = exposure_predicted_bad(labeled, pred[i], status[i])
                end
                # engineering carton check — INDEPENDENT of the image classifier:
                # it runs on every object exposure whether or not a 2D prediction
                # exists, because it is derived from targeting, not from pixels.
                is_obj = lowercase(strip(String(imtype[i]))) in ENGINEERING_CHECK_IMAGE_TYPES
                e = if !is_obj
                    (engineering = false, frac = NaN, carton = "", nsci = 0,
                        basis = "none")
                else
                    get!(engcache, Int(cfgid[i])) do
                        exposure_engineering_from_almanac(f, tele, mjd, cfgid[i],
                            imtype[i]; root = rawgrp)
                    end
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
            end
            out = create_group(tele_out, mjd)
            out["exposure"] = expnum
            out["predicted_bad"] = UInt8.(bad)
            out["exposure_class_pred"] = pred
            out["exposure_class_prob"] = prob
            out["exposure_class_status"] = status
            out["engineering"] = UInt8.(eng)
            out["engineering_frac"] = engfrac
            out["engineering_carton"] = engcarton
            out["engineering_basis"] = engbasis
            out["exposure_flags"] = exposure_flag_bits.(bad, eng)
            global nbad += sum(bad)
            global neng += sum(eng)
            global ntot += n
        end
    end
end
println("decorated $(parg["almanac_file"]): $ntot exposures, $nbad predicted_bad, " *
        "$neng engineering (of $nobj object exposures; $nobj_nocfg had no readable " *
        "configuration/carton table and are therefore NOT engineering by construction)")
if !isempty(eng_carton_counts)
    println("engineering exposures by dominant carton:")
    for (k, v) in sort(collect(eng_carton_counts), by = last, rev = true)
        println("  $(k): $(v)")
    end
end
