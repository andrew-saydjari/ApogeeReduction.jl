using Pkg;
Pkg.instantiate();
using HDF5, ArgParse, DataFrames
using ApogeeReduction: exposure_class_label, exposure_predicted_bad, initalize_git

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
## The `exposure_class` group carries git branch/commit/clean and the results
## file path as attributes, so the model + policy version is pinned to the
## pipeline git hash.

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
ntot = 0
h5open(parg["almanac_file"], "r+") do f
    rawgrp = haskey(f, "raw") ? "raw" : ""
    haskey(f, "exposure_class") && delete_object(f, "exposure_class")
    g = create_group(f, "exposure_class")
    attrs(g)["git_branch"] = git_branch
    attrs(g)["git_commit"] = string(git_commit)
    attrs(g)["git_clean"] = string(git_clean)
    attrs(g)["results_file"] = abspath(parg["results_file"])
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
            n = length(expnum)
            pred = fill("", n)
            prob = fill(NaN, n)
            status = fill("unclassified", n)
            bad = falses(n)
            for i in 1:n
                v = get(verdict, (tele, parse(Int, mjd), expnum[i]), nothing)
                isnothing(v) && continue
                pred[i], prob[i], status[i] = v
                labeled = exposure_class_label(imtype[i], lq[i], lt[i], lu[i])
                bad[i] = exposure_predicted_bad(labeled, pred[i], status[i])
            end
            out = create_group(tele_out, mjd)
            out["exposure"] = expnum
            out["predicted_bad"] = UInt8.(bad)
            out["exposure_class_pred"] = pred
            out["exposure_class_prob"] = prob
            out["exposure_class_status"] = status
            global nbad += sum(bad)
            global ntot += n
        end
    end
end
println("decorated $(parg["almanac_file"]): $ntot exposures, $nbad predicted_bad")
