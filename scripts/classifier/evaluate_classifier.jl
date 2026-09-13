# Evaluate the deployed (v6, predictions carried in the manifest from the
# almanac decoration) vs retrained (v7, cross-validated on train_pool rows /
# held-out final model on eval_only rows) exposure-type classifier on each
# weakness class the DR21 pass-1 sweep surfaced:
#   - twilight census truth classes (nightsky / twilight / daytime)
#   - special_cal_obs.txt sequences (uD/uF/sat/spQF/b2b/twL)
#   - engineering-carton exposures (no content truth; movement only)
#   - FPS-era config_id=-1 object frames
#   - lco 57800 stale frames (adversarial; includes an exact-duplicate
#     feature-vector check, which is what actually catches copied frames)
#
# Fairness note: v6 predictions on train_pool rows may include exposures v6
# itself trained on (favoring v6); v7 numbers on train_pool rows are always
# from the fold that did not train on the row.
using Pkg;
Pkg.instantiate();
using HDF5, ArgParse, DataFrames, CSV, StatsBase, Printf

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--manifest_final"
        required = true
        help = "final manifest CSV from train_exposure_classifier.jl"
        arg_type = String
        "--features"
        required = true
        arg_type = String
        "--outdir"
        required = true
        help = "directory for the evaluation CSV tables"
        arg_type = String
    end
    return parse_args(s)
end
parg = parse_commandline()
mkpath(parg["outdir"])

man = CSV.read(parg["manifest_final"], DataFrame)
man.v6_pred = coalesce.(man.v6_pred, "")
man.v7_pred = coalesce.(man.v7_pred, "")
ok = man.feat_ok .== true
println("manifest: ", nrow(man), ", features ok: ", sum(ok))

astwi(p) = startswith.(p, "twilightflat")
short(c) = c == "" ? "(none)" :
           replace(c, "arclamp_q0t1u0" => "thar", "arclamp_q0t0u1" => "une",
    "arclamp_q0t0u0" => "fpi", "twilightflat_faint" => "twi_faint",
    "twilightflat_q0t0u0" => "twilight", r"_q[\d?].*" => "")

# ---------- 1. twilight census ----------
println("\n================ TWILIGHT CENSUS ================")
res = DataFrame(truth = String[], n = Int[], v6_twi = Int[], v7_twi = Int[],
    v6_correct = Int[], v7_correct = Int[])
for t in ("nightsky", "twilight", "daytime")
    sub = man[coalesce.(man.census_truth .== t, false) .& ok, :]
    n = nrow(sub)
    n == 0 && continue
    corr6 = t == "nightsky" ? sum(sub.v6_pred .== "twilightflat_faint") :
            sum(astwi(sub.v6_pred))
    corr7 = t == "nightsky" ? sum(sub.v7_pred .== "twilightflat_faint") :
            sum(astwi(sub.v7_pred))
    push!(res, (t, n, sum(astwi(sub.v6_pred)), sum(astwi(sub.v7_pred)), corr6, corr7))
end
res.v6_pct = round.(100 .* res.v6_correct ./ res.n, digits = 1)
res.v7_pct = round.(100 .* res.v7_correct ./ res.n, digits = 1)
show(res, allrows = true, allcols = true)
CSV.write(joinpath(parg["outdir"], "twilight_v6v7.csv"), res)

for t in ("nightsky", "twilight", "daytime")
    sub = man[coalesce.(man.census_truth .== t, false) .& ok, :]
    nrow(sub) == 0 && continue
    println("\n--- $t: prediction movement (v6 -> v7) ---")
    mv = sort(combine(groupby(sub, [:v6_pred, :v7_pred]), nrow => :count),
        order(:count, rev = true))
    mv.v6_pred = short.(mv.v6_pred)
    mv.v7_pred = short.(mv.v7_pred)
    show(mv, allrows = true)
    println()
end

# ---------- 2. special_cal sequences ----------
println("\n================ SPECIAL CAL ================")
expect = Dict("uD" => ["dark_q0t0u0", "dark_persist"],
    "uF" => ["internalflat_q0t0u0"], "sat" => ["internalflat_q0t0u0"],
    "spQF" => ["quartzflat_q1t0u0"], "b2b" => ["quartzflat_q1t0u0", "domeflat_q0t0u0"],
    "rdL" => String[], "twL" => ["twilightflat_q0t0u0", "twilightflat_faint"])
sc = man[startswith.(man.truth_source, "special_cal_") .& ok, :]
sc.cal_type = replace.(sc.truth_source, "special_cal_" => "")
sres = DataFrame(cal_type = String[], n = Int[], v6_expected = Int[], v7_expected = Int[])
for g in groupby(sc, :cal_type)
    e = get(expect, g.cal_type[1], String[])
    isempty(e) && continue
    push!(sres, (g.cal_type[1], nrow(g), sum(in.(g.v6_pred, Ref(e))),
        sum(in.(g.v7_pred, Ref(e)))))
end
# note: uF includes the documented lamp-off/broken sequence (apo 58305) whose
# true content IS dark — expectation counts are a floor, not a target of 100%
show(sort(sres, :cal_type), allrows = true)
CSV.write(joinpath(parg["outdir"], "specialcal_v6v7.csv"), sres)
CSV.write(joinpath(parg["outdir"], "specialcal_movement.csv"),
    sort(combine(groupby(sc, [:cal_type, :labeled, :v6_pred, :v7_pred]), nrow => :count),
        [:cal_type, order(:count, rev = true)]))

# ---------- 3. engineering ----------
println("\n\n================ ENGINEERING (eval-only; no content truth) ================")
eng = man[(man.truth_source .== "engineering") .& ok, :]
emv = sort(combine(groupby(eng, [:v6_pred, :v7_pred]), nrow => :count),
    order(:count, rev = true))
emv.v6_pred = short.(emv.v6_pred)
emv.v7_pred = short.(emv.v7_pred)
show(emv[1:min(15, nrow(emv)), :], allrows = true)
CSV.write(joinpath(parg["outdir"], "engineering_movement.csv"), emv)
agree = sum(eng.v6_pred .== eng.v7_pred)
println("\nengineering: n=", nrow(eng), " v6==v7 agreement: ", agree, " (",
    round(100agree / max(nrow(eng), 1), digits = 1), "%)")

# ---------- 4. FPS config_id = -1 ----------
println("\n================ FPS config_id=-1 object frames ================")
cfg = man[(man.truth_source .== "fps_configless") .& ok, :]
cmv = sort(combine(groupby(cfg, [:v6_pred, :v7_pred]), nrow => :count),
    order(:count, rev = true))
cmv.v6_pred = short.(cmv.v6_pred)
cmv.v7_pred = short.(cmv.v7_pred)
show(cmv, allrows = true)
CSV.write(joinpath(parg["outdir"], "fps_configless_movement.csv"), cmv)
println("\nn with features: ", nrow(cfg), " of ",
    sum(man.truth_source .== "fps_configless"))

# ---------- 5. lco 57800 stale frames ----------
println("\n================ LCO 57800 STALE FRAMES ================")
fk = man[man.truth_source .== "lco57800_stale", :]
show(fk[:, [:mjd, :expnum, :labeled, :v6_pred, :v6_prob, :v7_pred, :v7_prob, :feat_ok]],
    allrows = true, allcols = true)
F = h5open(parg["features"], "r") do f
    read(f["features"])
end
idx = findall(man.truth_source .== "lco57800_stale")
println("\nfeature-vector duplicate groups among the stale frames (exact match):")
groups = Dict{Vector{Float64}, Vector{Int}}()
for j in idx
    push!(get!(groups, vec(F[:, :, j]), Int[]), man.expnum[j])
end
for v in values(groups)
    println("  identical group: expnums ", sort(v))
end

# ---------- 6. dark_persist (drives dark-frame masking) ----------
# `exposure_predicted_bad` masks a commanded dark whenever the prediction is not
# dark_q0t0u0, so the dark_persist class is the mechanism that keeps
# persistence-contaminated darks out of dark cal stacks. Its labels are
# constructed, not commanded, and the yield is ~2% of prev-bright darks — so it
# is the class most sensitive to how the training pool was sampled. Compare the
# two models on that population directly.
println("\n================ DARK_PERSIST (dark-frame masking) ================")
dk = man[(man.labeled .== "dark_q0t0u0") .& ok, :]
dres = DataFrame(model = String[], n = Int[], pred_dark = Int[], pred_persist = Int[],
    pred_other = Int[])
for (nm, p) in (("v6", dk.v6_pred), ("v7", dk.v7_pred))
    push!(dres, (nm, nrow(dk), sum(p .== "dark_q0t0u0"), sum(p .== "dark_persist"),
        sum((p .!= "dark_q0t0u0") .& (p .!= "dark_persist"))))
end
show(dres, allrows = true, allcols = true)
CSV.write(joinpath(parg["outdir"], "darkpersist_v6v7.csv"), dres)
np = sum(man.content_label .== "dark_persist")
nprec = sum((man.content_label .== "dark_persist") .& (man.v7_pred .== "dark_persist"))
println("\nv7 dark_persist content labels: ", np, "; held-out recall: ", nprec, "/", np)
