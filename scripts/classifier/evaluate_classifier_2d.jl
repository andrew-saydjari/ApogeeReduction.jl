# Evaluate a 2D content forest the way the swap from v6 needs to be judged:
#
#  1. COMMANDED OBJECTS: apply the approved object-frame rule
#     (ApogeeReduction.object_frame_rule) to the out-of-fold scores and the rule
#     features; write every acted-on exposure with the test that fired, counts by
#     reason, and the released-row impact (if the manifest has rel_nrows /
#     rel_ngood). Optionally diff against a reference acted list (--reference_acted).
#  2. FINAL MODEL vs OUT-OF-FOLD (--model + --features): score every commanded
#     object with the single deployed model and compare the rule outcome with the
#     out-of-fold outcome. Rows the final model trained on, and rows on nights
#     that contributed training rows, are scored in-sample by the final model.
#  3. CALIBRATION MASKING, v6 vs 2D: `exposure_predicted_bad` under the deployed
#     v6 verdicts (manifest v6_pred / v6_status) versus under the 2D forest
#     (translated v6 policy: exposure_check_category + exposure_predicted_bad),
#     per commanded class, both directions, with the list of changes.
#
# Inputs: --oof (train_classifier_2d.jl --out_oof; columns p_<class>, oh_rel,
#   n_G_line_z, n_G_tr_dyn, template_night, in_train), --manifest (label manifest
#   with lamp flags and v6_pred / v6_status; rel_nrows / rel_ngood optional).
# Run: JULIA_NUM_THREADS=n julia --project=<repo> scripts/classifier/evaluate_classifier_2d.jl \
#        --oof oof.csv --manifest M.csv --outdir eval/ [--model m.jld2 --features F.csv]
#        [--reference_acted acted.csv]
using ArgParse, CSV, DataFrames, DecisionTree, JLD2, Statistics, Printf
using ApogeeReduction
using ApogeeReduction: object_frame_rule, exposure_class_label, exposure_check_category,
                       exposure_predicted_bad, CLASSIFIER2D_CLASSES, CLASSIFIER_UNKNOWN_TAU,
                       CLASSIFIER_FLAG_TAU

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--oof"
        required = true
        arg_type = String
        "--manifest"
        required = true
        arg_type = String
        "--outdir"
        required = true
        arg_type = String
        "--prob_prefix"
        default = "p_"
        help = "column prefix of the class probabilities in --oof (prototype tables: p_full_)"
        arg_type = String
        "--model"
        default = ""
        arg_type = String
        "--features"
        default = ""
        help = "feature table (assemble_features_2d.jl) for --model scoring"
        arg_type = String
        "--reference_acted"
        default = ""
        help = "CSV (tele, mjd, expnum) of a reference acted-on list to diff against"
        arg_type = String
    end
    return parse_args(s)
end
parg = parse_commandline()
mkpath(parg["outdir"])
pp = parg["prob_prefix"]

O = CSV.read(parg["oof"], DataFrame; types = Dict(:tele => String))
M = CSV.read(parg["manifest"], DataFrame; types = Dict(:tele => String, :v6_pred => String,
    :v6_status => String, :label => String, :role => String))
mcols = intersect(["tele", "mjd", "expnum", "lamp_quartz", "lamp_thar", "lamp_une", "v6_pred",
        "v6_status", "rel_nrows", "rel_ngood"], names(M))
O = leftjoin(select(O, Not(intersect(names(O), setdiff(mcols, ["tele", "mjd", "expnum"])))),
    M[:, mcols], on = [:tele, :mjd, :expnum], order = :left)
for c in ("rel_nrows", "rel_ngood")
    c in names(O) ? (O[!, c] = coalesce.(O[!, c], 0)) : (O[!, c] = zeros(Int, nrow(O)))
end
O.labeled = [exposure_class_label(r.image_type, r.lamp_quartz, r.lamp_thar, r.lamp_une)
             for r in eachrow(O)]
O.p_sky_oof = Float64.(O[!, pp * "sky"])

function apply_rule!(T, pcol, suffix)
    rr = [object_frame_rule(1 - T[i, pcol], T.oh_rel[i], T.n_G_line_z[i], T.n_G_tr_dyn[i])
          for i in 1:nrow(T)]
    T[!, "act" * suffix] = getindex.(rr, :act)
    T[!, "reason" * suffix] = getindex.(rr, :reason)
    T
end

# ---------------- 1. commanded objects, out-of-fold ----------------
obj = O[O.image_type .== "object", :]
apply_rule!(obj, :p_sky_oof, "")
A = obj[obj.act, :]
sel = [:tele, :mjd, :expnum, :p_sky_oof, :oh_rel, :n_G_line_z, :n_G_tr_dyn, :reason,
    :template_night, :rel_nrows, :rel_ngood]
CSV.write(joinpath(parg["outdir"], "object_rule_acted_oof.csv"), sort(A[:, sel], [:tele, :mjd, :expnum]))
R = combine(groupby(obj, :reason), nrow => :exposures)
CSV.write(joinpath(parg["outdir"], "object_rule_reasons_oof.csv"), R)
@printf("commanded objects with features: %d; acted on (out-of-fold): %d\n", nrow(obj), nrow(A))
@printf("  in 0.4.0 (rel_nrows > 0): %d exposures, %d released rows, %d rows passing cuts\n",
    sum(A.rel_nrows .> 0), sum(A.rel_nrows), sum(A.rel_ngood))
show(R, allrows = true)
println()
if !isempty(parg["reference_acted"])
    Ref_ = CSV.read(parg["reference_acted"], DataFrame; select = [:tele, :mjd, :expnum],
        types = Dict(:tele => String))
    a = Set(zip(A.tele, A.mjd, A.expnum))
    b = Set(zip(Ref_.tele, Ref_.mjd, Ref_.expnum))
    @printf("reference acted list: %d; identical: %s (only here %d, only in reference %d)\n",
        length(b), a == b, length(setdiff(a, b)), length(setdiff(b, a)))
end

# ---------------- 2. final model vs out-of-fold ----------------
if !isempty(parg["model"])
    md = JLD2.load(parg["model"])
    cols = md["features"]
    F = CSV.read(parg["features"], DataFrame; types = Dict(:tele => String),
        select = (i, nm) -> String(nm) in vcat(["tele", "mjd", "expnum"], cols))
    if "n_read" in cols && !("n_read" in names(F))
        F = leftjoin(F, M[:, [:tele, :mjd, :expnum, :n_read]], on = [:tele, :mjd, :expnum])
    end
    F = innerjoin(F, obj[:, [:tele, :mjd, :expnum]], on = [:tele, :mjd, :expnum])
    X = Matrix{Float64}(undef, nrow(F), length(cols))
    for (j, c) in enumerate(cols)
        v = F[!, c]
        for i in 1:nrow(F)
            x = v[i]
            X[i, j] = (x === missing || !isfinite(x)) ? -9999.0 : Float64(x)
        end
    end
    P = apply_forest_proba(md["model"], X, md["classes"])
    F.p_sky_final = P[:, findfirst(==("sky"), md["classes"])]
    obj = leftjoin(obj, F[:, [:tele, :mjd, :expnum, :p_sky_final]], on = [:tele, :mjd, :expnum],
        order = :left)
    apply_rule!(obj, :p_sky_final, "_final")
    trn = Set((r.tele, r.mjd) for r in eachrow(O) if r.in_train)
    obj.night_has_training_rows = [(r.tele, r.mjd) in trn for r in eachrow(obj)]
    ch = obj[obj.act .!= obj.act_final, :]
    ch.direction = ifelse.(ch.act_final, "added_by_final", "removed_by_final")
    CSV.write(joinpath(parg["outdir"], "object_rule_final_vs_oof_changes.csv"),
        sort(ch[:, vcat(sel, [:p_sky_final, :reason_final, :in_train, :night_has_training_rows,
                :direction])], [:direction, :tele, :mjd, :expnum]))
    @printf("final model: acted %d (out-of-fold %d); added %d, removed %d\n", sum(obj.act_final),
        sum(obj.act), sum(ch.direction .== "added_by_final"), sum(ch.direction .== "removed_by_final"))
    @printf("  of the changes: %d rows were final-model training rows, %d on nights with training rows, %d template-night rows\n",
        sum(ch.in_train), sum(ch.night_has_training_rows), sum(ch.template_night))
    @printf("  |P(sky) final - out-of-fold|: median %.4f, 99th pct %.4f; on training-row nights %.4f vs other nights %.4f (median)\n",
        median(abs.(obj.p_sky_final .- obj.p_sky_oof)),
        quantile(abs.(obj.p_sky_final .- obj.p_sky_oof), 0.99),
        median(abs.(obj.p_sky_final[obj.night_has_training_rows] .- obj.p_sky_oof[obj.night_has_training_rows])),
        median(abs.(obj.p_sky_final[.!obj.night_has_training_rows] .- obj.p_sky_oof[.!obj.night_has_training_rows])))
end

# ---------------- 3. calibration masking: v6 verdict vs 2D forest ----------------
cal = O[(O.image_type .!= "object") .& .!ismissing.(O.v6_status), :]
P2 = Matrix{Float64}(cal[:, [pp * c for c in CLASSIFIER2D_CLASSES]])
cal.pred_2d = [CLASSIFIER2D_CLASSES[argmax(P2[i, :])] for i in 1:nrow(cal)]
cal.prob_2d = vec(maximum(P2, dims = 2))
cal.status_2d = [exposure_check_category(cal.labeled[i],
                     (pred = cal.pred_2d[i], prob = cal.prob_2d[i],
                         decision = cal.prob_2d[i] < CLASSIFIER_UNKNOWN_TAU ? "unknown" : "ok"),
                     CLASSIFIER_FLAG_TAU) for i in 1:nrow(cal)]
cal.bad_v6 = [exposure_predicted_bad(r.labeled, r.v6_pred, r.v6_status) for r in eachrow(cal)]
cal.bad_2d = [exposure_predicted_bad(r.labeled, r.pred_2d, r.status_2d) for r in eachrow(cal)]
S = combine(groupby(cal, :labeled), nrow => :exposures, :bad_v6 => sum => :bad_v6,
    :bad_2d => sum => :bad_2d,
    [:bad_v6, :bad_2d] => ((a, b) -> sum(.!a .& b)) => :ok_to_bad,
    [:bad_v6, :bad_2d] => ((a, b) -> sum(a .& .!b)) => :bad_to_ok)
sort!(S, :exposures, rev = true)
CSV.write(joinpath(parg["outdir"], "calibration_masking_v6_vs_2d_summary.csv"), S)
C = cal[cal.bad_v6 .!= cal.bad_2d, [:tele, :mjd, :expnum, :image_type, :labeled, :v6_pred,
    :v6_status, :bad_v6, :pred_2d, :prob_2d, :status_2d, :bad_2d, :template_night]]
C.direction = ifelse.(C.bad_2d, "ok_to_bad", "bad_to_ok")
CSV.write(joinpath(parg["outdir"], "calibration_masking_v6_vs_2d_changes.csv"),
    sort(C, [:labeled, :direction, :tele, :mjd, :expnum]))
println("\ncalibration masking, v6 verdict vs 2D forest (translated v6 policy):")
show(S, allrows = true)
println("\nwrote tables to ", parg["outdir"])
