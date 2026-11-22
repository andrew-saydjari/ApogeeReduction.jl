# Stage 4 of the 2D content-forest rebuild: train the forest, score every
# exposure out-of-fold by night, and write the deployable artifact.
#
# Recipe (identical to the prototype that produced the approved scores,
# 2026_11_19/classifier_2d s20_train_eval.jl, feature set "full"):
#   * features: classifier2d_feature_names() (175; order matters, see there),
#     every non-finite value -> -9999;
#   * classes: CLASSIFIER2D_CLASSES (sky, twilight, dark, dark_persist, ledflat,
#     quartzflat, domeflat, thar, une, fpi);
#   * training rows: manifest role "train" (strong) or "train_weak" with a label in
#     the classes; per class at most --cap (4000) rows, strong first, then weak,
#     each drawn with MersenneTwister(--seed_rows = 20261119) class by class in
#     class order, over the exposures sorted by (tele, mjd, expnum);
#   * folds: the unique nights (tele, mjd) in that sorted order, shuffled with
#     MersenneTwister(--seed_folds = 5), assigned round-robin to --nfold (5) folds;
#     template nights are always train-only and are never scored out-of-fold;
#   * forest: DecisionTree.build_forest(y, X, -1, 300, 0.7, -1; rng = 1)
#     (sqrt(n) features per split, 300 trees, 70 % bagging);
#   * out-of-fold probabilities for EVERY exposure with features (labeled or not)
#     from the fold model that did not train on its night; template-night rows
#     get the final model's (in-sample) probabilities, marked template_night.
#   * the final (deployed) model is trained on all training rows.
#
# Inputs:
#   --features   training feature table (assemble_features_2d.jl; or any CSV with
#                tele, mjd, expnum, image_type and the feature columns). n_read is
#                taken from the manifest when the table has no n_read column.
#   --manifest   label manifest CSV: tele, mjd, expnum, image_type, n_read, label, role
#                (role: train | train_weak | eval_only | unlabeled)
#   --template_nights  CSV tele, mjd (build_line_templates_2d.jl)
# Outputs: --out_model (JLD2: model, classes, features, unknown_tau, flag_tau,
#   provenance), --out_oof (CSV), --out_confusion (CSV, out-of-fold, non-template
#   training rows). The artifact is a calibration input: it is NOT committed;
#   record its sha256 in metadata/exposure_classifier_2d.sha256 and install it in
#   cal_ref/exposure_classifier/.
#
# Run: JULIA_NUM_THREADS=n julia --project=<repo> scripts/classifier/train_classifier_2d.jl \
#        --features F.csv --manifest M.csv --template_nights N.csv \
#        --out_model exposure_classifier_2d_v1.jld2 --out_oof oof.csv --out_confusion conf.csv
using ArgParse, CSV, DataFrames, DecisionTree, JLD2, Random, Statistics, Printf, Dates, SHA
using ApogeeReduction
using ApogeeReduction: classifier2d_feature_names, CLASSIFIER2D_CLASSES,
                       CLASSIFIER_UNKNOWN_TAU, CLASSIFIER_FLAG_TAU

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--features"
        required = true
        arg_type = String
        "--manifest"
        required = true
        arg_type = String
        "--template_nights"
        required = true
        arg_type = String
        "--out_model"
        required = true
        arg_type = String
        "--out_oof"
        default = ""
        arg_type = String
        "--out_confusion"
        default = ""
        arg_type = String
        "--final_only"
        action = :store_true
        help = "skip the cross-validation; train and write the final model only"
        "--cap"
        default = 4000
        arg_type = Int
        "--nfold"
        default = 5
        arg_type = Int
        "--ntrees"
        default = 300
        arg_type = Int
        "--seed_rows"
        default = 20261119
        arg_type = Int
        "--seed_folds"
        default = 5
        arg_type = Int
    end
    return parse_args(s)
end
parg = parse_commandline()
const CLASSES = CLASSIFIER2D_CLASSES
cols = classifier2d_feature_names()

t0 = time()
F = CSV.read(parg["features"], DataFrame; select = (i, nm) -> String(nm) in
                                                         vcat(["tele", "mjd", "expnum", "image_type"], cols),
    types = Dict(:tele => String))
L = CSV.read(parg["manifest"], DataFrame; select = [:tele, :mjd, :expnum, :n_read, :label, :role],
    types = Dict(:tele => String, :label => String, :role => String))
L.label = coalesce.(L.label, "")
L.role = coalesce.(L.role, "")
"n_read" in names(F) && select!(L, Not(:n_read))
D = innerjoin(F, L, on = [:tele, :mjd, :expnum])
sort!(D, [:tele, :mjd, :expnum])
missing_cols = setdiff(cols, names(D))
isempty(missing_cols) || error("feature table lacks: $(missing_cols)")
@printf("rows with features and a manifest entry: %d  (read in %.0f s)\n", nrow(D), time() - t0)
tn = CSV.read(parg["template_nights"], DataFrame; types = Dict(:tele => String))
tset = Set(zip(tn.tele, tn.mjd))
D.template_night = [(r.tele, r.mjd) in tset for r in eachrow(D)]

function mat(idx)
    X = Matrix{Float64}(undef, length(idx), length(cols))
    for (j, c) in enumerate(cols)
        v = D[!, c]
        for (k, i) in enumerate(idx)
            x = v[i]
            X[k, j] = (x === missing || !isfinite(x)) ? -9999.0 : Float64(x)
        end
    end
    X
end

rng = MersenneTwister(parg["seed_rows"])
isTrain = [(r.role == "train" || r.role == "train_weak") && r.label in CLASSES for r in eachrow(D)]
tr = Int[]
for c in CLASSES
    strong = findall(i -> isTrain[i] && D.label[i] == c && D.role[i] == "train", 1:nrow(D))
    weak = findall(i -> isTrain[i] && D.label[i] == c && D.role[i] == "train_weak", 1:nrow(D))
    s = vcat(strong[randperm(rng, length(strong))[1:min(parg["cap"], length(strong))]],
        weak[randperm(rng, length(weak))[1:min(max(0, parg["cap"] - length(strong)), length(weak))]])
    append!(tr, s)
    println(rpad(c, 13), " strong ", length(strong), " weak ", length(weak), " used ", length(s))
end
nightkey = [(r.tele, r.mjd) for r in eachrow(D)]
un = unique(nightkey)
perm = randperm(MersenneTwister(parg["seed_folds"]), length(un))
fold_of = Dict(un[perm[i]] => (i - 1) % parg["nfold"] + 1 for i in eachindex(un))
D.fold = [fold_of[k] for k in nightkey]

build(X, y) = build_forest(y, X, -1, parg["ntrees"], 0.7, -1; rng = 1)
P = fill(NaN, nrow(D), length(CLASSES))
if !parg["final_only"]
    for k in 1:parg["nfold"]
        t1 = time()
        trk = filter(i -> D.fold[i] != k || D.template_night[i], tr)
        model = build(mat(trk), String.(D.label[trk]))
        te = findall(i -> D.fold[i] == k && !D.template_night[i], 1:nrow(D))
        P[te, :] = apply_forest_proba(model, mat(te), CLASSES)
        @printf("  fold %d trained on %d, predicted %d  (%.0f s)\n", k, length(trk), length(te), time() - t1)
        flush(stdout)
    end
end
t1 = time()
final = build(mat(tr), String.(D.label[tr]))
@printf("final model trained on %d rows (%.0f s)\n", length(tr), time() - t1)
tmp = findall(D.template_night)
P[tmp, :] = apply_forest_proba(final, mat(tmp), CLASSES)

function git_info()
    try
        d = pkgdir(ApogeeReduction)
        c = strip(read(`git -C $d rev-parse HEAD`, String))
        b = strip(read(`git -C $d rev-parse --abbrev-ref HEAD`, String))
        clean = isempty(strip(read(`git -C $d status --porcelain`, String)))
        (c, b, clean)
    catch
        ("unknown", "unknown", false)
    end
end
gc, gb, gclean = git_info()
JLD2.jldsave(parg["out_model"]; model = final, classes = CLASSES, features = cols,
    unknown_tau = CLASSIFIER_UNKNOWN_TAU, flag_tau = CLASSIFIER_FLAG_TAU,
    git_commit = gc, git_branch = gb, git_clean = gclean,
    trained_utc = string(Dates.now(Dates.UTC)),
    manifest_sha256 = bytes2hex(open(sha256, parg["manifest"])),
    template_nights_sha256 = bytes2hex(open(sha256, parg["template_nights"])),
    n_train = length(tr))
println("wrote ", parg["out_model"], "  sha256 ", bytes2hex(open(sha256, parg["out_model"])))

if !parg["final_only"]
    ev = filter(i -> !D.template_night[i], tr)
    pred = [CLASSES[argmax(P[i, :])] for i in ev]
    @printf("CV accuracy %.4f\n", mean(pred .== D.label[ev]))
    cm = zeros(Int, length(CLASSES), length(CLASSES))
    for (j, i) in enumerate(ev)
        cm[findfirst(==(D.label[i]), CLASSES), findfirst(==(pred[j]), CLASSES)] += 1
    end
    CM = DataFrame(cm, CLASSES)
    insertcols!(CM, 1, :label => CLASSES)
    isempty(parg["out_confusion"]) || CSV.write(parg["out_confusion"], CM)
    show(CM, allrows = true)
    println()
end
if !isempty(parg["out_oof"])
    O = D[:, [:tele, :mjd, :expnum, :image_type, :label, :role, :fold, :template_night]]
    O.in_train = falses(nrow(D))
    O.in_train[tr] .= true
    for (j, c) in enumerate(CLASSES)
        O[!, "p_$(c)"] = P[:, j]
    end
    O.pred = [all(isnan, P[i, :]) ? "" : CLASSES[argmax(P[i, :])] for i in 1:nrow(D)]
    O.pmax = [all(isnan, P[i, :]) ? NaN : maximum(P[i, :]) for i in 1:nrow(D)]
    for c in ("oh_rel", "n_G_line_z", "n_G_tr_dyn", "cont_rel")
        O[!, c] = D[!, c]
    end
    O.scored_by = ifelse.(O.template_night, "final_in_sample", parg["final_only"] ? "none" : "out_of_fold")
    CSV.write(parg["out_oof"], O)
    println("wrote ", parg["out_oof"])
end
