# Rebuild the 2D content forest's LINE TEMPLATES from scratch
# (data/exposure_classifier_2d/line_templates.h5).
#
# A template is, per (kind, telescope, chip) with kind in oh/thar/une/fpi, the
# median over up to 300 labeled exposures of the normalized line excess
# ex / robsig(ex) of the chip's spectral profile
# (ApogeeReduction.classifier2d_line_excess, the same function the t_* features
# use), with 4 edge pixels zeroed on each side. Labels: sky -> oh, thar, une, fpi.
#
# Template exposures come ONLY from TEMPLATE NIGHTS, which the training script
# then treats as train-only (never scored out-of-fold), so the t_* features carry
# no evaluation leakage. Template nights are selected (seeded, class-aware):
#   * nights = unique (tele, mjd) of the exposures with features, in canonical
#     order (sorted by mjd, tele, expnum), shuffled with MersenneTwister(--seed=11);
#   * protected nights are never template nights: any night holding an eval_only
#     row or a commanded object whose label is not sky;
#   * walk the shuffled nights and take a night only if it adds training-labeled
#     (role train / train_weak) exposures of a template class (sky/thar/une/fpi,
#     per telescope) that still has fewer than 40; stop at 10 % of the eligible
#     nights or when every class has 40.
# Template exposures per (kind, telescope): the first 300 labeled train/train_weak
# exposures on template nights, in canonical order; a kind needs >= 5.
#
# Two-step use (avoids extracting the whole archive just for templates):
#   1. selection only:  --manifest M --universe U --out_template_nights N.csv
#                       --out_template_exposures X.csv
#      (U = CSV of tele, mjd, expnum of every exposure with features)
#   2. build:           --manifest M --template_nights N.csv --features F.h5 --outfile T.h5
#      (F.h5 from extract_features_2d.jl; must contain the template exposures)
# Or in one go with --universe omitted: the features file's ok rows are the universe.
using ArgParse, CSV, DataFrames, HDF5, Random, Statistics
using ApogeeReduction
using ApogeeReduction: classifier2d_line_excess, classifier2d_robsig

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--manifest"
        required = true
        help = "label manifest CSV: tele, mjd, expnum, image_type, label, role"
        arg_type = String
        "--universe"
        default = ""
        help = "CSV of tele, mjd, expnum of all exposures with features (default: ok rows of --features)"
        arg_type = String
        "--features"
        default = ""
        help = "extract_features_2d.jl output (PS profiles)"
        arg_type = String
        "--template_nights"
        default = ""
        help = "use these template nights (tele, mjd) instead of selecting them"
        arg_type = String
        "--out_template_nights"
        default = ""
        arg_type = String
        "--out_template_exposures"
        default = ""
        arg_type = String
        "--outfile"
        default = ""
        help = "templates HDF5"
        arg_type = String
        "--seed"
        default = 11
        arg_type = Int
        "--per_class"
        default = 40
        arg_type = Int
        "--max_per_template"
        default = 300
        arg_type = Int
    end
    return parse_args(s)
end
parg = parse_commandline()

const KINDS = (("oh", "sky"), ("thar", "thar"), ("une", "une"), ("fpi", "fpi"))
const TRAINROLES = ("train", "train_weak")

L = CSV.read(parg["manifest"], DataFrame; select = [:tele, :mjd, :expnum, :image_type, :label, :role],
    types = Dict(:tele => String, :image_type => String, :label => String, :role => String))
L.label = coalesce.(L.label, "")
L.role = coalesce.(L.role, "")

FH = nothing
U = if !isempty(parg["universe"])
    CSV.read(parg["universe"], DataFrame; select = [:tele, :mjd, :expnum], types = Dict(:tele => String))
else
    isempty(parg["features"]) && error("give --universe or --features")
    h5open(parg["features"], "r") do f
        ok = read(f["ok"]) .== 1
        DataFrame(tele = read(f["tele"])[ok], mjd = read(f["mjd"])[ok], expnum = read(f["expnum"])[ok])
    end
end
sort!(U, [:mjd, :tele, :expnum])                       # canonical order
F = leftjoin(U, L[:, [:tele, :mjd, :expnum, :label, :role]], on = [:tele, :mjd, :expnum],
    order = :left)
F.label = coalesce.(F.label, "")
F.role = coalesce.(F.role, "")
println("exposures with features: ", nrow(F))

tn = if !isempty(parg["template_nights"])
    CSV.read(parg["template_nights"], DataFrame; types = Dict(:tele => String))
else
    nights = unique(F[:, [:tele, :mjd]])
    rng = MersenneTwister(parg["seed"])
    nights = nights[randperm(rng, nrow(nights)), :]
    need = Dict((c, t) => parg["per_class"] for c in ("sky", "thar", "une", "fpi")
    for t in ("apo", "lco"))
    byn = Dict{Tuple{String, Int}, Vector{Int}}()
    for i in 1:nrow(F)
        F.role[i] in TRAINROLES || continue
        push!(get!(byn, (F.tele[i], F.mjd[i]), Int[]), i)
    end
    prot = Set((r.tele, r.mjd) for r in eachrow(L)
    if r.role == "eval_only" || (r.image_type == "object" && r.label != "sky"))
    nights = nights[[!((r.tele, r.mjd) in prot) for r in eachrow(nights)], :]
    println("protected nights excluded from template selection: ", length(prot))
    t = DataFrame(tele = String[], mjd = Int[])
    for r in eachrow(nights)
        nrow(t) >= 0.1 * nrow(nights) && break
        ii = get(byn, (r.tele, r.mjd), Int[])
        any(i -> get(need, (F.label[i], F.tele[i]), 0) > 0, ii) || continue
        push!(t, (r.tele, r.mjd))
        for i in ii
            k = (F.label[i], F.tele[i])
            haskey(need, k) && (need[k] -= 1)
        end
        all(<=(0), values(need)) && break
    end
    println("template nights: ", nrow(t), " of ", nrow(nights), "; remaining need: ", need)
    t
end
isempty(parg["out_template_nights"]) || CSV.write(parg["out_template_nights"], tn)
tset = Set(zip(tn.tele, tn.mjd))

sel = Dict{Tuple{String, String}, Vector{Int}}()
for (tmpl, lab) in KINDS, tele in ("apo", "lco")
    s = findall(i -> (F.tele[i], F.mjd[i]) in tset && F.label[i] == lab &&
                         F.role[i] in TRAINROLES && F.tele[i] == tele, 1:nrow(F))
    sel[(tmpl, tele)] = s[1:min(parg["max_per_template"], length(s))]
    println("$tmpl $tele: $(length(sel[(tmpl, tele)])) template exposures")
end
if !isempty(parg["out_template_exposures"])
    X = vcat([DataFrame(kind = tmpl, tele = F.tele[s], mjd = F.mjd[s], expnum = F.expnum[s])
              for ((tmpl, tele), s) in sel]...)
    X = leftjoin(X, L[:, [:tele, :mjd, :expnum, :image_type]], on = [:tele, :mjd, :expnum])
    CSV.write(parg["out_template_exposures"], sort(X, [:kind, :tele, :mjd, :expnum]))
end
isempty(parg["outfile"]) && (println("selection only (no --outfile)"); exit(0))
isempty(parg["features"]) && error("--outfile needs --features (PS profiles)")

# PS of the selected exposures, looked up by key in the features file
PSd = h5open(parg["features"], "r") do f
    keys_ = collect(zip(read(f["tele"]), read(f["mjd"]), read(f["expnum"])))
    ok = read(f["ok"]) .== 1
    idx = Dict(k => j for (j, k) in enumerate(keys_) if ok[j])
    need_rows = unique(vcat(values(sel)...))
    PS = f["PS"]
    Dict((F.tele[i], F.mjd[i], F.expnum[i]) => PS[:, :, idx[(F.tele[i], F.mjd[i], F.expnum[i])]]
    for i in need_rows)
end
out = Dict{String, Vector{Float64}}()
for (tmpl, lab) in KINDS, tele in ("apo", "lco")
    s = sel[(tmpl, tele)]
    length(s) < 5 && continue
    for (ci, chip) in enumerate(CHIP_LIST)
        E = hcat([(ex = classifier2d_line_excess(PSd[(F.tele[i], F.mjd[i], F.expnum[i])][:, ci])[1];
                   ex ./ classifier2d_robsig(ex)) for i in s]...)
        t = vec(median(E, dims = 2))
        t[1:4] .= 0
        t[(end - 3):end] .= 0
        out["$(tmpl)_$(tele)_$(chip)"] = t
    end
end
h5open(parg["outfile"], "w") do f
    for (k, v) in out
        f[k] = v
    end
end
println("wrote ", length(out), " templates to ", parg["outfile"])
