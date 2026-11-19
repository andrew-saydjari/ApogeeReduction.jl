# Build the exposure-type classifier training/eval manifest from a decorated
# almanac plus the truth sources surfaced by the DR21 pass-1 sweep:
#   - stratified almanac sample: ALL exposures from rare (obs, labeled-class)
#     strata, up to --nper random exposures from common strata (seeded rng)
#   - twilight census (sun-altitude ground truth; see the 2026_09_11
#     twilight_flats census): ALL in-range rows. Night-sky flats and true
#     twilights get census-truth content labels in training no matter what
#     image_type says — this is what fixes the object/dark-typed sky flats.
#   - metadata/special_cal_obs.txt sequences: ALL in-range rows
#   - engineering-carton exposures, FPS-era config_id=-1 object frames,
#     census daytime dome tests, and the lco 57800 stale frames: ALL included
#     but role=eval_only (commanded labels unreliable / adversarial; never
#     trained on)
#
# Output CSV columns: obs, mjd, expnum, image_type, labeled (image_type +
# lamp-flag class string), n_read, config_id, design_id, fps, engineering,
# truth_source, census_truth, census_signal, role (train_pool | eval_only),
# and the deployed-model columns v6_pred/v6_prob/v6_status carried from the
# almanac decoration for later old-vs-new comparison.
#
# The v7 manifest committed at metadata/exposure_classifier_training_manifest_v7.csv
# was built from /mnt/ceph/users/sdssv/work/asaydjari/pass1_dr21/almanac/
# allobs_57618_61230.h5 with the defaults below (see the 2026_09_13
# classifier_retrain report for the full provenance chain).
using Pkg;
Pkg.instantiate();
using HDF5, ArgParse, DataFrames, CSV, Random

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--almanac_file"
        required = true
        help = "decorated almanac (raw/ + exposure_class/ groups)"
        arg_type = String
        "--twilight_census"
        required = true
        help = "twilight census CSV with is_twilight/is_nightsky/is_daytime + signal"
        arg_type = String
        "--special_cal_obs"
        default = joinpath(dirname(dirname(@__DIR__)), "metadata", "special_cal_obs.txt")
        help = "special calibration sequence list"
        arg_type = String
        "--outfile"
        required = true
        help = "output manifest CSV"
        arg_type = String
        "--nper"
        default = 1000
        help = "max exposures sampled per common (obs, labeled) stratum"
        arg_type = Int
        "--seed"
        default = 2026
        help = "sampling rng seed"
        arg_type = Int
        "--mjd_min"
        default = 57618
        help = "restrict truth-set joins to mjd >= this"
        arg_type = Int
        "--dark_topup"
        default = 5000
        help = "extra commanded-dark exposures WITH A BRIGHT PREDECESSOR to add to the training pool beyond the stratified sample, split evenly apo/lco (0 disables). dark_persist labels are constructed, not commanded (labeled dark AND previous exposure bright AND flux above the 0.99 quantile of clean darks) and only ~2% of prev-bright darks qualify, so the plain --nper stratification starves the class: 978 prev-bright darks yielded 18 labels where v6 had 56, halving held-out recall. Since dark_persist is what makes exposure_predicted_bad keep persistence-contaminated darks out of dark cal stacks, that shortfall is a calibration problem, not a cosmetic one. Prev-bright darks are identifiable from the almanac alone, so this top-up costs no extra metadata, only feature extraction."
        arg_type = Int
    end
    return parse_args(s)
end
parg = parse_commandline()

# ---- flatten the almanac: raw metadata + exposure_class decoration ----
rows = DataFrame[]
h5open(parg["almanac_file"], "r") do f
    for tele in keys(f["raw"])
        for mjd in keys(f["raw/$tele"])
            haskey(f["raw/$tele/$mjd"], "exposures") || continue
            g = f["raw/$tele/$mjd/exposures"]
            n = length(read(g["exposure"]))
            d = DataFrame(obs = fill(tele, n), mjd = fill(parse(Int, mjd), n),
                expnum = read(g["exposure"]), image_type = read(g["image_type"]),
                lamp_quartz = string.(read(g["lamp_quartz"])),
                lamp_thar = string.(read(g["lamp_thar"])),
                lamp_une = string.(read(g["lamp_une"])),
                config_id = read(g["config_id"]), design_id = read(g["design_id"]),
                fps = read(g["fps"]), n_read = read(g["n_read"]))
            ec = f["exposure_class/$tele/$mjd"]
            @assert read(ec["exposure"]) == d.expnum "exposure_class misaligned at $tele $mjd"
            d.v6_pred = read(ec["exposure_class_pred"])
            d.v6_prob = read(ec["exposure_class_prob"])
            d.v6_status = read(ec["exposure_class_status"])
            d.engineering = read(ec["engineering"])
            push!(rows, d)
        end
    end
end
df = sort(vcat(rows...), [:obs, :mjd, :expnum])
println("almanac exposures: ", nrow(df))

lampstr(x) = (x in ("true", "True", "1")) ? "1" : ((x in ("false", "False", "0")) ? "0" : "?")
df.labeled = lowercase.(df.image_type) .* "_q" .* lampstr.(df.lamp_quartz) .*
             "t" .* lampstr.(df.lamp_thar) .* "u" .* lampstr.(df.lamp_une)
df.key = tuple.(df.obs, df.mjd, df.expnum)

# ---- truth-set memberships ----
tw = CSV.read(parg["twilight_census"], DataFrame)
tw = tw[tw.mjd .>= parg["mjd_min"], :]
tw.truth = ifelse.(tw.is_twilight, "twilight",
    ifelse.(tw.is_nightsky, "nightsky", ifelse.(tw.is_daytime, "daytime", "?")))
twmap = Dict(tuple(r.obs, r.mjd, r.exposure) => (r.truth, r.signal) for r in eachrow(tw))

# expid = (mjd - 55562) * 10000 + expnum  =>  expnum = expid mod 10000
spec = Dict{Tuple{String, Int, Int}, String}()
for line in eachline(parg["special_cal_obs"])
    (startswith(line, "#") || isempty(strip(line))) && continue
    p = split(line)
    r = split(p[3], "-")
    for e in mod(parse(Int, r[1]), 10000):mod(parse(Int, r[2]), 10000)
        spec[(String(p[1]), parse(Int, p[2]), e)] = String(p[4])
    end
end

df.truth_source = fill("almanac_label", nrow(df))
df.census_truth = fill("", nrow(df))
df.census_signal = fill(NaN, nrow(df))
df.role = fill("", nrow(df))
for (i, k) in enumerate(df.key)
    if haskey(twmap, k)
        t, s = twmap[k]
        df.truth_source[i] = "twilight_census_" * t
        df.census_truth[i] = t
        df.census_signal[i] = isfinite(s) ? s : NaN
    elseif haskey(spec, k)
        df.truth_source[i] = "special_cal_" * spec[k]
    end
end
isfake = (df.obs .== "lco") .& (df.mjd .== 57800)
df.truth_source[isfake] .= "lco57800_stale"
iseng = df.engineering .== 1
df.truth_source[iseng .& (df.truth_source .== "almanac_label")] .= "engineering"
iscfgless = (df.fps .== 1) .& (df.config_id .== -1) .& (df.image_type .== "object")
df.truth_source[iscfgless .& (df.truth_source .== "almanac_label")] .= "fps_configless"
evalonly = isfake .| iseng .| iscfgless .| (df.census_truth .== "daytime")
df.role[evalonly] .= "eval_only"

# ---- stratified training sample over the rest ----
rng = MersenneTwister(parg["seed"])
pool = df[.!evalonly, :]
picked = Set{Tuple{String, Int, Int}}()
for g in groupby(pool, [:obs, :labeled])
    n = nrow(g)
    idx = n <= parg["nper"] ? (1:n) : sort(randperm(rng, n)[1:parg["nper"]])
    for k in g.key[idx]
        push!(picked, k)
    end
end
truthin = .!evalonly .& (df.truth_source .!= "almanac_label")
for k in df.key[truthin]
    push!(picked, k)
end
df.role[in.(df.key, Ref(picked)) .& (df.role .== "")] .= "train_pool"

# ---- targeted top-up: prev-bright darks, for dark_persist label support ----
# (see --dark_topup help). Deliberately NOT a --nper increase: raising nper
# across the board would multiply the already-saturated common strata to buy the
# same handful of dark_persist rows.
if parg["dark_topup"] > 0
    PERSIST_SOURCES = ["internalflat", "quartzflat", "domeflat", "arclamp"]
    prevtype = Dict(df.key[j] => df.image_type[j] for j in 1:nrow(df))
    prev_bright(k) = get(prevtype, (k[1], k[2], k[3] - 1), "") in PERSIST_SOURCES
    iscand = (df.labeled .== "dark_q0t0u0") .& (df.role .== "") .& prev_bright.(df.key)
    rng2 = MersenneTwister(parg["seed"] + 1)
    nadd = 0
    for o in ("apo", "lco")
        idx = findall(iscand .& (df.obs .== o))
        n = min(length(idx), parg["dark_topup"] ÷ 2)
        sub = idx[sort(randperm(rng2, length(idx))[1:n])]
        df.role[sub] .= "train_pool"
        df.truth_source[sub] .= "almanac_label_darktopup"
        nadd += n
        println("dark top-up $o: $n of $(length(idx)) available prev-bright darks")
    end
    println("dark top-up total: ", nadd)
end

sel = df[df.role .!= "", :]
println("manifest rows: ", nrow(sel))
show(sort(combine(groupby(sel, [:role, :truth_source]), nrow => :count),
        [:role, :truth_source]), allrows = true)
println()

out = sel[:, [:obs, :mjd, :expnum, :image_type, :labeled, :n_read, :config_id,
    :design_id, :fps, :engineering, :truth_source, :census_truth, :census_signal,
    :role, :v6_pred, :v6_prob, :v6_status]]
sort!(out, [:obs, :mjd, :expnum])
CSV.write(parg["outfile"], out)
println("wrote ", parg["outfile"])
