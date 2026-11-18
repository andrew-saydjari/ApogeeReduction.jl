# Corpus-level FPI fiber-index homogeneity survey (task #54).
#
# `get_fpi_fiberIDs_from_almanac` (ar1D.jl) derives the FPI guide fibers per
# night from the almanac `bonus` rows, so a wrong or moved pair is no longer a
# hardcoded-constant bug -- but a night whose pair silently differs from the
# rest of the corpus is still worth a human look before delivery. This file
# implements the automated check that finds such nights. Design (AKS,
# 2026-09-08):
#
#   "survey the corpus of fpiguide fiber indices per telescope and test
#    whether it is homogeneous. Any night whose FPI fibers differ from that
#    telescope's dominant pair is flagged for review/exclusion. This
#    generalises past the three known nights and needs no constant."
#
# The check FLAGS inhomogeneity for review; it never hard-rejects. The known
# example is the LCO re-fibering epoch MJD 59810-59850 (FPI on fiber_id
# 142/153 instead of the usual 82/213), which AKS + ops ruled trustworthy and
# kept in DR21 -- a known-epochs file lets such reviewed epochs stay visible
# but visually distinct from anything new.
#
# Scope: FPS era only (mjd > get_fps_plate_divide(tele), i.e. apo > 59423,
# lco > 59808). In the plate era the FPI concept does not apply the same way
# (no confSummary, no `bonus` stubs), so plate-era nights are skipped by
# construction rather than reported as missing.

"""
    parse_known_epochs(path) -> Vector{NamedTuple}

Read a file of known/accepted FPI-deviation epochs, one per line:

    tele mjd_lo mjd_hi optional free-text label

Blank lines and `#` comments are ignored. The shipped corpus file is
`metadata/fpi_known_epochs.txt`. Pass `path = ""` for no known epochs.
"""
function parse_known_epochs(path)
    epochs = @NamedTuple{tele::String, mjd_lo::Int, mjd_hi::Int, label::String}[]
    isempty(path) && return epochs
    isfile(path) || error("known-epochs file not found: $(path)")
    for (lineno, line) in enumerate(eachline(path))
        s = strip(first(split(line, "#")))
        isempty(s) && continue
        toks = split(s)
        if length(toks) < 3
            error("known-epochs file $(path) line $(lineno): expected `tele mjd_lo mjd_hi [label]`, got: $(line)")
        end
        push!(epochs,
            (tele = lowercase(toks[1]), mjd_lo = parse(Int, toks[2]),
                mjd_hi = parse(Int, toks[3]), label = join(toks[4:end], " ")))
    end
    return epochs
end

"""
    known_epoch_label(known_epochs, tele, mjd) -> Union{Nothing, String}

The label of the first known epoch covering (`tele`, `mjd`), or `nothing` when
the night is not covered -- i.e. when a deviation there is NEW.
"""
function known_epoch_label(known_epochs, tele, mjd)
    for e in known_epochs
        if (e.tele == tele) && (e.mjd_lo <= mjd <= e.mjd_hi)
            return isempty(e.label) ? "known epoch $(e.mjd_lo)-$(e.mjd_hi)" : e.label
        end
    end
    return nothing
end

"""
    survey_fpi_homogeneity(almanac_file; teles = ["apo", "lco"]) -> Vector{NamedTuple}

Tabulate the FPI guide-fiber pair of every object-bearing FPS-era
configuration in a corpus almanac (raw/ layout), per telescope.

For each telescope this walks every FPS-era night, takes the `config_id` of
every `image_type == "object"` exposure, and reads the sorted `bonus`
fiber_ids of that configuration -- the same derivation
`get_fpi_fiberIDs_from_almanac` uses in production, except that ALL object
configurations are tabulated instead of the night's first. The most common
pair over configurations is the telescope's dominant pair; every
configuration whose pair differs is a deviation.

Returns one NamedTuple per telescope:

  `tele`                telescope name
  `dominant`            the dominant pair (sorted fiber_ids); empty if no data
  `n_configs`           object configurations with bonus rows
  `n_nights`            nights contributing at least one such configuration
  `pair_config_counts`  Dict pair => number of configurations
  `pair_night_mjds`     Dict pair => sorted mjds of nights showing that pair
  `deviations`          Vector of (mjd, config_id, pair) differing from dominant
  `missing_bonus`       Vector of (mjd, config_id) object configs with no bonus rows

Reading is sequential and light (one pass over the file); the almanac is
opened read-only.
"""
function survey_fpi_homogeneity(almanac_file; teles = ["apo", "lco"])
    isfile(almanac_file) || error("almanac file not found: $(almanac_file)")
    results = @NamedTuple{
        tele::String, dominant::Vector{Int}, n_configs::Int, n_nights::Int,
        pair_config_counts::Dict{Vector{Int}, Int},
        pair_night_mjds::Dict{Vector{Int}, Vector{Int}},
        deviations::Vector{@NamedTuple{mjd::Int, config_id::Int, pair::Vector{Int}}},
        missing_bonus::Vector{Tuple{Int, Int}}}[]
    h5open(almanac_file, "r") do f
        for tele in teles
            haskey(f, "raw/$(tele)") || continue
            mjdfps2plate = get_fps_plate_divide(tele)
            records = @NamedTuple{mjd::Int, config_id::Int, pair::Vector{Int}}[]
            missing_bonus = Tuple{Int, Int}[]
            mjd_strs = sort(
                keys(f["raw/$(tele)"]), by = s -> something(tryparse(Int, s), typemax(Int)))
            for mjd_str in mjd_strs
                mjd = tryparse(Int, mjd_str)
                isnothing(mjd) && continue
                # FPI rides the FPS fiber runs; the plate era is out of scope.
                (mjd > mjdfps2plate) || continue
                haskey(f, "raw/$(tele)/$(mjd_str)/exposures") || continue
                df_exp = read_almanac_exp_df(f, tele, mjd_str)
                (("config_id" in names(df_exp)) && ("image_type" in names(df_exp))) || continue
                config_ids = Int[]
                for row in eachrow(df_exp)
                    (row.image_type == "object") || continue
                    config_id = row.config_id isa Integer ? Int(row.config_id) :
                                something(tryparse(Int, string(row.config_id)), 0)
                    (config_id > 0) || continue
                    (config_id in config_ids) || push!(config_ids, config_id)
                end
                for config_id in config_ids
                    fibers_path = "raw/$(tele)/$(mjd_str)/fibers/$(config_id)"
                    if !haskey(f, fibers_path)
                        push!(missing_bonus, (mjd, config_id))
                        continue
                    end
                    df_fib = DataFrame(read(f[fibers_path]))
                    rename!(df_fib, lowercase.(names(df_fib)))
                    pair = if "category" in names(df_fib)
                        sort(unique(Int.(df_fib[df_fib[!, "category"] .== "bonus", "fiber_id"])))
                    else
                        Int[]
                    end
                    if isempty(pair)
                        push!(missing_bonus, (mjd, config_id))
                    else
                        push!(records, (mjd = mjd, config_id = config_id, pair = pair))
                    end
                end
            end
            pair_config_counts = Dict{Vector{Int}, Int}()
            pair_night_sets = Dict{Vector{Int}, Set{Int}}()
            for r in records
                pair_config_counts[r.pair] = get(pair_config_counts, r.pair, 0) + 1
                push!(get!(pair_night_sets, r.pair, Set{Int}()), r.mjd)
            end
            dominant = isempty(pair_config_counts) ? Int[] :
                       argmax(p -> pair_config_counts[p], collect(keys(pair_config_counts)))
            deviations = filter(r -> r.pair != dominant, records)
            push!(results,
                (tele = tele, dominant = dominant, n_configs = length(records),
                    n_nights = length(unique(r.mjd for r in records)),
                    pair_config_counts = pair_config_counts,
                    pair_night_mjds = Dict(p => sort(collect(s))
                    for (p, s) in pair_night_sets),
                    deviations = deviations, missing_bonus = missing_bonus))
        end
    end
    return results
end

_fmt_pair(pair) = isempty(pair) ? "(none)" : join(pair, "/")

"""
    report_fpi_homogeneity(io, results; known_epochs = []) -> NamedTuple

Format the output of [`survey_fpi_homogeneity`](@ref): a per-telescope summary
of every pair seen, then a FLAG section listing each deviating night with the
configurations and pair found there. Deviations inside a known epoch (see
[`parse_known_epochs`](@ref)) are annotated `[KNOWN]` with the epoch's label;
anything else is `[NEW]` and requires review before delivery.

Returns `(n_flagged_nights, n_known, n_new)`. Inhomogeneity is a flag, never
an error: callers should exit 0 regardless and route `n_new > 0` to a human.
"""
function report_fpi_homogeneity(io::IO, results;
        known_epochs = @NamedTuple{
            tele::String, mjd_lo::Int, mjd_hi::Int, label::String}[])
    println(io, "="^72)
    println(io, "FPI fiber-index homogeneity survey (task #54)")
    println(io, "="^72)
    n_known_total = 0
    n_new_total = 0
    n_flagged_total = 0
    for res in results
        mjdfps2plate = get_fps_plate_divide(res.tele)
        println(io, "")
        println(io, "--- $(res.tele) (FPS era, mjd > $(mjdfps2plate)) ---")
        if res.n_configs == 0
            println(io, "no object configurations with bonus rows found.")
            continue
        end
        println(io,
            "object configurations with FPI bonus rows: $(res.n_configs) on $(res.n_nights) nights")
        for pair in sort(collect(keys(res.pair_config_counts)),
            by = p -> -res.pair_config_counts[p])
            mjds = res.pair_night_mjds[pair]
            frac = res.pair_config_counts[pair] / res.n_configs
            marker = (pair == res.dominant) ? "  <- dominant" : ""
            println(io,
                "  pair $(_fmt_pair(pair)): $(res.pair_config_counts[pair]) configs " *
                "($(round(100 * frac, digits = 2))%) on $(length(mjds)) nights, " *
                "mjd $(first(mjds))..$(last(mjds))$(marker)")
        end
        if !isempty(res.missing_bonus)
            mjds = sort(unique(first.(res.missing_bonus)))
            println(io,
                "  note: $(length(res.missing_bonus)) object configs on " *
                "$(length(mjds)) nights carry no bonus rows " *
                "(mjd $(first(mjds))..$(last(mjds))); not counted above.")
        end
    end
    println(io, "")
    println(io, "="^72)
    if all(isempty(res.deviations) for res in results)
        println(io, "RESULT: HOMOGENEOUS -- every configuration matches its telescope's")
        println(io, "dominant FPI pair. Nothing to flag.")
    else
        println(io, "FLAG: nights whose FPI pair differs from the telescope's dominant pair")
        println(io, "(flagged for review/exclusion -- NOT auto-rejected)")
        for res in results
            isempty(res.deviations) && continue
            by_night = Dict{Int, Vector{@NamedTuple{mjd::Int, config_id::Int, pair::Vector{Int}}}}()
            for d in res.deviations
                push!(get!(by_night, d.mjd, []), d)
            end
            for mjd in sort(collect(keys(by_night)))
                n_flagged_total += 1
                devs = by_night[mjd]
                pairs = unique([d.pair for d in devs])
                label = known_epoch_label(known_epochs, res.tele, mjd)
                if isnothing(label)
                    n_new_total += 1
                    tag = "[NEW]  "
                    suffix = ""
                else
                    n_known_total += 1
                    tag = "[KNOWN]"
                    suffix = "  # $(label)"
                end
                println(io,
                    "  $(tag) $(res.tele) $(mjd): pair $(join(_fmt_pair.(pairs), ", ")) " *
                    "(dominant $(_fmt_pair(res.dominant))) on $(length(devs)) config" *
                    "$(length(devs) == 1 ? "" : "s") " *
                    "$(join([string(d.config_id) for d in devs], ","))$(suffix)")
            end
        end
        println(io, "")
        println(io,
            "FLAGGED nights: $(n_flagged_total) (KNOWN: $(n_known_total), NEW: $(n_new_total))")
        if n_new_total > 0
            println(io,
                "ACTION REQUIRED: $(n_new_total) night$(n_new_total == 1 ? "" : "s") deviate(s) " *
                "outside any known epoch -- review before delivery.")
        else
            println(io,
                "All flagged nights fall inside known/accepted epochs; no new deviations.")
        end
    end
    println(io, "="^72)
    return (n_flagged_nights = n_flagged_total, n_known = n_known_total, n_new = n_new_total)
end
