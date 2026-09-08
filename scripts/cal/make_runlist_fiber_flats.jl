using Pkg;
Pkg.instantiate();
using JLD2, ArgParse, DataFrames, HDF5
using ApogeeReduction: safe_jldsave, read_almanac_exp_df

## Parse command line arguments
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--tele"
        required = false
        help = "telescope name (apo or lco)"
        arg_type = String
        default = "both"
        "--almanac_file"
        required = true
        help = "path to the almanac file"
        arg_type = String
        default = ""
        "--output"
        required = true
        help = "path to output runlist file"
        arg_type = String
        default = ""
        "--flat_type"
        required = true
        help = "flat type, i.e. dome or quartz"
        arg_type = String
        default = "dome"
        "--mjds"
        required = false
        help = "comma-separated MJD list to restrict to (default \"0\" = all MJDs in the almanac file; a single value selects that day)"
        arg_type = String
        default = "0"
        "--use_exposure_class"
        required = false
        help = "drop exposures the exposure-type classifier marked predicted_bad (needs an almanac decorated by scripts/cal/decorate_almanac_exptype.jl; a no-op with a loud note when that group is absent)"
        arg_type = Bool
        default = true
    end
    return parse_args(s)
end

parg = parse_commandline()

"""
Per-(tele, mjd) `predicted_bad` mask written into the almanac by
`decorate_almanac_exptype.jl`, aligned to the `exposures` table by exposure
number. Returns `nothing` when the decoration is absent, which is the state of
every almanac built before that step joined the DAG — the caller must then keep
every exposure rather than silently dropping or silently keeping on a guess.
"""
function exposure_class_bad_set(f, tele, mjd)
    haskey(f, "exposure_class") || return nothing
    grp = "exposure_class/$(tele)/$(mjd)"
    haskey(f, grp) || return nothing
    g = f[grp]
    (haskey(g, "predicted_bad") && haskey(g, "exposure")) || return nothing
    pb = read(g["predicted_bad"])
    ex = read(g["exposure"])
    # 1 = judged bad. 0 = judged fine. Anything else (notably -1) = not judged,
    # and "not judged" is never grounds for exclusion.
    Set(Int(ex[i]) for i in eachindex(ex) if Int(pb[i]) == 1)
end

mjdexp_list = Int[]
expid_list = String[]
dfindx_list = Int[]
tele_list = String[]
n_dropped = 0
dropped_rows = String[]
n_nights_undecorated = 0
n_nights_total = 0
f = h5open(parg["almanac_file"])
tele2do = if parg["tele"] == "both"
    keys(f["raw"])
else
    [parg["tele"]]
end
for tele in tele2do
    mjd_list = keys(f["raw/$(tele)"])
    if parg["mjds"] != "0"
        # day-subset selection against a (possibly multi-day) almanac file
        mjds_sel = Set(String.(strip.(split(parg["mjds"], ","))))
        mjd_list = filter(in(mjds_sel), mjd_list)
    end
    for tstmjd in mjd_list
        tstmjd_int = parse(Int, tstmjd)
        df = read_almanac_exp_df(f, tele, tstmjd)
        good_exp = (df.image_type .== "$(parg["flat_type"])flat") .&
                   (df.lamp_une .== 0) .& (df.lamp_thar .== 0) .& (df.chip_flags .== 7) .& (df.flagged_bad .== 0)
        if parg["flat_type"] == "dome"
            good_exp .&= (df.n_read .> 3) .& (df.lamp_quartz .== 0)
        else
            #i.e. quartz
            good_exp .&= (df.n_read .>= 3) .& (df.lamp_quartz .== 1)
        end
        # Exposure-type classifier veto. This is the ONLY thing the predicted_bad
        # mask is allowed to exclude: bad flats from the trace/fluxing runlists.
        # It never removes an exposure from the reduction itself.
        global n_nights_total += 1
        badset = parg["use_exposure_class"] ?
                 exposure_class_bad_set(f, tele, tstmjd) : nothing
        if parg["use_exposure_class"] && isnothing(badset)
            global n_nights_undecorated += 1
        end

        dfindx_list_loc = findall(good_exp)
        for dfindx in dfindx_list_loc
            if !isnothing(badset) && (df.exposure[dfindx] in badset)
                global n_dropped += 1
                push!(dropped_rows,
                    "  DROPPED $(tele) $(tstmjd) exp $(df.exposure[dfindx]) " *
                    "($(parg["flat_type"])flat): exposure-type classifier predicted_bad=1")
                continue
            end
            push!(mjdexp_list, tstmjd_int)
            push!(expid_list, df.exposure_string[dfindx])
            push!(dfindx_list, dfindx)
            push!(tele_list, tele)
        end
    end
end

# Loud, never silent: a filter you cannot see in the log is the failure mode
# this whole path exists to fix.
println("make_runlist_fiber_flats: $(parg["flat_type"])flat runlist -> $(length(mjdexp_list)) exposures kept")
if !parg["use_exposure_class"]
    println("  exposure-class filter: DISABLED via --use_exposure_class false")
elseif n_nights_total == 0
    println("  exposure-class filter: no nights selected, nothing to filter")
elseif n_nights_undecorated == n_nights_total
    println("  exposure-class filter: NO-OP — no exposure_class group in $(parg["almanac_file"]) " *
            "for any of the $(n_nights_total) night(s). Run scripts/cal/decorate_almanac_exptype.jl " *
            "after pipeline.jl to populate it.")
else
    println("  exposure-class filter: ACTIVE on $(n_nights_total - n_nights_undecorated) of " *
            "$(n_nights_total) night(s) ($(n_nights_undecorated) undecorated, kept in full); " *
            "$(n_dropped) exposure(s) dropped")
    for row in dropped_rows
        println(row)
    end
end

safe_jldsave(parg["output"], Dict{String, Any}(); tele = tele_list, mjd = mjdexp_list, expid = expid_list, dfindx = dfindx_list)
