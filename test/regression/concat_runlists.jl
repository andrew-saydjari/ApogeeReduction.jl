# concat_runlists.jl — concatenate per-telescope runlist files into ONE
# combined runlist (tele/mjd/expid/dfindx vectors, input order preserved).
#
# Why: the golden test days are (tele, mjd) PAIRS with different day sets per
# telescope (apo at 5 mjds, lco only at 60000), but the runlist makers'
# --mjds selection applies to every telescope in one call. submit_goldens.sh
# therefore builds one runlist per telescope (each with its own --mjds) and
# concatenates them here into the single combined runlist that the
# run_bulk.sh-style chain expects (per-tele pipeline stages self-filter on
# the tele column; make_relFlux iterates the unique teles itself).
#
# Usage:
#   julia --project=<AR checkout> concat_runlists.jl <output> <input1> [<input2> ...]
#
# Run with the AR project (needs JLD2 + ApogeeReduction.safe_jldsave so the
# output is byte-compatible with runlists written by the makers).
using Pkg;
Pkg.instantiate();
using JLD2
using ApogeeReduction: safe_jldsave

if length(ARGS) < 2
    println(stderr, "usage: concat_runlists.jl <output> <input1> [<input2> ...]")
    exit(2)
end
outfname = ARGS[1]
tele = String[]
mjd = Int[]
expid = String[]
dfindx = Int[]
for fn in ARGS[2:end]
    append!(tele, load(fn, "tele"))
    append!(mjd, load(fn, "mjd"))
    append!(expid, load(fn, "expid"))
    append!(dfindx, load(fn, "dfindx"))
end
safe_jldsave(outfname, Dict{String, Any}(); tele, mjd, expid, dfindx)
println("wrote $(outfname): $(length(tele)) rows from $(length(ARGS) - 1) input(s)")
if length(tele) == 0
    println("WARNING: combined runlist is empty")
    exit(16)
end
