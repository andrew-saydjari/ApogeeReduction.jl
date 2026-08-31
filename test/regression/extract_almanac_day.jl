#!/usr/bin/env julia
# extract_almanac_day.jl — pull one (tele, mjd) out of a bulk almanac file into
# a single-day almanac file in the layout the current AR branch reads.
#
# Why this exists: the 2026_05_01 bulk almanac (allobs_57600_61160.h5) uses the
# post-PR#351 layout with a top-level `raw/` group (raw/<tele>/<mjd>/...), and
# also carries a top-level `exposure_class/` decoration group. The
# exposure-type-classifier branch of AR reads `<tele>/<mjd>/exposures` with NO
# `raw/` prefix (src/utils.jl read_almanac_exp_df), and make_runlist_all.jl
# treats every top-level key as a telescope name. Feeding it the bulk file
# directly therefore fails twice over. This script copies the day's group into
# a fresh file in the old layout, which is exactly what a daily-run almanac
# (allobs_<tele>_<mjd>.h5) looks like.
#
# Handles both source layouts (with and without the raw/ prefix).
#
# Usage:
#   julia --project=<this dir> extract_almanac_day.jl <src_almanac.h5> <tele> <mjd> <dest.h5>

using HDF5

function main()
    if length(ARGS) != 4
        println(stderr,
            "usage: extract_almanac_day.jl <src_almanac.h5> <tele> <mjd> <dest.h5>")
        exit(2)
    end
    src, tele, mjd, dest = ARGS
    isfile(src) || (println(stderr, "source almanac not found: $src"); exit(2))
    if isfile(dest)
        println("destination already exists, leaving untouched: $dest")
        exit(0)
    end
    h5open(src, "r") do f
        srcgroup = if haskey(f, "raw") && haskey(f["raw"], tele)
            "raw/$(tele)"
        elseif haskey(f, tele)
            tele
        else
            println(stderr, "telescope group '$tele' not found in $src " *
                            "(top-level keys: $(join(keys(f), ", ")))")
            exit(3)
        end
        haskey(f[srcgroup], mjd) ||
            (println(stderr, "mjd $mjd not found under $srcgroup in $src"); exit(3))
        mkpath(dirname(abspath(dest)))
        tmp = dest * ".part"
        h5open(tmp, "w") do fout
            g = create_group(fout, tele)
            copy_object(f["$(srcgroup)/$(mjd)"], g, mjd)
        end
        mv(tmp, dest; force = true)
        println("wrote $dest ($(tele)/$(mjd) from $(srcgroup)/$(mjd) of $src)")
    end
end

main()
