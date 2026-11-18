# Mode-consistency calibration ("final pass", stored as a GLOBAL calibration file; per-night processing stays
# independent): measures, from the FINISHED per-exposure products of a bulk run, the global offsets between the
# wavelength-solution modes and writes data/wavecal_mode_offsets_v<N>.csv. Read-only on the run.
#
# Rows applied by the pipeline (`refit_mode_offsets`, arc route: lambda = arc model x (1 + (theta - O)/c)):
#   mode arc_sky, era plate : O = survey median over plate nights of the night median (arc-route object exposures)
#                             of theta = pooled sky - arc fit at its fixed state, per chip and term (zp, stretch, tilt)
#   mode arc_sky, era fps   : zp      = survey median of the night median of (pooled sky - FPI-route model) over
#                                       FPI-route object exposures of clean-anchor nights (verified strict / nomove
#                                       pair, or the --nights_csv subset): FPS arc-route exposures then land on the
#                                       level of the FPI-fiber route (the arc-anchored comb), so both FPS modes agree
#                             stretch, tilt = as for plate, from FPS arc-route object exposures
# Every row also carries its night-to-night robust sigma and counts: the offsets are a MONITORED output (track them
# when the line fits move from Gaussians to the true LSF). The FPS zp is NOT taken from FPS arc-route exposures:
# FPS science sits at dither positions 0.1-0.3 px from the arcs, so their sky - arc zero point carries km/s of
# position (measured on run1's strict nights with fps_route=sky: chip medians -1.1..+1.6 km/s, robust sigma 0.3-1 km/s).
# The FPS stretch / tilt rows come from the few FPS arc-route nights and carry that contamination at a smaller level.
# The cause of the per-chip / per-telescope differences is AKS's hypothesis (LSF: Gaussian line fits), not tested.
#
# usage: julia --project=. scripts/cal/build_mode_offsets.jl --root RUNROOT [--fpi_root RUNROOT2]
#            [--nights_csv clean.csv] [--epochs tele_epoch_mjd.csv] [--merge_ms 10] --out FILE
# (--fpi_root: run whose FPI-route exposures give the FPS zp, default = --root; the FPI route does not depend on
#  the arc-route settings, so an earlier run with the same comb products serves.)
using ApogeeReduction, HDF5, Statistics, Dates, Printf
const AR = ApogeeReduction

robsig(x) = (y = filter(isfinite, x); length(y) < 3 ? NaN : (q = quantile(y, [0.16, 0.84]); (q[2] - q[1]) / 2))
nanmed(x) = (y = filter(isfinite, x); isempty(y) ? NaN : median(y))

function main(args)
    root, out, ncsv, froot, epf, mrg = "", "wavecal_mode_offsets.csv", "", "", "", 10.0
    i = 1
    while i <= length(args)
        args[i] == "--root" && (root = args[i + 1])
        args[i] == "--out" && (out = args[i + 1])
        args[i] == "--nights_csv" && (ncsv = args[i + 1])
        args[i] == "--fpi_root" && (froot = args[i + 1])
        args[i] == "--epochs" && (epf = args[i + 1])
        args[i] == "--merge_ms" && (mrg = parse(Float64, args[i + 1]))
        i += 2
    end
    isdir(joinpath(root, "exp")) || error("--root $(root): no exp/ directory")
    isempty(froot) && (froot = root)
    clean = Set{Tuple{String, Int}}()
    if !isempty(ncsv)
        for l in readlines(ncsv)[2:end]
            p = split(l, ",")
            push!(clean, (String(p[1]), parse(Int, p[2])))
        end
    end
    # per night medians: arc-route theta (3x3), FPI-route sky - model (3)
    acc = Dict{Tuple{String, String, String}, Vector{Float64}}()   # (tele, era, key) -> night values
    nexp = Dict{Tuple{String, String, String}, Int}()
    nmjd = Dict{Tuple{String, String, String}, Vector{Int}}()      # matching night MJDs (epoch grouping)
    cur = Ref(0)
    push_n!(k, v, n) = (isfinite(v) && (push!(get!(acc, k, Float64[]), v); push!(get!(nmjd, k, Int[]), cur[]);
                                        nexp[k] = get(nexp, k, 0) + n))
    rx = r"^waveCalRefit_(apo|lco)_(\d+)\.h5$"
    for (pass, rr) in ((:arc, root), (:fpi, froot)), fn in readdir(joinpath(rr, "exp"))
        m = match(rx, fn); m === nothing && continue
        tele, mjd = String(m.captures[1]), parse(Int, m.captures[2])
        era = AR.refit_era(tele, mjd)
        cur[] = mjd
        (pass == :fpi && era != "fps") && continue
        h5open(joinpath(rr, "exp", fn), "r") do f
            wt = read(f["wavecal_type"]); kd = read(f["exp_kind"])
            th = read(f["sky_terms_ms"]); vr = read(f["sky_minus_model_ms"])
            ka = [k for k in eachindex(wt) if kd[k] == "object" && startswith(wt[k], "arc_") && occursin("skyshift", wt[k])]
            for (q, qn) in enumerate(("zp", "stretch", "tilt")), (c, cn) in enumerate(("R", "G", "B"))
                (isempty(ka) || pass == :fpi) && continue
                v = nanmed([th[q, c, k] for k in ka])
                push_n!((tele, era, "arcroute_$(qn)_$(cn)"), v, count(k -> isfinite(th[q, c, k]), ka))
            end
            (era == "fps" && pass == :fpi) || return
            cf = joinpath(rr, "comb", "waveCalFPIArc_$(tele)_$(mjd).h5")
            okn = (tele, mjd) in clean
            if isempty(clean) && isfile(cf)
                kind = h5open(g -> haskey(g, "anchor_kind") ? String(read(g["anchor_kind"])) : "", cf, "r")
                okn = startswith(kind, "verified_pair:strict") || startswith(kind, "verified_pair:nomove")
            end
            okn || return
            kf = [k for k in eachindex(wt) if kd[k] == "object" && endswith(wt[k], "_fpifib")]
            for (c, cn) in enumerate(("R", "G", "B"))
                isempty(kf) && continue
                push_n!((tele, era, "fpifib_zp_$(cn)"), nanmed([vr[c, k] for k in kf]), count(k -> isfinite(vr[c, k]), kf))
            end
        end
    end
    rows = String[]
    row(tele, era, mode, cn, qn, key) = begin
        v = get(acc, (tele, era, key), Float64[])
        length(v) < 5 && return
        push!(rows, @sprintf("%s,%s,%s,%s,%s,%.1f,%.1f,%d,%d,%s,,,%.1f", tele, era, mode, cn, qn, median(v), robsig(v),
            length(v), get(nexp, (tele, era, key), 0), key, 1.2533 * robsig(v) / sqrt(length(v))))
    end
    for tele in ("apo", "lco"), cn in ("R", "G", "B")
        for qn in ("zp", "stretch", "tilt")
            row(tele, "plate", "arc_sky", cn, qn, "arcroute_$(qn)_$(cn)")
        end
        row(tele, "fps", "arc_sky", cn, "zp", "fpifib_zp_$(cn)")
        for qn in ("stretch", "tilt")
            row(tele, "fps", "arc_sky", cn, qn, "arcroute_$(qn)_$(cn)")
        end
    end
    # instrument / LSF epochs (--epochs tele,epoch,mjd_start,mjd_end; e.g. the trace-template epochs): zero points
    # per epoch x era, adjacent epochs of one telescope and era MERGED when every chip's zp differs by less than
    # max(3 sigma, --merge_ms): the coarsest epoch set the data support. Stretch / tilt stay per era (they are
    # measured on the arc route only). Rows carry mjd_lo, mjd_hi and win over the era rows inside their range.
    nepoch = 0
    if !isempty(epf)
        EL = [l for l in readlines(epf) if !isempty(strip(l)) && !startswith(strip(l), "#")]
        EP = [split(l, ",") for l in EL[2:end]]
        divide = Dict(t => AR.get_fps_plate_divide(t) for t in ("apo", "lco"))
        for tele in ("apo", "lco"), era in ("plate", "fps")
            key(cn) = era == "plate" ? "arcroute_zp_$(cn)" : "fpifib_zp_$(cn)"
            groups = Any[]
            for e in EP
                e[1] == tele || continue
                lo, hi = parse(Int, e[3]), parse(Int, e[4])
                lo, hi = era == "plate" ? (lo, min(hi, divide[tele])) : (max(lo, divide[tele] + 1), hi)
                lo > hi && continue
                g = Dict(cn => [v for (v, m) in zip(get(acc, (tele, era, key(cn)), Float64[]), get(nmjd, (tele, era, key(cn)), Int[]))
                                if lo <= m <= hi] for cn in ("R", "G", "B"))
                length(g["G"]) < 5 && continue
                push!(groups, (name = String(e[2]), lo = lo, hi = hi, v = g))
            end
            merged = Any[]
            for g in groups
                if !isempty(merged)
                    p = merged[end]
                    same = all(cn -> begin
                        a, b = p.v[cn], g.v[cn]
                        err = sqrt((1.2533 * robsig(a))^2 / length(a) + (1.2533 * robsig(b))^2 / length(b))
                        abs(median(a) - median(b)) < max(3err, mrg)
                    end, ("R", "G", "B"))
                    if same
                        merged[end] = (name = p.name * "+" * g.name, lo = p.lo, hi = g.hi,
                            v = Dict(cn => vcat(p.v[cn], g.v[cn]) for cn in ("R", "G", "B")))
                        continue
                    end
                end
                push!(merged, g)
            end
            for g in merged, cn in ("R", "G", "B")
                v = g.v[cn]
                push!(rows, @sprintf("%s,%s,%s,%s,%s,%.1f,%.1f,%d,%d,%s,%d,%d,%.1f", tele, g.name, "arc_sky", cn, "zp", median(v),
                    robsig(v), length(v), 0, key(cn), g.lo, g.hi, 1.2533 * robsig(v) / sqrt(length(v))))
                nepoch += 1
            end
        end
    end
    open(out, "w") do io
        println(io, "# wavecal mode offsets (m/s), built $(now()) by scripts/cal/build_mode_offsets.jl")
        println(io, "# code $(AR.git_branch) $(AR.git_commit) clean=$(AR.git_clean); run root $(root); FPI-route root $(froot)",
            isempty(ncsv) ? "" : "; clean FPS nights from $(ncsv)")
        println(io, "# value = survey median of per-night medians; robsig = night-to-night robust sigma (monitored)")
        isempty(epf) || println(io, "# epochs from $(epf); adjacent epochs merged when every chip's zp differs by < max(3 sigma, $(mrg) m/s)")
        println(io, "tele,era,mode,chip,term,value_ms,robsig_ms,n_nights,n_exp,source,mjd_lo,mjd_hi,err_ms")
        foreach(r -> println(io, r), rows)
    end
    println("$(length(rows)) rows ($(nepoch) epoch rows) -> $(out)")
end

main(ARGS)
