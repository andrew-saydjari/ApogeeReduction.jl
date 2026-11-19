# Stored FPI cavity model (the "cavity-first" calibration input; RefitConfig.cavity_ref_file, read per night by
# `comb_cavity_reference` / `comb_mode_reference`). Built from FINISHED arc-anchored comb products only; no night's
# processing reads another night: like the chip-scale file, this table is rebuilt from scratch before each release
# (and may be refreshed for the daily runs; cavity_first/CAVITY_FIRST.md: the level drifts <= 2.6 m/s/yr, a 60-day-old
# model is off by 3.5-4.7 m/s median, a year-old one by 3.7-5.5 m/s).
#
# Per comb night, the anchor level: verified strict candidates (vp_class == "strict", vp_verified) of a verified_pair
# product, median fixed-mode velocity; products without candidate tables (e.g. run1) contribute only nights listed as
# strict in --class_csv. Nights farther than --clip m/s from the telescope median are dropped (wrong pairs, moved pairs).
# Model (per telescope): eps FIXED at the median eps of those anchors (the comb cannot measure eps or L by itself, and
# the anchors' eps scatter is arc stretch noise, 9-24 m/s of R-B level); the level v(t) =
#   --model global : one constant (median),
#   --model smooth : running median over +-(--window/2) days, tabulated every --step days (rows mjd_lo..mjd_hi);
# L(t) then follows from (v, eps). Columns as data/comb_mode_velocity_ref_v1.csv: tele,mjd_lo,mjd_hi,vmode_ms,L_A,eps,
# n_nights,robsig_ms (n_nights / robsig_ms of the anchors inside the row's window).
#
# --from_monitor: also use stored-cavity products (comb_level = stored_cavity), whose per-night arc-pair levels are in
# the cavity monitor's candidate table (cavmon_vp_*; strict candidates, NOT pattern-verified: the monitor fits no comb).
#
# usage: julia --project=. scripts/cal/build_cavity_model.jl --out FILE [--model smooth|global] [--window 120] [--step 30]
#            [--class_csv F] [--clip 60] [--from_monitor] ROOT [ROOT ...]
using ApogeeReduction, HDF5, Statistics, Dates, Printf
const AR = ApogeeReduction
robsig(x) = (y = filter(isfinite, x); length(y) < 3 ? NaN : 1.4826 * median(abs.(y .- median(y))))

function anchor_levels(roots, strict; from_monitor = false)
    vals = Dict{Tuple{String, Int}, Tuple{Float64, Float64}}()   # (tele, mjd) -> (vmode, eps)
    rx = r"^waveCalFPIArc_(apo|lco)_(\d+)\.h5$"
    for root in roots, fn in readdir(joinpath(root, "comb"))
        m = match(rx, fn); m === nothing && continue
        key = (String(m.captures[1]), parse(Int, m.captures[2]))
        haskey(vals, key) && continue
        v = try
            h5open(joinpath(root, "comb", fn), "r") do f
                if from_monitor && haskey(f, "cavmon_vp_class") && !isempty(read(f["cavmon_vp_class"]))
                    c = String.(read(f["cavmon_vp_class"]))
                    s = findall(==("strict"), c)
                    isempty(s) && return nothing
                    return (median(Float64.(read(f["cavmon_vp_vmode_ms"]))[s]), median(Float64.(read(f["cavmon_vp_eps"]))[s]))
                end
                if haskey(f, "vp_class") && !isempty(read(f["vp_class"]))
                    c = String.(read(f["vp_class"])); ok = Bool.(read(f["vp_verified"]))
                    s = findall(i -> c[i] == "strict" && ok[i], eachindex(c))
                    isempty(s) && return nothing
                    return (median(Float64.(read(f["vp_vmode_ms"]))[s]), median(Float64.(read(f["vp_eps"]))[s]))
                end
                key in strict || return nothing
                L = Float64(read(f["fpi_cavity_size"])); e = Float64(read(f["fpi_m0_offset"]))
                return (AR.refit_comb_mode_velocity(L, e, key[1]), e)
            end
        catch
            nothing
        end
        (v !== nothing && all(isfinite, v)) && (vals[key] = v)
    end
    vals
end

"L such that (L, eps) has fixed-mode velocity v (inverse of AR.refit_comb_mode_velocity at fixed eps)."
function L_for(v, eps, tele)
    L0, e0 = AR.REFIT_MODE_REF[tele]
    lref = 2L0 / (AR.REFIT_MODE_M0 + e0)
    return lref * (v * (AR.REFIT_MODE_M0 + e0) / AR.REFIT_C_MPS + eps + AR.REFIT_MODE_M0) / 2
end

function main(args)
    out, model, window, step, ccsv, clip, roots = "cavity_model.csv", "smooth", 120, 30, "", 60.0, String[]
    from_monitor = false
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--out"; out = args[i + 1]; i += 2
        elseif a == "--model"; model = args[i + 1]; i += 2
        elseif a == "--window"; window = parse(Int, args[i + 1]); i += 2
        elseif a == "--step"; step = parse(Int, args[i + 1]); i += 2
        elseif a == "--class_csv"; ccsv = args[i + 1]; i += 2
        elseif a == "--clip"; clip = parse(Float64, args[i + 1]); i += 2
        elseif a == "--from_monitor"; from_monitor = true; i += 1
        else; push!(roots, a); i += 1
        end
    end
    model in ("smooth", "global") || error("--model smooth | global")
    isempty(roots) && error("give at least one finished run root")
    strict = Set{Tuple{String, Int}}()
    if !isempty(ccsv)
        L = readlines(ccsv); h = split(L[1], ",")
        it, im, ic = findfirst(==("tele"), h), findfirst(==("mjd"), h), findfirst(==("cls"), h)
        for l in L[2:end]
            p = split(l, ",")
            p[ic] == "strict" && push!(strict, (String(p[it]), parse(Int, p[im])))
        end
    end
    vals = anchor_levels(roots, strict; from_monitor)
    rows = String[]
    for tele in ("apo", "lco")
        ks = sort([k for k in keys(vals) if k[1] == tele], by = last)
        length(ks) < 3 && continue
        v = [vals[k][1] for k in ks]; e = [vals[k][2] for k in ks]; t = [k[2] for k in ks]
        keep = abs.(v .- median(v)) .<= clip
        v, e, t = v[keep], e[keep], t[keep]
        epsg = median(e)
        if model == "global"
            push!(rows, @sprintf("%s,%d,%d,%.2f,%.1f,%.6f,%d,%.2f", tele, 0, 99999, median(v), L_for(median(v), epsg, tele), epsg, length(v), robsig(v)))
            continue
        end
        grid = collect((minimum(t) - step ÷ 2):step:(maximum(t) + step ÷ 2))
        last_v = median(v)
        for (j, g) in enumerate(grid)
            w = findall(abs.(t .- g) .<= window / 2)
            vv = length(w) >= 3 ? median(v[w]) : last_v
            last_v = vv
            lo = j == 1 ? 0 : g - step ÷ 2
            hi = j == length(grid) ? 99999 : g + step - step ÷ 2 - 1
            push!(rows, @sprintf("%s,%d,%d,%.2f,%.1f,%.6f,%d,%.2f", tele, lo, hi, vv, L_for(vv, epsg, tele), epsg, length(w), robsig(v[w])))
        end
    end
    open(out, "w") do io
        commit = try
            readchomp(`git -C $(dirname(@__DIR__)) rev-parse HEAD`)
        catch
            "unknown"
        end
        println(io, "# stored FPI cavity model ($(model)", model == "smooth" ? ", window $(window) d, step $(step) d" : "", "), built $(now()) by scripts/cal/build_cavity_model.jl")
        println(io, "# code $(commit); roots: ", join(roots, " "), isempty(ccsv) ? "" : "; strict classes from $(ccsv)", from_monitor ? "; cavity-monitor strict candidates (unverified)" : "", "; clip $(clip) m/s; eps fixed per telescope")
        println(io, "tele,mjd_lo,mjd_hi,vmode_ms,L_A,eps,n_nights,robsig_ms")
        foreach(r -> println(io, r), rows)
    end
    println("wrote $(out): $(length(rows)) rows")
end
main(ARGS)
