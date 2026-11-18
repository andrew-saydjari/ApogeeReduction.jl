# Global stored comb cavity reference (per-epoch table) from FINISHED arc-anchored comb products, for the per-night
# cavity check and the no-adjacent-pair stored-cavity route of verified_pair (`comb_mode_reference`,
# `comb_cavity_reference`; RefitConfig.cavity_ref_file). Read-only on the runs; no night's processing reads another
# night: this table is a calibration input like the chip-scale file, rebuilt before each release.
#
# Per comb night, the strict-pair value:
#   * products with a candidate table (verified_pair): median over candidates with vp_class == "strict" and
#     vp_verified of vp_vmode_ms (fixed-mode velocity, m/s) and of vp_L;
#   * products without one (e.g. run1, anchor rule nearest_science): the anchor's own (L, eps) only on nights listed
#     as strict in --class_csv (tele,mjd,cls; e.g. anchor_verified/interim_classes.csv).
# Per telescope, epochs of --bin nights (MJD bins aligned on 59000): median vmode and L of the strict nights in the
# epoch (>= --min_nights, else the epoch is omitted and lookups use the nearest epoch); eps is stored consistent with
# (L, vmode) through refit_comb_mode_velocity. Columns: tele,mjd_lo,mjd_hi,vmode_ms,L_A,eps,n_nights,robsig_ms.
#
# usage: julia --project=. scripts/cal/build_comb_mode_reference.jl --out FILE [--class_csv F] [--bin 60]
#            [--min_nights 3] ROOT [ROOT ...]
using ApogeeReduction, HDF5, Statistics, Dates, Printf
const AR = ApogeeReduction
robsig(x) = (y = filter(isfinite, x); length(y) < 3 ? NaN : (q = quantile(y, [0.16, 0.84]); (q[2] - q[1]) / 2))

function main(args)
    out, ccsv, bin, minn, roots = "comb_mode_velocity_ref.csv", "", 60, 3, String[]
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--out"; out = args[i + 1]; i += 2
        elseif a == "--class_csv"; ccsv = args[i + 1]; i += 2
        elseif a == "--bin"; bin = parse(Int, args[i + 1]); i += 2
        elseif a == "--min_nights"; minn = parse(Int, args[i + 1]); i += 2
        else; push!(roots, a); i += 1
        end
    end
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
    vals = Dict{Tuple{String, Int}, Tuple{Float64, Float64}}()   # (tele, mjd) -> (vmode, L)
    rx = r"^waveCalFPIArc_(apo|lco)_(\d+)\.h5$"
    for root in roots, fn in readdir(joinpath(root, "comb"))
        m = match(rx, fn); m === nothing && continue
        key = (String(m.captures[1]), parse(Int, m.captures[2]))
        haskey(vals, key) && continue
        v = try
            h5open(joinpath(root, "comb", fn), "r") do f
                if haskey(f, "vp_class") && !isempty(read(f["vp_class"]))
                    c = String.(read(f["vp_class"])); ok = Bool.(read(f["vp_verified"]))
                    s = findall(i -> c[i] == "strict" && ok[i], eachindex(c))
                    isempty(s) && return nothing
                    vm = Float64.(read(f["vp_vmode_ms"]))[s]; LL = Float64.(read(f["vp_L"]))[s]
                    return (median(vm), median(LL))
                end
                key in strict || return nothing
                L = Float64(read(f["fpi_cavity_size"])); e = Float64(read(f["fpi_m0_offset"]))
                return (AR.refit_comb_mode_velocity(L, e, key[1]), L)
            end
        catch
            nothing
        end
        (v !== nothing && all(isfinite, v)) && (vals[key] = v)
    end
    rows = String[]
    for tele in ("apo", "lco")
        ks = sort([k for k in keys(vals) if k[1] == tele], by = last)
        isempty(ks) && continue
        bins = Dict{Int, Vector{Tuple{String, Int}}}()
        foreach(k -> push!(get!(bins, fld(k[2] - 59000, bin), Tuple{String, Int}[]), k), ks)
        for b in sort(collect(keys(bins)))
            kk = bins[b]
            length(kk) < minn && continue
            vm = [vals[k][1] for k in kk]; LL = [vals[k][2] for k in kk]
            v, L = median(vm), median(LL)
            push!(rows, @sprintf("%s,%d,%d,%.2f,%.1f,%.6f,%d,%.2f", tele, 59000 + b * bin, 59000 + (b + 1) * bin - 1,
                v, L, AR.refit_eps_for_vmode(L, v, tele), length(kk), robsig(vm)))
        end
    end
    open(out, "w") do io
        println(io, "# comb fixed-mode velocity reference (per-epoch, strict pairs), built $(now()) by scripts/cal/build_comb_mode_reference.jl")
        println(io, "# code $(AR.git_branch) $(AR.git_commit) clean=$(AR.git_clean); roots: $(join(roots, " "))",
            isempty(ccsv) ? "" : "; strict classes from $(ccsv)", "; bin $(bin) nights, min $(minn)")
        println(io, "tele,mjd_lo,mjd_hi,vmode_ms,L_A,eps,n_nights,robsig_ms")
        foreach(r -> println(io, r), rows)
    end
    println("$(length(rows)) epochs from $(length(vals)) strict nights -> $(out)")
end

main(ARGS)
