# Post-hoc, read-only cross-night QA of finished arc-anchored comb products. NOT part of the pipeline:
# every night is processed independently (AKS 2026-11-18: "NO, we do not do this. This violates our design
# principal for parallelism."), so nothing here changes, re-anchors or flags any product. It reads the comb
# products (`<root>/comb/waveCalFPIArc_<tele>_<mjd>.h5`) of one or more finished runs and writes a CSV:
#   fsr_continuity      lambda at fiber 150, chip G, pixel 1024 jumps by ~ a nonzero integer FSR (within 0.2)
#                       from the previous comb night
#   cavity_jump_8sigma  L against the trailing `window` nights (>= min_prev): |L - median| > 8 max(1.4826 MAD, 1 Å)
#   vmode_vs_reference  the fixed-mode comb velocity against the stored global reference (data/
#                       comb_mode_velocity_ref_v1.csv): |delta| > --vmode_ms (what the per-night check flags)
# usage: julia --project=. scripts/arcfirst/crossnight_comb_qa.jl --tele apo --out flags.csv ROOT [ROOT ...]
using ApogeeReduction, HDF5
const _AR = ApogeeReduction

"""
    crossnight_comb_flags(mjds, lam, L; window = 15, min_prev = 5, fsr_tol = 0.2, nsig = 8.0)
        -> Vector{(mjd, check, value)}

Pure kernel: per comb night in MJD order, FSR continuity against the previous comb night and the cavity L
against the trailing `window` nights. Each night is judged against EARLIER nights only.
"""
function crossnight_comb_flags(mjds::AbstractVector{<:Integer}, lam::AbstractVector{<:Real},
        L::AbstractVector{<:Real}; window = 15, min_prev = 5, fsr_tol = 0.2, nsig = 8.0)
    o = sortperm(mjds)
    M, Λ, LL = mjds[o], lam[o], L[o]
    flags = Tuple{Int, String, Float64}[]
    for i in 2:length(M)
        fsr = Λ[i]^2 / (2 * LL[i])
        d = (Λ[i] - Λ[i - 1]) / fsr
        k = round(d)
        (k != 0 && abs(d - k) < fsr_tol) && push!(flags, (M[i], "fsr_continuity", d))
        lo = max(1, i - window)
        if i - lo >= min_prev
            prev = LL[lo:(i - 1)]
            m = _qa_median(prev)
            s = max(1.4826 * _qa_median(abs.(prev .- m)), 1.0)
            abs(LL[i] - m) > nsig * s && push!(flags, (M[i], "cavity_jump_8sigma", (LL[i] - m) / s))
        end
    end
    return flags
end
_qa_median(x) = (y = sort(collect(Float64, x)); n = length(y); isodd(n) ? y[(n + 1) ÷ 2] : (y[n ÷ 2] + y[n ÷ 2 + 1]) / 2)

"Comb nights of `roots` for `tele`: Dict mjd => (lambda_G_f150_p1024, L, eps)."
function crossnight_read_combs(roots, tele)
    rx = Regex("^waveCalFPIArc_$(tele)_(\\d+)\\.h5\$")
    seen = Dict{Int, NTuple{3, Float64}}()
    for root in roots
        d = joinpath(root, "comb")
        isdir(d) || continue
        for fn in readdir(d)
            m = match(rx, fn)
            m === nothing && continue
            mj = parse(Int, m.captures[1])
            haskey(seen, mj) && continue
            v = try
                h5open(joinpath(d, fn), "r") do f
                    lin = read(f["linParams"]); nl = read(f["nlParams"])
                    (_AR.refit_prod_wave(lin, nl, 150, 2, [1024.0])[1], Float64(read(f["fpi_cavity_size"])),
                        Float64(read(f["fpi_m0_offset"])))
                end
            catch
                continue
            end
            all(isfinite, v) && (seen[mj] = v)
        end
    end
    return seen
end

function crossnight_qa_main(args)
    tele, out, vms, roots = "apo", "crossnight_flags.csv", 30.0, String[]
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--tele"; tele = args[i + 1]; i += 2
        elseif a == "--out"; out = args[i + 1]; i += 2
        elseif a == "--vmode_ms"; vms = parse(Float64, args[i + 1]); i += 2
        else; push!(roots, a); i += 1
        end
    end
    isempty(roots) && error("give at least one finished run root")
    seen = crossnight_read_combs(roots, tele)
    mj = sort(collect(keys(seen)))
    fl = crossnight_comb_flags(mj, [seen[m][1] for m in mj], [seen[m][2] for m in mj])
    for m in mj
        ref = _AR.comb_mode_reference("default", tele, m)
        v = _AR.refit_comb_mode_velocity(seen[m][2], seen[m][3], tele)
        (isfinite(ref) && abs(v - ref) > vms) && push!(fl, (m, "vmode_vs_reference", v - ref))
    end
    open(out, "w") do io
        println(io, "# post-hoc cross-night comb QA (read-only; changes nothing); roots: ", join(roots, " "))
        println(io, "tele,mjd,check,value,n_comb_nights_seen")
        foreach(x -> println(io, "$(tele),$(x[1]),$(x[2]),$(x[3]),$(length(mj))"), sort(fl))
    end
    println("$(length(fl)) flags over $(length(mj)) comb nights -> $(out)")
    return fl
end

abspath(PROGRAM_FILE) == @__FILE__ && crossnight_qa_main(ARGS)
