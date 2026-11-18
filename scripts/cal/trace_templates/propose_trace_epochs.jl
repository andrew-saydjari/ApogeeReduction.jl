# propose_trace_epochs.jl — propose instrument epochs (trace-template boundaries) from the data.
#
#   julia --project=<repo> propose_trace_epochs.jl --root <reduction outdir> --outdir <dir> \
#       [--tele apo,lco] [--mjd_start 57618] [--mjd_end 99999] [--k 10] [--thr 8] [--floor 0.15]
#
# Input: the night traces of a reduction (<root>/apred/<mjd>/traceMain_<tele>_<mjd>_<chip>.h5,
# opened read-only), e.g. a run with --trace_template_cal none (README, "Bootstrap"); it needs no
# template. Nights whose trace is not from one of their own flats (trace_level != 0, or the old
# quartz_fallback) are left out.
#
# Method (TRACE_ROBUST_V2.md Sec. 3.1, the script that proposed the first epoch table, with the
# lit-fiber-only pattern the design asks for): per (telescope, chip, night) a rigid fit of the
# centers at x = 512/1024/1536 against the archive-median geometry gives shift d0 (x = 1024),
# rotation (d1536 - d512)/1024, bow d1024 - (d512 + d1536)/2, scale s and the per-fiber residual
# pattern r_f (fibers below 5 % of the median height on that night are left out of the pattern).
# At each night, with k valid nights on each side:
#   z = |median(after) - median(before)| / sqrt(NMAD_a^2 + NMAD_b^2) * sqrt(k)  per parameter,
#   step = the change expressed as displacement at the chip edge (px);
# a candidate needs z >= thr AND step >= floor; greedy selection with separation k; each boundary
# is localized by the best two-level (L1) step fit of the driving parameter within +-k nights.
#
# Outputs in --outdir:
#   epoch_candidates.csv        one row per proposed boundary with the evidence (step sizes, z)
#   epoch_series.h5             the per-night registration series (for figures / review)
#   trace_epochs_draft.csv      tele,epoch,mjd_start,mjd_end from the candidates: a DRAFT. A person
#                               confirms or edits it (README, "Epochs") before templates are built.
using ArgParse, HDF5, CSV, DataFrames, Statistics
using ApogeeReduction: read_metadata

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--root"
        arg_type = String
        required = true
        "--outdir"
        arg_type = String
        required = true
        "--tele"
        arg_type = String
        default = "apo,lco"
        "--mjd_start"
        help = "first night considered (DR21 starts at 57618, after the Q3 amplifier fix)"
        arg_type = Int
        default = 57618
        "--mjd_end"
        arg_type = Int
        default = 99999
        "--k"
        help = "valid nights on each side of a boundary"
        arg_type = Int
        default = 10
        "--thr"
        help = "significance threshold"
        arg_type = Float64
        default = 8.0
        "--floor"
        help = "physical floor: smallest step (px at the chip edge) that makes a boundary"
        arg_type = Float64
        default = 0.15
        "--lit_frac"
        help = "fibers below this fraction of the night's median height are left out of the pattern"
        arg_type = Float64
        default = 0.05
    end
    return parse_args(s)
end
parg = parse_commandline()
k, thr, floor_px = parg["k"], parg["thr"], parg["floor"]
mkpath(parg["outdir"])
const XS = (512, 1024, 1536)
nmad(v) = 1.4826 * median(abs.(v .- median(v)))
fmed(v) = (w = filter(isfinite, v); isempty(w) ? NaN : median(w))

# ------------------------------------------------------------------ read the night traces
function read_series(tele, chip)
    apred = joinpath(parg["root"], "apred")
    M = Int[]; C = Array{Float64}[]
    for d in sort(readdir(apred))
        m = tryparse(Int, d)
        (isnothing(m) || !(parg["mjd_start"] <= m <= parg["mjd_end"])) && continue
        p = joinpath(apred, d, "traceMain_$(tele)_$(m)_$(chip).h5")
        isfile(p) || continue
        md = try
            read_metadata(p)
        catch
            continue
        end
        string(get(md, "trace_type", "")) == "quartz_fallback" && continue
        Int(get(md, "trace_level", 0)) == 0 || continue
        c = h5open(p, "r") do f
            ds = f["trace_params"]
            hcat([vec(ds[x, :, 2]) for x in XS]...)', vec(ds[1024, :, 1])
        end
        cen, h = c
        lit = h .> parg["lit_frac"] * fmed(h)
        push!(M, m); push!(C, vcat(cen, Float64.(lit)'))     # 4 x 300: centers at XS, lit mask
    end
    return M, C
end

function rigid(y, ref)
    g = isfinite.(y) .& isfinite.(ref)
    sum(g) < 50 && return (NaN, NaN, fill(NaN, length(y)))
    X = hcat(ones(sum(g)), ref[g] .- 1024)
    β = X \ (y[g] .- ref[g])
    for _ in 1:2
        rr = (y[g] .- ref[g]) .- X * β; σ = 1.4826 * median(abs.(rr)) + 0.01
        keep = abs.(rr) .< 5σ; β = X[keep, :] \ (y[g][keep] .- ref[g][keep])
    end
    (β[1], β[2], y .- (ref .+ β[1] .+ β[2] .* (ref .- 1024)))
end

allrows = DataFrame(); draft = DataFrame()
h5open(joinpath(parg["outdir"], "epoch_series.h5"), "w") do fout
for tele in split(parg["tele"], ",")
    per = Dict{String, Any}()
    for chip in ("R", "G", "B")
        M, Cs = read_series(tele, chip)
        n = length(M)
        n >= 2k + 1 || error("$(tele) $(chip): only $(n) usable nights")
        cen = cat([c[1:3, :] for c in Cs]...; dims = 3)        # 3 x 300 x n
        lit = cat([c[4, :] .> 0.5 for c in Cs]...; dims = 2)  # 300 x n
        ref = [vec(mapslices(fmed, cen[j, :, :], dims = 2)) for j in 1:3]
        P = zeros(n, 4); R = fill(NaN, n, 300)
        for i in 1:n
            d = [rigid(cen[j, :, i], ref[j]) for j in 1:3]
            P[i, 1] = d[2][1]; P[i, 2] = (d[3][1] - d[1][1]) / 1024; P[i, 3] = d[2][1] - (d[1][1] + d[3][1]) / 2
            P[i, 4] = d[2][2]; R[i, :] .= d[2][3]
            R[i, .!lit[:, i]] .= NaN                           # lit fibers only in the pattern
        end
        per[chip] = (M = M, P = P, R = R)
        fout["$(tele)/$(chip)/mjd"] = M; fout["$(tele)/$(chip)/P"] = P
    end
    Mc = sort(intersect([per[c].M for c in ("R", "G", "B")]...))
    idx = Dict(c => [findfirst(==(m), per[c].M) for m in Mc] for c in ("R", "G", "B"))
    n = length(Mc)
    names = ["shift", "rotation", "bow", "scale"]
    Z = zeros(n, 5); S = zeros(n, 5)
    edge = [1.0, 1024.0, 1.0, 1024.0]
    for c in ("R", "G", "B")
        P = per[c].P[idx[c], :]; R = per[c].R[idx[c], :]
        sp = zeros(n)
        for i in (k + 1):(n - k + 1)
            b = (i - k):(i - 1); a = i:(i + k - 1)
            for j in 1:4
                st = abs(median(P[a, j]) - median(P[b, j]))
                sc = sqrt(nmad(P[a, j])^2 + nmad(P[b, j])^2) / sqrt(k) + 1e-12
                S[i, j] = max(S[i, j], st * edge[j]); Z[i, j] = max(Z[i, j], st / sc)
            end
            ra = vec(mapslices(fmed, R[a, :], dims = 1)); rb = vec(mapslices(fmed, R[b, :], dims = 1))
            dd = filter(isfinite, ra .- rb)
            sp[i] = isempty(dd) ? 0.0 : sqrt(mean(dd .^ 2))
        end
        S[:, 5] .= max.(S[:, 5], sp)
        Z[:, 5] .= max.(Z[:, 5], sp ./ (median(sp[(k + 1):(n - k + 1)]) + 1e-9))
    end
    eff(i) = [(Z[i, j] >= thr && S[i, j] >= floor_px) ? S[i, j] : 0.0 for j in 1:5]
    zmax = [maximum(eff(i)) for i in 1:n]
    chosen = Int[]
    for i in sortperm(zmax, rev = true)
        zmax[i] <= 0 && break
        any(abs(i - j) < k for j in chosen) && continue
        push!(chosen, i)
    end
    function localize(i)
        jdrv = argmax(eff(i))
        jdrv == 5 && return i
        P = per["G"].P[idx["G"], jdrv]
        lo, hi = max(1, i - k), min(n, i + k)
        best, bi = Inf, i
        for t in (lo + 1):hi
            b = P[lo:(t - 1)]; a = P[t:hi]
            c = sum(abs.(b .- median(b))) + sum(abs.(a .- median(a)))
            c < best && ((best, bi) = (c, t))
        end
        bi
    end
    chosen = sort(unique(localize.(chosen)))
    for i in chosen
        e = eff(i); j = argmax(e)
        push!(allrows, (tele = tele, mjd_first_after = Mc[i], mjd_last_before = i > 1 ? Mc[i - 1] : -1,
            gap_days = i > 1 ? Mc[i] - Mc[i - 1] : 0, step_px = round(maximum(e), digits = 3),
            driver = j <= 4 ? names[j] : "pattern",
            S_shift = round(S[i, 1], digits = 3), S_rot_edge = round(S[i, 2], digits = 3), S_bow = round(S[i, 3], digits = 3),
            S_scale_edge = round(S[i, 4], digits = 3), S_pattern = round(S[i, 5], digits = 3),
            z_shift = round(Z[i, 1], digits = 1), z_rot = round(Z[i, 2], digits = 1), z_bow = round(Z[i, 3], digits = 1),
            z_scale = round(Z[i, 4], digits = 1), z_pattern = round(Z[i, 5], digits = 2), k = k, thr = thr, floor = floor_px))
    end
    starts = vcat(parg["mjd_start"] > Mc[1] ? parg["mjd_start"] : Mc[1], [Mc[i] for i in chosen])
    for (q, s0) in enumerate(starts)
        push!(draft, (tele = tele, epoch = "$(tele)E$(q)", mjd_start = s0,
            mjd_end = q < length(starts) ? starts[q + 1] - 1 : Mc[end]))
    end
    fout["$(tele)/mjd"] = Mc; fout["$(tele)/Z"] = Z; fout["$(tele)/S"] = S
end
end
CSV.write(joinpath(parg["outdir"], "epoch_candidates.csv"), allrows)
CSV.write(joinpath(parg["outdir"], "trace_epochs_draft.csv"), draft)
show(allrows, allrows = true, allcols = true); println()
println("DRAFT epoch table (confirm before building): ", joinpath(parg["outdir"], "trace_epochs_draft.csv"))
