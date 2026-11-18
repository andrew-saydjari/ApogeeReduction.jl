# From-scratch rebuild of the plate-era comb SHAPE template (design P2, stage 9 file) consumed by
# `refit_night_assemble` (`RefitConfig.shape_template_file`). The refit runs without it (arc
# shape, flagged WAVE_SHAPE_ARC).
#
# Definition (same as ARCFPI_VALIDATION build_template.jl, whose output template_shape_v1.h5 is the
# refit's default input): for every FPS comb night solution (production form linParams/nlParams)
# and fiber/chip, the exact pixel-displacement field q(p) - p against a fixed per-fiber reference
# comb solution (lambda_night(p) = lambda_ref(q)), on the 97-point grid p = 65:20:1985, is fitted
# with Legendre P0..P6 in t = (p - 1025)/960. P2..P6 do not change under a rigid detector shift,
# so they carry the comb SHAPE independently of where the spectra sit. The template is the
# per-(coefficient, chip, fiber) median over a block's nights (and its NMAD).
#
# Sources (`--source`):
#   pass1  pass-1 waveCalFPI_<tele>_<mjd>_arclamp.h5 of nights whose pass-1 FPI gate passed
#          (read-only); the shape does not depend on the comb's absolute level.
#   refit  waveCalFPIArc_<tele>_<mjd>.h5 of a refit output root (arc-anchored combs, no mode slip).
# Blocks: all; fps1/fps2 (APO 60930, LCO 60700); APO "home" (>= 60880, spectra where they sat
# in the plate era). Reference nights: APO 60440, LCO 60626 (as the study; arguments).
#
# usage: julia --project=. scripts/refit/build_comb_shape_template.jl --out FILE
#            [--source pass1|refit] [--pass1 DIR] [--refit_root DIR] [--ref_apo 60440] [--ref_lco 60626]
using ApogeeReduction, ArgParse, HDF5, Statistics, Dates, Printf
const AR = ApogeeReduction

s = ArgParseSettings()
@add_arg_table! s begin
    "--out"; arg_type = String; required = true
    "--source"; arg_type = String; default = "pass1"
    "--pass1"; arg_type = String; default = "/mnt/ceph/users/sdssv/work/asaydjari/pass1_dr21"
    "--refit_root"; arg_type = String; default = ""
    "--ref_apo"; arg_type = Int; default = 60440
    "--ref_lco"; arg_type = Int; default = 60626
    "--max_nights"; arg_type = Int; default = 0
end
const O = parse_args(s)
const FPS_START = Dict("apo" => 59558, "lco" => 59809)
const ERA_SPLIT = Dict("apo" => 60930, "lco" => 60700)
const PXG = collect(65.0:20.0:1985.0)
const TG = (PXG .- 1025) ./ 960
const LEG = hcat([AR.refit_legendre7(t) for t in TG]...)'    # 97 x 7

function night_solution(tele, mjd)
    if O["source"] == "pass1"
        na = AR.refit_night_ave(joinpath(O["pass1"], "wavecal"), tele, mjd)
        (na === nothing || !na.gate) && return nothing
        fn = joinpath(O["pass1"], "apred", string(mjd), "waveCalFPI_$(tele)_$(mjd)_arclamp.h5")
        isfile(fn) || return nothing
        return h5open(f -> (read(f["linParams"]), read(f["nlParams"])), fn, "r")
    else
        fn = AR.refit_path(O["refit_root"], "comb", tele, mjd)
        isfile(fn) || return nothing
        C = AR.load_refit_comb(fn)
        C.n_slip == 0 || return nothing
        return (C.lin, C.nl)
    end
end

function nights_of(tele)
    if O["source"] == "pass1"
        rx = Regex("^wavecalNightAve_$(tele)_(\\d+)\\.h5\$")
        ms = [parse(Int, m.captures[1]) for m in (match(rx, f) for f in readdir(joinpath(O["pass1"], "wavecal"))) if m !== nothing]
    else
        rx = Regex("^waveCalFPIArc_$(tele)_(\\d+)\\.h5\$")
        ms = [parse(Int, m.captures[1]) for m in (match(rx, f) for f in readdir(joinpath(O["refit_root"], "comb"))) if m !== nothing]
    end
    ms = sort(filter(>=(FPS_START[tele]), ms))
    O["max_nights"] > 0 && (ms = ms[1:min(end, O["max_nights"])])
    return ms
end

"Legendre coefficients (7 x 3 x 300) of the displacement of solution (lin, nl) against the reference."
function displacement_coefs(Tp, lin, nl)
    X = fill(NaN, 7, 3, AR.N_FIBERS)
    for f in 1:AR.N_FIBERS, c in 1:3
        (all(isfinite, view(lin, f, :)) && all(isfinite, view(Tp.ref_lin, f, :))) || continue
        lam = AR.refit_prod_wave(lin, nl, f, c, PXG)
        q = [AR.refit_ref_invert(Tp, f, c, lam[i], PXG[i]) for i in eachindex(PXG)]
        d = q .- PXG
        all(isfinite, d) || continue
        X[:, c, f] .= LEG \ d
    end
    return X
end

nanmed(x) = (y = filter(isfinite, x); isempty(y) ? NaN : median(y))
nmad(x) = (y = filter(isfinite, x); length(y) < 3 ? NaN : 1.4826 * median(abs.(y .- median(y))))

tmp = O["out"] * ".tmp$(getpid())"
h5open(tmp, "w") do fo
    prov = create_group(fo, "provenance")
    prov["script"] = read(@__FILE__, String)
    prov["script_fingerprint"] = AR.arclamp_file_fingerprint(@__FILE__)
    prov["code_commit"] = AR.git_commit
    prov["created"] = string(now())
    prov["source"] = O["source"] == "pass1" ? "pass1 waveCalFPI under $(O["pass1"]) (FPI gate passed)" :
                     "refit waveCalFPIArc under $(O["refit_root"]) (no mode slip)"
    prov["definition"] = "median over nights of Legendre P0..P6 (t=(p-1025)/960, p=65:20:1985) of q(p)-p, lambda_night(p)=lambda_ref(q)"
    for tele in ("apo", "lco")
        refm = tele == "apo" ? O["ref_apo"] : O["ref_lco"]
        R = night_solution(tele, refm)
        R === nothing && error("no reference solution for $(tele) $(refm)")
        Tp = AR.RefitTemplate(tele, "ref", zeros(7, 3, AR.N_FIBERS), R[1], R[2], "")
        C = Array{Float64, 3}[]
        mj = Int[]
        for m in nights_of(tele)
            S = night_solution(tele, m)
            S === nothing && continue
            push!(C, displacement_coefs(Tp, S[1], S[2]))
            push!(mj, m)
        end
        A = cat(C...; dims = 4)
        g = create_group(fo, tele)
        g["mjd"] = mj
        g["ref_lin"] = R[1]
        g["ref_nl"] = R[2]
        g["ref_mjd"] = refm
        blocks = ["all" => trues(length(mj)), "fps1" => mj .< ERA_SPLIT[tele], "fps2" => mj .>= ERA_SPLIT[tele]]
        tele == "apo" && push!(blocks, "home" => mj .>= 60880)
        for (nm, sel) in blocks
            any(sel) || continue
            g["T_$(nm)"] = [nanmed(A[k, c, f, sel]) for k in 1:7, c in 1:3, f in 1:AR.N_FIBERS]
            g["Tnmad_$(nm)"] = [nmad(A[k, c, f, sel]) for k in 1:7, c in 1:3, f in 1:AR.N_FIBERS]
            g["n_$(nm)"] = count(sel)
            @printf("%s %s: %d nights\n", tele, nm, count(sel))
        end
    end
end
mv(tmp, O["out"]; force = true)
println("wrote $(O["out"])")
