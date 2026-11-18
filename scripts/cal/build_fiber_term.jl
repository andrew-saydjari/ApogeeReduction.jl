# Stored per-fiber sky term (for RefitConfig.sky_fiber_term = "stored"), per telescope and era, built ONLY from sky
# fibers and faint-target fibers (target S/N < --snr_max, or no target row = sky / unplugged), because stellar lines
# bias the sky-line centroids on bright-target fibers (wave_stability Task 4). Read-only on the inputs.
# Value per (tele, era, fiber, chip) = median over object exposures of a finished run of the fiber-relative per-fiber
# sky zero point (fiber_sky_zp_ms, fiber median already removed per exposure), using only the exposures in which
# that fiber was sky/faint; columns tele,era,fiber,chip,zp_ms,robsig_ms,n_exp.
# usage: julia --project=. scripts/cal/build_fiber_term.jl --root RUN --out FILE [--allexp H5] [--snr_max 20]
#        [--min_exp 20]
using ApogeeReduction, HDF5, Statistics, Dates, Printf
const AR = ApogeeReduction
robsig(x) = (y = filter(isfinite, x); length(y) < 3 ? NaN : (q = quantile(y, [0.16, 0.84]); (q[2] - q[1]) / 2))

function main(args)
    root, out = "", "fiber_term.csv"
    ae = "/mnt/ceph/users/sdssv/work/scifest/0.4.0/public/allExposure_pass1_dr21.h5"
    smax, minn = 20.0, 20
    i = 1
    while i <= length(args)
        a, b = args[i], args[i + 1]
        a == "--root" && (root = b); a == "--out" && (out = b); a == "--allexp" && (ae = b)
        a == "--snr_max" && (smax = parse(Float64, b)); a == "--min_exp" && (minn = parse(Int, b))
        i += 2
    end
    D = h5open(ae, "r") do f; Dict(c => read(f[c]) for c in ("tele", "mjd", "expnum", "adjfiberindx", "snr")); end
    snr = Dict{NTuple{4, Int}, Float32}()
    for k in eachindex(D["mjd"])
        t = Int(D["tele"][k]); fb = Int(D["adjfiberindx"][k]) - 300 * (t - 1)
        snr[(t, Int(D["mjd"][k]), Int(D["expnum"][k]) % 10000, fb)] = Float32(D["snr"][k])
    end
    D = nothing; GC.gc()
    acc = Dict{Tuple{String, String, Int, Int}, Vector{Float32}}()
    rx = r"^waveCalRefit_(apo|lco)_(\d+)\.h5$"
    for fn in readdir(joinpath(root, "exp"))
        m = match(rx, fn); m === nothing && continue
        tele, mjd = String(m.captures[1]), parse(Int, m.captures[2])
        era = AR.refit_era(tele, mjd); ti = tele == "apo" ? 1 : 2
        h5open(joinpath(root, "exp", fn), "r") do f
            fz = read(f["fiber_sky_zp_ms"]); en = Int.(read(f["exp_expnum"])); kd = String.(read(f["exp_kind"]))
            for k in eachindex(en)
                kd[k] == "object" || continue
                for fb in 1:AR.N_FIBERS
                    s = get(snr, (ti, mjd, en[k] % 10000, fb), NaN32)
                    (isnan(s) || s < smax) || continue
                    for c in 1:3
                        z = fz[fb, c, k]
                        isfinite(z) && push!(get!(acc, (tele, era, fb, c), Float32[]), z)
                    end
                end
            end
        end
    end
    rows = String[]
    for key in sort(collect(keys(acc)))
        v = acc[key]; length(v) < minn && continue
        push!(rows, @sprintf("%s,%s,%d,%s,%.2f,%.2f,%d", key[1], key[2], key[3], ("R", "G", "B")[key[4]], median(v), robsig(v), length(v)))
    end
    open(out, "w") do io
        println(io, "# stored per-fiber sky term (m/s, fiber-relative), sky + faint-target fibers only (target S/N < $(smax));",
            " built $(now()) by scripts/cal/build_fiber_term.jl")
        println(io, "# code $(AR.git_branch) $(AR.git_commit) clean=$(AR.git_clean); run root $(root); S/N from $(ae)")
        println(io, "tele,era,fiber,chip,zp_ms,robsig_ms,n_exp")
        foreach(r -> println(io, r), rows)
    end
    println("$(length(rows)) rows -> $(out)")
end

main(ARGS)
