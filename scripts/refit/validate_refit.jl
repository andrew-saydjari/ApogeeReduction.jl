# Post-run validation of an arc + FPI refit output root (run after refit_wavecal_archive.jl phase 3).
#
#   1. counts per wavecal_type and per bitmsk_wavecal bit (by telescope and era); stage failures;
#   2. comparison with pass 1 per exposure: zero point, shape and across-fiber scatter per chip
#      (distributions by telescope/era), and the fiber profile (median dv per fiber and chip);
#   3. nights where the arcs changed the comb's integer order, mode slips, cross-night flags;
#   4. the first-order RV repeat test (method of arc_verdict2/rv_analysis.jl): per stable-star row
#      the prediction is the equal-weight 3-chip mean of the flat pixel mean of
#      c (lambda_new - lambda_pass1)/lambda_pass1 at the row's fiber (`dv_pass1_ms`); per night the
#      median of (dv + prediction); reported per telescope and era: the night bias (median predicted
#      shift) AND the change in night zero-point scatter (robust sigma, clipped rms, and rms after
#      removing the bias), with bootstrap-over-nights 95 % intervals, plus the row-level change.
# Writes <outroot>/validation/*.csv and VALIDATION_SUMMARY.txt; figures (dark theme) to --plotdir.
#
# usage: julia --project=. scripts/refit/validate_refit.jl --outroot DIR [--plotdir DIR] [--nboot 200]
using ApogeeReduction, ArgParse, HDF5, CSV, DataFrames, Statistics, Random, Printf, Dates
using CairoMakie
const AR = ApogeeReduction

s = ArgParseSettings()
@add_arg_table! s begin
    "--outroot"; arg_type = String; required = true
    "--plotdir"; arg_type = String; default = ""
    "--rvrows"; arg_type = String; default = "/mnt/ceph/users/sdssv/work/asaydjari/2026_11_17/wavecal_consistency/rv_rows.h5"
    "--nboot"; arg_type = Int; default = 200
    "--noplots"; action = :store_true
end
const O = parse_args(s)
const OUT = O["outroot"]
const VD = joinpath(OUT, "validation")
mkpath(VD)
era_of(tele, mjd) = mjd < (tele == "apo" ? 59558 : 59809) ? "plate" :
                    (tele == "apo" ? (mjd < 60930 ? "fps1" : "fps2") : (mjd < 60700 ? "fps1" : "fps2"))
nanmed(x) = (y = filter(isfinite, x); isempty(y) ? NaN : median(y))
robsig(x) = (y = filter(isfinite, x); length(y) < 3 ? NaN : (q = quantile(y, [0.16, 0.84]); (q[2] - q[1]) / 2))
rmsclip(x) = (y = filter(isfinite, x); isempty(y) ? NaN : sqrt(mean(clamp.(y, -1000, 1000) .^ 2)))
q(x, p) = (y = filter(isfinite, x); isempty(y) ? NaN : quantile(y, p))

# ------------------------------------------------------------------ read every night product
nights = CSV.read(joinpath(OUT, "nights.csv"), DataFrame; types = Dict(:tele => String, :mjd => Int))
E = DataFrame(tele = String[], mjd = Int[], era = String[], exp = Int[], kind = String[], wtype = String[],
    bits = Int64[], p1type = String[], xcheck = Float32[], cmax = Float64[],
    zp_R = Float64[], zp_G = Float64[], zp_B = Float64[], sh_R = Float64[], sh_G = Float64[], sh_B = Float64[],
    rs_R = Float64[], rs_G = Float64[], rs_B = Float64[])
DV = Dict{Tuple{String, Int, Int}, Matrix{Float32}}()
fails = DataFrame(tele = String[], mjd = Int[], stage = String[], status = String[], reason = String[])
for r in eachrow(nights)
    t, m = r.tele, r.mjd
    for st in AR.REFIT_STAGES
        fn = AR.refit_status_path(OUT, st, t, m)
        isfile(fn) || continue
        l = split(strip(read(fn, String)), '\t')
        (startswith(l[1], "fail") || l[1] == "deferred") &&
            push!(fails, (t, m, st, String(l[1]), length(l) > 1 ? String(l[2]) : ""))
    end
    ef = AR.refit_path(OUT, "assemble", t, m)
    isfile(ef) || continue
    h5open(ef, "r") do f
        en = read(f["exp_expnum"]); kd = read(f["exp_kind"]); wt = read(f["wavecal_type"]); bm = read(f["bitmsk_wavecal"])
        p1 = read(f["pass1_wavecal_type"]); xc = read(f["xcheck_resid_ms"]); cm = read(f["coeff_max_resid_A"])
        z = read(f["pass1_zp_ms"]); sh = read(f["pass1_shape_ms"]); rs = read(f["pass1_fiber_rsig_ms"]); dv = read(f["dv_pass1_ms"])
        for k in eachindex(en)
            push!(E, (t, m, era_of(t, m), en[k], kd[k], wt[k], bm[k], p1[k], xc[k], cm[k], z[1, k], z[2, k], z[3, k],
                sh[1, k], sh[2, k], sh[3, k], rs[1, k], rs[2, k], rs[3, k]))
            kd[k] == "object" && (DV[(t, m, en[k] % 10000)] = dv[:, :, k])
        end
    end
end
CSV.write(joinpath(VD, "exposures.csv"), E)
CSV.write(joinpath(VD, "failures.csv"), fails)
println("exposures: $(nrow(E)) on $(length(unique(zip(E.tele, E.mjd)))) nights; failures/deferred: $(nrow(fails))")

# ------------------------------------------------------------------ 1. counts
cnt = combine(groupby(E, [:tele, :era, :kind, :wtype]), nrow => :n)
CSV.write(joinpath(VD, "counts_wavecal_type.csv"), cnt)
bitrows = DataFrame(tele = String[], era = String[], bit = String[], n = Int[], frac = Float64[])
for g in groupby(E, [:tele, :era])
    for (i, nm) in enumerate(AR.WAVE_BIT_NAMES)
        n = count(b -> (b & (1 << (i - 1))) != 0, g.bits)
        push!(bitrows, (g.tele[1], g.era[1], nm, n, n / nrow(g)))
    end
end
CSV.write(joinpath(VD, "counts_bits.csv"), bitrows)

# ------------------------------------------------------------------ 2. vs pass 1
cmp = DataFrame()
for g in groupby(E[(E.kind .== "object") .& .!startswith.(E.wtype, "sky_") .& .!startswith.(E.wtype, "none_"), :], [:tele, :era, :wtype])
    for (c, ch) in enumerate("RGB")
        z = g[!, Symbol("zp_$(ch)")]; sh = g[!, Symbol("sh_$(ch)")]; rs = g[!, Symbol("rs_$(ch)")]
        push!(cmp, (tele = g.tele[1], era = g.era[1], wtype = g.wtype[1], chip = string(ch), n = count(isfinite, z),
            zp_median = nanmed(z), zp_robsig = robsig(z), zp_p05 = q(z, 0.05), zp_p95 = q(z, 0.95),
            shape_median = nanmed(sh), shape_p95 = q(sh, 0.95), fiber_rsig_median = nanmed(rs), fiber_rsig_p95 = q(rs, 0.95)))
    end
end
CSV.write(joinpath(VD, "vs_pass1_by_type.csv"), cmp)
# fiber profile per telescope/era/chip: median over exposures of dv[c, f]
prof = DataFrame(tele = String[], era = String[], chip = Int[], fiber = Int[], dv_median = Float64[], n = Int[])
for (t, e) in unique(zip(E.tele, E.era))
    ks = [k for k in keys(DV) if k[1] == t && era_of(k[1], k[2]) == e]
    isempty(ks) && continue
    for c in 1:3, f in 1:AR.N_FIBERS
        v = Float64[DV[k][c, f] for k in ks]
        push!(prof, (t, e, c, f, nanmed(v), count(isfinite, v)))
    end
end
CSV.write(joinpath(VD, "vs_pass1_fiber_profile.csv"), prof)

# ------------------------------------------------------------------ 3. combs: order changes, slips, cross-night
ct = isfile(joinpath(OUT, "comb_table.csv")) ? CSV.read(joinpath(OUT, "comb_table.csv"), DataFrame) : DataFrame()
if nrow(ct) > 0
    CSV.write(joinpath(VD, "nights_order_changed.csv"), ct[coalesce.(ct.n_order_changed, 0) .> 0, :])
end

# ------------------------------------------------------------------ 4. RV repeat test
rvres = DataFrame()
if isfile(O["rvrows"])
    R = h5open(O["rvrows"], "r") do f; Dict(k => read(f[k]) for k in ("dv", "is_apo", "mjd", "expnum", "fiber", "clean")); end
    NR = length(R["dv"])
    tele_r = [x == 1 ? "apo" : "lco" for x in R["is_apo"]]
    mjd_r = Int.(R["mjd"]); exp_r = Int.(R["expnum"]) .% 10000; fib_r = Int.(R["fiber"])
    pred = fill(NaN, NR); maxchip = fill(Inf, NR)
    for k in 1:NR
        d = get(DV, (tele_r[k], mjd_r[k], exp_r[k]), nothing)
        d === nothing && continue
        fb = fib_r[k]; 1 <= fb <= AR.N_FIBERS || continue
        v = d[:, fb]; all(isfinite, v) || continue
        pred[k] = mean(v); maxchip[k] = maximum(abs.(v))
    end
    dvm = Float64.(R["dv"]) .* 1000
    ok = (R["clean"] .== 1) .& isfinite.(pred) .& (maxchip .<= 1000)
    era_r = [era_of(tele_r[k], mjd_r[k]) for k in 1:NR]
    println("RV rows with a refit prediction: $(count(ok)) of $(NR) ($(count((R["clean"] .== 1) .& isfinite.(pred) .& (maxchip .> 1000))) excluded at > 1 km/s on a chip)")
    subsets = [(t, e) for t in ("apo", "lco") for e in ("plate", "fps", "fps1", "fps2")]
    for (t, e) in subsets
        rows = findall(k -> ok[k] && tele_r[k] == t && (e == "fps" ? era_r[k] != "plate" : era_r[k] == e), 1:NR)
        length(rows) < 100 && continue
        g = Dict{Int, Vector{Int}}()
        for k in rows; push!(get!(g, mjd_r[k], Int[]), k); end
        ks = [m for (m, v) in g if length(v) >= 10]
        nn = length(ks); nn < 5 && continue
        z0 = [median(dvm[g[m]]) for m in ks]
        z1 = [median(dvm[g[m]] .+ pred[g[m]]) for m in ks]
        sh = z1 .- z0
        bias = median(sh)
        rng = MersenneTwister(7)
        Bs = [rand(rng, 1:nn, nn) for _ in 1:O["nboot"]]
        ci(f) = (v = [f(b) for b in Bs]; (quantile(v, 0.025), quantile(v, 0.975)))
        bci = ci(b -> median(sh[b]))
        dsig = robsig(z1) - robsig(z0); dsci = ci(b -> robsig(z1[b]) - robsig(z0[b]))
        drms = rmsclip(z1) - rmsclip(z0); drci = ci(b -> rmsclip(z1[b]) - rmsclip(z0[b]))
        z1b = z1 .- bias
        drmsb = rmsclip(z1b) - rmsclip(z0); drbci = ci(b -> rmsclip(z1b[b]) - rmsclip(z0[b]))
        x0 = dvm[rows]; x1 = dvm[rows] .+ pred[rows]
        push!(rvres, (tele = t, era = e, n_nights = nn, n_rows = length(rows),
            night_bias = bias, night_bias_lo = bci[1], night_bias_hi = bci[2],
            zp_robsig_pass1 = robsig(z0), zp_robsig_refit = robsig(z1), d_zp_robsig = dsig, d_zp_robsig_lo = dsci[1], d_zp_robsig_hi = dsci[2],
            zp_rms_pass1 = rmsclip(z0), d_zp_rms = drms, d_zp_rms_lo = drci[1], d_zp_rms_hi = drci[2],
            d_zp_rms_biasremoved = drmsb, d_zp_rms_biasremoved_lo = drbci[1], d_zp_rms_biasremoved_hi = drbci[2],
            row_robsig_pass1 = robsig(x0), d_row_robsig = robsig(x1) - robsig(x0)))
    end
    CSV.write(joinpath(VD, "rv_repeat_test.csv"), rvres)
else
    println("no RV rows file $(O["rvrows"]): RV test skipped")
end

# ------------------------------------------------------------------ summary
open(joinpath(VD, "VALIDATION_SUMMARY.txt"), "w") do io
    println(io, "refit validation $(now()); outroot $(OUT)")
    println(io, "exposures $(nrow(E)); stage failures/deferred $(nrow(fails))")
    println(io, "\n== counts per wavecal_type (all kinds)")
    for g in groupby(combine(groupby(E, [:tele, :era, :wtype]), nrow => :n), [:tele, :era])
        println(io, g.tele[1], " ", g.era[1], ": ", join(["$(r.wtype)=$(r.n)" for r in eachrow(g)], " "))
    end
    println(io, "\n== bits set (fraction of exposures), nonzero only")
    for r in eachrow(bitrows)
        r.n > 0 && println(io, @sprintf("%s %-5s %-28s %7d  %.3f", r.tele, r.era, r.bit, r.n, r.frac))
    end
    println(io, "\n== vs pass 1, object exposures with a new solution (m/s): zp median / robsig / p05..p95; shape median; fiber rsig median")
    for r in eachrow(cmp)
        println(io, @sprintf("%s %-5s %-14s %s n=%6d  zp %7.1f / %6.1f / %7.1f..%7.1f  shape %6.1f  fib %6.1f",
            r.tele, r.era, r.wtype, r.chip, r.n, r.zp_median, r.zp_robsig, r.zp_p05, r.zp_p95, r.shape_median, r.fiber_rsig_median))
    end
    if nrow(ct) > 0
        println(io, "\n== combs: nights $(nrow(ct)); arc-anchored $(count(==("arc"), ct.anchor)); arcwin $(count(==("arcwin"), ct.anchor))")
        println(io, "nights where the arcs changed the integer order on >= 1 fiber-chip: $(count(coalesce.(ct.n_order_changed, 0) .> 0))")
        println(io, "nights with a residual mode slip vs the arcs: $(count(coalesce.(ct.n_slip, 0) .> 0))")
    end
    if isfile(joinpath(OUT, "crossnight_flags.csv"))
        cf = CSV.read(joinpath(OUT, "crossnight_flags.csv"), DataFrame)
        println(io, "cross-night comb flags: $(nrow(cf))")
    end
    if nrow(rvres) > 0
        println(io, "\n== RV repeat test (first order; m/s; [95 % bootstrap over nights])")
        for r in eachrow(rvres)
            println(io, @sprintf("%s %-5s nights %4d rows %7d | bias %+6.1f [%+6.1f,%+6.1f] | night zp robsig %5.1f -> %5.1f (d %+5.1f [%+5.1f,%+5.1f]) | d rms %+5.1f, bias-removed %+5.1f [%+5.1f,%+5.1f] | d row robsig %+5.1f",
                r.tele, r.era, r.n_nights, r.n_rows, r.night_bias, r.night_bias_lo, r.night_bias_hi, r.zp_robsig_pass1,
                r.zp_robsig_refit, r.d_zp_robsig, r.d_zp_robsig_lo, r.d_zp_robsig_hi, r.d_zp_rms, r.d_zp_rms_biasremoved,
                r.d_zp_rms_biasremoved_lo, r.d_zp_rms_biasremoved_hi, r.d_row_robsig))
        end
    end
end
print(read(joinpath(VD, "VALIDATION_SUMMARY.txt"), String))

# ------------------------------------------------------------------ figures (dark theme)
if !O["noplots"] && !isempty(O["plotdir"]) && nrow(E) > 0
    set_theme!(merge(theme_black(), theme_latexfonts()))
    CairoMakie.disable_mime!("svg", "pdf", "text/html")
    mkpath(O["plotdir"])
    obj = E[(E.kind .== "object") .& .!startswith.(E.wtype, "sky_") .& .!startswith.(E.wtype, "none_"), :]
    cols = [:deepskyblue, :lime, :orangered]
    fig = Figure(size = (1400, 800))
    for (i, t) in enumerate(("apo", "lco"))
        ax = Axis(fig[i, 1], xlabel = "MJD", ylabel = "refit - pass 1 zero point (m/s)", title = uppercase(t))
        for (c, ch) in enumerate("RGB")
            sel = obj.tele .== t
            scatter!(ax, obj.mjd[sel] .+ 0.0, clamp.(obj[sel, Symbol("zp_$(ch)")], -800, 800), markersize = 2, color = (cols[c], 0.5), label = string(ch))
        end
        axislegend(ax, position = :lt)
        ax2 = Axis(fig[i, 2], xlabel = "refit - pass 1 zero point (m/s)", ylabel = "exposures")
        for (c, ch) in enumerate("RGB")
            v = filter(isfinite, obj[obj.tele .== t, Symbol("zp_$(ch)")])
            isempty(v) || hist!(ax2, clamp.(v, -800, 800), bins = 160, color = (cols[c], 0.5), label = string(ch))
        end
    end
    save(joinpath(O["plotdir"], "refit_vs_pass1_zp.png"), fig)
    if nrow(prof) > 0
        fig = Figure(size = (1400, 700))
        for (i, t) in enumerate(("apo", "lco"))
            ax = Axis(fig[1, i], xlabel = "fiber", ylabel = "median refit - pass 1 (m/s)", title = uppercase(t))
            for (c, ch) in enumerate("RGB"), e in ("plate", "fps1", "fps2")
                p = prof[(prof.tele .== t) .& (prof.era .== e) .& (prof.chip .== c), :]
                nrow(p) == 0 && continue
                lines!(ax, p.fiber, p.dv_median, color = cols[c], linestyle = e == "plate" ? :dot : (e == "fps1" ? :solid : :dash),
                    label = "$(ch) $(e)")
            end
            axislegend(ax, position = :lb, nbanks = 3, labelsize = 9)
        end
        save(joinpath(O["plotdir"], "refit_vs_pass1_fiber_profile.png"), fig)
    end
    foreach(f -> chmod(joinpath(O["plotdir"], f), 0o644), readdir(O["plotdir"]))
    println("figures in $(O["plotdir"])")
end
