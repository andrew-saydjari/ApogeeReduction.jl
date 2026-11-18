# build_trace_template.jl — build ONE per-epoch trace template (one telescope, chip, epoch) from
# the 2D quartz flats alone. Needs no trace, no traceMain and no existing template file, so the
# template can be rebuilt from scratch (README, "Bootstrap").
#
#   julia --project=<repo> build_trace_template.jl --root <reduction outdir> --tele apo --chip R \
#       --epoch apoE5 --mjd_start 59534 --mjd_end 59815 [--nmax 0] --outdir <dir>
#
# Algorithm (TRACE_ROBUST_V2.md Sec. 7.2; the code that built the first packed calibration):
# 1. Candidates: the quartz flats of every night in [mjd_start, mjd_end]
#    (<root>/apred/<mjd>/ar2Dcal_<tele>_<mjd>_<exp>_<chip>_quartzflat.h5, or --flat_list), one per
#    night (the last exposure of the night), --nmax nights spread evenly over the epoch (0 = all).
# 2. ALL-LIT only, decided template-free: on the central 128 columns find the fiber peaks
#    (>= 15 % of the median peak), turn gaps into missing slots; a flat is all-lit if detected +
#    missing = 300 AND its set of missing slots equals the epoch's modal set (fibers dark on every
#    flat of the epoch are physically absent; a flat with any ADDITIONAL dark fiber is refused).
# 3. Bootstrap from the highest-S/N all-lit flat: number the slots 1..300 from the bottom of the
#    chip, refine every slot in the central band, walk band by band to the chip edges.
# 4. Drift within the epoch: every flat is registered to the current template (shift quadratic in
#    x, scale linear in x) and its free per-band centers are mapped back to the common frame BEFORE
#    the per-fiber median, so rigid drift never enters the template. Iterate (<= --maxit) until the
#    99th percentile of the template change (common shift removed) is < 0.003 px.
# 5. Pattern drift: per-fiber residual pattern of the last third of the epoch's flats minus the
#    first third (rms over lit fibers, px). Above --drift_warn the unit is flagged
#    (`drift_flag`): the design's remedy is a human decision (split the epoch / lower the floor).
# 6. Validation (--validate true): leave-nights-out (build on even-indexed nights, fit the odd
#    ones and vice versa) and no-template mode (fit each held-out flat with a template bootstrapped
#    from that flat alone and compare to the epoch-template fit).
# Writes <outdir>/template_<tele>_<chip>_<epoch>.h5 (+ _log.csv, _loo.csv, _notemplate.csv).
# Every input is opened read-only.
using ArgParse
include(joinpath(@__DIR__, "TemplateFit.jl"))
using .TemplateFit, HDF5, CSV, DataFrames, Statistics, Polynomials, Dates, SHA

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--root"
        help = "reduction output directory holding apred/<mjd>/ar2Dcal_*_quartzflat.h5 (read-only)"
        arg_type = String
        required = true
        "--tele"
        arg_type = String
        required = true
        "--chip"
        arg_type = String
        required = true
        "--epoch"
        help = "epoch name (the epochs CSV)"
        arg_type = String
        required = true
        "--mjd_start"
        arg_type = Int
        required = true
        "--mjd_end"
        arg_type = Int
        required = true
        "--nmax"
        help = "at most this many nights, spread evenly over the epoch (0 = all nights)"
        arg_type = Int
        default = 0
        "--flat_list"
        help = "optional CSV with columns tele,chip,kind,mjd,dfindx (e.g. the study's census_flats.csv); quartz rows are the candidates instead of the directory listing"
        arg_type = String
        default = ""
        "--outdir"
        arg_type = String
        required = true
        "--maxit"
        help = "maximum template iterations"
        arg_type = Int
        default = 4
        "--drift_warn"
        help = "pattern drift (px rms, last third - first third of the epoch) above which the unit is flagged (my choice from the design's ~0.01 px; needs AKS)"
        arg_type = Float64
        default = 0.01
        "--fit_rms_max"
        help = "after the first build, input flats whose template fit has fit_rms above this (px) are dropped and the template is rebuilt without them (in-focus flats only). Default = the trace gate's TRACE_GATE_FIT_RMS_MAX (src/trace_qa.jl). Inf = keep all (the first packed calibration was built that way). (my choice; needs AKS)"
        arg_type = Float64
        default = 0.02
        "--validate"
        help = "run the leave-nights-out and no-template validation"
        arg_type = Bool
        default = true
    end
    return parse_args(s)
end
parg = parse_commandline()
tele, chip, ename = parg["tele"], parg["chip"], parg["epoch"]
m0, m1 = parg["mjd_start"], parg["mjd_end"]
NMAX = parg["nmax"]; OUTD = parg["outdir"]; ROOT = parg["root"]
mkpath(OUTD)
const NB = 32
smed(v) = (w = filter(isfinite, v); isempty(w) ? NaN : median(w))
P = load_profiles(tele, chip)
logrows = DataFrame()

# ------------------------------------------------------------------ candidates
imgpath(m, e) = joinpath(ROOT, "apred", string(m), "ar2Dcal_$(tele)_$(m)_$(lpad(e, 4, '0'))_$(chip)_quartzflat.h5")
function candidate_flats()
    if parg["flat_list"] != ""
        F = CSV.read(parg["flat_list"], DataFrame)
        F = F[(F.tele .== tele) .& (F.chip .== chip) .& (F.kind .== "quartz") .& (m0 .<= F.mjd .<= m1), :]
        return [(r.mjd, r.dfindx) for r in eachrow(F)], "flat_list:" * parg["flat_list"]
    end
    out = Tuple{Int, Int}[]
    apred = joinpath(ROOT, "apred")
    for d in readdir(apred)
        m = tryparse(Int, d)
        (isnothing(m) || !(m0 <= m <= m1)) && continue
        pre = "ar2Dcal_$(tele)_$(m)_"; suf = "_$(chip)_quartzflat.h5"
        for f in readdir(joinpath(apred, d))
            (startswith(f, pre) && endswith(f, suf)) || continue
            e = tryparse(Int, f[(length(pre) + 1):(end - length(suf))])
            isnothing(e) || push!(out, (m, e))
        end
    end
    return out, "directory:" * apred
end
cands, cand_source = candidate_flats()
lastexp = Dict{Int, Int}()
for (m, e) in cands
    lastexp[m] = max(get(lastexp, m, e), e)
end
nights = sort(collect(keys(lastexp)))
use = nights
if NMAX > 0 && length(nights) > NMAX
    use = nights[unique(round.(Int, range(1, length(nights), length = NMAX)))]
end
println("$(tele) $(chip) $(ename) [$(m0),$(m1)]: $(length(nights)) nights with quartz flats, using $(length(use)) ($(cand_source))"); flush(stdout)

# ------------------------------------------------------------------ template-free peak census
function central_peaks(prof, perr, xc)
    cb = findall(x -> 960 <= x <= 1088, xc)
    p = vec(mean(prof[cb, :], dims = 1)); p[.!isfinite.(p)] .= 0
    w = exp.(-0.5 .* ((-3:3) ./ 0.8) .^ 2); w ./= sum(w)
    ps = [sum(w[k] * p[clamp(y + k - 4, 1, 2048)] for k in 1:7) for y in 1:2048]
    pk = [y for y in 4:2045 if ps[y] > ps[y - 1] && ps[y] >= ps[y + 1] && ps[y] > 0]
    h = ps[pk]
    hmed = median(h[h .> 0.3 * quantile(h, 0.9)])
    good = h .>= 0.15 * hmed
    pk = pk[good]; h = h[good]
    keep = trues(length(pk))          # merge peaks closer than 3 px (keep the higher)
    for i in 2:length(pk)
        if pk[i] - pk[i - 1] < 3
            h[i] > h[i - 1] ? (keep[i - 1] = false) : (keep[i] = false)
        end
    end
    pk = pk[keep]; h = h[keep]
    sp = median(diff(pk))
    slots = Int[1]; miss = Int[]
    for i in 2:length(pk)
        nstep = round(Int, (pk[i] - pk[i - 1]) / sp)
        for j in 1:(nstep - 1)
            push!(miss, slots[end] + j)
        end
        push!(slots, slots[end] + nstep)
    end
    return (peaks = pk, heights = h, slots = slots, missing = miss, nslots = slots[end], spacing = sp, hmed = hmed)
end

# ------------------------------------------------------------------ bootstrap from one flat
function bootstrap_template(prof, perr, xc, pc)
    pc.nslots == 300 || error("bootstrap needs 300 slots, got $(pc.nslots)")
    pos = zeros(300)
    pos[pc.slots] .= pc.peaks
    for m in pc.missing     # interpolate missing slots linearly between neighbours
        lo = maximum(filter(s -> s < m && !(s in pc.missing), 1:300)); hi = minimum(filter(s -> s > m && !(s in pc.missing), 1:300))
        pos[m] = pos[lo] + (pos[hi] - pos[lo]) * (m - lo) / (hi - lo)
    end
    nb = length(xc)
    T = zeros(nb, 300); W = fill(0.0, nb, 300); H = zeros(nb, 300)
    mid = argmin(abs.(xc .- 1024))
    for b in vcat(mid:nb, (mid - 1):-1:1)
        start = b == mid ? Float64.(pos) : (b > mid ? T[b - 1, :] : T[b + 1, :])
        w0 = b == mid ? fill(0.9, 300) : (b > mid ? W[b - 1, :] : W[b + 1, :])
        mu, w, h, _, _, _ = TemplateFit.refine_band(P, prof[b, :], perr[b, :], start, w0; niter = 8, maxmove = 2.0)
        dark = h .< 0.05 * median(h[h .> 0])      # dark slots keep the start position
        mu[dark] .= start[dark]; w[dark] .= w0[dark]
        T[b, :] .= mu; W[b, :] .= w; H[b, :] .= h
    end
    T, W, H
end

function smooth_rows(T, xc; deg = 4)        # per-fiber polynomial in u over bands
    u = (xc .- 1024) ./ 1024
    S = similar(T); C = zeros(deg + 1, size(T, 2))
    for f in 1:size(T, 2)
        g = isfinite.(T[:, f])
        p = Polynomials.fit(u[g], T[g, f], deg)
        S[:, f] .= p.(u); c = coeffs(p); C[1:length(c), f] .= c
    end
    S, C
end

# ------------------------------------------------------------------ load and census every candidate
flats = []
for m in use
    p = imgpath(m, lastexp[m])
    isfile(p) || continue
    d, iv, good = load_flat(p)
    xc, prof, perr = band_profiles(d, iv, good; nb = NB)
    pc = central_peaks(prof, perr, xc)
    push!(flats, (mjd = m, exp = lastexp[m], xc = xc, prof = prof, perr = perr, pc = pc))
end
missets = [Tuple(f.pc.missing) for f in flats]
modal = isempty(missets) ? () : first(sort(collect(Dict(s => count(==(s), missets) for s in unique(missets))), by = x -> -x[2]))[1]
alllit = [f.pc.nslots == 300 && Tuple(f.pc.missing) == modal for f in flats]
println("candidates $(length(flats)); modal missing slots $(modal); all-lit $(sum(alllit))"); flush(stdout)
for (f, a) in zip(flats, alllit)
    push!(logrows, (stage = "census", mjd = f.mjd, exp = f.exp, nslots = f.pc.nslots, nmissing = length(f.pc.missing), alllit = a,
        fit_rms = NaN, d_center = NaN, s_center = NaN, k = 0, resid_rms = NaN, note = join(f.pc.missing, ";")))
end
A = flats[alllit]
length(A) >= 3 || error("fewer than 3 all-lit flats; cannot build template")

# ------------------------------------------------------------------ bootstrap
snr = [f.pc.hmed / median(f.perr[16, :][isfinite.(f.perr[16, :])]) for f in A]
b0 = argmax(snr)
T0, W0, H0 = bootstrap_template(A[b0].prof, A[b0].perr, A[b0].xc, A[b0].pc)
T, _ = smooth_rows(T0, A[b0].xc); W, _ = smooth_rows(W0, A[b0].xc; deg = 2)
Href = vec(median(H0, dims = 1)); Href ./= median(Href[Href .> 0])
lit_epoch = trues(300); lit_epoch[collect(modal)] .= false
println("bootstrap from mjd $(A[b0].mjd) exp $(A[b0].exp)"); flush(stdout)
b0flat = "$(A[b0].mjd)_$(A[b0].exp)"; xc0 = A[b0].xc

u = (A[1].xc .- 1024) ./ 1024
function stack_fit(Fl, T, W, Href)
    Ts = fill(NaN, length(Fl), NB, 300); Ws = fill(NaN, length(Fl), NB, 300); Hs = fill(NaN, length(Fl), 300)
    stats = []
    for (i, f) in enumerate(Fl)
        res = fit_night(P, f.xc, f.prof, f.perr, T, W; Href = Href, dmax = 40.0)
        for b in 1:NB
            dd = res.d_poly(u[b]); ss = res.s_poly(u[b])
            Ts[i, b, :] .= (res.MU[b, :] .- dd .+ ss * 1024) ./ (1 + ss)   # back to the common frame
            Ws[i, b, :] .= res.WW[b, :]
        end
        Hs[i, :] .= vec(median(res.HH, dims = 1))
        Ts[i, :, .!res.lit] .= NaN
        push!(stats, (res.fit_rms, res.d_center, res.s_center, res.kbest))
    end
    Ts, Ws, Hs, stats
end
function combine_template(Ts, Ws, Hs)
    Tn = fill(NaN, NB, 300); Wn = fill(NaN, NB, 300)
    for b in 1:NB, f in 1:300
        v = filter(isfinite, Ts[:, b, f]); w = filter(isfinite, Ws[:, b, f])
        isempty(v) || (Tn[b, f] = median(v)); isempty(w) || (Wn[b, f] = median(w))
    end
    Hn = vec(mapslices(x -> median(filter(isfinite, x)), Hs, dims = 1)); Hn ./= median(Hn[Hn .> 0])
    Tn, Wn, Hn
end
function fill_unlit!(Tn, Wn, Tprev, Wprev)
    for f in 1:300
        if any(isnan, Tn[:, f])
            Tn[:, f] .= Tprev[:, f]; Wn[:, f] .= Wprev[:, f]
        end
    end
end
function build(Fl, T, W, Href; maxit = 4, tag = "")
    hist = Float64[]
    for it in 1:maxit
        Ts, Ws, Hs, st = stack_fit(Fl, T, W, Href)
        Tn, Wn, Hn = combine_template(Ts, Ws, Hs)
        fill_unlit!(Tn, Wn, T, W)
        Tn, _ = smooth_rows(Tn, Fl[1].xc); Wn, _ = smooth_rows(Wn, Fl[1].xc; deg = 2)
        dT = Tn[:, lit_epoch] .- T[:, lit_epoch]
        dT = dT .- median(dT, dims = 2)        # gauge: a common shift per band
        dmax = quantile(vec(abs.(dT)), 0.99)
        push!(hist, dmax)
        println("  $(tag) iter $(it): p99 template change (gauge removed) $(round(dmax, digits = 4)) px; median fit_rms $(round(median(first.(st)), digits = 4))"); flush(stdout)
        T, W, Href = Tn, Wn, Hn
        dmax < 0.003 && break
    end
    T, W, Href, hist
end
T, W, Href, hist = build(A, T, W, Href; maxit = parg["maxit"], tag = "all")
# in-focus flats only: drop inputs whose fit to the template is poor, rebuild without them
_, _, _, st0 = stack_fit(A, T, W, Href)
gated = [s[1] > parg["fit_rms_max"] for s in st0]
if any(gated)
    for (f, s) in zip(A[gated], st0[gated])
        push!(logrows, (stage = "gated_out", mjd = f.mjd, exp = f.exp, nslots = 300, nmissing = length(modal), alllit = true,
            fit_rms = s[1], d_center = s[2], s_center = s[3], k = s[4], resid_rms = NaN, note = "fit_rms > $(parg["fit_rms_max"])"))
    end
    println("dropping $(count(gated)) flat(s) with fit_rms > $(parg["fit_rms_max"]) px: ", join(["$(f.mjd)" for f in A[gated]], ", "))
    A = A[.!gated]
    length(A) >= 3 || error("fewer than 3 in-focus all-lit flats; cannot build template")
    T, W, Href, hist2 = build(A, T, W, Href; maxit = parg["maxit"], tag = "all(gated)")
    append!(hist, hist2)
end

# ------------------------------------------------------------------ drift and final residuals
Ts, Ws, Hs, st = stack_fit(A, T, W, Href)
resid = [sqrt(mean(filter(isfinite, vec(Ts[i, :, lit_epoch] .- T[:, lit_epoch])) .^ 2)) for i in 1:length(A)]
pm = [smed(Ts[i, :, f] .- T[:, f]) for i in 1:length(A), f in 1:300]
for (i, f) in enumerate(A)
    push!(logrows, (stage = "final", mjd = f.mjd, exp = f.exp, nslots = 300, nmissing = length(modal), alllit = true,
        fit_rms = st[i][1], d_center = st[i][2], s_center = st[i][3], k = st[i][4], resid_rms = resid[i], note = ""))
end
n = length(A); t1 = 1:max(1, n ÷ 3); t3 = (n - n ÷ 3 + 1):n
dpat = [smed(pm[t3, f]) - smed(pm[t1, f]) for f in findall(lit_epoch)]
drift_rms = sqrt(mean(filter(isfinite, dpat) .^ 2))
drift_flag = drift_rms > parg["drift_warn"]
println("drift: per-fiber pattern change first->last third rms $(round(drift_rms, digits = 4)) px", drift_flag ? "  ** above --drift_warn $(parg["drift_warn"]): consider splitting the epoch **" : "")

# ------------------------------------------------------------------ validation
loo = DataFrame(); nt = DataFrame()
if parg["validate"]
    ev = A[1:2:end]; od = A[2:2:end]
    for (trn, tst, nm_) in ((ev, od, "even->odd"), (od, ev, "odd->even"))
        Tt, Wt, Ht, _ = build(trn, smooth_rows(T0, xc0)[1], smooth_rows(W0, xc0; deg = 2)[1],
            (h = vec(median(H0, dims = 1)); h ./ median(h[h .> 0])); maxit = min(3, parg["maxit"]), tag = nm_)
        for f in tst
            res_h = fit_night(P, f.xc, f.prof, f.perr, Tt, Wt; Href = Ht)
            res_a = fit_night(P, f.xc, f.prof, f.perr, T, W; Href = Href)
            Ch, _ = predict_centers(res_h, x -> Tt[argmin(abs.(f.xc .- x)), :], x -> Wt[argmin(abs.(f.xc .- x)), :], f.xc)
            Ca, _ = predict_centers(res_a, x -> T[argmin(abs.(f.xc .- x)), :], x -> W[argmin(abs.(f.xc .- x)), :], f.xc)
            dc = vec((Ch .- Ca)[:, lit_epoch])
            push!(loo, (split = nm_, mjd = f.mjd, fit_rms_heldout = res_h.fit_rms, fit_rms_alltemplate = res_a.fit_rms,
                k = res_h.kbest, nmad_dc = 1.4826 * median(abs.(dc .- median(dc))), med_dc = median(dc), max_dc = maximum(abs.(dc))))
        end
        dT = vec((Tt .- T)[:, lit_epoch]); dT = dT .- median(dT)
        println("  LOO $(nm_): template vs all-flats template, rigid-removed: rms $(round(sqrt(mean(dT .^ 2)), digits = 4)) px, max $(round(maximum(abs.(dT)), digits = 3))")
    end
    for f in A[2:2:end]
        try
            Tb, Wb, _ = bootstrap_template(f.prof, f.perr, f.xc, f.pc)
            Tb, _ = smooth_rows(Tb, f.xc); Wb, _ = smooth_rows(Wb, f.xc; deg = 2)
            r0 = fit_night(P, f.xc, f.prof, f.perr, Tb, Wb; Href = nothing)
            ra = fit_night(P, f.xc, f.prof, f.perr, T, W; Href = Href)
            C0, _ = predict_centers(r0, x -> Tb[argmin(abs.(f.xc .- x)), :], x -> Wb[argmin(abs.(f.xc .- x)), :], f.xc)
            Ca, _ = predict_centers(ra, x -> T[argmin(abs.(f.xc .- x)), :], x -> W[argmin(abs.(f.xc .- x)), :], f.xc)
            dc = vec((C0 .- Ca)[:, lit_epoch])
            push!(nt, (mjd = f.mjd, fit_rms_notemplate = r0.fit_rms, fit_rms_template = ra.fit_rms, k_notemplate = r0.kbest,
                med_dc = median(dc), nmad_dc = 1.4826 * median(abs.(dc .- median(dc))), p99_dc = quantile(abs.(dc .- median(dc)), 0.99), max_dc = maximum(abs.(dc))))
        catch e
            push!(nt, (mjd = f.mjd, fit_rms_notemplate = NaN, fit_rms_template = NaN, k_notemplate = -99, med_dc = NaN, nmad_dc = NaN, p99_dc = NaN, max_dc = NaN))
        end
    end
end

# ------------------------------------------------------------------ write
base = "template_$(tele)_$(chip)_$(ename)"
_, Cc = smooth_rows(T, A[1].xc); _, Cw = smooth_rows(W, A[1].xc; deg = 2)
codesha = bytes2hex(sha256(read(@__FILE__))) * "/" * bytes2hex(sha256(read(joinpath(@__DIR__, "TemplateFit.jl"))))
tmp = joinpath(OUTD, base * ".h5.tmp$(getpid())")
h5open(tmp, "w") do f
    f["xc"] = A[1].xc; f["center"] = T; f["width"] = W; f["height_rel"] = Href
    f["center_poly_u"] = Cc; f["width_poly_u"] = Cw; f["lit_in_epoch"] = Int.(lit_epoch)
    f["input_mjd"] = [x.mjd for x in A]; f["input_exp"] = [x.exp for x in A]
    a = HDF5.attributes(f)
    a["tele"] = tele; a["chip"] = chip; a["epoch"] = ename; a["mjd_start"] = m0; a["mjd_end"] = m1
    a["n_candidates"] = length(flats); a["n_alllit"] = length(A); a["modal_missing_slots"] = join(modal, ";")
    a["bootstrap_flat"] = b0flat; a["iterations"] = length(hist); a["last_change_px"] = hist[end]
    a["fit_rms_max"] = parg["fit_rms_max"]; a["n_gated_out"] = count(gated)
    a["drift_pattern_rms_px"] = drift_rms; a["drift_flag"] = Int(drift_flag); a["drift_warn_px"] = parg["drift_warn"]
    a["code_sha256"] = codesha; a["candidate_source"] = cand_source; a["nmax"] = NMAX
    a["created"] = string(now()); a["u_definition"] = "u = (x - 1024)/1024; poly coeffs low->high"
end
CSV.write(joinpath(OUTD, base * "_log.csv"), logrows)
parg["validate"] && CSV.write(joinpath(OUTD, base * "_loo.csv"), loo)
parg["validate"] && CSV.write(joinpath(OUTD, base * "_notemplate.csv"), nt)
mv(tmp, joinpath(OUTD, base * ".h5"); force = true)    # last: the .h5 marks a finished unit
println("DONE $(base): n_alllit $(length(A)) of $(length(flats)); iterations $(length(hist)); drift $(round(drift_rms, digits = 4))")
