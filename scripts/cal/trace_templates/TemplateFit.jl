# TemplateFit.jl — the slit-model fitter used to BUILD the trace templates (scripts/cal/trace_templates).
#
# Vendored verbatim from the robust-trace study prototype (2026_11_17/trace_robust/RobustTrace.jl,
# the code that built the templates in the first packed calibration), except: the module name and
# the profile-table directory (the package's data/ folder instead of a fixed worktree path). It is
# kept separate from the pipeline's own slit model (src/) on purpose: the template build must be
# reproducible from the study's numbers until the pipeline fitter is validated to give the same
# templates; switching the builder to it is a follow-up (README, "Open items").
#
# (original header)
#
# Model:  y_f(x) = T_f(x) + d(x) + s(x) * (T_f(x) - 1024) + delta_f(x)
#   T_f(x)      epoch template of fiber centers (fixed per epoch, learned from the archive)
#   d(x), s(x)  nightly registration: shift (quadratic in x: offset, rotation, bow) and
#               scale (linear in x) -- 5 free numbers per night per chip
#   delta_f(x)  per-fiber residual, cubic in x, fitted only for fibers with light;
#               unlit fibers keep delta = 0 (template geometry), flagged
# Identification is never done by counting peaks: every fiber is placed by the registered
# template, and the registration is a global chi^2 over the whole slit with per-fiber free
# heights, so dead / parked / missing fibers simply fit with height ~ 0.
#
# Profile shape: the production Gauss-Hermite CDF tables (data/quartzTraceProfileParams_*),
# rebuilt here read-only, so centers are comparable with production.
module TemplateFit

using HDF5, Statistics, StatsBase, Polynomials, SpecialFunctions, LinearAlgebra

export load_flat, ProfileTables, load_profiles, band_profiles, register, fit_night,
       predict_centers

const BADBITS = 24566          # = bad_pix_bits in ApogeeReduction src/utils.jl
const NX = 2048
const NF = 300
# Gauss-Hermite profile tables (data/quartzTraceProfileParams_*): the package data folder
const DATA = Ref(normpath(joinpath(@__DIR__, "..", "..", "..", "data")))

# ------------------------------------------------------------------ I/O
function load_flat(path)
    h5open(path, "r") do f
        d = f["dimage"][1:NX, 1:NX]
        iv = f["ivarimage"][1:NX, 1:NX]
        m = f["pix_bitmask"][1:NX, 1:NX]
        good = ((m .& BADBITS) .== 0) .& (iv .> 0) .& isfinite.(d)
        return d, iv, good
    end
end

# ------------------------------------------------------------------ profiles (production GH tables)
_ghpoly(z, n) = n == 0 ? 0.0 : n == 1 ? 1.0 : n == 2 ? z : n == 3 ? z^2 - 1 :
                n == 4 ? z * (z^2 - 3) : n == 5 ? 3 + z^2 * (z^2 - 6) :
                n == 6 ? z * (15 + z^2 * (z^2 - 10)) :
                n == 7 ? -15 + z^2 * (45 + z^2 * (z^2 - 15)) : error("gh order")

struct ProfileTables
    cdf::Matrix{Float64}     # (nprof, nbins)
    zmin::Float64
    nsub::Int
    fmin::Int
    fmax::Int
end

function load_profiles(tele, chip)
    sjd, expid = tele == "apo" ? ("59549", "39870035") : ("60044", "44820015")
    p = joinpath(DATA[], "quartzTraceProfileParams_$(tele)_$(sjd)_$(expid)_$(chip).h5")
    h5open(p, "r") do f
        pf = Int.(read(f["fiber_index"])) .+ 1
        fibs = minimum(pf):maximum(pf)
        ng = Int(read(f["gh_order"])[1]) + 1
        H = zeros(length(fibs), ng)
        for j in 1:ng
            H[:, j] .= Polynomial(reverse(read(f["gh_$(j-1)_height_coeffs"]))).(fibs)
        end
        nsub = 100
        xb = range(-15.5, 15.5, step = 1 / nsub)
        C = zeros(length(fibs), length(xb))
        for (k, z) in enumerate(xb)
            g = -exp(-0.5 * z^2) / sqrt(2π)
            for i in eachindex(fibs)
                v = H[i, 1] * 0.5 * (1 + erf(z / sqrt(2)))
                for j in 2:ng
                    v += H[i, j] * _ghpoly(z, j - 1) * g   # gh_(j-1) term, = production int_gauss_hermite_term(z, j-1)
                end
                C[i, k] = v
            end
        end
        for i in axes(C, 1)
            C[i, :] .-= C[i, 1]
            C[i, :] ./= C[i, end]
        end
        return ProfileTables(C, first(xb), nsub, first(fibs), last(fibs))
    end
end

@inline function _cdf(P::ProfileTables, fib, z)
    k = clamp(round(Int, (z - P.zmin) * P.nsub) + 1, 1, size(P.cdf, 2))
    @inbounds P.cdf[clamp(fib, P.fmin, P.fmax) - P.fmin + 1, k]
end

# pixel-integrated unit profile of fiber `fib` at center mu, width w on pixels y0:y1
function unitprof!(out, P, fib, mu, w, y0, y1)
    prev = _cdf(P, fib, (y0 - 0.5 - mu) / w)
    @inbounds for (i, y) in enumerate(y0:y1)
        c = _cdf(P, fib, (y + 0.5 - mu) / w)
        out[i] = c - prev
        prev = c
    end
    out
end

# ------------------------------------------------------------------ band profiles
"""
Collapse the image into `nb` bands along the dispersion axis (x). Returns
(xc, prof[nb, 2048], err[nb, 2048]) using a masked mean with 4-sigma clipping about
the median. A smooth background (running 10th percentile, 61 px) is subtracted.
"""
function band_profiles(d, iv, good; nb = 32, bgsub = false, bw = 16)
    w = NX ÷ nb
    xc = Float64[]
    prof = fill(NaN, nb, NX)
    err = fill(Inf, nb, NX)
    for b in 1:nb
        c0 = (b - 1) * w + w ÷ 2
        cols = (c0 - bw ÷ 2 + 1):(c0 + bw ÷ 2)
        push!(xc, mean(cols))
        for y in 1:NX
            v = Float64[]
            s = Float64[]
            @inbounds for x in cols
                if good[x, y]
                    push!(v, d[x, y]); push!(s, 1 / sqrt(iv[x, y]))
                end
            end
            length(v) < 4 && continue
            m = median(v)
            mad = 1.4826 * median(abs.(v .- m)) + 1e-9
            k = abs.(v .- m) .< 4 * max(mad, median(s))
            sum(k) < 3 && continue
            prof[b, y] = mean(v[k])
            err[b, y] = median(s[k]) / sqrt(sum(k))
        end
        bgsub || continue
        p = prof[b, :]
        bg = similar(p)
        for y in 1:NX
            lo, hi = max(1, y - 30), min(NX, y + 30)
            seg = filter(isfinite, p[lo:hi])
            bg[y] = isempty(seg) ? 0.0 : quantile(seg, 0.1)
        end
        # smooth bg
        bgs = similar(bg)
        for y in 1:NX
            lo, hi = max(1, y - 30), min(NX, y + 30)
            bgs[y] = mean(bg[lo:hi])
        end
        prof[b, :] .-= bgs
    end
    return xc, prof, err
end

# ------------------------------------------------------------------ linear heights + chi2
"""
Given positions mu[f] and widths w[f] (fibers with mu outside the chip are skipped),
solve for heights by tridiagonal least squares against profile p with weights wt.
Returns (chi2, heights, nused).
"""
function heights_chi2(P, p, wt, mu, w; half = 6, fibs = 1:NF)
    n = length(mu)
    A = [zeros(2half + 1) for _ in 1:n]
    y0s = zeros(Int, n)
    use = falses(n)
    for f in 1:n
        yc = round(Int, mu[f])
        (yc - half < 1 || yc + half > NX) && continue
        use[f] = true
        y0s[f] = yc - half
        unitprof!(A[f], P, fibs[f], mu[f], w[f], yc - half, yc + half)
    end
    dg = zeros(n); off = zeros(max(n - 1, 0)); rhs = zeros(n)
    for f in 1:n
        use[f] || (dg[f] = 1.0; continue)
        @inbounds for i in 1:(2half + 1)
            y = y0s[f] + i - 1
            ww = wt[y]
            dg[f] += ww * A[f][i]^2
            rhs[f] += ww * A[f][i] * p[y]
        end
        if f < n && use[f + 1]
            for i in 1:(2half + 1)
                y = y0s[f] + i - 1
                j = y - y0s[f + 1] + 1
                (1 <= j <= 2half + 1) || continue
                off[f] += wt[y] * A[f][i] * A[f + 1][j]
            end
        end
    end
    dg .+= 1e-12
    h = SymTridiagonal(dg, off) \ rhs
    h[.!use] .= 0
    model = zeros(NX)
    for f in 1:n
        use[f] || continue
        for i in 1:(2half + 1)
            model[y0s[f] + i - 1] += h[f] * A[f][i]
        end
    end
    lo = max(1, floor(Int, minimum(mu[use]) - 10)); hi = min(NX, ceil(Int, maximum(mu[use]) + 10))
    # chi2 over the whole chip (light where the model has no fiber counts against it)
    chi2 = 0.0
    for y in 1:NX
        isfinite(p[y]) || continue
        chi2 += wt[y] * (p[y] - model[y])^2
    end
    return chi2, h, use, model
end

function band_weights(p, e; frac = 0.0)
    wt = zeros(NX)
    for y in 1:NX
        (isfinite(p[y]) && isfinite(e[y])) || continue
        wt[y] = 1 / (e[y]^2 + (frac * max(p[y], 0))^2)
    end
    wt
end

# ------------------------------------------------------------------ registration
"""
Register template centers T (vector over fibers, at one band) to profile p.
Grid over shift d in [-dmax, dmax] and scale s; returns best (d, s, chi2) and the
alias margin: chi2 at the best solution more than `alias_sep` px away, divided by best.
"""
function register(P, p, e, T, W; dmax = 40.0, dstep = 0.25, svals = -0.004:0.001:0.004,
        dcenter = 0.0, alias_sep = 4.0)
    wt = band_weights(p, e)
    pp = map(x -> isfinite(x) ? x : 0.0, p)
    best = (Inf, 0.0, 0.0)
    ds = collect((dcenter - dmax):dstep:(dcenter + dmax))
    chis = fill(Inf, length(ds))
    for s in (length(svals) > 5 ? (0.0,) : svals)
        for (k, d) in enumerate(ds)
            mu = T .+ d .+ s .* (T .- 1024)
            c, = heights_chi2(P, pp, wt, mu, W)
            if c < chis[k]
                chis[k] = c
            end
            if c < best[1]
                best = (c, d, s)
            end
        end
    end
    # refine d,s by local grid
    c0, d0, s0 = best
    for (dd, ss, ns) in ((0.1, 0.001, 4), (0.02, 0.0002, 5), (0.005, 0.00005, 4))
        for s in (s0 - ns * ss):ss:(s0 + ns * ss), d in (d0 - 5dd):dd:(d0 + 5dd)
            mu = T .+ d .+ s .* (T .- 1024)
            c, = heights_chi2(P, pp, wt, mu, W)
            if c < c0
                c0, d0, s0 = c, d, s
            end
        end
    end
    far = abs.(ds .- d0) .>= alias_sep
    alias = any(far) ? minimum(chis[far]) : Inf
    return (d = d0, s = s0, chi2 = c0, alias_chi2 = alias, margin = alias / c0,
        dgrid = ds, chigrid = chis)
end

# ------------------------------------------------------------------ local per-fiber refinement
"""
Gauss-Newton refinement of (height, center, width) for all fibers in one band with
neighbours held at their current model (Jacobi iterations). Centers are limited to
`maxmove` px from the start, widths to [0.6, 1.6] x start.
"""
function refine_band(P, p, e, mu0, w0; niter = 6, maxmove = 1.5, half = 4)
    n = length(mu0)
    wt = band_weights(p, e)
    pp = map(x -> isfinite(x) ? x : 0.0, p)
    mu = copy(mu0); w = copy(w0)
    _, h, use, model = heights_chi2(P, pp, wt, mu, w)
    buf = zeros(2half + 1); bp = similar(buf); bm = similar(buf); bw = similar(buf)
    for it in 1:niter
        newmu = copy(mu); neww = copy(w); newh = copy(h)
        for f in 1:n
            use[f] || continue
            yc = round(Int, mu[f]); y0 = yc - half; y1 = yc + half
            (y0 < 1 || y1 > NX) && continue
            unitprof!(buf, P, f, mu[f], w[f], y0, y1)
            # data minus neighbours
            r = pp[y0:y1] .- model[y0:y1] .+ h[f] .* buf
            ε = 0.01
            unitprof!(bp, P, f, mu[f] + ε, w[f], y0, y1)
            unitprof!(bm, P, f, mu[f] - ε, w[f], y0, y1)
            dmu = h[f] .* (bp .- bm) ./ (2ε)
            unitprof!(bw, P, f, mu[f], w[f] * (1 + ε), y0, y1)
            dw = h[f] .* (bw .- buf) ./ (ε * w[f])
            J = hcat(buf, dmu, dw)
            ww = wt[y0:y1]
            JTW = J' .* ww'
            M = JTW * J
            g = JTW * (r .- h[f] .* buf)
            dv = try
                (M + 1e-9I) \ g
            catch
                zeros(3)
            end
            all(isfinite, dv) || continue
            newh[f] = h[f] + dv[1]
            newmu[f] = clamp(mu[f] + clamp(dv[2], -0.5, 0.5), mu0[f] - maxmove, mu0[f] + maxmove)
            neww[f] = clamp(w[f] + clamp(dv[3], -0.2, 0.2), 0.6w0[f], 1.6w0[f])
        end
        mu, w = newmu, neww
        _, h, use, model = heights_chi2(P, pp, wt, mu, w)
    end
    # per-fiber local chi2 and S/N proxy
    snr = zeros(n)
    for f in 1:n
        use[f] || continue
        yc = round(Int, mu[f])
        y0 = max(1, yc - 2); y1 = min(NX, yc + 2)
        snr[f] = h[f] / sqrt(sum(e[y0:y1] .^ 2) + 1e-12)
    end
    return mu, w, h, use, snr, model
end

# ------------------------------------------------------------------ whole-night fit
"extend a template row (300 fibers) by `ne` slots at each end using the mean end spacing"
function extend_row(T, ne)
    s1 = mean(diff(T[1:6])); s2 = mean(diff(T[(end - 5):end]))
    vcat([T[1] - s1 * k for k in ne:-1:1], T, [T[end] + s2 * k for k in 1:ne])
end

safecorr(a, b) = (length(a) > 10 && std(a) > 0 && std(b) > 0) ? cor(a, b) : NaN

"""
fit_night(P, xc, prof, err, Tb, Wb; Href, dprior, ...)
 Tb, Wb: template centers / widths at band centers, (nb, 300)
 Href:   reference relative throughput per fiber (300; NaN = unknown) used to identify
         fibers (which trace is fiber 1) independently of peak counting
 dprior: optional expected central shift (e.g. from the previous good night)
Returns a NamedTuple with the registration, identification diagnostics, refined centers
per band, the smooth per-fiber model and QA metrics.
"""
function fit_night(P, xc, prof, err, Tb, Wb; Href = nothing, dprior = nothing, dmax = 40.0,
        lit_frac = 0.05, deg_delta = 3, NE = 3, kmax = 3)
    nb = length(xc)
    mid = clamp(searchsortedfirst(xc, 1024.0), 1, nb)
    cb = max(1, mid - 1):min(nb, mid)
    pc = vec(mean(prof[cb, :], dims = 1))
    ec = vec(sqrt.(mean(err[cb, :] .^ 2, dims = 1)) ./ sqrt(length(cb)))
    Tc = vec(mean(Tb[cb, :], dims = 1)); Wc = vec(mean(Wb[cb, :], dims = 1))
    # 1) coarse + fine registration of the template to the central profile
    reg = register(P, pc, ec, Tc, Wc; dmax = dmax, dcenter = isnothing(dprior) ? 0.0 : dprior)
    # 2) identification: heights in extended slots, compare label shifts k
    Te = extend_row(Tc, NE); We = vcat(fill(Wc[1], NE), Wc, fill(Wc[end], NE))
    fibe = vcat(fill(1, NE), 1:NF, fill(NF, NE))
    wt = band_weights(pc, ec); pp = map(x -> isfinite(x) ? x : 0.0, pc)
    mue = Te .+ reg.d .+ reg.s .* (Te .- 1024)
    _, he, usee, _ = heights_chi2(P, pp, wt, mue, We; fibs = fibe)
    hpos = max.(he, 0.0)
    href_med = median(hpos[(NE + 1):(NE + NF)])
    ks = -kmax:kmax
    tp = fill(NaN, length(ks)); chik = fill(Inf, length(ks)); edge = fill(NaN, length(ks))
    for (i, k) in enumerate(ks)
        sl = (NE + 1 + k):(NE + NF + k)
        if first(sl) < 1 || last(sl) > length(mue)
            continue
        end
        # chi2 with only the 300 slots of this labelling
        c, = heights_chi2(P, pp, wt, mue[sl], We[sl]; fibs = collect(1:NF))
        chik[i] = c
        # light left outside the labelled slots, relative to the median fiber
        out = setdiff(1:length(mue), sl)
        edge[i] = sum(hpos[out]) / href_med
        if !isnothing(Href)
            g = isfinite.(Href) .& (Href .> 0.02)
            a = log.(max.(hpos[sl][g], 1e-3 * href_med)); b = log.(Href[g])
            tp[i] = safecorr(a, b)
        end
    end
    i0 = findfirst(==(0), ks)
    kchi = ks[argmin(chik)]
    if !isnothing(Href) && any(isfinite, tp)
        o = sortperm(replace(tp, NaN => -Inf), rev = true)
        ktp = ks[o[1]]; tpmargin = tp[o[1]] - tp[o[2]]
    else
        ktp = kchi; tpmargin = NaN
    end
    # decision: throughput fingerprint if decisive, else chi2 (edge light), flag disagreement
    kbest = (isfinite(tpmargin) && tpmargin >= 0.05) ? ktp : kchi
    id_conflict = ktp != kchi
    sp = median(diff(Tc))
    # label-shifted template
    Tbe = hcat([extend_row(Tb[b, :], NE) for b in 1:nb]...)'   # nb x (300+2NE)
    Wbe = hcat([vcat(fill(Wb[b, 1], NE), Wb[b, :], fill(Wb[b, end], NE)) for b in 1:nb]...)'
    sl = (NE + 1 + kbest):(NE + NF + kbest)
    Tb2 = Matrix(Tbe[:, sl]); Wb2 = Wb   # widths follow the fiber, keep
    # 3) per-band registration around the central solution
    db = zeros(nb); sb = zeros(nb)
    for b in 1:nb
        r = register(P, prof[b, :], err[b, :], Tb2[b, :], Wb2[b, :]; dmax = 2.0, dstep = 0.25,
            svals = (reg.s - 0.002):0.001:(reg.s + 0.002), dcenter = reg.d)
        db[b] = r.d; sb[b] = r.s
    end
    u = (xc .- 1024) ./ 1024
    keep = trues(nb)
    pd = Polynomials.fit(u, db, 2); ps = Polynomials.fit(u, sb, 1)
    for _ in 1:2
        rd = db .- pd.(u)
        σ = 1.4826 * median(abs.(rd[keep])) + 0.02
        keep .= abs.(rd) .< 4σ
        pd = Polynomials.fit(u[keep], db[keep], 2); ps = Polynomials.fit(u[keep], sb[keep], 1)
    end
    # 4) refine every fiber in every band starting from the registered template
    MU = zeros(nb, NF); WW = zeros(nb, NF); HH = zeros(nb, NF); SN = zeros(nb, NF)
    MU0 = zeros(nb, NF)
    for b in 1:nb
        mu0 = Tb2[b, :] .+ pd(u[b]) .+ ps(u[b]) .* (Tb2[b, :] .- 1024)
        MU0[b, :] .= mu0
        mu, w, h, use, snr, _ = refine_band(P, prof[b, :], err[b, :], mu0, Wb2[b, :])
        MU[b, :] .= mu; WW[b, :] .= w; HH[b, :] .= h; SN[b, :] .= snr
    end
    hmed = vec(median(HH, dims = 1))
    ref = median(hmed[hmed .> 0])
    lit = hmed .> lit_frac * ref
    # 5) per-fiber residual delta_f(x), cubic, lit fibers only; unlit keep the template
    DEL = zeros(nb, NF)
    delta_poly = Vector{Polynomial{Float64}}(undef, NF)
    wpoly = Vector{Polynomial{Float64}}(undef, NF)
    for f in 1:NF
        r = MU[:, f] .- MU0[:, f]
        okb = (HH[:, f] .> lit_frac * ref) .& (abs.(r) .< 1.49)
        if lit[f] && sum(okb) >= deg_delta + 3
            k = copy(okb)
            pf = Polynomials.fit(u[k], r[k], deg_delta)
            for _ in 1:2
                rr = r .- pf.(u)
                σ = 1.4826 * median(abs.(rr[k])) + 0.005
                k .= okb .& (abs.(rr) .< 4σ)
                sum(k) >= deg_delta + 3 || break
                pf = Polynomials.fit(u[k], r[k], deg_delta)
            end
            delta_poly[f] = pf
            wpoly[f] = Polynomials.fit(u[k], WW[k, f], 2)
        else
            delta_poly[f] = Polynomial([0.0])
            wpoly[f] = Polynomials.fit(u, Wb2[:, f], 2)
        end
        DEL[:, f] .= delta_poly[f].(u)
    end
    # QA metrics
    rr = MU .- (MU0 .+ DEL)
    fit_rms = 1.4826 * median(abs.(rr[:, lit]))
    big_delta = sum(maximum(abs.(DEL), dims = 1) .> 1.0)
    return (reg = reg, d_center = pd(0.0), s_center = ps(0.0), kbest = kbest, ktp = ktp,
        kchi = kchi, tpmargin = tpmargin, tp = tp, chik = chik, edge = edge, ks = collect(ks),
        id_conflict = id_conflict, spacing = sp,
        d_poly = pd, s_poly = ps, db = db, sb = sb, Tb2 = Tb2, he = he, mue = mue,
        MU = MU, MU0 = MU0, WW = WW, HH = HH, SN = SN, lit = lit, nlit = sum(lit),
        delta_poly = delta_poly, wpoly = wpoly, DEL = DEL, xc = xc, fit_rms = fit_rms,
        big_delta = big_delta, label_shift = kbest)
end

"Evaluate final centers / widths at arbitrary columns x (vector) for all fibers: (nx, 300)."
function predict_centers(fitres, Tfun, Wfun, xs)
    u = (xs .- 1024) ./ 1024
    C = zeros(length(xs), NF); Wd = zeros(length(xs), NF)
    for (i, x) in enumerate(xs)
        T = extend_row(Tfun(x), 3)[(4 + fitres.label_shift):(3 + NF + fitres.label_shift)]
        C[i, :] .= T .+ fitres.d_poly(u[i]) .+ fitres.s_poly(u[i]) .* (T .- 1024) .+
                   [fitres.delta_poly[f](u[i]) for f in 1:NF]
        Wd[i, :] .= [fitres.wpoly[f](u[i]) for f in 1:NF]
    end
    C, Wd
end

end # module
