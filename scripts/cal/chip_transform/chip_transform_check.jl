# Stage 4 of the chip-transform rebuild: re-evaluate sky-line (and, if line tables are given,
# arc-lamp) wavelength solutions of a sample of nights with and without a calibration file,
# on HELD-OUT lines. Included by build_chip_transform.jl (needs chip_transform_stages.jl).
#
# Per tracer group (one object exposure for sky, one arc sequence for arcs) and fiber, the lines
# are fitted with the production model -- one quartic in the green-chip coordinate shared by the
# three chips, x' = a_c + b_c x on chips R and B -- under four chip-transform hypotheses:
#
#   scale1        b_c = 1, free R/B zero points            (production sky fit, cporder 0)
#   cal_scale     b_c from the file, free R/B zero points  (the change in this PR)
#   cal_scale_off b_c and a_c from the file, no zero points
#   free_scale    cal_scale + free R/B stretch             (cporder 1 freedom; production arc fit)
#
# a_c is the file's offset in every hypothesis; the free zero points absorb its error to first
# order, so the fits are linear and leave-one-line-out residuals are exact (e_i / (1 - h_ii)).
# All four fits use the same lines: the intersection of the lines kept by each model's robust
# fit (`clip` robust sigma).
#
# Per night and chip it reports: held-out (leave-one-line-out) line RMS; the change of the mean
# per-chip zero point against scale1 (the bias, m/s); and, where the night has a comb solution,
# the shape of (solution - comb) after removing a per-chip zero point and stretch (m/s). It also
# reports the chip scale the free_scale fit lands on (b_free), the tracer's own estimate of b.

const CHECK_MODELS = (:scale1, :cal_scale, :cal_scale_off, :free_scale)
const CHECK_PIX = collect(64:16:1984)          # evaluation grid (pixels), avoids the chip ends

function check_design(model, x, c, nl4)
    n = length(x)
    xp = similar(x, Float64)
    for i in 1:n
        if c[i] == 2
            xp[i] = x[i]
        else
            j = c[i] == 1 ? 0 : 2
            b = model == :scale1 ? 1.0 : nl4[j + 2]
            xp[i] = nl4[j + 1] + b * x[i]
        end
    end
    A = positional_poly_mat(xp, porder = 4)
    model == :cal_scale_off && return A
    E = zeros(n, model == :free_scale ? 4 : 2)
    for i in 1:n
        c[i] == 1 && (E[i, 1] = 1; model == :free_scale && (E[i, 3] = x[i]))
        c[i] == 3 && (E[i, 2] = 1; model == :free_scale && (E[i, 4] = x[i]))
    end
    return hcat(A, E)
end

function check_rfit(A, y; clip, iters = 3)
    keep = trues(length(y)); p = nothing
    for _ in 1:iters
        count(keep) < size(A, 2) + 3 && return nothing, keep
        p = A[keep, :] \ y[keep]
        r = y .- A * p; s = robust_sigma(r[keep])
        (isfinite(s) && s > 0) || break
        nk = abs.(r .- median(r[keep])) .< clip * s
        nk == keep && break
        keep = nk
    end
    return p, keep
end

"""Fit all models to one (group, fiber); returns nothing or a Dict model => (loo, chips, grid, p)."""
function check_fit_group(x, c, y, nl4; clip, min_per_chip)
    all(count(==(cc), c) >= min_per_chip for cc in 1:3) || return nothing
    keep = trues(length(y))
    for md in CHECK_MODELS
        p, k = check_rfit(check_design(md, x, c, nl4), y; clip = clip)
        p === nothing && return nothing
        keep .&= k
    end
    all(count(==(cc), c[keep]) >= min_per_chip for cc in 1:3) || return nothing
    xk = x[keep]; ck = c[keep]; yk = y[keep]
    out = Dict{Symbol, Any}()
    xg = ((CHECK_PIX .- (N_XPIX ÷ 2)) ./ N_XPIX)
    for md in CHECK_MODELS
        A = check_design(md, xk, ck, nl4)
        size(A, 1) > size(A, 2) + 2 || return nothing
        p = A \ yk
        r = yk .- A * p
        H = try
            A * (pinv(A' * A) * A')
        catch
            return nothing
        end
        h = diag(H)
        g = h .< 0.95
        loo = fill(NaN, length(r)); loo[g] .= r[g] ./ (1 .- h[g])
        G = hcat([check_design(md, xg, fill(cc, length(xg)), nl4) * p for cc in 1:3]...)
        out[md] = (loo = loo, chips = ck, grid = G, p = p)
    end
    return out
end

"""dλ/dx' of the fitted quartic at the chip centre x' = a_c (for b_free)."""
function quartic_slope(p, xp)
    return p[2] + 2p[3] * xp + 3p[4] * xp^2 + 4p[5] * xp^3
end

"""Remove a per-chip zero point and stretch (linear in pixel) from a velocity pattern."""
function remove_zp_stretch(v)
    u = (CHECK_PIX .- 1024.5) ./ 1024
    X = hcat(ones(length(u)), u)
    g = isfinite.(v)
    count(g) > 10 || return fill(NaN, length(v))
    b = X[g, :] \ v[g]
    return v .- X * b
end

"""
Sky lines of one exposure, exactly as `ingest_skyLines_exp` delivers them to the production fit
(sky_line_mat_clean x and wavelength, chips stacked, rows whose fiber-to-fiber IQR of x exceeds
0.002 masked), but opened strictly READ-ONLY: the package reader opens products "r+", which a
validation run over a finished reduction must never do.
"""
function read_sky_lines_readonly(fname)
    X = Matrix{Float64}[]; Y = Matrix{Float64}[]; C = Matrix{Int}[]
    for (ci, chip) in enumerate(CHIP_LIST)
        fn = replace(fname, "_$(FIRST_CHIP)_" => "_$(chip)_")
        isfile(fn) || continue
        M = h5open(f -> read(f["sky_line_mat_clean"]), fn, "r")
        push!(X, (M[:, 1, :] .- (N_XPIX ÷ 2)) ./ N_XPIX); push!(Y, M[:, 2, :]); push!(C, fill(ci, size(M, 1), size(M, 3)))
    end
    X = vcat(X...); Y = vcat(Y...); C = vcat(C...)
    if size(X, 1) > 1
        bad = dropdims(ApogeeReduction.nanzeroiqr(X, 2) .> 0.002, dims = 2)
        X[bad, :] .= NaN; Y[bad, :] .= NaN
    end
    return X, Y, C
end

"""Sky line groups of a night: up to `max_exp` object exposures evenly spaced through the night."""
function check_sky_groups(root, tele, mjd; max_exp)
    d = night_dir(root, mjd); isdir(d) || return []
    r = Regex("^skyLinePeaks_$(tele)_$(mjd)_(\\d+)_$(FIRST_CHIP)_object\\.h5\$")
    fns = sort([joinpath(d, fn) for fn in readdir(d) if occursin(r, fn)])
    fns = [f for f in fns if all(isfile(replace(f, "_$(FIRST_CHIP)_" => "_$(ch)_")) for ch in CHIP_LIST)]
    length(fns) > max_exp && (fns = fns[round.(Int, range(1, length(fns), length = max_exp))])
    groups = []
    for (gi, fn) in enumerate(fns)
        X, Y, C = read_sky_lines_readonly(fn)
        push!(groups, (id = gi, X = X, Y = Y, C = C))
    end
    return groups
end

"""Arc line groups from a line table `arclines_<tele>_<mjd>.h5` (datasets fiber, chip, group, pix,
wave: identified lamp lines of every arc sequence)."""
function check_arc_groups(arcdir, tele, mjd)
    isempty(arcdir) && return []
    p = joinpath(arcdir, "arclines_$(tele)_$(mjd).h5"); isfile(p) || return []
    T = h5open(f -> Dict(k => read(f[k]) for k in ("fiber", "chip", "group", "pix", "wave")), p, "r")
    groups = []
    for g in sort(unique(T["group"]))
        rows = findall(T["group"] .== g)
        nmax = maximum(count(==(f), T["fiber"][rows]) for f in 1:N_FIBERS)
        X = fill(NaN, nmax, N_FIBERS); Y = fill(NaN, nmax, N_FIBERS); C = zeros(Int, nmax, N_FIBERS)
        cnt = zeros(Int, N_FIBERS)
        for i in rows
            f = T["fiber"][i]; cnt[f] += 1
            X[cnt[f], f] = (T["pix"][i] - (N_XPIX ÷ 2)) / N_XPIX; Y[cnt[f], f] = T["wave"][i]; C[cnt[f], f] = T["chip"][i]
        end
        push!(groups, (id = g, X = X, Y = Y, C = C))
    end
    return groups
end

"""
Evaluate one night. Writes `workdir/stage4/<tag>/<tele>_<mjd>.h5` (resumable) holding per tracer,
model and chip the pooled held-out residuals' sum of squares and count, the per-(group, fiber)
zero-point change against scale1, the median shape pattern against the comb (if any), and the
per-fiber median b_free.
"""
function check_night(root, tele, mjd, cal, workdir, tag; arcdir = "", max_exp = 4, clip = 4.0)
    out = joinpath(workdir, "stage4", tag, "$(tele)_$(mjd).h5")
    isfile(out) && return :done
    blk = chip_transform_block(cal, tele, mjd)
    blk === nothing && return :no_block
    nl4 = [hcat(blk.offset[f, 1], blk.scale[f, 1], blk.offset[f, 2], blk.scale[f, 2]) for f in 1:N_FIBERS]
    comb = nothing
    cp = comb_path(root, tele, mjd)
    if isfile(cp)
        qa = read_fpi_qa(root, tele, mjd)
        (qa.has_qa && !qa.gate_pass) || (comb = h5open(f -> read(f["chipWaveSoln"]), cp, "r")[CHECK_PIX, :, :])
    end
    res = Dict{String, Any}()
    for (tr, groups, minc) in (("sky", check_sky_groups(root, tele, mjd; max_exp = max_exp), 3),
                               ("arc", check_arc_groups(arcdir, tele, mjd), 5))
        isempty(groups) && continue
        ss = Dict(md => zeros(3) for md in CHECK_MODELS); nn = Dict(md => zeros(Int, 3) for md in CHECK_MODELS)
        dzp = Dict(md => Vector{Vector{Float64}}() for md in CHECK_MODELS)
        pats = Dict(md => Vector{Matrix{Float64}}() for md in CHECK_MODELS)
        bfree = [Float64[] for _ in 1:N_FIBERS, _ in 1:2]
        nfit = 0
        for g in groups, f in 1:N_FIBERS
            x = g.X[:, f]; y = g.Y[:, f]; c = g.C[:, f]
            k = isfinite.(x) .& isfinite.(y) .& (c .> 0)
            count(k) < 12 && continue   # production sky fits have ~15-25 lines per fiber
            F = check_fit_group(x[k], c[k], y[k], vec(nl4[f]); clip = clip, min_per_chip = minc)
            F === nothing && continue
            nfit += 1
            G0 = F[:scale1].grid
            for md in CHECK_MODELS
                r = F[md]
                for cc in 1:3
                    v = r.loo[(r.chips .== cc) .& isfinite.(r.loo)]
                    ss[md][cc] += sum(abs2, v) * 1e6; nn[md][cc] += length(v)
                end
                push!(dzp[md], [mean(C_MPS .* (r.grid[:, cc] .- G0[:, cc]) ./ G0[:, cc]) for cc in 1:3])
                if comb !== nothing && all(isfinite, comb[:, f, :])
                    push!(pats[md], hcat([remove_zp_stretch(C_MPS .* (r.grid[:, cc] .- comb[:, f, cc]) ./ comb[:, f, cc]) for cc in 1:3]...))
                end
            end
            ## the chip scale the free fit lands on: b_ref + stretch / (dλ/dx' at the chip centre)
            p = F[:free_scale].p
            for (j, cc, col) in ((1, 1, 8), (2, 3, 9))
                b0 = nl4[f][2j]; a0 = nl4[f][2j - 1]
                push!(bfree[f, j], b0 + p[col] / quartic_slope(p, a0))
            end
        end
        nfit == 0 && continue
        res["$(tr)/n_fits"] = nfit
        res["$(tr)/n_groups"] = length(groups)
        for md in CHECK_MODELS
            res["$(tr)/$(md)/loo_ss"] = ss[md]; res["$(tr)/$(md)/loo_n"] = nn[md]
            res["$(tr)/$(md)/dzp"] = permutedims(hcat(dzp[md]...))      # n_fits x 3, m/s
            if !isempty(pats[md])
                P = cat(pats[md]...; dims = 3)
                res["$(tr)/$(md)/shape_pattern"] = [median(filter(isfinite, P[i, cc, :])) for i in 1:length(CHECK_PIX), cc in 1:3]
            end
        end
        res["$(tr)/b_free"] = [isempty(bfree[f, j]) ? NaN : median(bfree[f, j]) for f in 1:N_FIBERS, j in 1:2]
    end
    mkpath(dirname(out)); tmp = out * ".tmp$(getpid())"
    h5open(tmp, "w") do fo
        for (k, v) in res
            fo[k] = v
        end
        HDF5.attributes(fo)["calibration_version"] = cal.calibration_version
        HDF5.attributes(fo)["block_mjd_start"] = blk.mjd_start
        HDF5.attributes(fo)["has_comb"] = comb === nothing ? 0 : 1
    end
    mv(tmp, out; force = true)
    return isempty(res) ? :no_lines : :ok
end

# ------------------------------------------------------------------------------------------
# stage-4 summary
# ------------------------------------------------------------------------------------------

function bootstrap_ci(v; nboot = 2000, rng = MersenneTwister(1), stat = median)
    v = filter(isfinite, v); length(v) < 3 && return (NaN, NaN)
    s = [stat(v[rand(rng, 1:length(v), length(v))]) for _ in 1:nboot]
    return (quantile(s, 0.025), quantile(s, 0.975))
end

"""
Read the per-night stage-4 files of `tag` and return a long per-night table:
tele, mjd, category, tracer, model, chip, loo_rms_mA, dzp_ms (median over fibers/groups vs
scale1), shape_ms (RMS of the median pattern vs the comb), b_free_median.
"""
function check_collect(workdir, tag, nights::DataFrame)
    rows = NamedTuple[]
    for r in eachrow(nights)
        p = joinpath(workdir, "stage4", tag, "$(r.tele)_$(r.mjd).h5"); isfile(p) || continue
        h5open(p, "r") do f
            for tr in ("sky", "arc")
                haskey(f, tr) || continue
                bf = read(f["$(tr)/b_free"])
                for md in CHECK_MODELS
                    g = f["$(tr)/$(md)"]
                    ss = read(g["loo_ss"]); nn = read(g["loo_n"]); dz = read(g["dzp"])
                    pat = haskey(g, "shape_pattern") ? read(g["shape_pattern"]) : nothing
                    for cc in 1:3
                        push!(rows, (tele = r.tele, mjd = r.mjd, category = r.category, tracer = tr, model = String(md), chip = cc,
                            loo_rms_mA = nn[cc] > 0 ? sqrt(ss[cc] / nn[cc]) : NaN, n_lines = nn[cc],
                            dzp_ms = median(dz[:, cc]),
                            shape_ms = pat === nothing ? NaN : sqrt(mean(filter(isfinite, pat[:, cc]) .^ 2)),
                            b_free = cc == 2 ? NaN : medinit(filter(isfinite, bf[:, cc == 1 ? 1 : 2]), NaN)))
                    end
                end
            end
        end
    end
    return DataFrame(rows)
end

"""Paired summary of model `md` against `ref` per (tele, category, tracer, chip), with bootstrap
95 % intervals over nights and the fraction of nights improved."""
function check_summary(T::DataFrame; ref = "scale1")
    out = NamedTuple[]
    for k in groupby(T, [:tele, :category, :tracer, :chip])
        R = k[k.model .== ref, :]
        for md in unique(k.model)
            M = k[k.model .== md, :]
            J = innerjoin(select(R, :mjd, :loo_rms_mA => :loo0, :shape_ms => :shape0), select(M, :mjd, :loo_rms_mA, :shape_ms, :dzp_ms, :b_free), on = :mjd)
            d = J.loo_rms_mA .- J.loo0
            ci = bootstrap_ci(d)
            zci = bootstrap_ci(J.dzp_ms)
            push!(out, (tele = k.tele[1], category = k.category[1], tracer = k.tracer[1], chip = k.chip[1], model = md, ref = ref,
                n_nights = nrow(J), loo_rms_mA = median(J.loo_rms_mA), dloo_mA = median(d), dloo_lo = ci[1], dloo_hi = ci[2],
                frac_improved = mean(d .< 0), dzp_ms = median(J.dzp_ms), dzp_lo = zci[1], dzp_hi = zci[2], dzp_nmad = robust_sigma(J.dzp_ms),
                shape_ms = medinit(filter(isfinite, J.shape_ms), NaN), shape_ref_ms = medinit(filter(isfinite, J.shape0), NaN),
                b_free = medinit(filter(isfinite, J.b_free), NaN)))
        end
    end
    return DataFrame(out)
end

medinit(v, init) = isempty(v) ? init : median(v)

"""Flat 3-chip zero-point change per night (mean of the chips), summarised per group."""
function check_flat_rv(T::DataFrame)
    out = NamedTuple[]
    for k in groupby(T, [:tele, :category, :tracer, :model])
        n = combine(groupby(k, :mjd), :dzp_ms => (v -> mean(v)) => :flat)
        ci = bootstrap_ci(n.flat)
        push!(out, (tele = k.tele[1], category = k.category[1], tracer = k.tracer[1], model = k.model[1], n_nights = nrow(n),
            flat_rv_ms = median(n.flat), lo = ci[1], hi = ci[2], nmad = robust_sigma(n.flat)))
    end
    return DataFrame(out)
end
