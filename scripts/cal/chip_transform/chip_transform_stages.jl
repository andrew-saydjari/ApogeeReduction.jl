# Stage functions of the chip-transform calibration rebuild (see README.md in this folder).
# Included by build_chip_transform.jl on the main process and on every worker.
#
# Every number that ends up in the calibration file is computed here from the comb fits of a
# reduction; every tunable of the procedure is a command-line parameter of the driver, recorded
# in the file's provenance. Nothing learned is written into source code.

using ApogeeReduction
using ApogeeReduction: N_FIBERS, N_XPIX, CHIP_LIST, FIRST_CHIP, positional_poly_mat, transform_x_chips,
                       comb_exp_get_and_save_fpi_wavecal, fpi_light_check,
                       fpi_resid_stats, FPI_COMB_FIT_KWARGS, ChipTransformBlock, ChipTransformCal,
                       load_chip_transform_cal, chip_transform_block, write_chip_transform_cal,
                       robust_sigma, detect_blocks, robust_block_median, fiber_poly_fit, fiber_poly_eval
using HDF5, CSV, DataFrames, Statistics, Printf, Dates, Random, LinearAlgebra

const C_MPS = 299792458.0
const XIMPORT_GRID = ((1:N_XPIX) .- (N_XPIX ÷ 2)) ./ N_XPIX

expstr4(e) = lpad(string(e % 10000), 4, '0')
night_dir(root, mjd) = joinpath(root, "apred", string(mjd))
comb_path(root, tele, mjd) = joinpath(night_dir(root, mjd), "waveCalFPI_$(tele)_$(mjd)_arclamp.h5")
nightave_path(root, tele, mjd) = joinpath(root, "wavecal", "wavecalNightAve_$(tele)_$(mjd).h5")

"""All (tele, mjd) with a nightly wavelength file in <root>/wavecal (one directory listing)."""
function list_nights(root, tele; mjd_min = typemin(Int), mjd_max = typemax(Int))
    out = Int[]
    for fn in readdir(joinpath(root, "wavecal"))
        m = match(Regex("^wavecalNightAve_$(tele)_(\\d+)\\.h5\$"), fn)
        m === nothing && continue
        mjd = parse(Int, m.captures[1])
        mjd_min <= mjd <= mjd_max && push!(out, mjd)
    end
    return sort(out)
end

"""fpiPeaks files of chip FIRST_CHIP for the night (the night's FPI calibration exposures)."""
function fpi_peak_files(root, tele, mjd)
    d = night_dir(root, mjd); isdir(d) || return String[]
    r = Regex("^fpiPeaks_$(tele)_$(mjd)_\\d+_$(FIRST_CHIP)_arclamp\\.h5\$")
    return sort([joinpath(d, fn) for fn in readdir(d) if occursin(r, fn)])
end

function read_fpi_qa(root, tele, mjd)
    p = nightave_path(root, tele, mjd)
    isfile(p) || return (has_qa = false, gate_pass = false, failed_guard = "no_nightave")
    h5open(p, "r") do f
        haskey(f, "fpi_qa") || return (has_qa = false, gate_pass = false, failed_guard = "no_fpi_qa")
        g = f["fpi_qa"]
        gp = haskey(g, "gate_pass") ? Bool(read(g["gate_pass"])) : false
        fg = haskey(g, "failed_guard") ? string(read(g["failed_guard"])) : "unknown"
        return (has_qa = true, gate_pass = gp, failed_guard = fg)
    end
end

# ------------------------------------------------------------------------------------------
# stage 1: collect (or run) the comb fits
# ------------------------------------------------------------------------------------------

"""
Run the production FPI comb fit for a night whose comb file is missing, writing ONLY into
`workdir`: the fpiPeaks inputs are COPIED there (the package readers open their inputs
read-write, so they must never be pointed at a finished reduction), the fit's output lands next
to the copies, and the reduction root is only ever read. Same entry point, settings and acceptance gate as
`fpi_medwavecal_skyline_dither`, seeded from the night's stored sky solution.
"""
function run_missing_comb(root, tele, mjd, fpi_files, workdir)
    wd = joinpath(workdir, "stage1_combfits", string(mjd)); mkpath(wd)
    links = String[]
    for fn in fpi_files, chip in CHIP_LIST
        src = replace(fn, "_$(FIRST_CHIP)_" => "_$(chip)_")
        isfile(src) || continue
        dst = joinpath(wd, basename(src))
        isfile(dst) || cp(src, dst)
        chip == FIRST_CHIP && push!(links, dst)
    end
    isempty(links) && return nothing, "no_fpiPeaks"
    light = fpi_light_check(links)
    light.pass || return nothing, "fpi_light"
    p = nightave_path(root, tele, mjd)
    isfile(p) || return nothing, "no_nightave"
    linP, nlP = h5open(p, "r") do f
        read(f["sky/nightAve_linParams"]), read(f["sky/nightAve_nlParams"])
    end
    linP = reshape(linP, size(linP, 1), size(linP, 2)); nlP = reshape(nlP, size(nlP, 1), size(nlP, 2))
    out = comb_exp_get_and_save_fpi_wavecal(links, linP, nlP; FPI_COMB_FIT_KWARGS..., verbose = false,
        checkpoint_mode = "commit_same")
    out === nothing && return nothing, "fpi_solution"
    r = fpi_resid_stats(out[1])
    r.n_used == 0 && return nothing, "fpi_no_peaks_used"
    r.pass || return nothing, "fpi_resid_rms"
    return out[1], "ran"
end

"""
Stage-1 record of one night: per-fiber comb chip parameters (a_R, b_R, a_B, b_B), the number of
comb peaks used per fiber and chip, and the per-fiber residual RMS. Cached per night in
`workdir/stage1/<tele>/<mjd>.h5` (resumable).
"""
function stage1_night(root, tele, mjd, workdir; run_missing = false)
    cache = joinpath(workdir, "stage1", tele, "$(mjd).h5")
    isfile(cache) && return cache
    status = "used"; src = comb_path(root, tele, mjd)
    fpi_files = String[]
    if isfile(src)
        qa = read_fpi_qa(root, tele, mjd)
        if qa.has_qa && !qa.gate_pass
            status = "rejected_gate_" * qa.failed_guard
        end
    else
        fpi_files = fpi_peak_files(root, tele, mjd)
        qa = isempty(fpi_files) ? nothing : read_fpi_qa(root, tele, mjd)
        if isempty(fpi_files)
            status = "no_fpi_exposures"
        elseif qa.has_qa && !qa.gate_pass
            # the reduction ran the comb fit (or its light check) and the gate refused it
            status = "rejected_gate_" * qa.failed_guard
        elseif run_missing
            p, st = run_missing_comb(root, tele, mjd, fpi_files, workdir)
            if p === nothing
                status = "rejected_rerun_" * st
            else
                src = p; status = "used_rerun"
            end
        else
            status = "missing_comb"
        end
    end
    nl = fill(NaN, N_FIBERS, 4); npk = zeros(Int, N_FIBERS, 3); rrms = fill(NaN, N_FIBERS)
    if (startswith(status, "used") || startswith(status, "rejected_gate")) && isfile(src)
        h5open(src, "r") do f
            x = read(f["nlParams"])
            size(x, 2) == 4 && (nl .= x)
            if haskey(f, "resid_vec") && haskey(f, "resid_used_in_fit") && haskey(f, "resid_chipInts")
                R = read(f["resid_vec"]); U = read(f["resid_used_in_fit"]) .!= 0; C = read(f["resid_chipInts"])
                for fb in 1:N_FIBERS
                    u = U[:, fb] .& isfinite.(R[fb, :])
                    for c in 1:3
                        npk[fb, c] = count(u .& (C[:, fb] .== c))
                    end
                    any(u) && (rrms[fb] = sqrt(mean(R[fb, u] .^ 2)))
                end
            end
        end
    end
    mkpath(dirname(cache)); tmp = cache * ".tmp$(getpid())"
    h5open(tmp, "w") do f
        f["nlParams"] = nl; f["n_peaks"] = npk; f["resid_rms"] = rrms
        HDF5.attributes(f)["status"] = status; HDF5.attributes(f)["source"] = (status == "no_fpi_exposures" ? "" : src)
        HDF5.attributes(f)["n_fpi_exposure_files"] = length(fpi_files)
    end
    mv(tmp, cache; force = true)
    return cache
end

"""Gather the per-night stage-1 caches of one telescope into one table + arrays."""
function stage1_gather(workdir, tele, mjds)
    n = length(mjds)
    NL = fill(NaN, N_FIBERS, 4, n); NP = zeros(Int, N_FIBERS, 3, n); RR = fill(NaN, N_FIBERS, n)
    status = fill("", n); src = fill("", n)
    for (i, mjd) in enumerate(mjds)
        h5open(joinpath(workdir, "stage1", tele, "$(mjd).h5"), "r") do f
            NL[:, :, i] .= read(f["nlParams"]); NP[:, :, i] .= read(f["n_peaks"]); RR[:, i] .= read(f["resid_rms"])
            status[i] = read(HDF5.attributes(f)["status"]); src[i] = read(HDF5.attributes(f)["source"])
        end
    end
    return (; mjds, NL, NP, RR, status, src)
end

# ------------------------------------------------------------------------------------------
# stage 2: derive the chip transform per fiber, in blocks of time
# ------------------------------------------------------------------------------------------

"""
Per-night fiber rejection (returns a fibers x nights mask of usable fibers, and appends to
`rejected`): non-finite parameters; fewer than `min_peaks` comb peaks used on any chip; a scale
more than `nsigma` robust sigma off that night's quadratic-in-fiber trend. A night is rejected
if fewer than `min_fibers` fibers survive.
"""
function stage2_fiber_mask(S, tele; min_peaks, nsigma, min_fibers, rejected)
    n = length(S.mjds); ok = falses(N_FIBERS, n); night_ok = falses(n)
    fib = collect(1:N_FIBERS)
    for i in 1:n
        if !startswith(S.status[i], "used")
            S.status[i] == "no_fpi_exposures" ||
                push!(rejected, (tele = tele, mjd = S.mjds[i], fiber = 0, reason = S.status[i], value = NaN))
            continue
        end
        m = vec(all(isfinite, S.NL[:, :, i], dims = 2))
        lowpk = vec(any(S.NP[:, :, i] .< min_peaks, dims = 2)) .& m
        for f in findall(lowpk)
            push!(rejected, (tele = tele, mjd = S.mjds[i], fiber = f, reason = "few_comb_peaks", value = minimum(S.NP[f, :, i])))
        end
        m .&= .!lowpk
        for (j, nm) in ((2, "scale_R"), (4, "scale_B"))
            count(m) < min_fibers && break
            c, _ = fiber_poly_fit(fib[m], S.NL[m, j, i], 2; nsigma = nsigma)
            r = S.NL[:, j, i] .- fiber_poly_eval(c, fib)
            s = robust_sigma(r[m])
            bad = m .& (abs.(r) .> nsigma * s)
            for f in findall(bad)
                push!(rejected, (tele = tele, mjd = S.mjds[i], fiber = f, reason = "outlier_$(nm)", value = r[f] / s))
            end
            m .&= .!bad
        end
        if count(m) < min_fibers
            push!(rejected, (tele = tele, mjd = S.mjds[i], fiber = 0, reason = "too_few_fibers", value = count(m)))
            continue
        end
        ok[:, i] .= m; night_ok[i] = true
    end
    return ok, night_ok
end

"""
Derive the blocks and per-fiber transforms of one telescope. Returns a NamedTuple with the
blocks (`ChipTransformBlock`s whose MJD ranges are the DATA ranges; stage 3 sets coverage), the
night series, the chosen fiber model and its cross-validation scores.
"""
function stage2_derive(S, tele; P, rejected)
    ok, night_ok = stage2_fiber_mask(S, tele; min_peaks = P["fiber_min_peaks"], nsigma = P["clip_nsigma"],
        min_fibers = P["night_min_fibers"], rejected = rejected)
    n = length(S.mjds)
    ## night series: fiber medians of (b_R, b_B, a_R, a_B)
    ser = fill(NaN, n, 4)
    for i in 1:n, (q, j) in enumerate((2, 4, 1, 3))
        night_ok[i] && (ser[i, q] = median(S.NL[ok[:, i], j, i]))
    end
    used = findall(night_ok)
    isempty(used) && error("$(tele): no usable comb nights")
    ## block boundaries (indices into `used`)
    mode = P["block_mode"]
    starts_scale = detect_blocks(S.mjds[used], ser[used, 1:2]; min_nights = P["min_block_nights"],
        nsigma = P["block_nsigma"], min_step = [P["block_min_step"], P["block_min_step"]])
    ## offsets are NOT used for the blocks; their change points are reported for QA (hardware events)
    starts_offset = detect_blocks(S.mjds[used], ser[used, 3:4]; min_nights = P["min_block_nights"],
        nsigma = P["block_nsigma"], min_step = [P["offset_event_min_step"], P["offset_event_min_step"]])
    supplied = Int[]
    if mode in ("supplied", "supplied+detect")
        for l in eachline(P["boundaries"])
            s = split(strip(l)); (isempty(s) || startswith(s[1], "#")) && continue
            s[1] == tele && push!(supplied, parse(Int, s[2]))
        end
    end
    bmjd = if mode == "none"
        Int[]
    elseif mode == "detect"
        S.mjds[used][starts_scale]
    elseif mode == "supplied"
        supplied
    elseif mode == "supplied+detect"
        sort(unique(vcat(supplied, S.mjds[used][starts_scale])))
    else
        error("unknown block mode $(mode)")
    end
    edges = vcat(typemin(Int), sort(bmjd), typemax(Int))
    blocks = ChipTransformBlock[]; cv = NamedTuple[]
    fib = collect(1:N_FIBERS)
    for k in 1:(length(edges) - 1)
        inb = [i for i in used if edges[k] <= S.mjds[i] < edges[k + 1]]
        isempty(inb) && continue
        ## reject nights whose median scale is an outlier within the block
        for q in 1:2
            v = ser[inb, q]; m = median(v); s = robust_sigma(v)
            (isfinite(s) && s > 0) || continue
            bad = abs.(v .- m) .> P["clip_nsigma"] * s
            for i in inb[bad]
                push!(rejected, (tele = tele, mjd = S.mjds[i], fiber = 0, reason = "night_outlier_in_block_q$(q)", value = (ser[i, q] - m) / s))
            end
            inb = inb[.!bad]
        end
        length(inb) < P["block_min_used_nights"] && (@warn "$(tele) block $(k) has only $(length(inb)) nights"; length(inb) == 0 && continue)
        V = Dict(j => [ok[f, i] ? S.NL[f, j, i] : NaN for f in 1:N_FIBERS, i in inb] for j in 1:4)
        med = Dict{Int, Vector{Float64}}(); se = Dict{Int, Vector{Float64}}()
        for j in 1:4
            med[j], se[j], _ = robust_block_median(V[j]; nsigma = P["clip_nsigma"])
        end
        ## fiber model, chosen by cross-validation over nights (odd vs even nights of the block)
        cvrow = Dict{Symbol, Any}(:tele => tele, :block => k, :n_nights => length(inb))
        model = P["fiber_model"]
        for j in (2, 4)
            if length(inb) >= 4
                ho = isodd.(1:length(inb))
                mA, _, _ = robust_block_median(V[j][:, ho]; nsigma = P["clip_nsigma"])
                mB, _, _ = robust_block_median(V[j][:, .!ho]; nsigma = P["clip_nsigma"])
                g = isfinite.(mA) .& isfinite.(mB)
                cA, _ = fiber_poly_fit(fib[g], mA[g], P["fiber_poly_degree"]; nsigma = P["clip_nsigma"])
                cB, _ = fiber_poly_fit(fib[g], mB[g], P["fiber_poly_degree"]; nsigma = P["clip_nsigma"])
                pA = fiber_poly_eval(cA, fib[g]); pB = fiber_poly_eval(cB, fib[g])
                ## Split-half estimate of each model's error against the (unknown) true per-fiber
                ## scale, for the FULL block:
                ##   half-sample noise variance  s2h = mean((mA - mB)^2) / 2
                ##   per-fiber model (full block): s2h / 2
                ##   smooth model: mean((mB - pA)^2) - s2h  (what the curve misses; its own fit noise
                ##                 over ~300 fibers is negligible)
                s2h = mean((mA[g] .- mB[g]) .^ 2) / 2
                e_pf = sqrt(s2h / 2)
                e_po = sqrt(max(mean(vcat(mB[g] .- pA, mA[g] .- pB) .^ 2) - s2h, 0.0))
                nm = j == 2 ? "R" : "B"
                cvrow[Symbol("err_perfiber_$(nm)")] = e_pf
                cvrow[Symbol("err_poly_$(nm)")] = e_po
                cvrow[Symbol("halves_corr_of_fiber_deviation_$(nm)")] = cor(mA[g] .- pA, mB[g] .- pB)
            end
        end
        if model == "auto"
            ## the per-fiber model is used only if it has the smaller estimated error on both chips
            pf = all(get(cvrow, Symbol("err_perfiber_$(nm)"), Inf) < get(cvrow, Symbol("err_poly_$(nm)"), -Inf) for nm in ("R", "B"))
            model = pf ? "per_fiber" : "poly"
        end
        cvrow[:fiber_model] = model
        push!(cv, NamedTuple(cvrow))
        scale = zeros(N_FIBERS, 2); serr = zeros(N_FIBERS, 2); offs = zeros(N_FIBERS, 2); filled = falses(N_FIBERS)
        for (col, j) in ((1, 2), (2, 4))
            g = isfinite.(med[j])
            c, _ = fiber_poly_fit(fib[g], med[j][g], P["fiber_poly_degree"]; nsigma = P["clip_nsigma"])
            smooth = fiber_poly_eval(c, fib)
            scale[:, col] .= model == "per_fiber" ? ifelse.(g, med[j], smooth) : smooth
            serr[:, col] .= ifelse.(isfinite.(se[j]), se[j], NaN)
            filled .|= .!g
            go = isfinite.(med[j - 1])
            co, _ = fiber_poly_fit(fib[go], med[j - 1][go], P["fiber_poly_degree"]; nsigma = P["clip_nsigma"])
            offs[:, col] .= ifelse.(go, med[j - 1], fiber_poly_eval(co, fib))
        end
        dm = S.mjds[inb]
        push!(blocks, ChipTransformBlock(minimum(dm), maximum(dm), minimum(dm), maximum(dm), length(inb),
            scale, serr, offs, collect(filled)))
    end
    return (; blocks, ser, night_ok, ok, used, bmjd, starts_offset_mjd = S.mjds[used][starts_offset], cv, model = [r.fiber_model for r in cv])
end

"""Coverage: block k serves [first night of k (or `cover_from` for k = 1), first night of k+1 - 1]
(the last block up to `cover_to` or its last night)."""
function set_coverage(blocks::Vector{ChipTransformBlock}; cover_from = nothing, cover_to = nothing)
    out = ChipTransformBlock[]
    for (k, b) in enumerate(blocks)
        s = k == 1 ? (isnothing(cover_from) ? b.data_mjd_min : min(cover_from, b.data_mjd_min)) : b.data_mjd_min
        e = k == length(blocks) ? (isnothing(cover_to) ? b.data_mjd_max : max(cover_to, b.data_mjd_max)) :
            blocks[k + 1].data_mjd_min - 1
        push!(out, ChipTransformBlock(s, e, b.data_mjd_min, b.data_mjd_max, b.n_nights, b.scale, b.scale_err, b.offset, b.fiber_filled))
    end
    return out
end
