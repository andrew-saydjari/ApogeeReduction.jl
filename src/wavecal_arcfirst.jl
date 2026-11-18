# In-order ("arc-first") nightly wavelength chain, `pipeline_2d_1d.jl --wavecal_order arc_first`.
#
# AKS 2026-11-18: "I just don't think it is necessary to do sky, fpi, arc, then refit sky and FPI
# again. It seems like we should plan for the final state to just do them in order". Each stage
# runs ONCE, in this order, and reads only the products of earlier stages:
#
#   (a) peaks          sky / FPI / arc peak finding (unchanged pipeline functions)
#   (b) arc            refit_night_arc               -> <wroot>/arc/arclampWaveCalNight_*.h5
#   (c) comb           arcfirst_night_comb           -> <wroot>/comb/waveCalFPIArc_*.h5
#                      (FPI gate G1 on the fpiPeaks; provisional production comb fit -> measured comb
#                      exposure positions; refit_comb_anchor_stage (#420 anchor rule + arc seed))
#   (d,e) exposures    refit_night_fpifib            -> <wroot>/fpifib/fpiGuidePos_*.h5
#                      refit_night_assemble(inorder) -> <wroot>/exp/waveCalRefit_*.h5 and, per exposure,
#                                                       apred/<mjd>/waveCalExp_<tele>_<mjd>_<exp>_<kind>.h5
#   (f) resampling     reinterp_spectra(...; wavecal_order = "arc_first") reads waveCalExp (chipWaveSoln)
#
# A night without a usable arc fit runs the legacy sky(/FPI) chain instead (pipeline), and its
# exposures carry WAVE_NOT_ARCFIRST. The archive refit (scripts/refit/) calls the same stage
# functions with inputs gathered from a previous run; the in-order inputs are gathered here.

const ARCFIRST_ORDERS = ("legacy", "arc_first")
const ARCFIRST_NOARC_FALLBACKS = ("legacy",)
const ARCFIRST_FIBER_FILLS = ("fill_only", "none")

"Validate a `--wavecal_order` value (error naming the allowed values)."
function arcfirst_check_order(s::AbstractString)
    s in ARCFIRST_ORDERS || throw(ArgumentError("--wavecal_order \"$(s)\" not in $(ARCFIRST_ORDERS)"))
    return String(s)
end

"""Stage table of the in-order chain: (stage, products read, products written). `test/wavecal_arcfirst.jl`
checks that no stage reads a product of the same or a later stage."""
const ARCFIRST_STAGES = [
    (stage = "peaks", reads = ["ar1Dcal"], writes = ["skyLinePeaks", "fpiPeaks", "arclampPeaks"]),
    (stage = "arc", reads = ["ar1Dcal", "almanac", "arc_linelist", "chip_transform"], writes = ["arclampWaveCalNight"]),
    (stage = "comb", reads = ["fpiPeaks", "ar1Dcal", "almanac", "arclampWaveCalNight"], writes = ["waveCalFPIArc"]),
    (stage = "fpifib", reads = ["ar1Dcal", "waveCalFPIArc"], writes = ["fpiGuidePos"]),
    (stage = "assemble", reads = ["ar1Dcal", "skyLinePeaks", "arclampWaveCalNight", "waveCalFPIArc",
        "fpiGuidePos", "shape_template", "arc_linelist"], writes = ["waveCalRefit", "waveCalExp"]),
    (stage = "resample", reads = ["ar1Dcal", "waveCalExp"], writes = ["ar1Duni"]),
]

"Root of the in-order night products of a pipeline output directory."
arcfirst_root(outdir::AbstractString) = joinpath(outdir, "wavecal_arcfirst")

"Per-exposure product read by the resampling (design 2.1 `waveCalExp`)."
arcfirst_exp_name(apred_dir, tele, mjd, e, kind) =
    joinpath(apred_dir, string(mjd), "waveCalExp_$(tele)_$(mjd)_$(_arc_expstr(e))_$(kind).h5")

"""Night-level route of the in-order chain from the arc stage status: `:arcfirst` when the arc
fit is usable, else `:legacy` (the legacy sky/FPI chain runs for the night, WAVE_NOT_ARCFIRST)."""
arcfirst_night_route(arc_usable::Bool, noarc_fallback::AbstractString = "legacy") =
    arc_usable ? :arcfirst : (noarc_fallback == "legacy" ? :legacy :
                              throw(ArgumentError("noarc_fallback \"$(noarc_fallback)\" not implemented")))

"""Comb decision of the in-order chain (status string written to `status/comb_*.txt`):
no FPI exposures -> "nocomb"; FPI gate G1 (light in the FPI fibers) fails -> "nocomb_gated";
provisional comb fit fails or gate G2 (residual RMS) fails -> "nocomb_gated"; no arc sequence
matches a comb exposure by measured position -> "nocomb_unmatched"; else "ok_arc". Without a
comb the night's FPS exposures take #420's arc route flagged WAVE_FPS_NO_COMB (+ WAVE_FPI_GATED_OFF
when FPI frames existed)."""
function arcfirst_comb_decision(n_fpi::Integer, g1_pass::Bool, prov_ok::Bool, g2_pass::Bool, anchored::Bool)
    n_fpi == 0 && return "nocomb"
    g1_pass || return "nocomb_gated"
    (prov_ok && g2_pass) || return "nocomb_gated"
    anchored || return "nocomb_unmatched"
    return "ok_arc"
end

"""Legacy-chain `wavecal_type` of an exposure resampled in an arc_first run without an arc-first
solution -> (new-scheme name, bits incl. WAVE_NOT_ARCFIRST)."""
function arcfirst_legacy_type(wt::AbstractString)
    nm, b = refit_pass1_type(wt)
    return nm, b | WAVE_NOT_ARCFIRST
end

"""
    arcfirst_comb_refs_provisional(fpi_fnames, A, mdir, tele, mjd; kw...) -> (combs, info)

In-order replacement of the previous run's per-exposure comb solutions (`refit_comb_refs_pass1`):
the PRODUCTION comb fit (`comb_fpi_fit_core` with `FPI_COMB_FIT_KWARGS`, hard-coded initial cavity,
cavity refit by the production loop) on tonight's fpiPeaks, with integer orders seeded from the arc
fit at its reference sequence instead of the sky night average. Its per-exposure grids (frame +
the fit's own per-exposure dither transform) give each comb exposure's MEASURED position in the same
absolute frame the previous run's comb solutions had, so the #420 pairing (arc sequence state vs
comb exposure, < `pos_match_px`) is unchanged. `info` carries the provisional cavity and the gate-G2
residual statistics. Returns `(nothing, info)` if the fit fails.
"""
function arcfirst_comb_refs_provisional(fpi_fnames, A::RefitArc, mdir, tele, mjd)
    st = refit_arc_state(A, NaN, A.clus_u[1]; mode = :nearest)
    lin0, nl0, _ = refit_arc_lin_nl(A, st)
    R = comb_fpi_fit_core(fpi_fnames, lin0, nl0; FPI_COMB_FIT_KWARGS..., verbose = false)
    R === nothing && return nothing, (; ok = false, rms = NaN, n_used = 0, L = NaN, eps = NaN)
    r = filter(isfinite, R.resid_vec'[R.good_peaks])
    rms = isempty(r) ? NaN : sqrt(sum(abs2, r) / length(r))
    info = (; ok = true, rms, n_used = length(r), L = R.cavity_size, eps = R.m_offset,
        g2_pass = !isempty(r) && isfinite(rms) && rms <= FPI_RESID_RMS_MAX)
    px = 880:1170
    ximport = (px .- (N_XPIX ÷ 2)) ./ N_XPIX
    combs = map(enumerate(fpi_fnames)) do (k, fn)
        e = parse(Int, split(basename(fn), "_")[4])
        t = arclamp_mjd_mid(ar1dcal_name(mdir, tele, mjd, e, FIRST_CHIP, "arclamp"))
        W = fill(NaN, length(px), length(ARCGATE_FIBERS))
        for (jf, f) in enumerate(ARCGATE_FIBERS)
            cpp = copy(R.chipPolyParams0[f, :, :])
            params2ChipPolyParams!(cpp, R.nlParams[f, :], R.cporder)
            xt = transform_x_chips(ximport, cpp[2, :])
            xt = transform_x_chips(xt, R.ditherPolyParams[k, f, :])
            W[:, jf] .= positional_poly_mat(xt, porder = R.wporder) * R.linParams[f, :]
        end
        (; ref = ArcObjRef(e, t, "fpi", first(px), ARCGATE_FIBERS, W), t)
    end
    return combs, info
end

"""
    arcfirst_night_comb(tele, mjd; apred_dir, outroot, cfg, fpi_fnames, sci_exps, fpi_fiberIndxs)

Stage (c) in order: FPI gate G1 on tonight's fpiPeaks (`fpi_light_check`), provisional production
comb fit (measured comb-exposure positions, gate G2 on its residual RMS), then the shared #420
anchor stage (`refit_comb_anchor_stage`) with science exposure TIMES only (`anchor_samepos` must
not be "pass1"). Returns the status string (also written to `status/comb_*.txt`).
"""
function arcfirst_night_comb(tele, mjd; apred_dir, outroot, cfg::RefitConfig, fpi_fnames,
        sci_exps, fpi_fiberIndxs = Int[], exclude_seqs = Int[], cache = Dict{Symbol, Any}())
    out = refit_path(outroot, "comb", tele, mjd)
    isfile(out) && return "done"
    cfg.anchor_samepos == "pass1" &&
        throw(ArgumentError("anchor_samepos = \"pass1\" reads a previous run's per-exposure solutions; not allowed in order"))
    mdir = joinpath(apred_dir, string(mjd))
    fn = filter(isfile, collect(String, fpi_fnames))
    isempty(fn) && return refit_write_status(outroot, "comb", tele, mjd,
        arcfirst_comb_decision(0, false, false, false, false), "no FPI exposures")
    light = fpi_light_check(fn; fpi_fiberIndxs = fpi_fiberIndxs)
    light.pass || return refit_write_status(outroot, "comb", tele, mjd,
        arcfirst_comb_decision(length(fn), false, false, false, false),
        "G1 fpi_light lit=$(round(light.lit, sigdigits = 3)) < $(FPI_LIT_FRACTION_MIN)")
    arcf = refit_path(outroot, "arc", tele, mjd)
    A = isfile(arcf) ? load_refit_arc(arcf) : nothing
    (A === nothing || !refit_arc_usable(A)) &&
        return refit_write_status(outroot, "comb", tele, mjd, "nocomb_noarc", "no usable arc fit")
    # the provisional fit does not depend on the anchor: computed once per night (re-anchoring reuses it)
    combs, info = get!(() -> arcfirst_comb_refs_provisional(fn, A, mdir, tele, mjd), cache, :provisional)
    if combs === nothing || !info.g2_pass
        return refit_write_status(outroot, "comb", tele, mjd,
            arcfirst_comb_decision(length(fn), true, combs !== nothing, false, false),
            combs === nothing ? "provisional comb fit failed" :
            "G2 resid rms=$(round(info.rms, sigdigits = 3)) A (n=$(info.n_used)) > $(FPI_RESID_RMS_MAX)")
    end
    scis = [(; ref = nothing, t = arclamp_mjd_mid(ar1dcal_name(mdir, tele, mjd, e, FIRST_CHIP, "object")))
            for e in sci_exps]
    prov = Dict("chip_transform_fingerprint" => refit_fingerprint(cfg.chip_transform_file),
        "n_science_refs" => string(length(scis)), "inorder" => "true",
        "comb_pos_ref" => "provisional production comb fit, integers from the arc reference sequence",
        "provisional_cavity_size" => string(info.L), "provisional_m0_offset" => string(info.eps),
        "g1_lit_fraction" => string(light.lit), "g2_resid_rms_A" => string(info.rms))
    prov["excluded_anchor_seqs"] = join(exclude_seqs, ",")
    st, reason = refit_comb_anchor_stage(tele, mjd; outroot, fnames = fn, combs, scis, cfg,
        pass1_comb = nothing, provenance = prov, exclude_seqs)
    st === nothing || return st
    return refit_write_status(outroot, "comb", tele, mjd,
        arcfirst_comb_decision(length(fn), true, true, true, false), reason)
end

"""Write one exposure's in-order solution (`waveCalExp`, read by `reinterp_spectra` in arc_first
mode): `chipWaveSoln` (N_XPIX x N_FIBERS x 3, the format of the legacy per-exposure products),
`wavecal_type`, `bitmsk_wavecal`, `wavecal_xcheck_resid` (m/s). Fibers still non-finite (no route
gave them a solution) copy the nearest finite fiber of the same chip and set WAVE_FIBER_FILLED.
Returns the bits it added."""
function arcfirst_write_exp(fname, G::AbstractArray, wtype, bits, xres; tele, mjd, expnum, kind, provenance = Dict{String, Any}())
    add = 0
    G = copy(G)
    for c in 1:N_CHIPS
        good = [all(isfinite, view(G, :, f, c)) for f in 1:N_FIBERS]
        all(good) && continue
        any(good) || continue
        gi = findall(good)
        for f in findall(.!good)
            G[:, f, c] .= G[:, gi[argmin(abs.(gi .- f))], c]
        end
        add |= WAVE_FIBER_FILLED
    end
    meta = Dict{String, Any}("product" => "waveCalExp", "tele" => String(tele), "mjd" => Int(mjd),
        "expnum" => Int(expnum), "kind" => String(kind), "wavecal_order" => "arc_first",
        "bit_names" => join(WAVE_BIT_NAMES, ","))
    for (k, v) in provenance
        meta[k] = v
    end
    mkpath(dirname(fname))
    tmp = fname * ".tmp$(getpid())"
    safe_jldsave(tmp, meta; chipWaveSoln = Float64.(G), wavecal_type = String(wtype),
        bitmsk_wavecal = Int64(bits | add), wavecal_xcheck_resid = Float32(xres))
    mv(tmp, fname; force = true)
    return add
end

"""
    arcfirst_check_night(outroot, tele, mjd; sky_check_ms = 1000.0) -> NamedTuple

END-OF-NIGHT verification of the comb anchor, using later-stage information ONLY as a check: per
chip, the median over the night's FPI-fiber exposures of the pooled sky-line offset against the
final solution (`sky_minus_model_ms`, before any sky tweak; the FPI route takes no sky level), and
the comb mode-slip count. A wrong anchor (arc block and paired comb exposure at different
mechanism positions) shifts every FPI exposure by the position difference, i.e. a dither step
(0.5 px = 2.1 / 2.6 / 3.2 km/s on R / G / B), against an arc-vs-sky ensemble offset of 40-140 m/s
per chip (design A5/D7). Fails if any chip median exceeds `sky_check_ms` or the comb slipped.
Unverifiable (no FPI-fiber exposure with sky) counts as pass, `verifiable = false`.
Returns (pass, verifiable, metric_ms (max |chip median|), chip_ms, n_exp, n_slip).
"""
function arcfirst_check_night(outroot, tele, mjd; sky_check_ms::Real = 1000.0)
    asm = refit_path(outroot, "assemble", tele, mjd)
    combf = refit_path(outroot, "comb", tele, mjd)
    isfile(asm) || return (; pass = true, verifiable = false, metric_ms = NaN, chip_ms = fill(NaN, 3), n_exp = 0, n_slip = 0)
    wt, v = h5open(f -> (String.(read(f["wavecal_type"])), read(f["sky_minus_model_ms"])), asm, "r")
    nslip = isfile(combf) ? h5open(f -> Int(read(f["n_mode_slip"])), combf, "r") : 0
    ks = [k for k in eachindex(wt) if endswith(wt[k], "_fpifib") && any(isfinite, view(v, :, k))]
    chip = [let x = filter(isfinite, v[c, ks]); isempty(x) ? NaN : median(x) end for c in 1:3]
    fin = filter(isfinite, chip)
    metric = isempty(fin) ? NaN : maximum(abs.(fin))
    verifiable = !isempty(fin)
    pass = nslip == 0 && (!verifiable || metric <= sky_check_ms)
    return (; pass, verifiable, metric_ms = metric, chip_ms = chip, n_exp = length(ks), n_slip = nslip)
end

"Write the night's end-of-night check record (status/check_<tele>_<mjd>.txt, one line per attempt)."
function arcfirst_write_check(outroot, tele, mjd, lines::Vector{String})
    fn = joinpath(outroot, "status", "check_$(tele)_$(mjd).txt")
    mkpath(dirname(fn))
    open(fn, "w") do io
        foreach(l -> println(io, l), lines)
    end
end

"Remove a night's in-order products (all stages) and its waveCalExp files."
function arcfirst_clear_night(outroot, apred_dir, tele, mjd, exps)
    for st in REFIT_STAGES
        rm(refit_path(outroot, st, tele, mjd); force = true)
        rm(refit_status_path(outroot, st, tele, mjd); force = true)
    end
    for (e, kind) in exps
        rm(arcfirst_exp_name(apred_dir, tele, mjd, e, kind); force = true)
    end
end

"""
    arcfirst_night(tele, mjd; outdir, almanac_file, exps, fpi_fnames, cfg, checkpoint_mode = "commit_same",
        fiber_fill = "fill_only", noarc_fallback = "legacy", ll = nothing) -> NamedTuple

Stages (b)-(e) of one night in order (see the file header). `exps`: this run's exposures of the
night as (expnum, kind) (from the ar1Dcal file names); `fpi_fnames`: tonight's FPI fpiPeaks files
(R-chip names). Returns `(route, arc, comb, fpifib, assemble, n_exp_written)`; `route == :legacy`
means the caller must run the legacy chain for the night (its exposures get WAVE_NOT_ARCFIRST at
resampling). Checkpointing: the night is redone from the arc stage unless its assembly product
passes `check_file(...; mode = checkpoint_mode)` (clobber: always redone).
"""
function arcfirst_night(tele, mjd; outdir, almanac_file, exps, fpi_fnames, cfg::RefitConfig,
        checkpoint_mode = "commit_same", fiber_fill = "fill_only", noarc_fallback = "legacy", ll = nothing,
        max_reanchor::Integer = 2, sky_check_ms::Real = 1000.0)
    fiber_fill in ARCFIRST_FIBER_FILLS || throw(ArgumentError("fiber_fill \"$(fiber_fill)\" not in $(ARCFIRST_FIBER_FILLS)"))
    refit_check_config(cfg)
    apred_dir = joinpath(outdir, "apred")
    outroot = arcfirst_root(outdir)
    exps = sort(unique(collect(Tuple{Int, String}, exps)))
    asm = refit_path(outroot, "assemble", tele, mjd)
    have_all = isfile(asm) && all(x -> isfile(arcfirst_exp_name(apred_dir, tele, mjd, x...)), exps)
    if have_all && check_file(asm, mode = checkpoint_mode)
        return (; route = :arcfirst, arc = "done", comb = "done", fpifib = "done", assemble = "done", n_exp_written = 0,
            n_attempts = 0)
    end
    arcfirst_clear_night(outroot, apred_dir, tele, mjd, exps)
    s_arc = refit_night_arc(tele, mjd; apred_dir, almanac_file, outroot, cfg)
    arcf = refit_path(outroot, "arc", tele, mjd)
    usable = isfile(arcf) && refit_arc_usable(load_refit_arc(arcf))
    route = arcfirst_night_route(usable, noarc_fallback)
    route == :legacy && return (; route, arc = s_arc, comb = "", fpifib = "", assemble = "", n_exp_written = 0,
        n_attempts = 0)
    fpi_ids = isempty(almanac_file) ? Int[] : get_fpi_fiberIDs_from_almanac(almanac_file, tele, mjd)
    fpi_idx = try
        filter(i -> 1 <= i <= N_FIBERS, fiberID2fiberIndx.(fpi_ids))
    catch
        Int[]
    end
    sci = [e for (e, k) in exps if k == "object"]
    A = load_refit_arc(arcf)
    cache = Dict{Symbol, Any}()
    exclude = Int[]                  # used-sequence indices rejected by the end-of-night check
    tried = Tuple{Int, Float64}[]    # (anchor j_seq, check metric) per attempt
    checklog = String[]
    s_comb = s_fpifib = s_asm = ""
    nw = Ref(0)
    attempt = 0
    best_j = 0
    while true
        # stages (c)-(e), a pure function of (arc fit, tonight's peaks, excluded anchors)
        extra = attempt == 0 ? 0 : WAVE_REANCHORED
        final_fail = attempt > max_reanchor
        final_fail && (extra |= WAVE_ANCHOR_CHECK_FAILED)
        s_comb = arcfirst_night_comb(tele, mjd; apred_dir, outroot, cfg, fpi_fnames, sci_exps = sci,
            fpi_fiberIndxs = fpi_idx, exclude_seqs = exclude, cache)
        s_fpifib = refit_night_fpifib(tele, mjd; apred_dir, outroot, cfg, exps)
        cst = refit_read_status(outroot, "comb", tele, mjd)
        gated = cst == "nocomb_gated"
        nw[] = 0
        prov = Dict{String, Any}("arc_product" => arcf, "comb_status" => cst,
            "chip_transform_fingerprint" => refit_fingerprint(cfg.chip_transform_file),
            "fiber_fill" => String(fiber_fill), "anchor_attempt" => attempt,
            "excluded_anchor_seqs" => join(exclude, ","))
        writer = function (k, e, kind, G, wtype, bits, xres)
            nw[] += 1
            return arcfirst_write_exp(arcfirst_exp_name(apred_dir, tele, mjd, e, kind), G, wtype, bits, xres;
                tele, mjd, expnum = e, kind, provenance = prov)
        end
        s_asm = refit_night_assemble(tele, mjd; apred_dir, outroot, cfg, ll, exps, inorder = true,
            on_grid = writer, fiber_fill, fpi_gated_off_in = gated, extra_bits = extra)
        combf = refit_path(outroot, "comb", tele, mjd)
        isfile(combf) || break                     # no comb: nothing to re-anchor
        aseq = h5open(f -> Int(read(f["anchor_seq"])), combf, "r")
        j = something(findfirst(==(aseq), A.used), 0)
        chk = arcfirst_check_night(outroot, tele, mjd; sky_check_ms)
        push!(checklog, "attempt=$(attempt)\tanchor_seq=$(aseq)\tpass=$(chk.pass)\tverifiable=$(chk.verifiable)" *
            "\tmetric_ms=$(round(chk.metric_ms, digits = 1))\tchip_ms=$(join(round.(chk.chip_ms, digits = 1), ","))" *
            "\tn_exp=$(chk.n_exp)\tn_slip=$(chk.n_slip)")
        (chk.pass || final_fail || j == 0) && break
        push!(tried, (j, isfinite(chk.metric_ms) ? chk.metric_ms : Inf))
        attempt += 1
        if attempt > max_reanchor
            # cap reached: rerun once with the best tried candidate, flagged WAVE_ANCHOR_CHECK_FAILED
            best_j = tried[argmin(last.(tried))][1]
            exclude = [jj for (jj, _) in tried if jj != best_j]
        else
            push!(exclude, j)
        end
        for st in ("comb", "fpifib", "assemble")
            rm(refit_path(outroot, st, tele, mjd); force = true)
            rm(refit_status_path(outroot, st, tele, mjd); force = true)
        end
    end
    isempty(checklog) || arcfirst_write_check(outroot, tele, mjd, checklog)
    return (; route, arc = s_arc, comb = s_comb, fpifib = s_fpifib, assemble = s_asm, n_exp_written = nw[],
        n_attempts = attempt + 1)
end

"""Exposures (expnum, kind) of one night from ar1Dcal (or ar2D/ar1D) file names of FIRST_CHIP."""
function arcfirst_exps_from_names(fnames, tele, mjd)
    rx = Regex("_$(tele)_$(mjd)_(\\d{4})_$(FIRST_CHIP)_([a-z]+)\\.h5\$")
    out = Tuple{Int, String}[]
    for f in fnames
        m = match(rx, basename(f))
        m === nothing && continue
        push!(out, (parse(Int, m.captures[1]), String(m.captures[2])))
    end
    return sort(unique(out))
end
