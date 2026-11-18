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
    (stage = "anchorcheck", reads = ["waveCalFPIArc"], writes = ["anchorCheck"]),   # this + neighbouring nights' combs
    (stage = "fpifib", reads = ["ar1Dcal", "waveCalFPIArc"], writes = ["fpiGuidePos"]),
    (stage = "assemble", reads = ["ar1Dcal", "skyLinePeaks", "arclampWaveCalNight", "waveCalFPIArc",
        "fpiGuidePos", "anchorCheck", "shape_template", "arc_linelist"], writes = ["waveCalRefit", "waveCalExp"]),
    (stage = "crossnight", reads = ["waveCalFPIArc"], writes = ["crossnight_flags"]),   # this + EARLIER nights
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
        sci_exps, fpi_fiberIndxs = Int[], exclude_seqs = Int[], cache = Dict{Symbol, Any}(),
        almanac_file::AbstractString = "", force_pair = nothing, keep = nothing)
    out = refit_path(outroot, "comb", tele, mjd)
    (isfile(out) && force_pair === nothing) && return "done"
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
        pass1_comb = nothing, provenance = prov, exclude_seqs, almanac_file, force_pair, keep)
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

"In-order inputs of a night: FPI fiber indices (almanac) and the science exposure numbers."
function arcfirst_night_inputs(tele, mjd, almanac_file, exps)
    fpi_ids = isempty(almanac_file) ? Int[] : get_fpi_fiberIDs_from_almanac(almanac_file, tele, mjd)
    fpi_idx = try
        filter(i -> 1 <= i <= N_FIBERS, fiberID2fiberIndx.(fpi_ids))
    catch
        Int[]
    end
    return fpi_idx, [e for (e, k) in exps if k == "object"]
end

"""
    arcfirst_night_phase1(tele, mjd; outdir, almanac_file, exps, fpi_fnames, cfg, checkpoint_mode, noarc_fallback)

Phase 1 of a night (stages (b) arc and (c) comb), run for EVERY night before any night's phase 2
so that the cross-night anchor continuity step (phase 2) sees its neighbours' comb products.
Returns `(route, arc, comb, done)`; `done = true` when the night's products already pass the
checkpoint (phase 2 then does nothing); `route == :legacy`: no usable arc fit (legacy chain).
"""
function arcfirst_night_phase1(tele, mjd; outdir, almanac_file, exps, fpi_fnames, cfg::RefitConfig,
        checkpoint_mode = "commit_same", noarc_fallback = "legacy")
    refit_check_config(cfg)
    apred_dir = joinpath(outdir, "apred")
    outroot = arcfirst_root(outdir)
    exps = sort(unique(collect(Tuple{Int, String}, exps)))
    asm = refit_path(outroot, "assemble", tele, mjd)
    have_all = isfile(asm) && all(x -> isfile(arcfirst_exp_name(apred_dir, tele, mjd, x...)), exps)
    (have_all && check_file(asm, mode = checkpoint_mode)) &&
        return (; route = :arcfirst, arc = "done", comb = "done", done = true)
    arcfirst_clear_night(outroot, apred_dir, tele, mjd, exps)
    s_arc = refit_night_arc(tele, mjd; apred_dir, almanac_file, outroot, cfg)
    arcf = refit_path(outroot, "arc", tele, mjd)
    usable = isfile(arcf) && refit_arc_usable(load_refit_arc(arcf))
    route = arcfirst_night_route(usable, noarc_fallback)
    route == :legacy && return (; route, arc = s_arc, comb = "", done = false)
    fpi_idx, sci = arcfirst_night_inputs(tele, mjd, almanac_file, exps)
    s_comb = arcfirst_night_comb(tele, mjd; apred_dir, outroot, cfg, fpi_fnames, sci_exps = sci,
        fpi_fiberIndxs = fpi_idx, almanac_file)
    return (; route, arc = s_arc, comb = s_comb, done = false)
end

"""
    arcfirst_night_phase2(tele, mjd; outdir, almanac_file, exps, fpi_fnames, cfg, fiber_fill, ll,
        max_reanchor = 2, sky_check_ms = 1000.0, history_roots = String[]) -> NamedTuple

Phase 2 of a night: (c') the cross-night anchor continuity step of verified_pair
(`refit_night_anchorcheck`: neighbours = this run's comb products of all nights (bulk) plus
`history_roots` (daily: earlier runs, trailing nights); a swap refits the comb through the
in-order gather), (d) FPI fibers, (e) assembly + waveCalExp, then the end-of-night check with the
capped conditional re-anchor (a re-anchor removes the anchor-check product: the continuity step is
not rerun on a re-anchored comb).
"""
function arcfirst_night_phase2(tele, mjd; outdir, almanac_file, exps, fpi_fnames, cfg::RefitConfig,
        fiber_fill = "fill_only", ll = nothing, max_reanchor::Integer = 2, sky_check_ms::Real = 1000.0,
        history_roots = String[])
    fiber_fill in ARCFIRST_FIBER_FILLS || throw(ArgumentError("fiber_fill \"$(fiber_fill)\" not in $(ARCFIRST_FIBER_FILLS)"))
    apred_dir = joinpath(outdir, "apred")
    outroot = arcfirst_root(outdir)
    exps = sort(unique(collect(Tuple{Int, String}, exps)))
    arcf = refit_path(outroot, "arc", tele, mjd)
    A = load_refit_arc(arcf)
    fpi_idx, sci = arcfirst_night_inputs(tele, mjd, almanac_file, exps)
    combfn = (fp, kp) -> arcfirst_night_comb(tele, mjd; apred_dir, outroot, cfg, fpi_fnames, sci_exps = sci,
        fpi_fiberIndxs = fpi_idx, almanac_file, force_pair = fp, keep = kp)
    s_check = refit_night_anchorcheck(tele, mjd; apred_dir, outroot, cfg, almanac_file, history_roots,
        refit_comb = combfn)
    cache = Dict{Symbol, Any}()
    exclude = Int[]                  # used-sequence indices rejected by the end-of-night check
    tried = Tuple{Int, Float64}[]    # (anchor j_seq, check metric) per attempt
    checklog = String[]
    s_comb = refit_read_status(outroot, "comb", tele, mjd)
    s_fpifib = s_asm = ""
    nw = Ref(0)
    attempt = 0
    best_j = 0
    clear_downstream() = for st in ("comb", "anchorcheck", "fpifib", "assemble")
        rm(refit_path(outroot, st, tele, mjd); force = true)
        rm(refit_status_path(outroot, st, tele, mjd); force = true)
    end
    while true
        # stages (c)-(e), a pure function of (arc fit, tonight's peaks, excluded anchors); attempt 0
        # uses phase 1's comb
        extra = attempt == 0 ? 0 : WAVE_REANCHORED
        final_fail = attempt > max_reanchor
        final_fail && (extra |= WAVE_ANCHOR_CHECK_FAILED)
        attempt == 0 || (s_comb = arcfirst_night_comb(tele, mjd; apred_dir, outroot, cfg, fpi_fnames, sci_exps = sci,
            fpi_fiberIndxs = fpi_idx, exclude_seqs = exclude, cache, almanac_file))
        s_fpifib = refit_night_fpifib(tele, mjd; apred_dir, outroot, cfg, exps)
        cst = refit_read_status(outroot, "comb", tele, mjd)
        gated = cst == "nocomb_gated"
        nw[] = 0
        prov = Dict{String, Any}("arc_product" => arcf, "comb_status" => cst,
            "chip_transform_fingerprint" => refit_fingerprint(cfg.chip_transform_file),
            "fiber_fill" => String(fiber_fill), "anchor_attempt" => attempt,
            "excluded_anchor_seqs" => join(exclude, ","), "anchorcheck_status" => String(s_check))
        writer = function (k, e, kind, G, wtype, bits, xres)
            nw[] += 1
            return arcfirst_write_exp(arcfirst_exp_name(apred_dir, tele, mjd, e, kind), G, wtype, bits, xres;
                tele, mjd, expnum = e, kind, provenance = prov)
        end
        s_asm = refit_night_assemble(tele, mjd; apred_dir, outroot, cfg, ll, exps, inorder = true,
            on_grid = writer, fiber_fill, fpi_gated_off_in = gated, extra_bits = extra)
        combf = refit_path(outroot, "comb", tele, mjd)
        if !isfile(combf)
            # no comb on the first attempt: nothing to re-anchor. On a re-anchor: the candidates ran out
            # before the cap -> back to the best tried candidate, flagged WAVE_ANCHOR_CHECK_FAILED
            (attempt == 0 || final_fail || isempty(tried)) && break
            push!(checklog, "attempt=$(attempt)\tno further anchor candidate ($(cst)); back to the best tried one")
            best_j = tried[argmin(last.(tried))][1]
            exclude = [jj for (jj, _) in tried if jj != best_j]
            attempt = max_reanchor + 1
            clear_downstream()
            continue
        end
        aseq = h5open(f -> Int(read(f["anchor_seq"])), combf, "r")
        j = something(findfirst(==(aseq), A.used), 0)
        chk = arcfirst_check_night(outroot, tele, mjd; sky_check_ms)
        push!(checklog, "attempt=$(attempt)\tanchor_seq=$(aseq)\tpass=$(chk.pass)\tverifiable=$(chk.verifiable)" *
            "\tmetric_ms=$(round(chk.metric_ms, digits = 1))\tchip_ms=$(join(round.(chk.chip_ms, digits = 1), ","))" *
            "\tn_exp=$(chk.n_exp)\tn_slip=$(chk.n_slip)\tanchorcheck=$(s_check)")
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
        clear_downstream()
        s_check = "skipped_reanchor"
    end
    isempty(checklog) || arcfirst_write_check(outroot, tele, mjd, checklog)
    return (; route = :arcfirst, comb = s_comb, anchorcheck = s_check, fpifib = s_fpifib, assemble = s_asm,
        n_exp_written = nw[], n_attempts = attempt + 1)
end

"""
    arcfirst_night(tele, mjd; outdir, almanac_file, exps, fpi_fnames, cfg, ...) -> NamedTuple

One night end to end (phase 1 then phase 2); the continuity step then sees whatever neighbouring
comb products exist (`history_roots` and earlier nights of the same root). A run over many nights
(`pipeline_2d_1d.jl`) calls phase 1 for all nights first.
"""
function arcfirst_night(tele, mjd; outdir, almanac_file, exps, fpi_fnames, cfg::RefitConfig,
        checkpoint_mode = "commit_same", fiber_fill = "fill_only", noarc_fallback = "legacy", ll = nothing,
        max_reanchor::Integer = 2, sky_check_ms::Real = 1000.0, history_roots = String[])
    p1 = arcfirst_night_phase1(tele, mjd; outdir, almanac_file, exps, fpi_fnames, cfg, checkpoint_mode, noarc_fallback)
    p1.done && return (; route = :arcfirst, arc = "done", comb = "done", anchorcheck = "done", fpifib = "done",
        assemble = "done", n_exp_written = 0, n_attempts = 0)
    p1.route == :legacy && return (; route = :legacy, arc = p1.arc, comb = "", anchorcheck = "", fpifib = "",
        assemble = "", n_exp_written = 0, n_attempts = 0)
    p2 = arcfirst_night_phase2(tele, mjd; outdir, almanac_file, exps, fpi_fnames, cfg, fiber_fill, ll,
        max_reanchor, sky_check_ms, history_roots)
    return merge((; arc = p1.arc), p2)
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

# ------------------------------------------------------------------------------------------
# cross-night comb step (after stage (c) of every night): bulk = a pass over the run's comb
# products; daily = the same kernel over trailing nights only (history root)
# ------------------------------------------------------------------------------------------
"""
    crossnight_comb_flags(mjds, lam, L; window = 15, min_prev = 5, fsr_tol = 0.2, nsig = 8.0)
        -> Vector{(mjd, check, value)}

Pure kernel of the cross-night comb checks (#420 phase 3; shared with the refit driver): per comb
night in MJD order, (i) FSR continuity: λ at a fixed fiber/pixel (`lam`) jumping by ~ a nonzero
integer FSR (λ²/2L) from the previous comb night; (ii) cavity L against the TRAILING `window`
nights (>= `min_prev` of them): |L - median| > `nsig` x max(1.4826 MAD, 1 Å). Uses only earlier
nights for each night, so the daily (trailing history) and a bulk run give the same flags.
A cross-night verification rule for the anchor (e.g. comb continuity at a fixed mode, branch
anchor-verified) drops in next to these checks.
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
            m = median(prev)
            s = max(1.4826 * median(abs.(prev .- m)), 1.0)
            abs(LL[i] - m) > nsig * s && push!(flags, (M[i], "cavity_jump_8sigma", (LL[i] - m) / s))
        end
    end
    return flags
end

"""
    arcfirst_crossnight(outroot, tele, target_mjds; history_roots = String[]) -> Vector

The cross-night step of the in-order chain for `target_mjds`: reads the arc-anchored comb products
(`<root>/comb/waveCalFPIArc_<tele>_<mjd>.h5`) of this run and of `history_roots` (e.g. the
previous nights' daily output roots), runs `crossnight_comb_flags` over all of them (each night is
judged against EARLIER nights only) and writes the flags of the target nights to
`<outroot>/status/crossnight_<tele>.csv`. Bulk: called once after every night's chain; daily: the
history supplies the trailing nights. Flags only (QA); a flagged night is not re-anchored here.
"""
function arcfirst_crossnight(outroot, tele, target_mjds; history_roots = String[])
    rx = Regex("^waveCalFPIArc_$(tele)_(\\d+)\\.h5\$")
    seen = Dict{Int, Tuple{Float64, Float64}}()
    for root in vcat(String[outroot], collect(String, history_roots))
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
                    (refit_prod_wave(lin, nl, 150, 2, [1024.0])[1], Float64(read(f["fpi_cavity_size"])))
                end
            catch
                continue
            end
            all(isfinite, v) && (seen[mj] = v)
        end
    end
    mj = sort(collect(keys(seen)))
    fl = crossnight_comb_flags(mj, [seen[m][1] for m in mj], [seen[m][2] for m in mj])
    fl = filter(x -> x[1] in target_mjds, fl)
    fn = joinpath(outroot, "status", "crossnight_$(tele).csv")
    mkpath(dirname(fn))
    open(fn, "w") do io
        println(io, "tele,mjd,check,value,n_comb_nights_seen")
        foreach(x -> println(io, "$(tele),$(x[1]),$(x[2]),$(x[3]),$(length(mj))"), fl)
    end
    return fl
end
