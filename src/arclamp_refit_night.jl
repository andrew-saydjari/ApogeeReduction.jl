# Per-night orchestration of the arc + FPI refit (stages A-D) on existing pass-1 products.
# Every input is opened read-only; every output goes under an output root:
#   <outroot>/arc/arclampWaveCalNight_<tele>_<mjd>.h5      stage A (stage-1 product, refit settings)
#   <outroot>/sky/skyNightRefit_<tele>_<mjd>.h5            sky stage (nights without a usable arc fit only;
#                                                          src/arclamp_refit_noarc.jl)
#   <outroot>/comb/waveCalFPIArc_<tele>_<mjd>.h5           stage B (stage 3)
#   <outroot>/fpifib/fpiGuidePos_<tele>_<mjd>.h5           stage C (stage 4)
#   <outroot>/exp/waveCalRefit_<tele>_<mjd>.h5             stage D (stages 5 + 6, per exposure)
#   <outroot>/status/<stage>_<tele>_<mjd>.txt              one line: status [reason]
# A stage is skipped when its product or a terminal status exists (resumable); `deferred` (stage B
# waiting for neighbour nights; only with the opt-in arcwin) and `fail` are not terminal for the stage
# that wrote them when the driver is asked to retry. Every night is processed independently of every
# other night (AKS 2026-11-18); global inputs are stored calibration files only.

const REFIT_STAGES = ("arc", "sky", "comb", "fpifib", "assemble")
# routes whose exposures are shifted to their own pooled sky lines minus the global arc-sky offset
const REFIT_SKYSHIFT_ROUTES = (:arc, :noarc_tmpl, :noarc_sky, :noarc_nightave)

refit_path(outroot, stage, tele, mjd) = joinpath(outroot,
    stage == "arc" ? "arc" : stage == "sky" ? "sky" : stage == "comb" ? "comb" : stage == "fpifib" ? "fpifib" : "exp",
    (stage == "arc" ? "arclampWaveCalNight" : stage == "sky" ? "skyNightRefit" : stage == "comb" ? "waveCalFPIArc" :
     stage == "fpifib" ? "fpiGuidePos" : "waveCalRefit") * "_$(tele)_$(mjd).h5")
refit_status_path(outroot, stage, tele, mjd) = joinpath(outroot, "status", "$(stage)_$(tele)_$(mjd).txt")

function refit_write_status(outroot, stage, tele, mjd, status::AbstractString, reason = "")
    fn = refit_status_path(outroot, stage, tele, mjd)
    mkpath(dirname(fn))
    tmp = fn * ".tmp$(getpid())"
    open(tmp, "w") do io
        println(io, status, isempty(reason) ? "" : "\t" * replace(String(reason), '\n' => ' '))
    end
    mv(tmp, fn; force = true)
    return status
end
function refit_read_status(outroot, stage, tele, mjd)
    fn = refit_status_path(outroot, stage, tele, mjd)
    isfile(fn) || return ""
    s = strip(read(fn, String))
    return String(split(s, '\t')[1])
end

# ------------------------------------------------------------------------------------------
# night inputs (pass 1, read-only)
# ------------------------------------------------------------------------------------------
"Exposures of a night with ar1Dunical products: sorted vector of (expnum, kind)."
function refit_night_exposures(mdir::AbstractString, tele::AbstractString, mjd::Integer)
    isdir(mdir) || return Tuple{Int, String}[]
    rx = Regex("^ar1Dunical_$(tele)_$(mjd)_(\\d{4})_([a-z]+)\\.h5\$")
    out = Tuple{Int, String}[]
    for fn in readdir(mdir)
        m = match(rx, fn)
        m === nothing && continue
        push!(out, (parse(Int, m.captures[1]), String(m.captures[2])))
    end
    return sort(out)
end

"Pass-1 night-average file: best type, the FPI files used (gated), FPI gate, sky seed."
function refit_night_ave(wavecal_dir, tele, mjd)
    fn = joinpath(wavecal_dir, "wavecalNightAve_$(tele)_$(mjd).h5")
    isfile(fn) || return nothing
    h5open(fn, "r") do f
        bwt = haskey(f, "best_wave_type") ? String(read(f["best_wave_type"])) : ""
        used = (haskey(f, "fpi") && haskey(f["fpi"], "used_fnames")) ? String.(read(f["fpi/used_fnames"])) : String[]
        gate = (haskey(f, "fpi_qa") && haskey(f["fpi_qa"], "gate_pass")) ? Bool(read(f["fpi_qa/gate_pass"])) : (bwt == "fpi")
        nfpi = (haskey(f, "fpi_qa") && haskey(f["fpi_qa"], "n_fpi_exposures")) ? Int(read(f["fpi_qa/n_fpi_exposures"])) : length(used)
        slin = haskey(f, "sky/nightAve_linParams") ? read(f["sky/nightAve_linParams"]) : nothing
        snl = haskey(f, "sky/nightAve_nlParams") ? read(f["sky/nightAve_nlParams"]) : nothing
        (; fname = fn, bwt, used, gate, nfpi, sky_lin = slin, sky_nl = snl)
    end
end

"Read one metadata scalar of an HDF5 product (default if absent)."
function refit_meta(fn::AbstractString, key::AbstractString, default)
    isfile(fn) || return default
    h5open(fn, "r") do f
        (haskey(f, "metadata") && haskey(f["metadata"], key)) ? read(f["metadata/$(key)"]) : default
    end
end

ar1dcal_name(mdir, tele, mjd, e, chip, kind) =
    joinpath(mdir, "ar1Dcal_$(tele)_$(mjd)_$(_arc_expstr(e))_$(chip)_$(kind).h5")

# ------------------------------------------------------------------------------------------
# stage A: arc night fit (stage-1 code, refit settings)
# ------------------------------------------------------------------------------------------
function refit_night_arc(tele, mjd; apred_dir, almanac_file, outroot, cfg::RefitConfig)
    out = refit_path(outroot, "arc", tele, mjd)
    isfile(out) && return "done"
    mkpath(dirname(out))
    ll, neon_used = refit_linelist(cfg, tele)
    cs = refit_chip_scales(cfg.chip_transform_file, tele, mjd)
    tmp = out * ".tmp$(getpid())"
    r = get_and_save_arclamp_night_wavecal(tele, mjd; apred_dir = apred_dir,
        almanac_file = almanac_file, outname = tmp, linelist = ll, chip_scales = cs,
        chip_scale_source = cs === nothing ? "free" : "fixed:" * basename(cfg.chip_transform_file),
        nsig = cfg.nsig, checkpoint_mode = "clobber")
    if r === nothing
        return refit_write_status(outroot, "arc", tele, mjd, "noarc", "no usable arc sequence")
    end
    aq = refit_quartz_pattern(tmp, almanac_file, tele, mjd)
    h5open(tmp, "r+") do f
        f["metadata/refit_neon_mode"] = neon_used
        f["metadata/refit_chip_transform_fingerprint"] = refit_fingerprint(cfg.chip_transform_file)
        f["refit_seq_after_quartz_same_dither"] = Int8.(aq.after)
        f["metadata/refit_quartz_pattern_n"] = count(aq.after)
        f["metadata/refit_quartz_pattern_gated"] = aq.n_gated
    end
    mv(tmp, out; force = true)
    refit_arc_usable(load_refit_arc(out)) ||
        return refit_write_status(outroot, "arc", tele, mjd, "unusable", "no sequence/fibers in the final fit")
    return refit_write_status(outroot, "arc", tele, mjd, "ok")
end

"""ARCFPI_VALIDATION Q7 pattern: arc sequences whose first exposure directly follows a quartz flat
with the same commanded dither (no commanded move). Such evening blocks sit > 0.02 px off in 78 %
of cases at LCO from MJD 60470 (23 % APO FPS). The commanded label is used here ONLY to name the
pattern; positions stay measured, and the #417 arc-only rejection decides exclusion. Returns
(after::Vector{Bool}, n_gated = number of such sequences that were rejected (ARCSEQ_REJECTED))."""
function refit_quartz_pattern(arc_fname, almanac_file, tele, mjd)
    the, une, sfl = h5open(arc_fname, "r") do f
        Int.(read(f["seq_thar_exp"])), Int.(read(f["seq_une_exp"])), Int.(read(f["seq_flags"]))
    end
    after = falses(length(the))
    isfile(almanac_file) || return (; after, n_gated = 0)
    ex, it, di = h5open(almanac_file, "r") do f
        g = f["raw/$(tele)/$(mjd)/exposures"]
        Int.(read(g["exposure"])), String.(read(g["image_type"])), Float64.(read(g["dithered_pixels"]))
    end
    for s in eachindex(the)
        e0 = min(the[s], une[s])
        prev = findall(<(e0), ex)
        isempty(prev) && continue
        j = prev[argmax(ex[prev])]
        i0 = findfirst(==(e0), ex)
        i0 === nothing && continue
        (occursin("quartz", lowercase(it[j])) && abs(di[j] - di[i0]) < 0.01) && (after[s] = true)
    end
    ng = count(s -> after[s] && (sfl[s] & ARCSEQ_REJECTED) != 0, eachindex(after))
    return (; after, n_gated = ng)
end

# ------------------------------------------------------------------------------------------
# stage B: comb (arc-anchored; arcwin for nights without an arc anchor)
# ------------------------------------------------------------------------------------------
"Nearest arc-anchored comb of the same telescope within ±window nights: (L, eps, mjd) or nothing."
function refit_neighbour_cavity(outroot, tele, mjd, window)
    d = joinpath(outroot, "comb")
    isdir(d) || return nothing
    best = nothing
    bd = typemax(Int)
    rx = Regex("^waveCalFPIArc_$(tele)_(\\d+)\\.h5\$")
    for fn in readdir(d)
        m = match(rx, fn)
        m === nothing && continue
        mj = parse(Int, m.captures[1])
        (mj == mjd || abs(mj - mjd) > window) && continue
        abs(mj - mjd) < bd || continue
        C = try
            load_refit_comb(joinpath(d, fn))
        catch
            continue
        end
        (C.anchor == COMB_ANCHOR_ARC && C.n_slip == 0) || continue
        bd = abs(mj - mjd)
        best = (C.L, C.eps, mj)
    end
    return best
end

"""
    refit_comb_anchor_stage(tele, mjd; outroot, fnames, combs, scis, cfg, pass1_comb = nothing,
        provenance = Dict{String, String}()) -> (status or nothing, reason)

The arc-anchored comb of one night from EXPLICIT inputs; shared by the archive refit
(`refit_night_comb`, inputs gathered from a previous run's products) and the in-order nightly
pipeline (`arcfirst_night_comb`, inputs measured in this run). Inputs:
  * `fnames`: the night's FPI comb peak files (R-chip names), already accepted by the caller's
    FPI gate;
  * `combs`: per file `(ref::ArcObjRef or nothing, t)`, the comb exposure's MEASURED position
    reference (chip-G slice of a per-exposure comb solution) and its mid-exposure MJD;
  * `scis`: per science exposure `(ref, t)` (`ref` only for `cfg.anchor_samepos == "pass1"`);
  * `pass1_comb`: optional previous comb (`refit_read_pass1_comb`) for the order-change QA only.
Reads the night's arc product (`refit_path(outroot, "arc", ...)`, the earlier stage). On success
writes the comb product and returns `("ok_arc", "")`; otherwise `(nothing, reason)` (no status
written: the caller decides the fallback).
"""
function refit_comb_anchor_stage(tele, mjd; outroot, fnames, combs, scis, cfg::RefitConfig,
        pass1_comb = nothing, provenance = Dict{String, String}(), exclude_seqs = Int[],
        almanac_file::AbstractString = "", force_pair = nothing, keep = nothing)
    out = refit_path(outroot, "comb", tele, mjd)
    mkpath(dirname(out))
    arcf = refit_path(outroot, "arc", tele, mjd)
    A = isfile(arcf) ? load_refit_arc(arcf) : nothing
    (A !== nothing && !refit_arc_usable(A)) && (A = nothing)
    A === nothing && return nothing, (isfile(arcf) ? "arc fit unusable" : "no arc product")
    sfl = h5open(f -> Int.(read(f["seq_flags"])), arcf, "r")
    # verified_pair: sequence adjacency of every (used arc sequence, FPI file) from the almanac
    # (exposure order and commanded moves only; an in-order input)
    adjacency = nothing; lampoff = nothing
    if cfg.anchor_block == "verified_pair"
        ex, code, dith = refit_exposure_sequence(almanac_file, tele, mjd)
        the, une, lo = h5open(f -> (Int.(read(f["seq_thar_exp"])), Int.(read(f["seq_une_exp"])),
            Float64.(read(f["seq_lamp_offset_px"]))), arcf, "r")
        fe = [parse(Int, split(basename(fn), "_")[4]) for fn in fnames]
        adjacency = [refit_pair_adjacency((the[s], une[s]), e, ex, code, dith) for s in A.used, e in fe]
        lampoff = lo[A.used]
        provenance["anchor_adjacency_counts"] = join(["$(c)=$(count(==(c), adjacency))" for c in ("strict", "nomove", "other", "unknown")], ",")
    end
    # experiment (comb_level = "stored_cavity"): the comb's level from the stored global cavity on EVERY night (no
    # nightly arc -> comb level transfer; arcs only seed the integer orders and calibrate the table globally)
    if cfg.comb_level == "stored_cavity" && cfg.anchor_block == "verified_pair"
        sc = comb_cavity_reference(cfg.cavity_ref_file, String(tele), mjd)
        # without a stored cavity (no file / no row / no L_A,eps columns) the night falls back to the pair route
        # below; recorded so a product never claims a stored level it did not use
        provenance["comb_level_effective"] = sc === nothing ? "pair (no stored cavity in cavity_ref_file)" : "stored_cavity"
        if sc !== nothing
            sscore = h5open(f -> haskey(f, "seq_reject_score") ? Float64.(read(f["seq_reject_score"]))[A.used] : nothing, arcf, "r")
            # frame + integer seed from the first sequence-adjacent pair when there is one (QA at the same position)
            pc = refit_adjacent_candidates(adjacency; seq_flags = sfl[A.used], seq_score = sscore, exclude_seqs,
                anchor_adjacency = cfg.anchor_adjacency, anchor_max_candidates = 1)
            Qs = refit_stored_cavity_comb(fnames, A, tele, sc...; seq_flags = sfl[A.used], seq_score = sscore,
                exclude_seqs, pass1_comb, comb_arc_qa_ms = cfg.comb_arc_qa_ms,
                pair = isempty(pc) ? nothing : (pc[1][1], pc[1][2]))
            Qs.ok || (provenance["comb_level_effective"] = "pair (stored-cavity comb fit failed: $(Qs.reason))")
            if Qs.ok
                # monitoring only (changes no solution): tonight's arc-pair cavity level vs the stored cavity
                mon = cfg.cavity_monitor_ms > 0 ?
                      refit_cavity_monitor(fnames, A, tele, sc...; adjacency, seq_flags = sfl[A.used],
                    seq_score = sscore, exclude_seqs, anchor_adjacency = cfg.anchor_adjacency,
                    anchor_max_candidates = cfg.anchor_max_candidates, tol_ms = cfg.cavity_monitor_ms) :
                      (; cavmon_decision = "off", cavmon_bits = 0)
                Qs = merge(Qs, mon, (; anchor_kind = "stored_cavity_level", anchor_bits = mon.cavmon_bits))
                keep === nothing || (Qs = merge(Qs, keep))
                ispath(out) && rm(out)
                write_refit_comb(out, Qs, fnames; tele, mjd, cfg, arc_fname = arcf, provenance = provenance)
                return refit_write_status(outroot, "comb", tele, mjd, "ok_arc"), ""
            end
        end
    end
    Q = refit_arc_anchored_comb(fnames, A, combs, tele; seq_flags = sfl, pass1_comb = pass1_comb,
        pos_match_px = cfg.pos_match_px, comb_arc_qa_ms = cfg.comb_arc_qa_ms, scis,
        anchor_block = cfg.anchor_block, anchor_prefer = cfg.anchor_prefer, anchor_samepos = cfg.anchor_samepos,
        anchor_pairing = cfg.anchor_pairing, adjacency, lamp_off_px = lampoff,
        anchor_max_candidates = cfg.anchor_max_candidates, anchor_pattern_tol_px = cfg.anchor_pattern_tol_px,
        anchor_lamp_tol_px = cfg.anchor_lamp_tol_px, anchor_adjacency = cfg.anchor_adjacency,
        anchor_block_h = cfg.anchor_block_h, exclude_seqs, force = force_pair,
        cavity_ref = cfg.cavity_ref_check == "on" ? comb_mode_reference(cfg.cavity_ref_file, String(tele), mjd) : NaN,
        cavity_ref_ms = cfg.cavity_ref_ms,
        stored_cavity = cfg.noadjacent_route == "stored_cavity" ? comb_cavity_reference(cfg.cavity_ref_file, String(tele), mjd) : nothing,
        seq_score = h5open(f -> haskey(f, "seq_reject_score") ? Float64.(read(f["seq_reject_score"])) : nothing, arcf, "r"))
    Q.ok || return nothing, "arc_unmatched: " * Q.reason
    keep === nothing || (Q = merge(Q, keep))
    ispath(out) && rm(out)        # a forced refit replaces the product (or a symlink to a reused one)
    write_refit_comb(out, Q, fnames; tele, mjd, cfg, arc_fname = arcf, provenance = provenance)
    return refit_write_status(outroot, "comb", tele, mjd, "ok_arc"), ""
end

"""
    refit_object_refs_pass1(apred_mjd_dir, tele, mjd)

REFIT ONLY (previous-run products; formerly #417's gate-B `arclamp_object_refs`): for every object
exposure with a per-exposure wavelength solution (`waveCalNightfpiDither_*_object.h5`, else
`waveCalNightskyDither_*_object.h5`), a chip-G slice of it (pixels 880-1170, fibers 100-200) and its
mid-exposure MJD. Opened read-only. Feeds `anchor_samepos = "pass1"` and the science times.
"""
function refit_object_refs_pass1(apred_mjd_dir::AbstractString, tele::AbstractString, mjd::Integer)
    isdir(apred_mjd_dir) || return ArcObjRef[]
    fs = readdir(apred_mjd_dir)
    found = Dict{Int, Tuple{String, String}}()
    for (kind, pre) in (("sky", "waveCalNightskyDither"), ("fpi", "waveCalNightfpiDither"))
        rx = Regex("^$(pre)_$(tele)_$(mjd)_(\\d{4})_object\\.h5\$")
        for fn in fs
            m = match(rx, fn)
            m === nothing && continue
            found[parse(Int, m.captures[1])] = (kind, joinpath(apred_mjd_dir, fn))   # fpi wins
        end
    end
    px0, px1 = 880, 1170
    refs = ArcObjRef[]
    for e in sort(collect(keys(found)))
        kind, fn = found[e]
        W = try
            h5open(fn, "r") do f
                Float64.(f["chipWaveSoln"][px0:px1, ARCGATE_FIBERS, 2])
            end
        catch
            continue
        end
        t = arclamp_mjd_mid(joinpath(apred_mjd_dir,
            "ar1Dcal_$(tele)_$(mjd)_$(_arc_expstr(e))_$(FIRST_CHIP)_object.h5"))
        push!(refs, ArcObjRef(e, t, kind, px0, ARCGATE_FIBERS, W))
    end
    return refs
end

"""Per comb exposure of a previous run (REFIT ONLY): (ref, t) with `ref` a chip-G slice of that
run's per-exposure comb solution `waveCalNightfpiDither_*_arclamp.h5` (a stored solution of a
previous run; the in-order pipeline measures this with `arcfirst_comb_refs_provisional`)."""
function refit_comb_refs_pass1(mdir, tele, mjd, fnames)
    map(fnames) do fn
        e = parse(Int, split(basename(fn), "_")[4])
        t = arclamp_mjd_mid(ar1dcal_name(mdir, tele, mjd, e, "R", "arclamp"))
        rf = joinpath(mdir, "waveCalNightfpiDither_$(tele)_$(mjd)_$(_arc_expstr(e))_arclamp.h5")
        ref = isfile(rf) ? ArcObjRef(e, t, "fpi", 880, ARCGATE_FIBERS,
            h5open(f -> Float64.(f["chipWaveSoln"][880:1170, ARCGATE_FIBERS, 2]), rf, "r")) : nothing
        (; ref, t)
    end
end

function refit_night_comb(tele, mjd; apred_dir, wavecal_dir, outroot, cfg::RefitConfig,
        allow_arcwin::Bool = false, almanac_file::AbstractString = "", force_pair = nothing, keep = nothing)
    out = refit_path(outroot, "comb", tele, mjd)
    (isfile(out) && force_pair === nothing) && return "done"
    mkpath(dirname(out))
    mdir = joinpath(apred_dir, string(mjd))
    # previous-run inputs (refit only): the pass-1 night-average file names the FPI files used and
    # carries the pass-1 FPI gate verdict and the sky seed (arcwin)
    NA = refit_night_ave(wavecal_dir, tele, mjd)
    NA === nothing && return refit_write_status(outroot, "comb", tele, mjd, "nocomb", "no pass-1 night-average file")
    isempty(NA.used) && return refit_write_status(outroot, "comb", tele, mjd, "nocomb", "no FPI exposures used by pass 1")
    NA.gate || return refit_write_status(outroot, "comb", tele, mjd, "nocomb_gated", "pass-1 FPI gate failed")
    fnames = [joinpath(mdir, basename(u)) for u in NA.used]
    all(isfile, fnames) || return refit_write_status(outroot, "comb", tele, mjd, "nocomb", "FPI peak files missing")
    p1 = try
        refit_read_pass1_comb(joinpath(mdir, "waveCalFPI_$(tele)_$(mjd)_arclamp.h5"))
    catch
        nothing
    end
    # per comb exposure: a slice of its pass-1 per-exposure comb solution (measured position) and its time
    combs = refit_comb_refs_pass1(mdir, tele, mjd, fnames)
    # science exposures: object-exposure times; their pass-1 per-exposure solutions (a later-stage
    # product) only feed the optional same-position input (cfg.anchor_samepos)
    scis = cfg.anchor_block == "verified_pair" ? NamedTuple[] :
           [(; ref = cfg.anchor_samepos == "pass1" ? r : nothing, t = r.mjd_mid) for r in refit_object_refs_pass1(mdir, tele, mjd)]
    arcf = refit_path(outroot, "arc", tele, mjd)
    prov = Dict("chip_transform_fingerprint" => refit_fingerprint(cfg.chip_transform_file),
        "n_science_refs" => string(length(scis)))
    # no usable arc fit: within-night ladder, rung 1 (comb integers from the night's own sky solution, cavity
    # fixed at the stored global cavity); never a neighbour night
    if refit_night_is_noarc(outroot, tele, mjd, cfg)
        cfg.noarc_route == "ladder" || return refit_write_status(outroot, "comb", tele, mjd, "nocomb_noarc",
            "no usable arc fit; noarc_route = $(cfg.noarc_route) (no comb)")
        skyf = refit_path(outroot, "sky", tele, mjd)
        isfile(skyf) || return refit_write_status(outroot, "comb", tele, mjd, "nocomb_noarc",
            "no usable arc fit and no night sky solution (sky stage: $(refit_read_status(outroot, "sky", tele, mjd)))")
        sc = comb_cavity_reference(cfg.cavity_ref_file, String(tele), mjd)
        sc === nothing && return refit_write_status(outroot, "comb", tele, mjd, "nocomb_noarc",
            "no usable arc fit and no stored cavity in cavity_ref_file")
        S = load_refit_sky(skyf)
        skyG = refit_eval_grid((f, c, px) -> refit_sky_eval(S.lin, S.nl, (0.0, 1.0), f, c, px), REFIT_PXG)
        Q = refit_skyseed_comb(fnames, S.lin, S.nl, sc...; skyG, pass1_comb = p1)
        Q.ok || return refit_write_status(outroot, "comb", tele, mjd, "nocomb_noarc", "sky-seeded comb: " * Q.reason)
        prov["comb_level_effective"] = "stored_cavity (no arcs: integers from the night sky solution)"
        prov["qa_reference"] = "night-average sky solution (chip_arc_minus_comb_ms = sky - comb, at different frames)"
        keep === nothing || (Q = merge(Q, keep))
        ispath(out) && rm(out)
        write_refit_comb(out, Q, fnames; tele, mjd, cfg, arc_fname = "", provenance = prov)
        return refit_write_status(outroot, "comb", tele, mjd, "ok_noarc_skyseed")
    end
    st, reason = refit_comb_anchor_stage(tele, mjd; outroot, fnames, combs, scis, cfg, pass1_comb = p1,
        provenance = prov, almanac_file, force_pair, keep)
    st === nothing || return st
    # arcwin (a neighbouring night's cavity) is a cross-night input: opt-in only (arcwin_nights > 0)
    cfg.arcwin_nights > 0 || return refit_write_status(outroot, "comb", tele, mjd, "nocomb_noanchor",
        reason * "; no arc anchor tonight (arcwin off)")
    allow_arcwin || return refit_write_status(outroot, "comb", tele, mjd, "deferred", reason)
    nb = refit_neighbour_cavity(outroot, tele, mjd, cfg.arcwin_nights)
    nb === nothing && return refit_write_status(outroot, "comb", tele, mjd, "nocomb_noanchor",
        reason * "; no arc-anchored neighbour within $(cfg.arcwin_nights) nights")
    (NA.sky_lin === nothing || NA.sky_nl === nothing) &&
        return refit_write_status(outroot, "comb", tele, mjd, "nocomb_noseed", reason * "; no sky seed")
    lin0 = ndims(NA.sky_lin) == 3 ? NA.sky_lin[:, :, 1] : NA.sky_lin
    nl0 = ndims(NA.sky_nl) == 3 ? NA.sky_nl[:, :, 1] : NA.sky_nl
    Q = refit_arcwin_comb(fnames, lin0, nl0, nb[1], nb[2]; pass1_comb = p1,
        comb_arc_qa_ms = cfg.comb_arc_qa_ms, neighbour_mjd = nb[3])
    Q.ok || return refit_write_status(outroot, "comb", tele, mjd, "fail", Q.reason)
    prov["arcwin_reason"] = reason
    write_refit_comb(out, Q, fnames; tele, mjd, cfg, arc_fname = isfile(arcf) ? arcf : "", provenance = prov)
    return refit_write_status(outroot, "comb", tele, mjd, "ok_arcwin", reason)
end

# ------------------------------------------------------------------------------------------
# stage C: FPI fibers of every exposure
# ------------------------------------------------------------------------------------------
"The night's FPI-fiber indices (fiberTypeList == \"fpiguide\" on the object exposures)."
function refit_night_fpi_fibers(mdir, tele, mjd, exps)
    fibs = Int[]
    for (e, kind) in exps
        kind == "object" || continue
        fn = ar1dcal_name(mdir, tele, mjd, e, "R", kind)
        isfile(fn) || continue
        ft = h5open(f -> String.(read(f["fiberTypeList"])), fn, "r")
        append!(fibs, findall(==("fpiguide"), ft))
        length(unique(fibs)) >= 2 && break
    end
    return sort(unique(fibs))
end

function refit_night_fpifib(tele, mjd; apred_dir, outroot, cfg::RefitConfig, exps = nothing)
    out = refit_path(outroot, "fpifib", tele, mjd)
    isfile(out) && return "done"
    combf = refit_path(outroot, "comb", tele, mjd)
    isfile(combf) || return refit_write_status(outroot, "fpifib", tele, mjd, "nocomb")
    mdir = joinpath(apred_dir, string(mjd))
    # exposure list: given by the caller (in-order pipeline: this run's list), else (refit) the
    # night's ar1Dunical files of the previous run
    exps === nothing && (exps = refit_night_exposures(mdir, tele, mjd))
    fibs = refit_night_fpi_fibers(mdir, tele, mjd, exps)
    isempty(fibs) && return refit_write_status(outroot, "fpifib", tele, mjd, "nofibers")
    C = load_refit_comb(combf)
    n = length(exps)
    a = fill(NaN, 2, 3, n); b = fill(NaN, 2, 3, n); ea = fill(NaN, 2, 3, n); eb = fill(NaN, 2, 3, n)
    np = zeros(Int, 2, 3, n); rms = fill(NaN, 2, 3, n); mfx = fill(NaN, 2, 3, n)
    for (k, (e, kind)) in enumerate(exps)
        fns = [ar1dcal_name(mdir, tele, mjd, e, ch, kind) for ch in CHIP_LIST]
        all(isfile, fns) || continue
        P = try
            refit_fpifib_exposure(fns, fibs, C; nsig = cfg.nsig)
        catch
            continue
        end
        a[:, :, k] .= P.a; b[:, :, k] .= P.b; ea[:, :, k] .= P.ea; eb[:, :, k] .= P.eb
        np[:, :, k] .= P.n; rms[:, :, k] .= P.rms; mfx[:, :, k] .= P.mfrac_max
    end
    mkpath(dirname(out))
    tmp = out * ".tmp$(getpid())"
    meta = Dict{String, Any}("product" => "fpiGuidePos", "stage" => 4, "tele" => String(tele),
        "mjd" => Int(mjd), "comb_product" => combf, "finder_nsig" => cfg.nsig,
        "frame" => "night comb solution at its first FPI exposure")
    safe_jldsave(tmp, meta; exp_expnum = [x[1] for x in exps], exp_kind = [x[2] for x in exps],
        fibers = vcat(fibs[1:min(end, 2)], zeros(Int, max(0, 2 - length(fibs)))),
        a_px = a, b_px = b, ea_px = ea, eb_px = eb, n_peaks = np, rms_px = rms, mfrac_max = mfx)
    mv(tmp, out; force = true)
    return refit_write_status(outroot, "fpifib", tele, mjd, "ok")
end

function load_refit_fpifib(fname)
    h5open(fname, "r") do f
        (; exp = Int.(read(f["exp_expnum"])), fibers = Int.(read(f["fibers"])), a = read(f["a_px"]),
            b = read(f["b_px"]), ea = read(f["ea_px"]), eb = read(f["eb_px"]), n = read(f["n_peaks"]),
            rms = read(f["rms_px"]), mfrac = read(f["mfrac_max"]))
    end
end

# ------------------------------------------------------------------------------------------
# stage D: assembly per exposure (stages 5 + 6)
# ------------------------------------------------------------------------------------------
"Pixel-grid evaluation of a model callback (f, c, pxvec) -> λ vector over chips and fibers."
function refit_eval_grid(modelfun, px::AbstractVector{<:Real})
    G = fill(NaN, length(px), N_FIBERS, N_CHIPS)
    pxv = collect(Float64, px)
    for c in 1:N_CHIPS, f in 1:N_FIBERS
        G[:, f, c] .= modelfun(f, c, pxv)
    end
    return G
end

"""
    refit_night_assemble(tele, mjd; apred_dir, wavecal_dir, outroot, cfg, ll = nothing)

Stages 5 + 6 for every exposure of the night (see the file header for the decision chain).
Writes `waveCalRefit_<tele>_<mjd>.h5` with, per exposure: Legendre coefficients of the new
solution (`wave_coeffs`, evaluate with `refit_exposure_wave`), `wavecal_type`, `bitmsk_wavecal`,
the sky cross-check (`xcheck_resid_ms`), the per-chip pooled sky offsets, the FPI-fiber terms,
the comparison with pass 1 (`dv_pass1_ms` = first-order RV prediction per chip and fiber), and
night-level QA (per-chip arc - comb residual, mode slips, order changes).
"""
function refit_night_assemble(tele, mjd; apred_dir, wavecal_dir = "", outroot, cfg::RefitConfig, ll = nothing,
        exps = nothing, inorder::Bool = false, on_grid = nothing, fiber_fill::AbstractString = "none",
        extra_bits::Integer = 0, fpi_gated_off_in = nothing)
    out = refit_path(outroot, "assemble", tele, mjd)
    isfile(out) && return "done"
    mdir = joinpath(apred_dir, string(mjd))
    # previous-run inputs are read only when `inorder == false` (archive refit): the exposure list
    # and times from its ar1Dunical files, its adopted wavecal_type (pass-1 copy route and QA), its
    # night-average file (FPI gated-off flag) and its adopted solutions (pass-1 copy + comparison).
    # `inorder = true` (nightly pipeline): exposures from the caller, times from the ar1Dcal files,
    # FPI gated-off from the caller, no pass-1 copy (route :pass1 -> no grid; the caller falls back).
    exps === nothing && (exps = refit_night_exposures(mdir, tele, mjd))
    isempty(exps) && return refit_write_status(outroot, "assemble", tele, mjd, "noexp")
    arcf = refit_path(outroot, "arc", tele, mjd)
    combf = refit_path(outroot, "comb", tele, mjd)
    fpif = refit_path(outroot, "fpifib", tele, mjd)
    A = isfile(arcf) ? load_refit_arc(arcf) : nothing
    (A !== nothing && !refit_arc_usable(A)) && (A = nothing)
    cfg.ignore_arcs && (A = nothing)      # validation only: the no-arc ladder on a normal night
    # no usable arc fit (archive refit): the within-night ladder (src/arclamp_refit_noarc.jl); the in-order
    # pipeline sends such nights to the legacy chain instead (route :pass1 -> none_none -> caller)
    noarc = !inorder && A === nothing
    skyf = refit_path(outroot, "sky", tele, mjd)
    S = (noarc && isfile(skyf)) ? load_refit_sky(skyf) : nothing
    nfill = 0
    (A !== nothing && fiber_fill == "fill_only") && ((A, nfill) = refit_arc_fill_fibers(A))
    C = isfile(combf) ? load_refit_comb(combf) : nothing
    FP = isfile(fpif) ? load_refit_fpifib(fpif) : nothing
    NA = inorder ? nothing : refit_night_ave(wavecal_dir, tele, mjd)
    ll === nothing && (ll = refit_linelist(cfg, tele)[1])
    arcbits = A === nothing ? 0 : refit_arc_night_bits(A.night_flags, A.seq_flags)
    if A !== nothing && isfile(arcf)
        qp = refit_meta(arcf, "refit_quartz_pattern_gated", 0)
        (qp isa Integer && qp > 0) && (arcbits |= WAVE_ARC_QUARTZ_PATTERN)
    end
    fpi_gated_off = fpi_gated_off_in !== nothing ? Bool(fpi_gated_off_in) :
                    (C === nothing && NA !== nothing && NA.nfpi > 0)
    # era from the plate/FPS MJD boundary only (never from FPI light / comb presence)
    era = refit_era(tele, mjd)
    shape_choice = era == "plate" ? cfg.plate_shape : (tele == "apo" ? cfg.fps_shape_apo : cfg.fps_shape_lco)
    # FPS era without a night comb (no FPI exposures, pass-1 FPI gate failed, comb fit failed or no
    # anchor): every exposure takes the arc-route (or pass-1 copy) fallback, flagged on every exposure
    fps_no_comb = era == "fps" && C === nothing
    Tp = (A !== nothing && shape_choice == "template") ?
         load_refit_template(cfg.shape_template_file, tele, cfg.shape_template_block) : nothing
    # no-arc rung 2 (template shifted to the sky) vs rung 3 (sky-only chain)
    TpN = noarc ? load_refit_template(cfg.shape_template_file, tele, cfg.shape_template_block) : nothing
    TpN_ok = TpN === nothing ? falses(N_FIBERS) : refit_template_ok(TpN)
    noarc_tmpl = TpN !== nothing && (cfg.noarc_route == "template_sky" ||
                  (cfg.noarc_route == "ladder" && (era == "fps" || cfg.noarc_plate_shape == "template")))
    smodel(d) = (f, c, px) -> refit_sky_eval(S.lin, S.nl, view(d, f, :), f, c, px)
    d_nightave = repeat([0.0 1.0], N_FIBERS)
    # template fibers without a complete template take the exposure's sky-chain model (`fb`; WAVE_FIBER_FILLED)
    tmodel(dl, fb) = (f, c, px) -> TpN_ok[f] ? [refit_template_noarc(TpN, dl, f, c, p) for p in px] :
                                   (fb === nothing ? fill(NaN, length(px)) : fb(f, c, px))
    TP = Tp === nothing ? nothing : refit_fit_template(Tp, A, ll)
    (TP !== nothing && count(TP[2]) < 200) && (Tp = nothing; TP = nothing)
    combbits = 0
    if C !== nothing
        C.anchor == COMB_ANCHOR_ARCWIN && (combbits |= WAVE_ABS_ARC_WINDOW)
        C.n_slip > 0 && (combbits |= WAVE_COMB_MODE_SLIP)
        C.inconsistent && (combbits |= WAVE_COMB_ARC_INCONSISTENT)
        combbits |= C.anchor_bits
        C.anchor == COMB_ANCHOR_ARC && (combbits |= arcbits)   # an arcwin comb carries no arc flags
    end
    vpxG = tele == "apo" ? 5253.0 : 5281.0
    vpx = tele == "apo" ? (4231.0, 5253.0, 6328.0) : (4256.0, 5281.0, 6359.0)
    # GLOBAL arc-sky offset per chip (stored calibration, not tonight's data): the sky line list's level
    # against the arc fit at the same position; zeros (flagged) without the file
    Osky, omiss = refit_mode_offsets(cfg.mode_offset_file, String(tele), "arc_sky"; era, mjd)
    nE = length(exps)
    tmid = fill(NaN, nE)
    p1type = fill("", nE)
    route = fill(:pass1, nE)            # :fpifib, :arc, :pass1
    absname = fill("", nE)
    relname = fill("", nE)
    bits = zeros(Int, nE)
    models = Vector{Any}(nothing, nE)    # (f, c, pxvec) -> λ
    clus = zeros(Int, nE)
    hassky = falses(nE)
    pos_off_px = fill(NaN, nE)
    vraw = fill(NaN, 3, nE)              # pooled sky - applied relative model (m/s), before any tweak
    theta = fill(NaN, 3, 3, nE)          # sky per-exposure terms (zp, stretch, tilt) x chip vs the model
    nlines = zeros(Int, 3, nE)
    fz = fill(NaN32, N_FIBERS, 3, nE)    # per-fiber sky terms (fiber-relative)
    fs = fill(NaN32, N_FIBERS, 3, nE)
    tilt = fill(NaN, 3, nE)
    nfpi_used = zeros(Int, nE)
    arcmodel(st) = (f, c, px) -> A.okf[f] ?
        arcfit_eval(view(A.TH, f, :), A.g, A.lay, f, c, collect(Float64, px);
            grpv = st[1], d = st[2], t = st[3], tau = st[4]) : fill(NaN, length(px))
    # fiber_fill == "fill_only": a fiber without its own template fit (e.g. an arc fiber filled
    # across the slit) takes the arc shape instead of NaN (flagged WAVE_FIBER_FILLED)
    shapemodel(st) = Tp === nothing ? arcmodel(st) :
        (f, c, px) -> TP[2][f] ? [refit_template_lambda(Tp, view(TP[1], :, f), f, c,
            p + refit_state_shift_px(A, f, c, p, st)) for p in px] :
            (fiber_fill == "fill_only" ? arcmodel(st)(f, c, px) : fill(NaN, length(px)))
    tp_missing = (Tp !== nothing && fiber_fill == "fill_only") ? count(f -> !TP[2][f], 1:N_FIBERS) : 0
    # no-arc base model of one exposure: rung 2 (template shifted to the pooled sky), rung 3 (sky-chain dither
    # solution), else the night-average sky solution; `nothing` when no rung is possible
    function noarc_base(e, sky, skyP)
        js = S === nothing ? nothing : findfirst(==(e), S.expnum)
        fb = S === nothing ? nothing : (js !== nothing && S.ok[js]) ? smodel(S.dither[:, :, js]) : smodel(d_nightave)
        if sky !== nothing && noarc_tmpl
            m0 = tmodel(zeros(3), fb)
            v0 = refit_sky_velocities(skyP, (f, c, p) -> m0(f, c, [p])[1])
            z0, _, _ = refit_pooled_sky(skyP, v0)
            dl = refit_template_shift_px(TpN, z0)
            all(isfinite, dl) && return (; model = tmodel(dl, fb), route = :noarc_tmpl, relname = "templateskyshift",
                bits = WAVE_NOARC_TEMPLATE_SKY)
        end
        (sky !== nothing && js !== nothing && S.ok[js]) && return (; model = smodel(S.dither[:, :, js]),
            route = :noarc_sky, relname = "skydither", bits = WAVE_NOARC_SKY_ONLY)
        S === nothing && return nothing
        return (; model = smodel(d_nightave), route = :noarc_nightave, relname = "skynightave",
            bits = WAVE_NOARC_SKY_ONLY | WAVE_NIGHT_AVE | (sky === nothing ? WAVE_POS_UNMATCHED : 0))
    end
    # no-arc per-fiber zero points (m/s, fiber-relative per chip): neither the template (no per-fiber fit without
    # arcs) nor the sky chain (smoothed across fibers) carries a fiber's own offset, and offsets beyond the pooled
    # estimator's 600 m/s row clip never reach the static per-fiber sky term. Robust median over ALL of the night's
    # sky rows of (row velocity - the exposure's pooled chip value), >= 5 rows, |.| < 5 km/s; applied to every
    # no-arc model (the static term then corrects the remainder as on every route).
    fibz = zeros(N_FIBERS, 3)
    if noarc && (TpN !== nothing || S !== nothing)
        acc = [Float64[] for _ in 1:N_FIBERS, _ in 1:3]
        skydir0 = isempty(cfg.sky_peaks_dir) ? mdir : joinpath(cfg.sky_peaks_dir, string(mjd))
        for (e, kind) in exps
            kind == "object" || continue
            sk = refit_read_sky(skydir0, tele, mjd, e)
            (sk === nothing || isempty(sk.pix)) && continue
            nb = noarc_base(e, sk, sk)
            (nb === nothing || nb.route == :noarc_nightave) && continue
            v = refit_sky_velocities(sk, (f, c, p) -> nb.model(f, c, [p])[1])
            zp, _, _ = refit_pooled_sky(sk, v)
            for i in eachindex(v)
                c = sk.chip[i]
                (isfinite(v[i]) && isfinite(zp[c]) && abs(v[i] - zp[c]) < 5000) && push!(acc[sk.fiber[i], c], v[i] - zp[c])
            end
        end
        for c in 1:3
            for f in 1:N_FIBERS
                length(acc[f, c]) >= 5 && (fibz[f, c] = median(acc[f, c]))
            end
            ok = [length(acc[f, c]) >= 5 for f in 1:N_FIBERS]
            any(ok) && (fibz[ok, c] .-= median(fibz[ok, c]))
        end
    end
    noarc_fibcorr(m) = (f, c, px) -> m(f, c, px) .* (1 + fibz[f, c] / REFIT_C_MPS)
    # ---------------- pass A: relative model and sky measurements per exposure
    for (k, (e, kind)) in enumerate(exps)
        if inorder
            tm = arclamp_mjd_mid(ar1dcal_name(mdir, tele, mjd, e, FIRST_CHIP, kind))
        else
            uf = joinpath(mdir, "ar1Dunical_$(tele)_$(mjd)_$(_arc_expstr(e))_$(kind).h5")
            p1type[k] = String(refit_meta(uf, "wavecal_type", ""))
            tm = refit_meta(uf, "mjd_mid_exposure", NaN)
        end
        tmid[k] = tm isa Real ? Float64(tm) : NaN
        bits[k] |= extra_bits
        # sky-line peaks: this run's apred tree, or a separate peaks tree (refit: re-measured peaks, #407 list)
        skydir = isempty(cfg.sky_peaks_dir) ? mdir : joinpath(cfg.sky_peaks_dir, string(mjd))
        sky = kind == "object" ? refit_read_sky(skydir, tele, mjd, e) : nothing
        (sky !== nothing && isempty(sky.pix)) && (sky = nothing)
        # rows that set the POOLED per-exposure terms (and the position match): all fibers, or only sky
        # fibers and faint-target fibers (`sky_line_fibers = "faint"`: stellar lines bias the sky-line
        # centroids on bright-target fibers); the per-fiber terms always use every fiber
        skyP = sky
        if sky !== nothing && cfg.sky_line_fibers == "faint"
            fm = refit_faint_fibers(mdir, tele, mjd, e, cfg.sky_faint_snr)
            skyP = fm === nothing ? sky : refit_sky_subset(sky, fm)
            isempty(skyP.pix) && (skyP = nothing)
        end
        skyP === nothing && (sky = nothing)
        hassky[k] = sky !== nothing
        # (1) FPI fibers on a comb night (a no-arc night's comb only under the ladder: rung 1)
        if C !== nothing && FP !== nothing && cfg.fps_route == "fpifib" && (!noarc || cfg.noarc_route == "ladder")
            kk = findfirst(==(e), FP.exp)
            if kk !== nothing
                P = (; a = FP.a[:, :, kk], b = FP.b[:, :, kk], ea = FP.ea[:, :, kk], eb = FP.eb[:, :, kk],
                    n = FP.n[:, :, kk], rms = FP.rms[:, :, kk], fibers = filter(>(0), FP.fibers))
                dfun, nf, tl, fl = refit_fpifib_transform(P, cfg.exp_terms)
                if dfun !== nothing
                    lin, nl = C.lin, C.nl
                    models[k] = (f, c, px) -> refit_prod_wave(lin, nl, f, c, px .- [dfun(f, c, p) for p in px])
                    route[k] = :fpifib
                    absname[k] = C.anchor
                    relname[k] = "fpifib"
                    bits[k] |= combbits | fl
                    tilt[:, k] .= tl
                    nfpi_used[k] = nf
                end
            end
        end
        # (2) arc route: plate era, nights without comb, FPS exposures without usable FPI fibers
        if route[k] == :pass1 && A !== nothing
            best = 0
            bo = Inf
            if sky !== nothing   # measured position: the cluster whose (fixed-state) arc fit the sky matches best
                for cl in refit_arc_clusters(A)
                    st = refit_arc_state(A, NaN, cl; mode = :nearest)
                    st === nothing && continue
                    m = arcmodel(st)
                    v = refit_sky_velocities(skyP, (f, c, p) -> m(f, c, [p])[1])
                    zp, _, _ = refit_pooled_sky(skyP, v)
                    vg = isfinite(zp[2]) ? zp[2] : _arc_nanmean(zp)
                    isfinite(vg) || continue
                    if abs(vg / vpxG) < abs(bo)
                        bo = vg / vpxG
                        best = cl
                    end
                end
            end
            usesky = best != 0
            if best == 0
                # no sky: the cluster of the nearest used sequence in time (position NOT measured)
                j = isfinite(tmid[k]) ? argmin(abs.(A.tseq_u .- tmid[k])) : 1
                best = A.clus_u[j]
                bits[k] |= WAVE_POS_UNMATCHED
            elseif abs(bo) > cfg.pos_match_px
                bits[k] |= WAVE_POS_UNMATCHED
            end
            clus[k] = best
            pos_off_px[k] = isfinite(bo) ? bo : NaN
            # night solution = the arc fit at a FIXED state (no arc time model): with sky lines, the matched
            # cluster's reference sequence (the exposure is then shifted onto its own sky lines in pass B);
            # without sky lines (calibration frames), the state of the nearest arc sequence in time
            st = usesky ? refit_arc_state(A, NaN, best; mode = :nearest) : refit_arc_state(A, tmid[k], best; mode = :nearest)
            models[k] = shapemodel(st)
            Tp === nothing && (bits[k] |= WAVE_SHAPE_ARC)
            route[k] = :arc
            absname[k] = "arc"
            # relative name says which shape: arc shape ("arc...") or the comb template ("comb...")
            relname[k] = (Tp === nothing ? "arc" : "comb") * (usesky ? "skyshift" : "nearest")
            bits[k] |= arcbits
            usesky || (bits[k] |= WAVE_REL_ARC_NEAREST)
            (C !== nothing || (FP !== nothing && kind == "object")) && (bits[k] |= WAVE_FPIFIB_UNAVAILABLE)
        end
        # (2b) no usable arc fit: rung 2 (template shape shifted to this exposure's sky), rung 3 (sky-only chain
        # of this exposure), else the night-average sky solution (no own sky lines: position not measured)
        if route[k] == :pass1 && noarc
            nb = noarc_base(e, sky, skyP)
            if nb !== nothing
                models[k] = noarc_fibcorr(nb.model)
                route[k] = nb.route
                relname[k] = nb.relname
                bits[k] |= nb.bits
            end
            route[k] == :pass1 && (bits[k] |= WAVE_POS_UNMATCHED)    # nothing possible: no solution
            absname[k] = "noarc"
            (route[k] != :pass1 && (C !== nothing || (FP !== nothing && kind == "object"))) &&
                (bits[k] |= WAVE_FPIFIB_UNAVAILABLE)
        end
        fpi_gated_off && (bits[k] |= WAVE_FPI_GATED_OFF)
        fps_no_comb && (bits[k] |= WAVE_FPS_NO_COMB)
        # (3) sky measurements against the relative model
        if models[k] !== nothing && sky !== nothing
            m = models[k]
            v = refit_sky_velocities(sky, (f, c, p) -> m(f, c, [p])[1])
            vP = skyP === sky ? v : refit_sky_velocities(skyP, (f, c, p) -> m(f, c, [p])[1])
            zp, nlc, keepP = refit_pooled_sky(skyP, vP)
            vraw[:, k] .= zp
            nlines[:, k] .= nlc
            theta[:, :, k] .= refit_sky_terms(skyP, vP, keepP, zp, nlc,
                route[k] in REFIT_SKYSHIFT_ROUTES ? cfg.sky_exp_terms : "chipzp")
            if cfg.sky_fiber_term in ("static", "full")
                keep = skyP === sky ? keepP : refit_pooled_sky(sky, v)[3]
                z, s = refit_fiber_sky_terms(sky, v, keep)
                fz[:, :, k] .= Float32.(z)
                fs[:, :, k] .= Float32.(s)
            end
        end
    end
    fpianchor = [let s = [k for k in 1:nE if route[k] == :fpifib && isfinite(vraw[c, k])]
                     isempty(s) ? NaN : median(vraw[c, s])
                 end for c in 1:3]
    static = fill(NaN, N_FIBERS, 3)
    if cfg.sky_fiber_term == "static"
        for f in 1:N_FIBERS, c in 1:3
            static[f, c] = _arc_nanmedian(Float64.(fz[f, c, :]))
        end
    elseif cfg.sky_fiber_term == "stored"
        # GLOBAL stored per-fiber term (telescope, era): built from sky / faint-target fibers only
        static .= refit_fiber_term_table(cfg.fiber_term_file, String(tele), era)
    end
    # ---------------- pass B: final grids, coefficients, QA
    B = refit_coeff_basis(cfg.coeff_degree)
    deg1 = cfg.coeff_degree + 1
    coeffs = fill(NaN, deg1, N_FIBERS, N_CHIPS, nE)
    cmax = fill(NaN, nE)
    applied = zeros(3, 3, nE)
    xres = fill(NaN32, nE)
    xres_chip = fill(NaN32, 3, nE)
    dvp1 = fill(NaN32, N_CHIPS, N_FIBERS, nE)
    p1zp = fill(NaN, 3, nE); p1shape = fill(NaN, 3, nE); p1rs = fill(NaN, 3, nE)
    wtype = fill("", nE)
    px_all = collect(1.0:N_XPIX)
    for (k, (e, kind)) in enumerate(exps)
        src = inorder ? nothing : refit_pass1_source(apred_dir, wavecal_dir, tele, mjd, e, kind, p1type[k])
        G = nothing
        if route[k] == :pass1
            # in-order: no solution, the caller's fallback decides; archive refit: no rung of the no-arc ladder
            # was possible (no sky solution on the night): no solution (NaN coefficients). The pass-1 copy is
            # retired (WAVE_RESERVED_18).
            wtype[k] = "none_none"
        else
            wtype[k] = absname[k] * "_" * relname[k]
            G = refit_eval_grid(models[k], px_all)
            T = zeros(3, 3)
            if route[k] in REFIT_SKYSHIFT_ROUTES && hassky[k] && any(isfinite, view(theta, 1, :, k))
                # shift the arc fit onto this exposure's own sky lines: T = theta - O (global arc-sky offset);
                # a chip with < 3 sky lines takes the measured chips' zero point in pixels, stretch/tilt 0
                pxs = [(theta[1, c, k] - Osky[1, c]) / vpx[c] for c in 1:3 if isfinite(theta[1, c, k])]
                for c in 1:3
                    if isfinite(theta[1, c, k])
                        for q in 1:3
                            T[q, c] = isfinite(theta[q, c, k]) ? theta[q, c, k] - Osky[q, c] : 0.0
                        end
                    else
                        T[1, c] = median(pxs) * vpx[c]
                    end
                end
                omiss && (bits[k] |= WAVE_MODE_OFFSET_MISSING)
            elseif route[k] == :noarc_nightave
                # no own sky lines: the night-average sky solution (sky level) put on the same level as every
                # sky-shifted route, minus the global arc-sky offset
                for c in 1:3, q in 1:3
                    T[q, c] = -Osky[q, c]
                end
                omiss && (bits[k] |= WAVE_MODE_OFFSET_MISSING)
            elseif route[k] == :fpifib
                d0 = [vraw[c, k] - fpianchor[c] for c in 1:3]
                xres_chip[:, k] .= Float32.(d0)
                if cfg.sky_fiber_term == "pooled"
                    for c in 1:3
                        isfinite(d0[c]) && (T[1, c] = d0[c])
                    end
                end
            end
            applied[:, :, k] .= T
            refit_apply_sky_terms!(G, T)
            if cfg.sky_fiber_term in ("static", "stored", "full")
                u = (px_all .- 1024.5) ./ 1024
                for c in 1:3, f in 1:N_FIBERS
                    if cfg.sky_fiber_term in ("static", "stored") && isfinite(static[f, c])
                        G[:, f, c] .*= 1 + static[f, c] / REFIT_C_MPS
                    elseif cfg.sky_fiber_term == "full" && isfinite(fz[f, c, k])
                        G[:, f, c] .*= 1 .+ (fz[f, c, k] .+ fs[f, c, k] .* u) ./ REFIT_C_MPS
                    end
                end
            end
            xc = Float64.(xres_chip[:, k])
            any(isfinite, xc) && (xres[k] = Float32(_arc_nanmean(xc)))
            (isfinite(xres[k]) && abs(xres[k]) > cfg.xcheck_threshold_ms) && (bits[k] |= WAVE_LEVEL_INCONSISTENT)
        end
        G === nothing && continue
        # filled arc fibers matter only where the arc model is evaluated (arc route); on the FPI route
        # the arc only seeds the comb (seed fibers are filled anyway, refit_arc_lin_nl)
        if route[k] == :arc && (nfill > 0 || tp_missing > 0)
            bits[k] |= WAVE_FIBER_FILLED
        elseif route[k] == :noarc_tmpl && !all(TpN_ok)
            bits[k] |= WAVE_FIBER_FILLED
        end
        # the in-order caller writes the per-exposure product here and may add bits (e.g. a final
        # fill of non-finite fibers); its return value is OR'd into the exposure's bits
        on_grid === nothing || (bits[k] |= on_grid(k, e, kind, G, wtype[k], bits[k], xres[k]))
        coeffs[:, :, :, k], cmax[k] = refit_grid_to_coeffs(G, B)
        if src !== nothing
            Gp1 = try
                h5open(f -> Float64.(f[src[2]][REFIT_PXG, :, :]), src[1], "r")
            catch
                nothing
            end
            if Gp1 !== nothing
                d, z, s, r = refit_vs_pass1(G[REFIT_PXG, :, :], Gp1)
                dvp1[:, :, k] .= d
                p1zp[:, k] .= z; p1shape[:, k] .= s; p1rs[:, k] .= r
            end
        end
    end
    # ---------------- write
    mkpath(dirname(out))
    tmp = out * ".tmp$(getpid())"
    meta = Dict{String, Any}("product" => "waveCalRefit", "stage" => 6, "tele" => String(tele),
        "mjd" => Int(mjd), "format_version" => REFIT_FORMAT_VERSION,
        "coeff_basis" => "Legendre in t = (pix - 1024.5)/1023.5, degree $(cfg.coeff_degree); evaluate with refit_exposure_wave",
        "arc_product" => isfile(arcf) ? arcf : "", "comb_product" => isfile(combf) ? combf : "",
        "fpifib_product" => isfile(fpif) ? fpif : "", "comb_anchor" => C === nothing ? "none" : C.anchor,
        "template" => Tp === nothing ? "none" : "$(cfg.shape_template_file)#$(cfg.shape_template_block)",
        "shape_choice" => shape_choice, "era" => era,
        "era_boundary_mjd" => get_fps_plate_divide(String(tele)), "fps_no_comb" => fps_no_comb,
        "comb_anchor_kind" => C === nothing ? "none" : C.anchor_kind,
        "template_fingerprint" => Tp === nothing ? "none" : Tp.fingerprint,
        "linelist_fingerprint" => ll.fingerprint, "linelist_selection" => ll.selection,
        "chip_transform_fingerprint" => refit_fingerprint(cfg.chip_transform_file),
        "neon_offset_fingerprint" => refit_fingerprint(cfg.neon_offset_file),
        "arc_route" => "arc fit (fixed state) shifted to each exposure's pooled sky lines minus the global arc-sky offset",
        "mode_offset_file" => _refit_default_path(cfg.mode_offset_file, MODE_OFFSETS_DEFAULT),
        "mode_offset_fingerprint" => refit_fingerprint(_refit_default_path(cfg.mode_offset_file, MODE_OFFSETS_DEFAULT)),
        "mode_offset_missing" => omiss, "bit_names" => join(WAVE_BIT_NAMES, ","),
        "inorder" => inorder, "fiber_fill" => String(fiber_fill), "n_arc_fibers_filled" => nfill,
        # where each scale enters (ARCFPI_VALIDATION scales list)
        "scale_arc_chip" => A === nothing ? "none" : A.chip_scale_source,
        "scale_comb" => "comb chip transform (cporder 1) and quartic fitted per night; cavity (L, eps) per night",
        "scale_exposure_fpifib" => cfg.exp_terms * " (per-chip zero point and stretch from the two FPI fibers)",
        "scale_exposure_sky" => cfg.sky_exp_terms * " (per chip, pooled sky minus the global arc-sky offset, applied directly)",
        "scale_night_average" => noarc ? "no-arc night: sky-only chain night average for exposures without sky lines" :
                                 "none built",
        "noarc" => noarc, "noarc_route" => cfg.noarc_route, "noarc_template_rung" => noarc_tmpl,
        "sky_product" => (noarc && isfile(skyf)) ? skyf : "")
    for (kx, v) in refit_config_dict(cfg)
        meta["cfg_" * kx] = v
    end
    safe_jldsave(tmp, meta; exp_expnum = [x[1] for x in exps], exp_kind = [x[2] for x in exps],
        exp_mjd_mid = tmid, wavecal_type = wtype, bitmsk_wavecal = Int64.(bits),
        pass1_wavecal_type = p1type, wave_coeffs = coeffs, coeff_max_resid_A = cmax,
        xcheck_resid_ms = xres, xcheck_resid_chip_ms = xres_chip, sky_minus_model_ms = vraw,
        sky_terms_ms = theta, sky_terms_applied_ms = applied, arc_sky_offset_ms = Osky,
        sky_nlines = nlines, sky_anchor_fpi_ms = fpianchor,
        fiber_sky_zp_ms = fz, fiber_sky_stretch_ms = fs, pos_cluster = clus, pos_offset_px = pos_off_px,
        fpifib_tilt_px = tilt, fpifib_n_used = nfpi_used, dv_pass1_ms = dvp1, pass1_zp_ms = p1zp,
        pass1_shape_ms = p1shape, pass1_fiber_rsig_ms = p1rs,
        chip_arc_minus_comb_ms = C === nothing ? fill(NaN, 3) : C.chip_arc_minus_comb_ms,
        comb_n_mode_slip = C === nothing ? -1 : C.n_slip, noarc_fiber_zp_ms = fibz)
    mv(tmp, out; force = true)
    cnt = Dict{String, Int}()
    for w in wtype
        cnt[w] = get(cnt, w, 0) + 1
    end
    return refit_write_status(outroot, "assemble", tele, mjd, "ok",
        join(["$(k)=$(v)" for (k, v) in sort(collect(cnt))], ","))
end
