# Nights WITHOUT a usable arc fit (arc stage `noarc` / `unusable`): the within-night fallback ladder of the
# archive refit (`RefitConfig.noarc_route = "ladder"`). Every rung uses only the night's own data plus stored,
# survey-level calibration files (AKS 2026-11-19: "We should have a fallback plan, just not one that is bad and
# uses neighbors"; no step reads another night).
#
#   sky stage  (`refit_night_sky`): the legacy sky-only chain computed fresh on this run's sky-line peaks:
#              per object exposure `get_and_save_sky_wavecal` (chip scales from the chip-transform file) ->
#              `get_ave_night_wave_soln` (night average, measured night frame) -> `get_sky_dither_per_fiber`
#              per exposure. Product <outroot>/sky/skyNightRefit_<tele>_<mjd>.h5 (parameters only, ~1 MB).
#   rung 1     FPS night with an FPI comb (pass-1 gate passed): comb integer orders seeded by the night's own
#              sky solution, cavity (L, eps) fixed at the stored global cavity -> the normal FPI-fiber route.
#              (A seed within 1/2 FSR, ~32 km/s, is enough; a FIXED stored seed is not: cavity_first agent.)
#   rung 2     other exposures with sky lines: per-fiber shape from the stored comb shape template,
#              shifted to the exposure's own pooled sky lines minus the global arc-sky offset (as the arc route).
#   rung 3     sky only: the sky stage's per-exposure solution, shifted the same way (level = sky - O, as every
#              sky-shifted route); exposures without sky lines: the night-average sky solution minus O.
# Without the template file rung 2 becomes rung 3; without the cavity file rung 1 is skipped; without sky
# peaks nothing is possible (no solution, NaN coefficients). Bits 2^39-2^41.

const COMB_ANCHOR_SKYSEED = "skyseed"

"""True when the night has no usable arc fit in `outroot` (or `cfg.ignore_arcs`, validation only)."""
function refit_night_is_noarc(outroot, tele, mjd, cfg::RefitConfig)
    cfg.ignore_arcs && return true
    arcf = refit_path(outroot, "arc", tele, mjd)
    isfile(arcf) || return true
    return !refit_arc_usable(load_refit_arc(arcf))
end

"""
    refit_sky_chain(tele, mjd, skydir, objs; chip_transform_file = "", almanac_file = "", workroot = tempdir())

The legacy within-night sky-only wavelength chain on the sky-line peak files `skydir/skyLinePeaks_*` of the object
exposures `objs` (exposure numbers). The peak files are symlinked into a private scratch directory (the legacy
functions write next to their inputs), deleted on return. Steps (unchanged legacy code): per exposure
`get_and_save_sky_wavecal` (cporder 0 + the chip-transform file's chip scales), `get_ave_night_wave_soln`
(fit_dither, measured night frame; the commanded dither labels are only a cross-check there), and
`get_sky_dither_per_fiber` (dporder 1) of every exposure against the night average.
Returns `nothing` without any usable exposure, else (lin, nl, expnum, dither (N_FIBERS x 2 x n), ok, n_lines,
frame_shift_px).
"""
function refit_sky_chain(tele, mjd, skydir, objs; chip_transform_file::AbstractString = "",
        almanac_file::AbstractString = "", workroot::AbstractString = tempdir())
    isdir(skydir) || return nothing
    wd = mktempdir(workroot; prefix = "skychain_$(tele)_$(mjd)_")
    try
        ctc = isempty(chip_transform_file) ? nothing : load_chip_transform_cal(chip_transform_file)
        peaks = String[]; wfn = String[]; ee = Int[]; nlin = Int[]
        for e in objs
            src = [joinpath(skydir, "skyLinePeaks_$(tele)_$(mjd)_$(_arc_expstr(e))_$(ch)_object.h5") for ch in CHIP_LIST]
            all(isfile, src) || continue
            foreach(s -> symlink(s, joinpath(wd, basename(s))), src)
            pf = joinpath(wd, basename(src[1]))
            w = try
                get_and_save_sky_wavecal(pf; checkpoint_mode = "clobber", chip_transform_cal = ctc)
            catch err
                @warn "sky chain: exposure fit failed" tele mjd e exception = err
                nothing
            end
            w === nothing && continue
            ux = ingest_skyLines_exp(pf)[1]
            push!(peaks, pf); push!(wfn, w); push!(ee, e); push!(nlin, count(r -> any(isfinite, r), eachrow(ux)))
        end
        isempty(wfn) && return nothing
        dpx = get_exp_dithpix_dict(almanac_file, tele, mjd)
        lin, nl, _, anchor = get_ave_night_wave_soln(copy(wfn); fit_dither = true, exp_dithpix = dpx,
            dithpix_direction = get(NIGHT_FRAME_DITHPIX_DIRECTION, String(tele), 1))
        ndims(lin) == 3 && (lin = lin[:, :, 1])      # nanmedian over exposures keeps a trailing singleton
        ndims(nl) == 3 && (nl = nl[:, :, 1])
        dith = fill(NaN, N_FIBERS, 2, length(ee))
        ok = falses(length(ee))
        for (k, pf) in enumerate(peaks)
            r = try
                get_sky_dither_per_fiber(pf, lin, nl; dporder = 1, max_offset = 1.0)
            catch err
                @warn "sky chain: dither fit failed" tele mjd ee[k] exception = err
                nothing
            end
            r === nothing && continue
            dith[:, :, k] .= r[1]
            ok[k] = all(isfinite, r[1])
        end
        return (; lin, nl, expnum = ee, dither = dith, ok, n_lines = nlin,
            frame_shift_px = anchor === nothing ? NaN : Float64(anchor.shift_pix))
    finally
        rm(wd; recursive = true, force = true)
    end
end

"""Wavelength (Å) of the sky-chain solution (lin, nl; per-fiber dither `d = (offset, scale)`) for fiber f, chip c at
pixels `px` (same transform chain as `get_sky_dither_per_fiber`; `d = (0, 1)` = the night average)."""
function refit_sky_eval(lin::AbstractMatrix, nl::AbstractMatrix, d, f::Integer, c::Integer, px)
    wporder = size(lin, 2) - 1
    cporder = size(nl, 2) ÷ 2 - 1
    cpp = zeros(Float64, N_CHIPS, cporder + 1)
    cporder > 0 && (cpp[2, 2] = 1.0)
    params2ChipPolyParams!(cpp, nl[f, :], cporder)
    x = (collect(Float64, px) .- (N_XPIX ÷ 2)) ./ N_XPIX
    xt = transform_x_chips(x, cpp[c, :])
    xt = transform_x_chips(xt, [Float64(d[1]), Float64(d[2])])
    return positional_poly_mat(xt, porder = wporder) * lin[f, :]
end

"Stage: the night's sky-only chain (only on nights without a usable arc fit, or `cfg.ignore_arcs`)."
function refit_night_sky(tele, mjd; apred_dir, outroot, cfg::RefitConfig, almanac_file::AbstractString = "",
        workroot::AbstractString = tempdir())
    out = refit_path(outroot, "sky", tele, mjd)
    isfile(out) && return "done"
    refit_night_is_noarc(outroot, tele, mjd, cfg) ||
        return refit_write_status(outroot, "sky", tele, mjd, "notneeded")
    mdir = joinpath(apred_dir, string(mjd))
    skydir = isempty(cfg.sky_peaks_dir) ? mdir : joinpath(cfg.sky_peaks_dir, string(mjd))
    objs = [e for (e, k) in refit_night_exposures(mdir, tele, mjd) if k == "object"]
    isempty(objs) && return refit_write_status(outroot, "sky", tele, mjd, "nosky", "no object exposures")
    S = refit_sky_chain(tele, mjd, skydir, objs; chip_transform_file = cfg.chip_transform_file,
        almanac_file, workroot)
    S === nothing && return refit_write_status(outroot, "sky", tele, mjd, "nosky", "no object exposure with a sky-line fit")
    mkpath(dirname(out))
    tmp = out * ".tmp$(getpid())"
    meta = Dict{String, Any}("product" => "skyNightRefit", "tele" => String(tele), "mjd" => Int(mjd),
        "format_version" => REFIT_FORMAT_VERSION, "sky_peaks_dir" => skydir,
        "chip_transform_fingerprint" => refit_fingerprint(cfg.chip_transform_file),
        "method" => "legacy sky-only chain: get_and_save_sky_wavecal per exposure -> get_ave_night_wave_soln " *
                    "(fit_dither) -> get_sky_dither_per_fiber (dporder 1); evaluate with refit_sky_eval")
    for (k, v) in refit_config_dict(cfg)
        meta["cfg_" * k] = v
    end
    safe_jldsave(tmp, meta; night_linParams = S.lin, night_nlParams = S.nl, exp_expnum = S.expnum,
        ditherParams = S.dither, exp_ok = Vector{Bool}(S.ok), exp_n_lines = S.n_lines,
        frame_shift_px = S.frame_shift_px)
    mv(tmp, out; force = true)
    return refit_write_status(outroot, "sky", tele, mjd, "ok", "$(count(S.ok)) of $(length(objs)) object exposures")
end

function load_refit_sky(fname::AbstractString)
    h5open(fname, "r") do f
        (; lin = read(f["night_linParams"]), nl = read(f["night_nlParams"]), expnum = Int.(read(f["exp_expnum"])),
            dither = read(f["ditherParams"]), ok = Bool.(read(f["exp_ok"])))
    end
end

"""Rung 1: FPS comb of a night without arcs. Integer orders from the night's own sky solution (`lin0`, `nl0`),
cavity (L, eps) FIXED at the stored global cavity. QA (mode slips only) against the night-average sky grid `skyG`;
the sky solution sits at the night frame, not at the comb frame, so the arc - comb level QA is not meaningful
here (WAVE_COMB_ARC_INCONSISTENT is never set by this rung)."""
function refit_skyseed_comb(fpi_fnames, lin0, nl0, L::Real, eps::Real; skyG = nothing, pass1_comb = nothing,
        n_sigma = 4, max_iter = 2)
    R = comb_fpi_fit_core(fpi_fnames, lin0, nl0; cporder = 1, wporder = 4, dporder = 2,
        n_sigma = n_sigma, max_ang_sigma = 0.2, max_iter = max_iter, verbose = false,
        cavity_override = L, m_offset_override = eps, fix_LE = true)
    R === nothing && return (; ok = false, reason = "sky-seeded comb fit returned nothing")
    return refit_comb_qa(R, skyG; anchor = COMB_ANCHOR_SKYSEED, pass1_comb = pass1_comb, comb_arc_qa_ms = Inf,
        extra = (; arc_cluster = 0, arc_cluster_offset_px = NaN, cluster_offsets_px = Float64[],
            arc_state = Float64[], n_seed_fibers_filled = 0, frame_offset_verified = false, neighbour_mjd = -1,
            anchor_bits = WAVE_NOARC_COMB_SKYSEED, anchor_seq = 0, anchor_dt_h = NaN,
            fnames_ordered = collect(fpi_fnames), anchor_morning = false,
            anchor_kind = "noarc_skyseed:stored_cavity", anchor_rank = 0, anchor_dt_sci_h = NaN, anchor_quality = -1))
end

"""Rung 2 model: template wavelength of fiber f, chip c at pixel p with per-exposure chip shifts `delta` (px)
added to the template's own per-fiber displacement (coefficient P0); no per-fiber fit (no arcs)."""
refit_template_noarc(Tp::RefitTemplate, delta, f, c, p) =
    refit_template_lambda(Tp, (Tp.T[1, 1, f] + delta[1], Tp.T[1, 2, f] + delta[2], Tp.T[1, 3, f] + delta[3], 0.0),
        f, c, p)

"Fibers with a complete template (all P0..P6 finite on every chip)."
refit_template_ok(Tp::RefitTemplate) = [all(isfinite, view(Tp.T, :, :, f)) for f in 1:N_FIBERS]

"""Per-chip pixel shift (px) that moves the template model by the pooled sky offsets `zp` (m/s); a chip without a
pooled value takes the median shift of the others (NaN if none)."""
function refit_template_shift_px(Tp::RefitTemplate, zp::AbstractVector)
    d = fill(NaN, 3)
    for c in 1:3
        isfinite(zp[c]) || continue
        lam, disp = refit_ref_wave_disp(Tp, 150, c, 1024.0)
        d[c] = zp[c] / REFIT_C_MPS * lam / disp
    end
    fin = filter(isfinite, d)
    isempty(fin) && return d
    m = median(fin)
    return [isfinite(x) ? x : m for x in d]
end
