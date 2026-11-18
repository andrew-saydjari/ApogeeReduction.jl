# Arc-lamp + FPI wavelength REFIT, design stages 3-6 (ARCLAMP_PRODUCTION_DESIGN.md, revision 2,
# with the PI decisions of 2026-11-17 that override it):
#
#   stage 3  arc-anchored comb per night: the comb's integer orders are assigned from the arc
#            solution at the first comb exposure's MEASURED position and its cavity (L, eps) is
#            re-initialised from the arc per night ("ok, let's keep the cavity refit per night"),
#            cavity only, no per-chip zero points; the per-chip arc - comb residual is QA.
#   stage 4  FPI-fiber per-exposure transform (FPS era): the two FPI fibers of every exposure,
#            measured against the night comb, give the exposure's position and drift.
#   stage 5  pooled sky per-exposure transform (no FPI fibers): a robust, line-balanced pooled
#            OH-line offset per chip, anchored to the arcs (never an absolute sky level).
#   plate    comb SHAPE template (file + rebuild script, `scripts/refit/build_comb_shape_template.jl`);
#            without the file the arc shape is used and flagged (WAVE_SHAPE_ARC).
#   stage 6  assembly per exposure: absolute + relative + per-exposure terms, `wavecal_type =
#            <absolute>_<relative>`, `bitmsk_wavecal` (design section 2.4), QA (per-chip arc - comb
#            residual, cross-checks against the sky and against pass 1).
#
# Everything undecided is a field of `RefitConfig` with a documented default, so a variant is a
# rerun of the affected stage only. Inputs are opened read-only. Recorded dither positions are
# never used to place a frame: every position is measured (arc sequences, FPI fibers, sky lines).
# Nothing in the production pipeline calls this file.

using HDF5
using CSV
using DataFrames
using LinearAlgebra
using Statistics: median, mean, quantile

const REFIT_C_MPS = 299792458.0
const REFIT_PXG = 65:20:1985                       # 97-point flat pixel grid (HOLDOUT_RV / AV2)
const REFIT_FORMAT_VERSION = 1

# ------------------------------------------------------------------------------------------
# per-exposure wavelength bitmask (design section 2.4; bits 2^0-2^17 as in the design table)
# ------------------------------------------------------------------------------------------
const WAVE_SKY = 1 << 0
const WAVE_NIGHT_AVE = 1 << 1
const WAVE_ROUGH = 1 << 2
const WAVE_FPI_GATED_OFF = 1 << 3
const WAVE_NIGHTAVE_POS_MISMATCH = 1 << 4
const WAVE_ABS_ARC_WINDOW = 1 << 5
const WAVE_REL_ARC_TIME = 1 << 6
const WAVE_REL_ARC_NEAREST = 1 << 7
const WAVE_FPIFIB_UNAVAILABLE = 1 << 8
const WAVE_ARC_SEQ_GATED = 1 << 9
const WAVE_ARC_ALL_GATED = 1 << 10
const WAVE_ARC_UNGATED = 1 << 11
const WAVE_ARC_MULTIPOS = 1 << 12
const WAVE_COMB_MODE_SLIP = 1 << 13
const WAVE_LEVEL_INCONSISTENT = 1 << 14
const WAVE_HIGHORDER_EVENT = 1 << 15       # needs the cross-night pass; not set per night
const WAVE_CHIPSCALE_FREE = 1 << 16
const WAVE_FPIFIB_ROTATION = 1 << 17
# refit-only bits (proposed here; not yet in the design table)
const WAVE_PASS1_COPY = 1 << 18            # no new solution possible: the pass-1 adopted solution is copied
const WAVE_SHAPE_ARC = 1 << 19             # no comb shape (no comb, no template file): arc shape used
const WAVE_SKY_ANCHOR_NIGHT = 1 << 20      # pooled-sky tweak anchored on the whole night (no exposure near an arc)
const WAVE_POS_UNMATCHED = 1 << 21         # measured position matched no arc position cluster
const WAVE_FPIFIB_ONE_FIBER = 1 << 22      # only one FPI fiber usable (no across-slit term)
const WAVE_COMB_ARC_INCONSISTENT = 1 << 23 # |arc - comb| per chip above the QA limit (no slip)
const WAVE_ARC_QUARTZ_PATTERN = 1 << 24    # an arc sequence right after a same-dither quartz flat (no
                                           # commanded move) failed the 0.02 px gate (Q7 pattern)

const WAVE_BIT_NAMES = [
    "WAVE_SKY", "WAVE_NIGHT_AVE", "WAVE_ROUGH", "WAVE_FPI_GATED_OFF", "WAVE_NIGHTAVE_POS_MISMATCH",
    "WAVE_ABS_ARC_WINDOW", "WAVE_REL_ARC_TIME", "WAVE_REL_ARC_NEAREST", "WAVE_FPIFIB_UNAVAILABLE",
    "WAVE_ARC_SEQ_GATED", "WAVE_ARC_ALL_GATED", "WAVE_ARC_UNGATED", "WAVE_ARC_MULTIPOS",
    "WAVE_COMB_MODE_SLIP", "WAVE_LEVEL_INCONSISTENT", "WAVE_HIGHORDER_EVENT", "WAVE_CHIPSCALE_FREE",
    "WAVE_FPIFIB_ROTATION", "WAVE_PASS1_COPY", "WAVE_SHAPE_ARC", "WAVE_SKY_ANCHOR_NIGHT",
    "WAVE_POS_UNMATCHED", "WAVE_FPIFIB_ONE_FIBER", "WAVE_COMB_ARC_INCONSISTENT",
    "WAVE_ARC_QUARTZ_PATTERN"]

"Night-level arc flags of the stage-1 product mapped onto the per-exposure bits 2^9-2^12, 2^16."
function refit_arc_night_bits(night_flags::Integer, seq_flags::AbstractVector{<:Integer})
    b = 0
    any(s -> (s & ARCSEQ_GATE_FAIL) != 0, seq_flags) && (b |= WAVE_ARC_SEQ_GATED)
    (night_flags & ARCNIGHT_ALL_GATED) != 0 && (b |= WAVE_ARC_ALL_GATED)
    (night_flags & ARCNIGHT_NO_GATE) != 0 && (b |= WAVE_ARC_UNGATED)
    (night_flags & ARCNIGHT_MULTIPOS) != 0 && (b |= WAVE_ARC_MULTIPOS)
    (night_flags & ARCNIGHT_CHIPSCALE_FREE) != 0 && (b |= WAVE_CHIPSCALE_FREE)
    return b
end

"Pass-1 `wavecal_type` -> (new-scheme name, bits) for exposures that keep the pass-1 solution."
function refit_pass1_type(wt::AbstractString)
    wt == "waveCalNightfpiDither" && return ("sky_skydith", WAVE_SKY)
    wt == "waveCalNightskyDither" && return ("sky_skydith", WAVE_SKY)
    startswith(wt, "error_night_ave") &&
        return ("sky_nightpos", WAVE_SKY | WAVE_NIGHT_AVE | WAVE_NIGHTAVE_POS_MISMATCH)
    wt == "error_fixed_fallback" && return ("none_rough", WAVE_ROUGH)
    return ("none_rough", WAVE_ROUGH)
end

# ------------------------------------------------------------------------------------------
# configuration: every undecided choice is a switch with a documented default
# ------------------------------------------------------------------------------------------
"""
    RefitConfig

Switches of the refit (defaults = the design's recommendation where it has one):

  * `neon` ("none"): arc line list neon handling, design L4 (OPEN-UNDER-TEST). "none" = shape-clean
    Th/U/Ar (`use_default`, the stage-1 default); "lab" = plus lab neon (`use_with_neon`);
    "lab_offset" = lab neon plus a measured per-lamp offset from `neon_offset_file` (learned set
    `neon_offset_set`); without the file "lab_offset" reduces to "lab" (recorded).
  * `chip_transform_file`: comb chip-scale calibration (PR #415 format); "" = free chip scale (flagged).
  * `exp_terms` ("chipzp_stretch_tilt"): per-exposure position transform from the FPI fibers, design
    D2/D3 (OPEN-UNDER-TEST): "shift" (one pixel shift, all chips), "chipzp" (per-chip zero point),
    "chipzp_stretch" (+ per-chip stretch; design D2'), "chipzp_stretch_tilt" (+ linear in fiber
    through the two FPI fibers; design D3).
  * `sky_fiber_term` ("none"): per-fiber sky term, design D6 (OPEN-UNDER-TEST; the design's earlier
    recommendation is marked not established): "none"; "pooled" (fiber-common per-exposure per-chip
    sky tweak applied in every era, arc-anchored); "static" (per-fiber night-median sky offset,
    fiber median removed); "full" (per-fiber per-exposure sky zero point + stretch, fiber median
    removed).
  * `prefps_time` ("regularised"): time dependence where no FPI fibers are usable (plate era, FPS
    exposures without FPI light), ARCFPI_VALIDATION Q4: "regularised" = option (c), the arc time
    model of the night's REFERENCE position cluster plus `reg_weight_<tele>` x the pooled-sky
    departure (APO 0.825, LCO 0.65, chosen on held-out nights); "sky" = option (a), a fixed arc
    reference state plus the full pooled-sky departure; "arc" = option (b), the arc time model of
    the exposure's own cluster, no sky. In (a) and (c) the step BETWEEN position clusters comes
    from the sky (median sky departure of the cluster's exposures), never from the arcs (the
    arc-measured step is biased on APO chip G: +59..+72 m/s zero point, -127..-133 m/s stretch).
  * `sky_exp_terms` ("chipzp_stretch_tilt"): terms of the pooled-sky per-exposure transform,
    per chip: "chipzp", "chipzp_stretch", "chipzp_stretch_tilt" (the APO G/B tilt moves in the
    plate era too: 93/45 m/s within a position against noise of 8).
  * `reg_weight_apo` (0.825), `reg_weight_lco` (0.65): weights of option (c).
  * `shape_template_file`, `shape_template_block` ("all"): comb shape template for nights without
    a comb (plate era); "" = arc shape (flagged WAVE_SHAPE_ARC).
  * `arcwin_nights` (15): neighbour-night window for the comb cavity when a night has no arc anchor.
  * `nsig` (8): arc and FPI-fiber peak detection threshold (stage-1 finder).
  * `sky_anchor_window_h` (1.5): exposures within this time of a used arc sequence anchor the
    pooled sky tweak to the arcs.
  * `pos_match_px` (0.25): position-cluster match radius (stage-1 `group_px`).
  * `xcheck_threshold_ms` (60): `WAVE_LEVEL_INCONSISTENT` threshold on the sky cross-check.
  * `comb_arc_qa_ms` (150): `WAVE_COMB_ARC_INCONSISTENT` threshold on |arc - comb| per chip.
  * `coeff_degree` (9): Legendre degree of the stored per-fiber, per-chip solution.
"""
Base.@kwdef struct RefitConfig
    neon::String = "none"
    neon_offset_file::String = ""
    neon_offset_set::String = "all"
    chip_transform_file::String = ""
    exp_terms::String = "chipzp_stretch_tilt"
    sky_fiber_term::String = "none"
    prefps_time::String = "regularised"
    sky_exp_terms::String = "chipzp_stretch_tilt"
    reg_weight_apo::Float64 = 0.825
    reg_weight_lco::Float64 = 0.65
    shape_template_file::String = ""
    shape_template_block::String = "all"
    arcwin_nights::Int = 15
    nsig::Float64 = 8.0
    sky_anchor_window_h::Float64 = 1.5
    pos_match_px::Float64 = 0.25
    xcheck_threshold_ms::Float64 = 60.0
    comb_arc_qa_ms::Float64 = 150.0
    coeff_degree::Int = 9
end

const REFIT_SWITCH_VALUES = Dict(
    :neon => ("none", "lab", "lab_offset"),
    :exp_terms => ("shift", "chipzp", "chipzp_stretch", "chipzp_stretch_tilt"),
    :sky_fiber_term => ("none", "pooled", "static", "full"),
    :prefps_time => ("regularised", "sky", "arc"),
    :sky_exp_terms => ("chipzp", "chipzp_stretch", "chipzp_stretch_tilt"))

"Check the switch values; throws an ArgumentError naming the bad field."
function refit_check_config(cfg::RefitConfig)
    for (k, vals) in REFIT_SWITCH_VALUES
        v = getfield(cfg, k)
        v in vals || throw(ArgumentError("RefitConfig.$(k) = \"$(v)\" not in $(vals)"))
    end
    cfg.coeff_degree >= 4 || throw(ArgumentError("coeff_degree must be >= 4"))
    return true
end

"All switch values as strings (recorded in every product)."
refit_config_dict(cfg::RefitConfig) =
    Dict{String, String}(string(k) => string(getfield(cfg, k)) for k in fieldnames(RefitConfig))

"Which stages a changed switch invalidates (for the plan; the driver uses the same table)."
const REFIT_SWITCH_STAGES = Dict(
    "neon" => "arc,comb,fpifib,assemble", "neon_offset_file" => "arc,comb,fpifib,assemble",
    "neon_offset_set" => "arc,comb,fpifib,assemble", "chip_transform_file" => "arc,comb,fpifib,assemble",
    "nsig" => "arc,comb,fpifib,assemble", "arcwin_nights" => "comb,fpifib,assemble",
    "pos_match_px" => "comb,fpifib,assemble", "comb_arc_qa_ms" => "comb,fpifib,assemble",
    "exp_terms" => "assemble", "sky_fiber_term" => "assemble", "prefps_time" => "assemble",
    "sky_exp_terms" => "assemble", "reg_weight_apo" => "assemble", "reg_weight_lco" => "assemble",
    "shape_template_file" => "assemble", "shape_template_block" => "assemble",
    "sky_anchor_window_h" => "assemble", "xcheck_threshold_ms" => "assemble",
    "coeff_degree" => "assemble")

"FNV-1a fingerprint of a file, or \"none\" for an empty path / missing file."
refit_fingerprint(path::AbstractString) =
    (isempty(path) || !isfile(path)) ? "none" : arclamp_file_fingerprint(path)

# ------------------------------------------------------------------------------------------
# calibration inputs (all optional; the code runs without them and records/flags that)
# ------------------------------------------------------------------------------------------
"""Comb chip scales (N_FIBERS x 2: b_R, b_B) for (tele, mjd) from a PR #415 chip-transform file,
or `nothing` (no file / no block). Minimal reader of the #415 format (`/<tele>/mjd_start`,
`mjd_end`, `scale`); to be replaced by #415's `chip_scales` once it is merged."""
function refit_chip_scales(path::AbstractString, tele::AbstractString, mjd::Integer)
    (isempty(path) || !isfile(path)) && return nothing
    h5open(path, "r") do f
        haskey(f, tele) || return nothing
        ms = read(f["$(tele)/mjd_start"])
        me = read(f["$(tele)/mjd_end"])
        S = read(f["$(tele)/scale"])
        k = findfirst(i -> ms[i] <= mjd <= me[i], eachindex(ms))
        k === nothing && return nothing
        return Matrix{Float64}(S[:, :, k])
    end
end

"""Per-lamp neon offsets (Dict "ThAr"/"UNe" => Å) for a telescope from the neon calibration file
(`learned/<tele>/<set>/lamp_offset_A`, ARCFPI_VALIDATION build_neon_offsets.jl format), or `nothing`."""
function refit_neon_offsets(path::AbstractString, tele::AbstractString, set::AbstractString)
    (isempty(path) || !isfile(path)) && return nothing
    h5open(path, "r") do f
        k = "learned/$(tele)/$(set)/lamp_offset_A"
        haskey(f, k) || return nothing
        v = read(f[k])
        return Dict("ThAr" => Float64(v[1]), "UNe" => Float64(v[2]))
    end
end

"The arc line list for a configuration and telescope; returns (list, neon_mode_used)."
function refit_linelist(cfg::RefitConfig, tele::AbstractString)
    cfg.neon == "none" && return load_arclamp_linelist(selection = "use_default"), "none"
    if cfg.neon == "lab_offset"
        off = refit_neon_offsets(cfg.neon_offset_file, tele, cfg.neon_offset_set)
        off === nothing && return load_arclamp_linelist(selection = "use_with_neon"), "lab(no offset file)"
        return load_arclamp_linelist(selection = "use_with_neon", neon_offset_A = off), "lab_offset"
    end
    return load_arclamp_linelist(selection = "use_with_neon"), "lab"
end

# ------------------------------------------------------------------------------------------
# stage-1 arc night product: loader and arc states
# ------------------------------------------------------------------------------------------
"""Stored stage-1 arc night fit, with the final-fit layout rebuilt. Sequence arrays with suffix
`_u` are over the sequences used in the final fit (in almanac order)."""
struct RefitArc
    fname::String
    TH::Matrix{Float64}
    okf::Vector{Bool}
    g::Vector{Float64}
    lay::ArcGlobalLayout
    used::Vector{Int}
    grp_u::Vector{Int}
    tau_u::Vector{Float64}
    tseq_u::Vector{Float64}
    clus_u::Vector{Int}
    tmean::Float64
    night_flags::Int
    seq_flags::Vector{Int}
    chip_scale_source::String
    lines::NamedTuple
end

function load_refit_arc(fname::AbstractString)
    h5open(fname, "r") do f
        TH = read(f["fiber_params"])
        okf = Bool.(read(f["fiber_ok"]))
        g = Float64.(read(f["global_params"]))
        used = findall(Bool.(read(f["seq_used"])))
        has2 = Bool(read(f["has_pos2"]))
        grp = Int.(read(f["seq_pos_group"]))
        tau = Float64.(read(f["seq_tau"]))
        grp_u = has2 ? grp[used] : ones(Int, length(used))
        lay = ArcGlobalLayout(length(used), grp_u, tau[used])
        lay.n == length(g) || error("$(fname): stored global terms do not match the rebuilt layout")
        css = haskey(f, "metadata/chip_scale_source") ? String(read(f["metadata/chip_scale_source"])) : "unknown"
        lines = (seq = Int.(read(f["lines_seq"])), fiber = Int.(read(f["lines_fiber"])),
            chip = Int.(read(f["lines_chip"])), id = Int.(read(f["lines_id"])),
            pix = Float64.(read(f["lines_pix"])), pixvar = Float64.(read(f["lines_pixvar"])),
            keep = Bool.(read(f["lines_keep"])))
        RefitArc(String(fname), TH, okf, g, lay, used, grp_u, tau[used],
            Float64.(read(f["seq_mjd_mid"]))[used], Int.(read(f["seq_pos_cluster"]))[used],
            Float64(read(f["tau_mjd_mean"])), Int(read(f["night_flags"])),
            Int.(read(f["seq_flags"])), css, lines)
    end
end

"An arc night fit is usable: at least one sequence in the final fit and >= 100 fitted fibers."
refit_arc_usable(A::RefitArc) = !isempty(A.used) && count(A.okf) >= 100

"Position clusters present among the used sequences."
refit_arc_clusters(A::RefitArc) = sort(unique(A.clus_u))

"Linear interpolation in time of (d, t) between sequences (constant outside their time range)."
function refit_interp_dt(ts::AbstractVector, D::AbstractVector, T::AbstractVector, t::Real)
    o = sortperm(ts)
    ts = ts[o]
    D = D[o]
    T = T[o]
    n = length(ts)
    (n == 1 || t <= ts[1]) && return D[1], T[1]
    t >= ts[n] && return D[n], T[n]
    k = searchsortedlast(ts, t)
    w = (t - ts[k]) / (ts[k + 1] - ts[k])
    return (1 - w) * D[k] + w * D[k + 1], (1 - w) * T[k] + w * T[k + 1]
end

"""Arc state (grpv, d, t, tau) at time `t` (MJD) for position cluster `cl`.
`mode = :time`: (d, t) interpolated in time between the cluster's sequences, time terms at the
exposure time (clamped to the sequence range) = the arc time model. `mode = :nearest`: the state
of the cluster's sequence nearest in time. Returns `nothing` if the cluster has no sequence."""
function refit_arc_state(A::RefitArc, t::Real, cl::Integer; mode::Symbol = :time)
    js = [j for j in eachindex(A.used) if A.clus_u[j] == cl]
    isempty(js) && return nothing
    if mode == :nearest || !isfinite(t)
        j = isfinite(t) ? js[argmin(abs.(A.tseq_u[js] .- t))] : js[1]
        d, tt = arcfit_seq_state(A.g, A.lay, j)
        return (A.grp_u[j], d, tt, A.tau_u[j])
    end
    Ds = [arcfit_seq_state(A.g, A.lay, j)[1] for j in js]
    Ts = [arcfit_seq_state(A.g, A.lay, j)[2] for j in js]
    d, tt = refit_interp_dt(A.tseq_u[js], Ds, Ts, t)
    lo, hi = extrema(A.tau_u)
    return (A.grp_u[js[1]], d, tt, clamp((t - A.tmean) * 2.0, lo, hi))
end

"Shift of the arc state on chip c of fiber f as A + B x (x = (pix-1024)/2048), in z units."
function refit_arc_shift_AB(A::RefitArc, f::Integer, c::Integer, st)
    grpv, d, t, tau = st
    lay = A.lay
    a = d + t * arcfit_phi(f)
    b = 0.0
    if lay.has2 && grpv == 2
        ia = lay.idx_alpha[c]
        ib = lay.idx_beta[c]
        ia > 0 && (a += A.g[ia])
        ib > 0 && (b += A.g[ib])
    end
    if lay.timeok
        iu = lay.idx_u[c]
        iv = lay.idx_v[c]
        iu > 0 && (a += A.g[iu] * tau)
        iv > 0 && (b += A.g[iv] * tau)
    end
    return a, b
end

"Arc model wavelengths (px, fibers, chips) at a state on pixel vector `px`."
function refit_arc_grid(A::RefitArc, st, px::AbstractVector{<:Real})
    grpv, d, t, tau = st
    G = fill(NaN, length(px), N_FIBERS, N_CHIPS)
    pxv = collect(Float64, px)
    for f in 1:N_FIBERS
        A.okf[f] || continue
        th = view(A.TH, f, :)
        for c in 1:N_CHIPS
            G[:, f, c] .= arcfit_eval(th, A.g, A.lay, f, c, pxv; grpv = grpv, d = d, t = t, tau = tau)
        end
    end
    return G
end

"""
    refit_arc_lin_nl(A, st) -> (linParams N_FIBERS x 5, nlParams N_FIBERS x 4, n_filled)

The arc solution at state `st` rewritten EXACTLY in the production comb form (quartic in the
green-chip coordinate, chip G the identity, chips R/B linear: `transform_x_chips` with cporder
1), so it can seed `comb_fpi_fit_core`. With xs_c = O_c + S_c x the arc's argument on chip c and
lambda = P(xs), the production form is Q(u) = P(O_G + S_G u) with o_c = (O_c - O_G)/S_G and
s_c = S_c/S_G. Fibers without an arc fit take the nearest good fiber's parameters (seed only).
"""
function refit_arc_lin_nl(A::RefitArc, st)
    lin = fill(NaN, N_FIBERS, 5)
    nl = fill(NaN, N_FIBERS, 4)
    for f in 1:N_FIBERS
        A.okf[f] || continue
        th = A.TH[f, :]
        all(isfinite, th) || continue
        aR, bR = refit_arc_shift_AB(A, f, 1, st)
        aG, bG = refit_arc_shift_AB(A, f, 2, st)
        aB, bB = refit_arc_shift_AB(A, f, 3, st)
        O = (th[1] + aR, aG, th[3] + aB)
        S = (th[2] + bR, 1.0 + bG, th[4] + bB)
        q = Polynomial(th[5:9])(Polynomial([O[2], S[2]])).coeffs
        lin[f, :] .= 0.0
        lin[f, 1:min(5, length(q))] .= q[1:min(5, length(q))]
        nl[f, :] .= ((O[1] - O[2]) / S[2], S[1] / S[2], (O[3] - O[2]) / S[2], S[3] / S[2])
    end
    good = findall(f -> all(isfinite, view(lin, f, :)), 1:N_FIBERS)
    isempty(good) && return lin, nl, N_FIBERS
    nfill = 0
    for f in 1:N_FIBERS
        f in good && continue
        j = good[argmin(abs.(good .- f))]
        lin[f, :] .= lin[j, :]
        nl[f, :] .= nl[j, :]
        nfill += 1
    end
    return lin, nl, nfill
end

"Production-form (cporder 1 or 0) solution of fiber f, chip c at pixels px."
function refit_prod_wave(lin::AbstractMatrix, nl::AbstractMatrix, f::Integer, c::Integer, px)
    x = (px .- (N_XPIX ÷ 2)) ./ N_XPIX
    cp = size(nl, 2) ÷ 2 - 1
    if c == 2
        xt = x
    else
        o = c == 1 ? 0 : cp + 1
        xt = cp == 0 ? nl[f, o + 1] .+ x : nl[f, o + 1] .+ nl[f, o + 2] .* x
    end
    return positional_poly_mat(xt, porder = size(lin, 2) - 1) * lin[f, :]
end

# ------------------------------------------------------------------------------------------
# stage 3: arc-anchored comb
# ------------------------------------------------------------------------------------------
const COMB_ANCHOR_ARC = "arc"
const COMB_ANCHOR_ARCWIN = "arcwin"

"""Per fiber-chip comb - arc difference in units of the local FSR (λ²/2L), median over the pixel
grid, and the flat-mean velocity difference arc - comb (m/s). Returns (dfsr, dv) each N_FIBERS x 3."""
function refit_comb_vs_arc(lin, nl, arcG::AbstractArray, L::Real; px = REFIT_PXG)
    dfsr = fill(NaN, N_FIBERS, N_CHIPS)
    dv = fill(NaN, N_FIBERS, N_CHIPS)
    pxv = collect(Float64, px)
    for f in 1:N_FIBERS, c in 1:N_CHIPS
        a = view(arcG, :, f, c)
        all(isfinite, a) || continue
        lc = refit_prod_wave(lin, nl, f, c, pxv)
        all(isfinite, lc) || continue
        fsr = lc .^ 2 ./ (2L)
        dfsr[f, c] = median((lc .- a) ./ fsr)
        dv[f, c] = mean((a .- lc) ./ lc) * REFIT_C_MPS
    end
    return dfsr, dv
end

"""Mode-slip check C3 (i): fiber-chips whose comb - arc difference is within `tol` FSR of a
NONZERO integer. Returns (count, mask N_FIBERS x 3)."""
function refit_mode_slips(dfsr::AbstractMatrix; tol = 0.2)
    m = falses(size(dfsr))
    for i in eachindex(dfsr)
        x = dfsr[i]
        isfinite(x) || continue
        k = round(x)
        (k != 0 && abs(x - k) < tol) && (m[i] = true)
    end
    return count(m), m
end

"""Integer-order change per fiber-chip between two peak-integer tables of the SAME peak ingest
(median over the fiber-chip's peaks). Returns an Int16 matrix N_FIBERS x 3 (0 = unchanged)."""
function refit_order_change(pk_new::AbstractMatrix, pk_old::AbstractMatrix, chipInt::AbstractMatrix,
        ok::AbstractMatrix)
    out = zeros(Int16, N_FIBERS, N_CHIPS)
    size(pk_new) == size(pk_old) || return fill(Int16(-999), N_FIBERS, N_CHIPS)
    for f in 1:N_FIBERS, c in 1:N_CHIPS
        m = [i for i in axes(pk_new, 1) if chipInt[i, f] == c && ok[i, f]]
        isempty(m) && continue
        out[f, c] = Int16(round(median(pk_new[m, f] .- pk_old[m, f])))
    end
    return out
end

"""
    refit_arc_anchored_comb(fpi_fnames, A, frame_ref, t_frame; pass1_comb = nothing, ...)

Stage 3 for a night with an arc anchor. `fpi_fnames`: the night's FPI peak files (R chip names,
read-only), first = the frame exposure. `frame_ref`: an `ArcObjRef` slice of a MEASURED solution
at the frame exposure (pass-1 per-exposure comb solution), used only to pick the arc position
cluster at the frame exposure (|offset| < `pos_match_px`); `t_frame`: its mid-exposure MJD.
The arc time model of that cluster at `t_frame` seeds the comb: integer orders from the arc,
(L, eps) by least squares to the arc, then the production comb loop (L, eps refit per night).
Returns a NamedTuple (`ok = false` with a `reason` if it cannot anchor).
"""
function refit_arc_anchored_comb(fpi_fnames::AbstractVector{<:AbstractString}, A::RefitArc,
        frame_ref, t_frame::Real; pass1_comb = nothing, pos_match_px = 0.25, comb_arc_qa_ms = 150.0,
        n_sigma = 4, max_iter = 2)
    cls = refit_arc_clusters(A)
    offs = fill(NaN, length(cls))
    sts = Any[nothing for _ in cls]
    for (k, cl) in enumerate(cls)
        st = refit_arc_state(A, t_frame, cl; mode = :time)
        st === nothing && continue
        sts[k] = st
        frame_ref === nothing && continue
        offs[k] = arclamp_offset_to_ref(A.TH, A.okf, A.g, A.lay, st[1], st[2], st[3], st[4], frame_ref)
    end
    kbest = 0
    if frame_ref === nothing
        length(cls) == 1 && (kbest = 1)       # one cluster: unverified, recorded
    elseif any(isfinite, offs)
        kb = argmin([isfinite(o) ? abs(o) : Inf for o in offs])
        abs(offs[kb]) < pos_match_px && (kbest = kb)
    end
    kbest == 0 && return (; ok = false, reason = "no arc position cluster at the frame exposure",
        cluster_offsets_px = offs)
    st = sts[kbest]
    lin0, nl0, nfill = refit_arc_lin_nl(A, st)
    R = comb_fpi_fit_core(fpi_fnames, lin0, nl0; cporder = 1, wporder = 4, dporder = 2,
        n_sigma = n_sigma, max_ang_sigma = 0.2, max_iter = max_iter, verbose = false,
        seed_LE_from_initial = true)
    R === nothing && return (; ok = false, reason = "comb fit returned nothing", cluster_offsets_px = offs)
    arcG = refit_arc_grid(A, st, REFIT_PXG)
    return refit_comb_qa(R, arcG; anchor = COMB_ANCHOR_ARC, pass1_comb = pass1_comb,
        comb_arc_qa_ms = comb_arc_qa_ms, extra = (; arc_cluster = cls[kbest],
            arc_cluster_offset_px = offs[kbest], cluster_offsets_px = offs, arc_state = collect(Float64, st),
            n_seed_fibers_filled = nfill, frame_offset_verified = frame_ref !== nothing))
end

"""Stage 3 for a night WITHOUT an arc anchor: integer orders from the pass-1 seed (`lin0`, `nl0`,
the night-average sky solution), cavity (L, eps) FIXED at a neighbouring night's arc-anchored
values (`arcwin`). `arcG` (optional) is an arc grid at the frame for QA only."""
function refit_arcwin_comb(fpi_fnames, lin0, nl0, L::Real, eps::Real; pass1_comb = nothing,
        comb_arc_qa_ms = 150.0, n_sigma = 4, max_iter = 2, neighbour_mjd = -1)
    R = comb_fpi_fit_core(fpi_fnames, lin0, nl0; cporder = 1, wporder = 4, dporder = 2,
        n_sigma = n_sigma, max_ang_sigma = 0.2, max_iter = max_iter, verbose = false,
        cavity_override = L, m_offset_override = eps, fix_LE = true)
    R === nothing && return (; ok = false, reason = "comb fit returned nothing")
    return refit_comb_qa(R, nothing; anchor = COMB_ANCHOR_ARCWIN, pass1_comb = pass1_comb,
        comb_arc_qa_ms = comb_arc_qa_ms, extra = (; arc_cluster = 0, arc_cluster_offset_px = NaN,
            cluster_offsets_px = Float64[], arc_state = Float64[], n_seed_fibers_filled = 0,
            frame_offset_verified = false, neighbour_mjd = neighbour_mjd))
end

"QA of a comb fit result (C3 (i), arc - comb per chip, order changes vs pass 1)."
function refit_comb_qa(R, arcG; anchor, pass1_comb, comb_arc_qa_ms, extra)
    gp = R.good_peaks
    resid = filter(isfinite, R.resid_vec'[gp])
    rms_mA = isempty(resid) ? NaN : 1000 * sqrt(mean(resid .^ 2))
    if arcG !== nothing
        dfsr, dv = refit_comb_vs_arc(R.linParams, R.nlParams, arcG, R.cavity_size)
        nslip, slipmask = refit_mode_slips(dfsr)
        chip_dv = [_arc_nanmedian(dv[:, c]) for c in 1:N_CHIPS]
        chip_dv_rsig = [_arc_robstd(dv[:, c]) for c in 1:N_CHIPS]
        n_far = count(x -> isfinite(x) && abs(x) > 0.3 && abs(x - round(x)) >= 0.2, dfsr)
    else
        dfsr = fill(NaN, N_FIBERS, N_CHIPS)
        dv = fill(NaN, N_FIBERS, N_CHIPS)
        nslip, slipmask = 0, falses(N_FIBERS, N_CHIPS)
        chip_dv = fill(NaN, 3)
        chip_dv_rsig = fill(NaN, 3)
        n_far = 0
    end
    ordch = fill(Int16(-999), N_FIBERS, N_CHIPS)
    p1_L = NaN
    p1_eps = NaN
    if pass1_comb !== nothing
        ordch = refit_order_change(R.peak_ints, pass1_comb.peak_ints, R.fpi_line_chipInt,
            gp .& pass1_comb.used)
        p1_L = pass1_comb.L
        p1_eps = pass1_comb.eps
    end
    n_ord = count(x -> x != 0 && x != -999, ordch)
    inconsistent = any(x -> isfinite(x) && abs(x) > comb_arc_qa_ms, chip_dv)
    return merge((; ok = true, reason = "", anchor, R, rms_mA, dfsr, dv_arc_minus_comb = dv,
        chip_arc_minus_comb_ms = chip_dv, chip_arc_minus_comb_rsig_ms = chip_dv_rsig, n_slip = nslip,
        slipmask, n_far_fiberchips = n_far, order_change_vs_pass1 = ordch,
        n_fiberchip_order_changed = n_ord, pass1_L = p1_L, pass1_eps = p1_eps,
        comb_arc_inconsistent = inconsistent), extra)
end

"Pass-1 comb (waveCalFPI) integers, use flags and cavity, opened read-only; `nothing` if absent."
function refit_read_pass1_comb(fname::AbstractString)
    isfile(fname) || return nothing
    h5open(fname, "r") do f
        (; peak_ints = read(f["resid_peak_ints"]), used = Bool.(read(f["resid_used_in_fit"])),
            L = Float64(read(f["fpi_cavity_size"])), eps = Float64(read(f["fpi_m0_offset"])),
            lin = read(f["linParams"]), nl = read(f["nlParams"]))
    end
end

"Write the stage-3 night product `waveCalFPIArc_<tele>_<mjd>.h5` (atomic)."
function write_refit_comb(outname::AbstractString, Q, fpi_fnames; tele, mjd, cfg::RefitConfig,
        arc_fname = "", provenance = Dict{String, String}())
    R = Q.R
    meta = Dict{String, Any}("product" => "waveCalFPIArc", "stage" => 3, "tele" => String(tele),
        "mjd" => Int(mjd), "anchor" => Q.anchor, "arc_product" => String(arc_fname),
        "format_version" => REFIT_FORMAT_VERSION)
    for (k, v) in refit_config_dict(cfg)
        meta["cfg_" * k] = v
    end
    for (k, v) in provenance
        meta[k] = v
    end
    tmp = outname * ".tmp$(getpid())"
    safe_jldsave(tmp, meta; linParams = R.linParams, nlParams = R.nlParams,
        ditherParams = R.ditherParams, ditherPolyParams = R.ditherPolyParams,
        exposure_names = [basename(x) for x in fpi_fnames], fpi_cavity_size = R.cavity_size,
        fpi_m0_offset = R.m_offset, fpi_m0 = R.med_m0, seed_cavity_size = R.seed_cavity_size,
        seed_m0_offset = R.seed_m_offset, init_cavity_size = R.init_cavity_size,
        init_m0_offset = R.init_m_offset, first_dither_offset = R.first_dither_offset,
        comb_resid_rms_mA = Q.rms_mA, n_good_peaks = count(R.good_peaks),
        dfsr_comb_minus_arc = Q.dfsr, dv_arc_minus_comb_ms = Q.dv_arc_minus_comb,
        chip_arc_minus_comb_ms = Q.chip_arc_minus_comb_ms,
        chip_arc_minus_comb_rsig_ms = Q.chip_arc_minus_comb_rsig_ms, n_mode_slip = Q.n_slip,
        mode_slip_mask = Matrix{Bool}(Q.slipmask), n_far_fiberchips = Q.n_far_fiberchips,
        order_change_vs_pass1 = Q.order_change_vs_pass1,
        n_fiberchip_order_changed = Q.n_fiberchip_order_changed, pass1_cavity_size = Q.pass1_L,
        pass1_m0_offset = Q.pass1_eps, comb_arc_inconsistent = Q.comb_arc_inconsistent,
        arc_cluster = Q.arc_cluster, arc_cluster_offset_px = Q.arc_cluster_offset_px,
        cluster_offsets_px = Float64.(Q.cluster_offsets_px), arc_state = Float64.(Q.arc_state),
        n_seed_fibers_filled = Q.n_seed_fibers_filled,
        frame_offset_verified = Q.frame_offset_verified,
        neighbour_mjd = hasproperty(Q, :neighbour_mjd) ? Q.neighbour_mjd : -1)
    mv(tmp, outname; force = true)
    return outname
end

"Stored stage-3 comb: (lin, nl, L, eps, anchor, flags...)."
function load_refit_comb(fname::AbstractString)
    h5open(fname, "r") do f
        (; lin = read(f["linParams"]), nl = read(f["nlParams"]), L = Float64(read(f["fpi_cavity_size"])),
            eps = Float64(read(f["fpi_m0_offset"])), anchor = String(read(f["metadata/anchor"])),
            n_slip = Int(read(f["n_mode_slip"])), inconsistent = Bool(read(f["comb_arc_inconsistent"])),
            chip_arc_minus_comb_ms = Float64.(read(f["chip_arc_minus_comb_ms"])),
            exposure_names = String.(read(f["exposure_names"])))
    end
end

# ------------------------------------------------------------------------------------------
# stage 4: FPI fibers in every exposure (position and drift against the night comb)
# ------------------------------------------------------------------------------------------
"""Weighted line fit d = a + b t (t = (x - 1024)/1024) with iterative clipping (error floor
0.003 px), as in the FPI_FIBER_DITHER study. Returns (a, b, ea, eb, chi2r, n, rms)."""
function refit_fit_ab(x, d, err; clip = 4.0, floor = 0.003, iters = 5, minn = 8)
    t = (x .- 1024.0) ./ 1024.0
    e = sqrt.(err .^ 2 .+ floor^2)
    ok = isfinite.(d) .& isfinite.(e)
    count(ok) < minn && return (NaN, NaN, NaN, NaN, NaN, 0, NaN)
    med = median(d[ok])
    s = 1.4826 * median(abs.(d[ok] .- med))
    s = (isfinite(s) && s > 0) ? s : 0.01
    ok .&= abs.(d .- med) .< max(8 * s, 0.05)
    a = b = ea = eb = chi = NaN
    for _ in 1:iters
        count(ok) < minn && return (NaN, NaN, NaN, NaN, NaN, 0, NaN)
        w = 1 ./ e[ok] .^ 2
        X = hcat(ones(count(ok)), t[ok])
        M = X' * (w .* X)
        p = M \ (X' * (w .* d[ok]))
        C = inv(M)
        a, b = p
        ea, eb = sqrt(C[1, 1]), sqrt(C[2, 2])
        r = (d .- (a .+ b .* t)) ./ e
        chi = sum(abs2, r[ok]) / max(1, count(ok) - 2)
        sc = sqrt(max(chi, 1.0))
        newok = isfinite.(d) .& isfinite.(e) .& (abs.(r) .< clip * sc) .&
                (abs.(d .- med) .< max(8 * s, 0.05))
        newok == ok && break
        ok = newok
    end
    rms = sqrt(mean(abs2, (d .- (a .+ b .* t))[ok]))
    return (a, b, ea, eb, chi, count(ok), rms)
end

"""Per-peak displacement (px) of comb peaks at pixels `x` on fiber f, chip c against the comb
night solution (lin, nl, L, eps): delta = (λ(x) - λ_m)/(dλ/dx), λ_m = 2L/(m + eps), m the nearest
integer. Returns (delta, m_frac)."""
function refit_comb_delta(x::AbstractVector, f::Integer, c::Integer, lin, nl, L, eps)
    n = length(x)
    dl = fill(NaN, n)
    mf = fill(NaN, n)
    n == 0 && return dl, mf
    l0 = refit_prod_wave(lin, nl, f, c, x)
    lp = refit_prod_wave(lin, nl, f, c, x .+ 0.5)
    lm = refit_prod_wave(lin, nl, f, c, x .- 0.5)
    for i in 1:n
        isfinite(l0[i]) || continue
        mc = 2L / l0[i] - eps
        m = round(mc)
        mf[i] = mc - m
        dl[i] = (l0[i] - 2L / (m + eps)) / (lp[i] - lm[i])
    end
    return dl, mf
end

"""
    refit_fpifib_exposure(fnames_by_chip, fibers, comb; nsig = 8.0) -> NamedTuple

Stage 4 for one exposure: the local-noise finder on the FPI fibers (`fibers`, up to 2) of the
three ar1Dcal files (read-only), per-peak displacement against the night comb, and per (fiber,
chip) a + b t fits. Arrays are (slot, chip): a, b, ea, eb, n, rms, mfrac_max.
"""
function refit_fpifib_exposure(fnames_by_chip::AbstractVector{<:AbstractString},
        fibers::AbstractVector{<:Integer}, comb; nsig = 8.0)
    ns = 2
    a = fill(NaN, ns, N_CHIPS)
    b = fill(NaN, ns, N_CHIPS)
    ea = fill(NaN, ns, N_CHIPS)
    eb = fill(NaN, ns, N_CHIPS)
    n = zeros(Int, ns, N_CHIPS)
    rms = fill(NaN, ns, N_CHIPS)
    mfx = fill(NaN, ns, N_CHIPS)
    for (c, fn) in enumerate(fnames_by_chip)
        isfile(fn) || continue
        h5open(fn, "r") do fh
            for (s, fb) in enumerate(fibers[1:min(end, ns)])
                fl = fh["flux_1d"][:, fb]
                iv = fh["ivar_1d"][:, fb]
                mk = fh["mask_1d"][:, fb]
                _, prm, cov, _ = get_arclamp_peaks_local(fl, iv, mk; nsig = nsig)
                size(prm, 1) == 0 && continue
                x = prm[:, 2]
                er = sqrt.(max.(cov[:, 2, 2], 0.0))
                dl, mf = refit_comb_delta(x, fb, c, comb.lin, comb.nl, comb.L, comb.eps)
                r = refit_fit_ab(x, dl, er)
                a[s, c], b[s, c], ea[s, c], eb[s, c] = r[1], r[2], r[3], r[4]
                n[s, c] = r[6]
                rms[s, c] = r[7]
                fm = filter(isfinite, mf)
                mfx[s, c] = isempty(fm) ? NaN : maximum(abs.(fm))
            end
        end
    end
    return (; a, b, ea, eb, n, rms, mfrac_max = mfx, fibers = collect(Int, fibers[1:min(end, ns)]))
end

"""A (fiber-slot, chip) FPI-fiber fit is usable: it looks like a comb (>= 80 clipped peaks on the
chip; a lit comb gives 120-210) with a comb-like residual (rms < 0.1 px; arc-lamp lines landing on
an FPI fiber give 3-5 px) and a centroid-level zero-point error (< 0.005 px)."""
refit_fpifib_usable(P, s, c) = P.n[s, c] >= 80 && isfinite(P.a[s, c]) && P.ea[s, c] < 0.005 &&
                               (!hasproperty(P, :rms) || P.rms[s, c] < 0.1)

"""Per-exposure displacement model from the FPI-fiber fits for a term set (`exp_terms`).
Returns a function (f, c, p) -> δ px, plus (n_fibers_used, tilt_px per chip, flags)."""
function refit_fpifib_transform(P, exp_terms::AbstractString)
    use = [all(refit_fpifib_usable(P, s, c) for c in 1:N_CHIPS) for s in 1:length(P.fibers)]
    nf = count(use)
    nf == 0 && return nothing, 0, fill(NaN, 3), WAVE_FPIFIB_UNAVAILABLE
    flags = 0
    ss = findall(use)
    za = [mean(P.a[ss, c]) for c in 1:N_CHIPS]
    zb = [mean(P.b[ss, c]) for c in 1:N_CHIPS]
    tilt = fill(NaN, 3)
    if exp_terms == "shift"
        sh = mean(za)
        return ((f, c, p) -> sh), nf, tilt, flags
    elseif exp_terms == "chipzp"
        return ((f, c, p) -> za[c]), nf, tilt, flags
    elseif exp_terms == "chipzp_stretch" || nf < 2
        exp_terms == "chipzp_stretch_tilt" && (flags |= WAVE_FPIFIB_ONE_FIBER)
        return ((f, c, p) -> za[c] + zb[c] * (p - 1024.0) / 1024.0), nf, tilt, flags
    end
    f1, f2 = P.fibers[1], P.fibers[2]
    for c in 1:N_CHIPS
        tilt[c] = P.a[2, c] - P.a[1, c]      # px between the two FPI fibers
        er = sqrt(P.ea[1, c]^2 + P.ea[2, c]^2)
        abs(tilt[c]) > 3 * max(er, 0.001) && (flags |= WAVE_FPIFIB_ROTATION)
    end
    a1, b1, a2, b2 = P.a[1, :], P.b[1, :], P.a[2, :], P.b[2, :]
    fn = function (f, c, p)
        w = (f - f1) / (f2 - f1)
        u = (p - 1024.0) / 1024.0
        return (1 - w) * (a1[c] + b1[c] * u) + w * (a2[c] + b2[c] * u)
    end
    return fn, nf, tilt, flags
end

# ------------------------------------------------------------------------------------------
# stage 5: sky lines (pooled, robust to misassigned OH lines) and per-fiber sky terms
# ------------------------------------------------------------------------------------------
"""OH-line rows of an object exposure from the pass-1 `skyLinePeaks_*` files (read-only), with
pass-1's ingest mask (a line whose across-fiber IQR in x exceeds 0.002 is dropped). Returns a
NamedTuple of vectors (fiber, chip, line, pix, lam) or `nothing` if a chip file is missing."""
function refit_read_sky(mdir::AbstractString, tele::AbstractString, mjd::Integer, e::Integer)
    rows = (fiber = Int[], chip = Int[], line = Int[], pix = Float64[], lam = Float64[])
    for (c, chip) in enumerate(CHIP_LIST)
        fn = joinpath(mdir, "skyLinePeaks_$(tele)_$(mjd)_$(_arc_expstr(e))_$(chip)_object.h5")
        isfile(fn) || return nothing
        M = h5open(f -> read(f["sky_line_mat_clean"]), fn, "r")
        size(M, 3) == N_FIBERS || return nothing
        for l in 1:size(M, 1)
            xx = filter(isfinite, (M[l, 1, :] .- 1024) ./ 2048)
            xx = filter(!=(0), xx)
            length(xx) < 10 && continue
            q = quantile(xx, [0.25, 0.75])
            q[2] - q[1] > 0.002 && continue
            for f in 1:N_FIBERS
                p = M[l, 1, f]
                w = M[l, 2, f]
                (isfinite(p) && isfinite(w) && p > 1) || continue
                push!(rows.fiber, f)
                push!(rows.chip, c)
                push!(rows.line, round(Int, M[l, 9, f]))
                push!(rows.pix, p)
                push!(rows.lam, w)
            end
        end
    end
    return rows
end

"Velocity (m/s) of OH rows against a model: c (λ_OH - λ_model(pix))/λ_model, model = callback."
function refit_sky_velocities(rows, modelfun)
    v = fill(NaN, length(rows.pix))
    for i in eachindex(v)
        lm = modelfun(rows.fiber[i], rows.chip[i], rows.pix[i])
        isfinite(lm) && (v[i] = (rows.lam[i] - lm) / lm * REFIT_C_MPS)
    end
    return v
end

"""Pooled per-chip sky offset (m/s), robust to misassigned OH lines (the ARCFPI_VALIDATION
estimator): per line, rows farther than max(5 NMAD, 600 m/s) from the line median are dropped (a
wrong OH component sits >= 0.1 Å = 1.8 km/s away); a line needs `minfib` surviving fibers and
contributes its median (line-balanced); the chip value is a clipped mean over lines (2 rounds of
2.5 sigma about the median, sigma floor 5 m/s); a chip needs `minlines` lines.
Returns (zp[3], nlines[3], keep::Vector{Bool})."""
function refit_pooled_sky(rows, v; minfib = 30, minlines = 3)
    zp = fill(NaN, 3)
    nl = zeros(Int, 3)
    keep = falses(length(v))
    for c in 1:3
        ii = findall(==(c), rows.chip)
        isempty(ii) && continue
        lv = Float64[]
        lrows = Vector{Int}[]
        for l in unique(rows.line[ii])
            jj = ii[rows.line[ii] .== l]
            vv = v[jj]
            ok = isfinite.(vv)
            count(ok) < minfib && continue
            m = median(vv[ok])
            s = 1.4826 * median(abs.(vv[ok] .- m))
            gsel = jj[ok .& (abs.(vv .- m) .< max(5s, 600.0))]
            length(gsel) < minfib && continue
            push!(lrows, gsel)
            push!(lv, median(v[gsel]))
        end
        length(lv) < minlines && continue
        use = trues(length(lv))
        for _ in 1:2
            m = median(lv[use])
            s = max(1.4826 * median(abs.(lv[use] .- m)), 5.0)
            nu = abs.(lv .- m) .< 2.5s
            count(nu) == 0 && break
            use = nu
        end
        nl[c] = count(use)
        nl[c] < minlines && continue
        zp[c] = mean(lv[use])
        for k in findall(use)          # rows of whole-line misassignments are not kept either
            keep[lrows[k]] .= true
        end
    end
    return zp, nl, keep
end

"""Per-fiber sky terms (m/s) on the rows kept by `refit_pooled_sky`: per fiber and chip a robust
zero point and stretch (v = z + s u, u = (pix - 1024.5)/1024), fiber median of z removed (fiber-
relative only). Returns (z, s) each N_FIBERS x 3 (NaN where < 3 lines)."""
function refit_fiber_sky_terms(rows, v, keep)
    z = fill(NaN, N_FIBERS, N_CHIPS)
    s = fill(NaN, N_FIBERS, N_CHIPS)
    for c in 1:N_CHIPS
        for f in 1:N_FIBERS
            ii = [i for i in eachindex(v) if keep[i] && rows.chip[i] == c && rows.fiber[i] == f]
            length(ii) < 3 && continue
            u = (rows.pix[ii] .- 1024.5) ./ 1024
            y = v[ii]
            if length(ii) >= 5 && (maximum(u) - minimum(u)) > 0.3
                X = hcat(ones(length(y)), u)
                w = ones(length(y))
                p = zeros(2)
                for _ in 1:3
                    p = (X .* w) \ (y .* w)
                    r = y .- X * p
                    sg = max(1.4826 * median(abs.(r)), 5.0)
                    w = sqrt.(min.(1.0, 1.5sg ./ max.(abs.(r), 1e-9)))
                end
                z[f, c], s[f, c] = p[1], p[2]
            else
                z[f, c] = median(y)
                s[f, c] = 0.0
            end
        end
        m = _arc_nanmedian(z[:, c])
        z[:, c] .-= m
    end
    return z, s
end

"""Per-chip pooled-sky per-exposure transform (m/s) on the rows kept by `refit_pooled_sky` (the QA
that removes wrong-peak matches): theta[:, c] = (zp, stretch, tilt) of v = zp + st u + tl phi,
u = (pix - 1024.5)/1024, phi = (fiber - 150.5)/149.5. "chipzp": zp = the line-balanced pooled
value, st = tl = 0. Otherwise a robust (Huber IRLS) fit on the kept rows of the terms requested.
Chips with fewer than 3 lines are NaN."""
function refit_sky_terms(rows, v, keep, zp_pooled, nlines, terms::AbstractString; minlines = 3)
    th = fill(NaN, 3, 3)
    for c in 1:3
        nlines[c] >= minlines || continue
        if terms == "chipzp"
            th[:, c] .= (zp_pooled[c], 0.0, 0.0)
            continue
        end
        ii = [i for i in eachindex(v) if keep[i] && rows.chip[i] == c && isfinite(v[i])]
        length(ii) < 20 && continue
        u = (rows.pix[ii] .- 1024.5) ./ 1024
        ph = (rows.fiber[ii] .- 150.5) ./ 149.5
        y = v[ii]
        X = terms == "chipzp_stretch" ? hcat(ones(length(y)), u) : hcat(ones(length(y)), u, ph)
        w = ones(length(y))
        b = zeros(size(X, 2))
        for _ in 1:4
            b = (X .* w) \ (y .* w)
            r = y .- X * b
            sg = max(1.4826 * median(abs.(r)), 5.0)
            w = sqrt.(min.(1.0, 1.5sg ./ max.(abs.(r), 1e-9)))
        end
        th[1, c] = b[1]
        th[2, c] = b[2]
        th[3, c] = size(X, 2) == 3 ? b[3] : 0.0
    end
    return th
end

"Apply per-chip sky terms theta (3 x 3, m/s) multiplicatively to a grid (N_XPIX x N_FIBERS x 3)."
function refit_apply_sky_terms!(G::AbstractArray, theta::AbstractMatrix)
    u = ((1:size(G, 1)) .- 1024.5) ./ 1024
    for c in 1:3
        all(isfinite, view(theta, :, c)) || continue
        for f in 1:size(G, 2)
            ph = (f - 150.5) / 149.5
            G[:, f, c] .*= 1 .+ (theta[1, c] .+ theta[2, c] .* u .+ theta[3, c] * ph) ./ REFIT_C_MPS
        end
    end
    return G
end

# ------------------------------------------------------------------------------------------
# plate era: comb shape template (file written by scripts/refit/build_comb_shape_template.jl)
# ------------------------------------------------------------------------------------------
"Legendre polynomials P0..P6 at t."
@inline function refit_legendre7(t)
    P = zeros(7)
    P[1] = 1.0
    P[2] = t
    for k in 2:6
        P[k + 1] = ((2k - 1) * t * P[k] - (k - 1) * P[k - 1]) / k
    end
    return P
end

"""Comb shape template of one telescope and block: per (Legendre coefficient P0..P6, chip,
fiber) the median pixel-displacement field of the comb night solutions against a fixed reference
solution (`ref_lin`, `ref_nl`); t = (p - 1025)/960. P2..P6 carry the shape."""
struct RefitTemplate
    tele::String
    block::String
    T::Array{Float64, 3}
    ref_lin::Matrix{Float64}
    ref_nl::Matrix{Float64}
    fingerprint::String
end

function load_refit_template(path::AbstractString, tele::AbstractString, block::AbstractString)
    (isempty(path) || !isfile(path)) && return nothing
    h5open(path, "r") do f
        (haskey(f, tele) && haskey(f[tele], "T_$(block)")) || return nothing
        RefitTemplate(String(tele), String(block), read(f["$(tele)/T_$(block)"]),
            read(f["$(tele)/ref_lin"]), read(f["$(tele)/ref_nl"]), arclamp_file_fingerprint(path))
    end
end

"Reference solution wavelength and dλ/dpix of fiber f, chip c at (possibly fractional) pixel p."
function refit_ref_wave_disp(Tp::RefitTemplate, f, c, p)
    x = (p - 1024) / 2048
    cp = size(Tp.ref_nl, 2) ÷ 2 - 1
    if c == 2
        xt, dc = x, 1.0
    else
        o = c == 1 ? 0 : cp + 1
        xt, dc = cp == 0 ? (Tp.ref_nl[f, o + 1] + x, 1.0) :
                 (Tp.ref_nl[f, o + 1] + Tp.ref_nl[f, o + 2] * x, Tp.ref_nl[f, o + 2])
    end
    lin = view(Tp.ref_lin, f, :)
    lam = evalpoly(xt, lin)
    dl = 0.0
    for k in 2:length(lin)
        dl += (k - 1) * lin[k] * xt^(k - 2)
    end
    return lam, dl * dc / 2048
end

"Pixel where the reference solution of fiber f, chip c equals `lambda` (Newton from q0)."
function refit_ref_invert(Tp::RefitTemplate, f, c, lambda, q0)
    q = q0
    for _ in 1:30
        l, s = refit_ref_wave_disp(Tp, f, c, q)
        st = (l - lambda) / s
        q -= st
        abs(st) < 1e-9 && break
    end
    return q
end

"Template model wavelength of fiber f, chip c at reference-frame pixel pr; P = (zp_R, zp_G, zp_B, stretch)."
function refit_template_lambda(Tp::RefitTemplate, P, f, c, pr)
    t = (pr - 1025) / 960
    Lg = refit_legendre7(t)
    q = pr + P[c] + (Tp.T[2, c, f] + P[4]) * t
    for k in 3:7
        q += Tp.T[k, c, f] * Lg[k]
    end
    return refit_ref_wave_disp(Tp, f, c, q)[1]
end

"Pixel shift (px) of an arc state on chip c of fiber f at pixel p (reference-frame pixel = p + shift)."
function refit_state_shift_px(A::RefitArc, f, c, p, st)
    a, b = refit_arc_shift_AB(A, f, c, st)
    x = (p - 1024) / 2048
    bc = c == 1 ? A.TH[f, 2] : c == 3 ? A.TH[f, 4] : 1.0
    return N_XPIX * (a + b * x) / bc
end

"""Fit the template to the night's kept arc lines (lab wavelengths of list `ll`): per fiber 3 chip
zero points + 1 common stretch, the global position terms from the arc fit. Returns (P 4 x N, ok)."""
function refit_fit_template(Tp::RefitTemplate, A::RefitArc, ll::ArcLineList)
    P = fill(NaN, 4, N_FIBERS)
    ok = falses(N_FIBERS)
    L = A.lines
    jmap = Dict(s => j for (j, s) in enumerate(A.used))
    rowsf = [Int[] for _ in 1:N_FIBERS]
    for i in eachindex(L.id)
        L.keep[i] || continue
        haskey(jmap, L.seq[i]) || continue
        haskey(ll.labwave, L.id[i]) || continue
        push!(rowsf[L.fiber[i]], i)
    end
    for f in 1:N_FIBERS
        A.okf[f] || continue
        ii = rowsf[f]
        length(ii) < 30 && continue
        nr = length(ii)
        X = zeros(nr, 4)
        y = zeros(nr)
        w = zeros(nr)
        for (k, i) in enumerate(ii)
            c = L.chip[i]
            j = jmap[L.seq[i]]
            d, t = arcfit_seq_state(A.g, A.lay, j)
            st = (A.grp_u[j], d, t, A.tau_u[j])
            pr = L.pix[i] + refit_state_shift_px(A, f, c, L.pix[i], st)
            q = refit_ref_invert(Tp, f, c, ll.labwave[L.id[i]], pr)
            tt = (pr - 1025) / 960
            Lg = refit_legendre7(tt)
            yy = q - pr - Tp.T[2, c, f] * tt
            for kk in 3:7
                yy -= Tp.T[kk, c, f] * Lg[kk]
            end
            X[k, c] = 1.0
            X[k, 4] = tt
            y[k] = yy
            w[k] = 1.0 / (max(L.pixvar[i], 0.0) + 0.02^2)
        end
        all(any(X[:, c] .!= 0) for c in 1:3) || continue
        good = isfinite.(y)
        X, y, w = X[good, :], y[good], w[good]
        length(y) < 30 && continue
        keep = abs.(y .- median(y)) .< 3.0
        b = zeros(4)
        for _ in 1:4
            sw = sqrt.(w[keep])
            b = (X[keep, :] .* sw) \ (y[keep] .* sw)
            r = y .- X * b
            z = r .* sqrt.(w)
            s = 1.4826 * median(abs.(z[keep]))
            keep = (abs.(r) .< 1.3) .& (abs.(z) .< max(5s, 1e-6))
            count(keep) < 10 && break
        end
        P[:, f] .= b
        ok[f] = true
    end
    return P, ok
end

# ------------------------------------------------------------------------------------------
# stored solution: per fiber, per chip Legendre coefficients in t = (pix - 1024.5)/1023.5
# ------------------------------------------------------------------------------------------
"Legendre basis (N_XPIX x (deg+1)) on t = (pix - 1024.5)/1023.5."
function refit_coeff_basis(deg::Integer)
    t = ((1:N_XPIX) .- 1024.5) ./ 1023.5
    B = zeros(N_XPIX, deg + 1)
    B[:, 1] .= 1.0
    deg >= 1 && (B[:, 2] .= t)
    for k in 2:deg
        B[:, k + 1] .= ((2k - 1) .* t .* B[:, k] .- (k - 1) .* B[:, k - 1]) ./ k
    end
    return B
end

"""Fit the per-fiber, per-chip coefficients to a wavelength grid (N_XPIX x N_FIBERS x 3).
Returns (C (deg+1) x N_FIBERS x 3, max |residual| Å). Non-finite fiber-chips give NaN."""
function refit_grid_to_coeffs(G::AbstractArray{<:Real, 3}, B::AbstractMatrix)
    deg1 = size(B, 2)
    C = fill(NaN, deg1, N_FIBERS, N_CHIPS)
    F = qr(B)
    mx = 0.0
    for c in 1:N_CHIPS, f in 1:N_FIBERS
        y = view(G, :, f, c)
        all(isfinite, y) || continue
        cf = F \ Vector{Float64}(y)
        C[:, f, c] .= cf
        mx = max(mx, maximum(abs.(B * cf .- y)))
    end
    return C, mx
end

"Evaluate stored coefficients to a wavelength grid (N_XPIX x N_FIBERS x 3), comparable with chipWaveSoln."
function refit_coeffs_to_grid(C::AbstractArray{<:Real, 3})
    B = refit_coeff_basis(size(C, 1) - 1)
    G = fill(NaN, N_XPIX, size(C, 2), size(C, 3))
    for c in axes(C, 3), f in axes(C, 2)
        all(isfinite, view(C, :, f, c)) || continue
        G[:, f, c] .= B * C[:, f, c]
    end
    return G
end

"""
    refit_exposure_wave(fname, expnum) -> Array (N_XPIX, N_FIBERS, 3)

The refit wavelength solution of exposure `expnum` from a night product
`waveCalRefit_<tele>_<mjd>.h5` (directly comparable with pass 1's `chipWaveSoln`).
"""
function refit_exposure_wave(fname::AbstractString, expnum::Integer)
    h5open(fname, "r") do f
        en = read(f["exp_expnum"])
        k = findfirst(==(expnum), en)
        k === nothing && throw(ArgumentError("exposure $(expnum) not in $(fname)"))
        refit_coeffs_to_grid(f["wave_coeffs"][:, :, :, k])
    end
end

# ------------------------------------------------------------------------------------------
# pass-1 comparison
# ------------------------------------------------------------------------------------------
"Pass-1 adopted solution file and dataset for an exposure (read-only), or `nothing`."
function refit_pass1_source(apred_dir, wavecal_dir, tele, mjd, e, kind, wt)
    mdir = joinpath(apred_dir, string(mjd))
    es = _arc_expstr(e)
    if wt == "waveCalNightfpiDither" || wt == "waveCalNightskyDither"
        fn = joinpath(mdir, "$(wt)_$(tele)_$(mjd)_$(es)_$(kind).h5")
        return isfile(fn) ? (fn, "chipWaveSoln") : nothing
    elseif startswith(wt, "error_night_ave")
        fn = joinpath(wavecal_dir, "wavecalNightAve_$(tele)_$(mjd).h5")
        grp = endswith(wt, "fpi") ? "fpi" : "sky"
        return isfile(fn) ? (fn, "$(grp)/nightAve_wave_soln") : nothing
    end
    return nothing
end

"""Comparison of a new grid on the flat pixel grid against pass 1: dv[c, f] = flat mean over
REFIT_PXG of c (λ_new - λ_p1)/λ_p1 (m/s; the first-order RV prediction of HOLDOUT_RV/AV2), the
per-chip fiber-median zero point, the fiber-median shape rms (dv about its flat mean) and the
across-fiber robust sigma of dv."""
function refit_vs_pass1(Gnew_pxg::AbstractArray, Gp1_pxg::AbstractArray)
    dv = fill(NaN32, N_CHIPS, N_FIBERS)
    shape = fill(NaN, N_CHIPS, N_FIBERS)
    for f in 1:N_FIBERS, c in 1:N_CHIPS
        a = view(Gnew_pxg, :, f, c)
        b = view(Gp1_pxg, :, f, c)
        (all(isfinite, a) && all(isfinite, b)) || continue
        v = (a .- b) ./ b .* REFIT_C_MPS
        m = mean(v)
        dv[c, f] = Float32(m)
        shape[c, f] = sqrt(mean((v .- m) .^ 2))
    end
    zp = [_arc_nanmedian(Float64.(dv[c, :])) for c in 1:N_CHIPS]
    sh = [_arc_nanmedian(shape[c, :]) for c in 1:N_CHIPS]
    rs = [_arc_robstd(Float64.(dv[c, :])) for c in 1:N_CHIPS]
    return dv, zp, sh, rs
end
