# 2D content classifier ("classifier_2d") and the approved commanded-object rule.
#
# Port of the phase-1 prototype (2026_11_19/classifier_2d, AKS 2026-11-19..22):
# a random forest over 173 image features that answers "what light is in this
# frame" with ten content classes, plus a thresholded rule that decides when a
# COMMANDED-OBJECT exposure should not be used for science.
#
# Status (AKS 2026-11-22): this forest REPLACES the v6 exposure-type classifier
# (pipeline.jl --exp_class_model), and for commanded objects the approved rule
# sets the existing EXPFLAG_PREDICTED_BAD bit (see `exposure_predicted_bad` in
# exposureClassifier.jl). No explicit twilight / sun-altitude rule (twilight is
# to be revisited later).
#
# Inputs are the same ar2D products the deployed check sees (`dimage` per chip),
# plus, for chip G only, the up-the-ramp products that live in the same file and
# are in memory in `process_2Dcal` (`last_unsaturated`, `chisqimage`, `CRimage`,
# `metadata/ndiff_used`). No fiber map and
# no wavelength solution.
#
# STORED REFERENCES (calibration inputs; each has a from-scratch rebuild script):
#   * line templates, 24 vectors of 2048 px: median normalized line excess of
#     labeled sky (OH) / ThAr / UNe / FPI exposures per (telescope, chip), from
#     template nights only. Feed the 12 `t_*` features.
#     File: data/exposure_classifier_2d/line_templates.h5
#     Rebuild: scripts/classifier/build_line_templates_2d.jl
#   * light reference, 12 numbers: per (telescope, era) median chip-G OH flux,
#     continuum and trace-height q90 of a reference set of release-NORMAL object
#     exposures. Feed `oh_rel`, `cont_rel`, `trq90_rel`.
#     File: data/exposure_classifier_2d/light_reference.csv
#     Rebuild: scripts/classifier/build_light_reference_2d.jl
#   * the forest itself (JLD2, ~290 MB, not in git; sha256 in
#     metadata/exposure_classifier_2d.sha256).
#     Rebuild: scripts/classifier/train_classifier_2d.jl
# When any of these is absent `load_exposure_classifier_2d` returns `nothing`
# and logs that the 2D content check is DISABLED; the reference-free features
# can still be computed.
#
# NUMERICAL FIDELITY. Everything below reproduces the prototype bit for bit,
# including two of its quirks that matter for exact reproduction:
#   1. the prototype stored the spectral (ps) and fiber (pf) profiles as Float32
#      and computed every profile feature and template correlation from the
#      Float32-rounded profiles. `classifier2d_chip_products` rounds the same way.
#   2. (v1 only) the approved prototype forest also used two READ-COUNT inputs that
#      are not image content: `n_read`, the commanded number of reads (it entered
#      through the prototype's "n_" column prefix, from the almanac), and
#      `n_ndiff`, the number of read differences the ramp fit used
#      (metadata/ndiff_used). AKS 2026-11-22: "retrain without the number of
#      reads" -> the v2 feature set (`classifier2d_feature_names()`, 173
#      features) drops both. `n_read` is no longer computed; `n_ndiff` stays in
#      the feature Dict as a diagnostic only. The telescope flag `d_is_apo` is a
#      non-image input too; it was in the approved v1 model (and in v6) and is
#      kept. A v1 artifact cannot be scored by this code (it lacks `n_read`) and
#      is rejected by `classify_exposure_2d`.

using FFTW: rfft
using LinearAlgebra: dot
using Statistics, StatsBase, CSV, DataFrames, HDF5, JLD2
using DecisionTree: apply_forest_proba

## ---------------------------------------------------------------------------
## feature names
## ---------------------------------------------------------------------------
const CLASSIFIER2D_DEPNAMES = vcat(["q$(q)" for q in CLASSIFIER_QUANTS],
    ["mad", "nanfrac", "ps_med", "spec_maxrat", "spec_npk5", "spec_npk20",
        "spec_conc20", "spec_acpk", "pf_med", "fib_contrast", "fib_npk2", "fib_acpk"])
const CLASSIFIER2D_PROFNAMES = ["oh", "cont", "line_cont", "nline5", "line_maxz",
    "line_z", "comb_pow", "comb_period", "comb_top5", "spacing_cv", "ntrace", "tr_q10",
    "tr_q50", "tr_q90", "tr_q99", "tr_floor", "tr_dyn", "lit_frac", "blk_min",
    "blk_ndark", "tr_width", "tr_contrast"]
const CLASSIFIER2D_EXTRANAMES = ["q9999", "pmax", "fpos100", "fpos1000", "fneg"]
const CLASSIFIER2D_QUARTNAMES = ["pfq_min_over_max", "pfq_slope"]
const CLASSIFIER2D_RAMPNAMES = ["sat_frac", "sat_frac_lit", "chi2r_med", "chi2r_q99",
    "chi2r_frac10", "cr_frac", "ndiff"]
const CLASSIFIER2D_TEMPLATE_KINDS = ("oh", "thar", "une", "fpi")
# the ten content classes of the forest, in the artifact's order
const CLASSIFIER2D_CLASSES = ["sky", "twilight", "dark", "dark_persist", "ledflat",
    "quartzflat", "domeflat", "thar", "une", "fpi"]

"Features of the v1 (approved prototype) forest that are read counts, not image content; removed in v2."
const CLASSIFIER2D_READCOUNT_FEATURES = ["n_read", "n_ndiff"]

"""
    classifier2d_feature_names(version = 2) -> Vector{String}

The features of the 2D content forest, in training column order (`d_*`
deployed block incl. the telescope flag `d_is_apo`; `n_*` per-chip profile and
image features, quarters, ramp, continuum colors; `t_*` template correlations;
`oh_rel`, `cont_rel`, `trq90_rel`).
- `version = 2` (deployed, AKS 2026-11-22): 173 features, without the read
  counts `CLASSIFIER2D_READCOUNT_FEATURES`.
- `version = 1`: the 175 features of the approved prototype forest, incl.
  `n_read` (after the `t_*` block) and `n_ndiff`; kept so the v1 retrain stays
  reproducible.
The order matters: the forest's random feature subsets index columns, so it is
part of what makes a retrain reproduce an artifact.
"""
function classifier2d_feature_names(version::Integer = 2)
    v1 = _classifier2d_feature_names_v1()
    version == 1 && return v1
    version == 2 && return filter(!in(CLASSIFIER2D_READCOUNT_FEATURES), v1)
    throw(ArgumentError("unknown classifier_2d feature-set version $version"))
end
function _classifier2d_feature_names_v1()
    d = vcat(["d_" * c * "_" * n for c in CHIP_LIST for n in CLASSIFIER2D_DEPNAMES],
        ["d_col_rg", "d_col_bg", "d_col_rb", "d_is_apo"])
    n = vcat(["n_" * c * "_" * nm for c in CHIP_LIST
              for nm in vcat(CLASSIFIER2D_PROFNAMES, CLASSIFIER2D_EXTRANAMES)],
        "n_" .* CLASSIFIER2D_QUARTNAMES, "n_" .* CLASSIFIER2D_RAMPNAMES,
        ["n_lcont_RG", "n_lcont_BG"])
    t = ["t_$(tm)_$(c)" for tm in CLASSIFIER2D_TEMPLATE_KINDS for c in CHIP_LIST]
    vcat(d, n, t, ["n_read", "oh_rel", "cont_rel", "trq90_rel"])
end

"""
First MJD of the FPS era per telescope as used by the prototype's `era` (and so
by the light reference). n.b. `get_fps_plate_divide` returns 59808 for LCO; the
prototype used `mjd >= 59809`. Kept separate so the reference keys stay exactly
the ones the approved scores were computed with.
"""
const CLASSIFIER2D_ERA2_START = Dict("apo" => 59423, "lco" => 59809)
classifier2d_era(tele, mjd) = mjd >= CLASSIFIER2D_ERA2_START[lowercase(String(tele))] ? 2 : 1

## ---------------------------------------------------------------------------
## profile helpers (verbatim prototype algorithms)
## ---------------------------------------------------------------------------
function classifier2d_runmed(x, w)
    n = length(x)
    o = similar(x, Float64)
    h = w ÷ 2
    for i in 1:n
        o[i] = median(@view x[max(1, i - h):min(n, i + h)])
    end
    o
end
function classifier2d_runmin(x, w)
    n = length(x)
    o = similar(x, Float64)
    h = w ÷ 2
    for i in 1:n
        o[i] = minimum(@view x[max(1, i - h):min(n, i + h)])
    end
    o
end
classifier2d_robsig(x) = 1.4826 * median(abs.(x .- median(x))) + 1e-9

"""
    classifier2d_line_excess(ps) -> (ex, c)

Line excess of a spectral profile: `ex = ps - 41-px running median`, with the 4
edge pixels on each side replaced by their neighbor first. Used for the template
correlations and to BUILD the templates (the same function on both sides).
"""
function classifier2d_line_excess(ps)
    p = Float64.(ps)
    p[1:4] .= p[5]
    p[(end - 3):end] .= p[end - 4]
    c = classifier2d_runmed(p, 41)
    ex = p .- c
    ex, c
end

"""
    classifier2d_tcorr(ex, t; L = 12)

Normalized cross-correlation of the line excess `ex` with template `t`,
maximized over integer lags |l| <= L (absorbs small dither / wavelength shifts).
"""
function classifier2d_tcorr(ex, t; L = 12)
    best = -1.0
    a = ex .- mean(ex)
    b = t .- mean(t)
    for l in (-L):L
        if l >= 0
            x = @view a[(1 + l):end]
            y = @view b[1:(end - l)]
        else
            x = @view a[1:(end + l)]
            y = @view b[(1 - l):end]
        end
        den = sqrt(sum(abs2, x) * sum(abs2, y)) + 1e-12
        c = dot(x, y) / den
        c > best && (best = c)
    end
    best
end

"""Comb periodicity of a high-passed spectral profile: (log10 peak/median power for
periods 4-80 px, the period of that peak, top-5 power fraction)."""
function classifier2d_comb_stats(ex)
    x = ex .- mean(ex)
    n = length(x)
    P = abs2.(abs.(rfft(x)))
    k = 2:length(P)
    per = n ./ (k .- 1)
    sel = (per .>= 4) .& (per .<= 80)
    Pb = P[k][sel]
    isempty(Pb) && return (0.0, 0.0, 0.0)
    m = argmax(Pb)
    (log10(Pb[m] / (median(Pb) + 1e-30)), per[sel][m],
        sum(sort(Pb, rev = true)[1:min(5, end)]) / (sum(Pb) + 1e-30))
end

"""
    classifier2d_profile_features(ps, pf) -> Vector{Float64} (CLASSIFIER2D_PROFNAMES order)

Per-chip profile features from the spectral profile `ps` (mean over the fiber
axis) and the fiber profile `pf` (mean over the spectral axis). Chip G's
`oh` / `cont` / `tr_*` are the EXPOSURE_LIGHT 2D metrics (2026_11_18,
s31_2dmetrics.jl):
- `oh`: sum of the positive excess above a 41-px running median, over pixels
  above 5 sigma (sigma = 1.4826 median |ps - runmed|);
- `line_z = log10(oh / (sigma sqrt(n_selected)))`: detection significance of the
  line flux (-3 when nothing is selected);
- `tr_*`: trace peaks of `pf` (local maxima over +-2 px above a 9-px running
  minimum); `tr_dyn = log10(q90 / q10)` of the trace heights, NaN when <= 50
  peaks. Strong, non-uniform star traces give a large `tr_dyn`.
"""
function classifier2d_profile_features(ps, pf)
    ps = Float64.(ps)
    pf = Float64.(pf)
    c = classifier2d_runmed(ps, 41)
    ex = ps .- c
    σ = 1.4826 * median(abs.(ps .- c)) + 1e-9
    sel = ex .> 5σ
    oh = sum(ex[sel])
    cont = median(c)
    line_cont = log10(oh / (abs(cont) * 2048 + 1e-9) + 1e-12)
    pk = [j for j in 2:(length(ex) - 1) if ex[j] > 5σ && ex[j] >= ex[j - 1] && ex[j] >= ex[j + 1]]
    nline5 = length(pk)
    line_maxz = log10(max(maximum(ex[5:(end - 4)]) / σ, 1e-3))
    line_z = sum(sel) > 0 ? log10(max(oh / (σ * sqrt(sum(sel))), 1e-3)) : -3.0
    exc = copy(ex)
    exc[1:4] .= 0
    exc[(end - 3):end] .= 0
    cp, cper, ctop = classifier2d_comb_stats(exc)
    spacing_cv = nline5 >= 6 ? (d = diff(pk); std(d) / mean(d)) : NaN
    # traces (EXPOSURE_LIGHT s31 algorithm)
    fl = classifier2d_runmin(pf, 9)
    tk = [j for j in 3:(length(pf) - 2) if pf[j] >= pf[j - 1] && pf[j] >= pf[j + 1] &&
                                         pf[j] >= pf[j - 2] && pf[j] >= pf[j + 2]]
    h = pf[tk] .- fl[tk]
    keep = h .> 0
    tk = tk[keep]
    h = h[keep]
    q = length(h) > 50 ? quantile(h, [0.1, 0.5, 0.9, 0.99]) : fill(NaN, 4)
    floor_ = median(fl)
    tr_dyn = length(h) > 50 ? log10(max(q[3], 1e-6) / max(q[1], 1e-6)) : NaN
    # lit fibers: trace peaks above 20 % of the q90 trace height and 5x the profile noise
    pn = 1.4826 * median(abs.(diff(pf))) / sqrt(2) + 1e-9
    lit = length(h) > 0 ? sum((h .> 0.2 * max(q[3], 0)) .& (h .> 5pn)) : 0
    lit_frac = lit / 300
    # 10 blocks along the slit: median trace height per block relative to the overall median
    blk = Float64[]
    for b in 0:9
        m = (tk .> b * 204.8) .& (tk .<= (b + 1) * 204.8)
        push!(blk, sum(m) >= 3 ? median(h[m]) : 0.0)
    end
    mb = median(blk) + 1e-9
    blk_min = minimum(blk) / mb
    blk_ndark = sum(blk ./ mb .< 0.2)
    # trace width: second moment within +-3 px of the 30 highest peaks
    ord = sortperm(h, rev = true)[1:min(30, length(h))]
    w = Float64[]
    for i in ord
        j = tk[i]
        (j <= 3 || j >= length(pf) - 3) && continue
        y = pf[(j - 3):(j + 3)] .- fl[j]
        y = max.(y, 0)
        s = sum(y)
        s <= 0 && continue
        x = -3:3
        push!(w, sqrt(sum(y .* x .^ 2) / s))
    end
    tr_width = isempty(w) ? NaN : median(w)
    tr_contrast = length(h) > 0 ? log10(max(median(h), 1e-6) / max(abs(floor_), 1e-3)) : NaN
    Float64[oh, cont, line_cont, nline5, line_maxz, line_z, cp, cper, ctop, spacing_cv,
        length(h), q..., floor_, tr_dyn, lit_frac, blk_min, blk_ndark, tr_width, tr_contrast]
end

"""Chip-G spectral-quarter illumination: (min/max of the 4 quarter medians, slope)."""
function classifier2d_quarter_features(PFQ)
    m = [median(@view PFQ[:, k]) for k in 1:4]
    mx = maximum(abs.(m)) + 1e-9
    [minimum(m) / mx, (m[4] - m[1]) / (abs(m[4]) + abs(m[1]) + 1e-9)]
end

## ---------------------------------------------------------------------------
## per-chip products from the in-memory ar2D arrays
## ---------------------------------------------------------------------------
"""
    classifier2d_chip_products(dimage; last_unsaturated = nothing, chisqimage = nothing,
                               CRimage = nothing, ndiff_used = -1, quarters = false)

Everything the 2D classifier needs from ONE chip's ar2D arrays (dim1 spectral,
dim2 fiber axis; rows 1:N_XPIX are the science region), as a NamedTuple:
- `dep`: the 21 deployed features, `exposure_class_features(dimage)` itself;
- `ps`, `pf`: spectral / fiber profiles, rounded to Float32 (prototype storage);
- `ext`: q99.99, max, fraction > 100, fraction > 1000, fraction < -20 of the
  2x2-subsampled image (NaN set to 0);
- `pfq`: fiber profiles of the 4 spectral quarters (Float32), when `quarters`;
- `ramp`: the 7 up-the-ramp features (`CLASSIFIER2D_RAMPNAMES`) when the ramp
  arrays are passed, else NaN. The prototype computes them for chip G only,
  sampling every 8th fiber-axis column.

Call it on the UNTOUCHED ar2D arrays (before dark subtraction), as
`process_2Dcal` does for `exposure_class_features`; it does not mutate them.
"""
function classifier2d_chip_products(dimage; last_unsaturated = nothing,
        chisqimage = nothing, CRimage = nothing, ndiff_used = -1, quarters::Bool = false)
    dep = exposure_class_features(dimage)
    sci = @view dimage[1:N_XPIX, :]
    nanmask = isnan.(sci)
    s = copy(sci)
    s[nanmask] .= 0.0
    ps = Float32.(vec(mean(s, dims = 2)))
    pf = Float32.(vec(mean(s, dims = 1)))
    sub = vec(@view s[1:2:end, 1:2:end])
    ext = [quantile(sub, 0.9999), maximum(sub), mean(sub .> 100), mean(sub .> 1000),
        mean(sub .< -20)]
    pfq = zeros(Float32, N_XPIX, 4)
    if quarters
        for k in 1:4
            pfq[:, k] = vec(mean(@view(s[((k - 1) * 512 + 1):(k * 512), :]), dims = 1))
        end
    end
    ramp = fill(NaN, length(CLASSIFIER2D_RAMPNAMES))
    if !(last_unsaturated === nothing || chisqimage === nothing || CRimage === nothing)
        lu = last_unsaturated[1:N_XPIX, 1:8:N_XPIX]
        chi = chisqimage[1:N_XPIX, 1:8:N_XPIX]
        cr = CRimage[1:N_XPIX, 1:8:N_XPIX]
        nr = Float64(ndiff_used)
        lit = @view(s[:, 1:8:N_XPIX]) .> 20
        # last_unsaturated = index of the last unsaturated read difference
        # (= ndiff_used when never saturated, -1 for bad pixels); chisqimage is
        # the total up-the-ramp chi2 -> divide by ndiff for a reduced chi2
        satm = (lu .< nr) .& (lu .>= 0)
        c2 = Float64.(chi) ./ max(nr, 1)
        c2 = c2[isfinite.(c2)]
        ramp .= [mean(satm), sum(lit) > 0 ? mean(satm[lit]) : 0.0,
            isempty(c2) ? NaN : median(c2), isempty(c2) ? NaN : quantile(c2, 0.99),
            isempty(c2) ? NaN : mean(c2 .> 10), mean(cr .> 0), nr]
    end
    (dep = dep, ps = ps, pf = pf, ext = ext, pfq = pfq, ramp = ramp)
end

"""
    read_classifier2d_chip_products(path, chip) -> NamedTuple

Read one ar2D chip file READ-ONLY and return `classifier2d_chip_products` for
it. Chip G also reads the ramp products and `metadata/ndiff_used`. The returned
tuple carries `nread_total` (from `metadata/nread_total`, -1 when absent).
"""
function read_classifier2d_chip_products(path, chip)
    h5open(path, "r") do f
        d = read(f["dimage"])
        nread = haskey(f, "metadata/nread_total") ? Int(read(f["metadata/nread_total"])) : -1
        p = if chip == "G"
            nd = haskey(f, "metadata/ndiff_used") ? read(f["metadata/ndiff_used"]) : -1
            classifier2d_chip_products(d; last_unsaturated = read(f["last_unsaturated"]),
                chisqimage = read(f["chisqimage"]), CRimage = read(f["CRimage"]),
                ndiff_used = nd, quarters = true)
        else
            classifier2d_chip_products(d)
        end
        merge(p, (nread_total = nread,))
    end
end

## ---------------------------------------------------------------------------
## references (templates + light reference): calibration inputs
## ---------------------------------------------------------------------------
"""
    load_classifier2d_references(templates_path, light_reference_path) -> NamedTuple or nothing

Load the two stored references. Returns `nothing` (and logs why) when either file
is absent or unreadable; callers then run without the 2D content check.
- templates: HDF5 with one 2048-vector per `<kind>_<tele>_<chip>`
  (kind in `CLASSIFIER2D_TEMPLATE_KINDS`);
- light reference: CSV with columns `tele, era, n_ref, oh, cont, tr_q90`.
"""
function load_classifier2d_references(templates_path, light_reference_path)
    for (what, p) in (("line templates", templates_path),
        ("light reference", light_reference_path))
        if p === nothing || isempty(p) || !isfile(p)
            @warn "2D content check DISABLED: $what file not found ($(repr(p))). " *
                  "Rebuild it with scripts/classifier/ (see src/exposureClassifier2D.jl)."
            return nothing
        end
    end
    try
        T = h5open(templates_path, "r") do f
            Dict{String, Vector{Float64}}(k => read(f[k]) for k in keys(f))
        end
        L = CSV.read(light_reference_path, DataFrame)
        ref = Dict{Tuple{String, Int}, NTuple{3, Float64}}()
        for r in eachrow(L)
            ref[(String(r.tele), Int(r.era))] = (Float64(r.oh), Float64(r.cont), Float64(r.tr_q90))
        end
        return (templates = T, light = ref, templates_path = String(templates_path),
            light_reference_path = String(light_reference_path))
    catch e
        @warn "2D content check DISABLED: could not read the references" templates_path light_reference_path exception = e
        return nothing
    end
end

## ---------------------------------------------------------------------------
## the feature dictionary
## ---------------------------------------------------------------------------
"""
    exposure_class_features_2d(prods, tele, mjd; refs = nothing) -> Dict{String, Float64}

All features of the 2D content forest for one exposure, keyed by the names the
model artifact stores (`d_*`, `n_*`, `t_*`, `oh_rel`, `cont_rel`, `trq90_rel`;
plus the diagnostic `n_ndiff`, which the v2 forest does not use). `prods` maps chip ("R", "G", "B") to `classifier2d_chip_products`
(chip G with the ramp products and quarters).

`d_*` follow `classify_exposure_type` exactly (non-finite -> -9999 inside this
block, as the prototype stored them). Other features keep NaN; the scorer maps
every non-finite value to -9999 before the forest, as training did.

Without `refs` the template correlations (`t_*`) and the relative light metrics
(`*_rel`) are NaN: the reference-free features are still all there, but the
approved rule's no-sky-light test cannot run (it needs `oh_rel`).
"""
function exposure_class_features_2d(prods::AbstractDict, tele, mjd; refs = nothing)
    tele = lowercase(String(tele))
    F = Dict{String, Float64}()
    # deployed block (identical to classify_exposure_type's design row)
    x = Float64[]
    for chip in CHIP_LIST
        append!(x, prods[chip].dep)
    end
    q90i = findfirst(==(0.9), CLASSIFIER_QUANTS)
    r, g, b = (prods[c].dep[q90i] for c in CHIP_LIST)
    append!(x, [(r - g) / (abs(r) + abs(g) + 1e-3), (b - g) / (abs(b) + abs(g) + 1e-3),
        (r - b) / (abs(r) + abs(b) + 1e-3)])
    push!(x, tele == "apo" ? 1.0 : 0.0)
    x[.!isfinite.(x)] .= -9999.0
    dnames = vcat([c * "_" * n for c in CHIP_LIST for n in CLASSIFIER2D_DEPNAMES],
        ["col_rg", "col_bg", "col_rb", "is_apo"])
    for (nm, v) in zip(dnames, x)
        F["d_" * nm] = v
    end
    # profile + full-image extras per chip
    for chip in CHIP_LIST
        p = prods[chip]
        for (nm, v) in zip(CLASSIFIER2D_PROFNAMES, classifier2d_profile_features(p.ps, p.pf))
            F["n_$(chip)_$(nm)"] = v
        end
        for (nm, v) in zip(CLASSIFIER2D_EXTRANAMES, p.ext)
            F["n_$(chip)_$(nm)"] = v
        end
    end
    pG = prods["G"]
    for (nm, v) in zip(CLASSIFIER2D_QUARTNAMES, classifier2d_quarter_features(pG.pfq))
        F["n_" * nm] = v
    end
    for (nm, v) in zip(CLASSIFIER2D_RAMPNAMES, pG.ramp)
        F["n_" * nm] = v
    end
    cR, cG, cB = F["n_R_cont"], F["n_G_cont"], F["n_B_cont"]
    F["n_lcont_RG"] = log10(max(cR, 1e-3) / max(cG, 1e-3))
    F["n_lcont_BG"] = log10(max(cB, 1e-3) / max(cG, 1e-3))
    classifier2d_add_reference_features!(F, Dict(c => prods[c].ps for c in CHIP_LIST),
        tele, mjd, refs)
    F
end

"""
    classifier2d_add_reference_features!(F, ps_by_chip, tele, mjd, refs) -> F

Fill the reference-dependent features into the feature Dict `F` (which must
already hold `n_G_oh`, `n_G_cont`, `n_G_tr_q90`): the 12 template correlations
`t_<kind>_<chip>` from the (Float32) spectral profiles `ps_by_chip`, and
`oh_rel`, `cont_rel`, `trq90_rel` from the light reference. With `refs ===
nothing` they are all NaN. Shared by the pipeline path and the training
assembly (scripts/classifier/assemble_features_2d.jl).
"""
function classifier2d_add_reference_features!(F::AbstractDict, ps_by_chip::AbstractDict,
        tele, mjd, refs)
    tele = lowercase(String(tele))
    for tm in CLASSIFIER2D_TEMPLATE_KINDS, chip in CHIP_LIST
        F["t_$(tm)_$(chip)"] = NaN
    end
    F["oh_rel"] = F["cont_rel"] = F["trq90_rel"] = NaN
    refs === nothing && return F
    for chip in CHIP_LIST
        ex = classifier2d_line_excess(ps_by_chip[chip])[1]
        exn = ex ./ classifier2d_robsig(ex)
        for tm in CLASSIFIER2D_TEMPLATE_KINDS
            key = "$(tm)_$(tele)_$(chip)"
            haskey(refs.templates, key) &&
                (F["t_$(tm)_$(chip)"] = classifier2d_tcorr(exn, refs.templates[key]))
        end
    end
    k = (tele, classifier2d_era(tele, mjd))
    if haskey(refs.light, k)
        oh0, cont0, q900 = refs.light[k]
        F["oh_rel"] = log10(max(F["n_G_oh"], 1e-6) / oh0)
        F["cont_rel"] = log10(max(F["n_G_cont"], 1e-6) / cont0)
        F["trq90_rel"] = log10(max(F["n_G_tr_q90"], 1e-6) / q900)
    end
    F
end

## ---------------------------------------------------------------------------
## the approved commanded-object rule
## ---------------------------------------------------------------------------
# AKS approved these four thresholds on 2026-11-21/22 (the "recommended"
# operating point of 2026_11_19/classifier_2d/tuning_object/TUNING_OBJECT.md);
# the tr_dyn guard threshold was approved "for now". Archive-wide (out-of-fold
# forest scores) the rule acts on 437 commanded-object exposures.
"forest P(non-sky) = 1 - P(sky) at or above this -> non-sky content (approved 2026-11-21/22)"
const OBJECT_RULE_P_NONSKY_MIN = 0.8
"no-sky-light test, level: chip-G `oh_rel` strictly below this (approved 2026-11-21/22)"
const OBJECT_RULE_OH_REL_MAX = -0.6
"no-sky-light test, detection: chip-G `line_z` strictly below this (approved 2026-11-21/22)"
const OBJECT_RULE_LINE_Z_MAX = 2.2
"starlight guard: chip-G `tr_dyn` at or above this exempts the frame (approved 2026-11-21/22, \"for now\")"
const OBJECT_RULE_TR_DYN_GUARD = 0.9

_rule_float(x) = (x === missing || x === nothing) ? NaN : Float64(x)

"""
    object_frame_rule(p_nonsky, oh_rel, line_z, tr_dyn; p_min, oh_max, z_max, guard_min)

The approved rule for a COMMANDED-OBJECT exposure (pure function; the caller is
responsible for applying it only to `image_type == "object"`):

    act = (P(non-sky) >= 0.8  OR  (oh_rel < -0.6 AND line_z < 2.2))  AND NOT (tr_dyn >= 0.9)

Returns a NamedTuple:
- `act::Bool` — the frame should not be used for science;
- `reason::String` — which test decided it:
  `"forest"`, `"no_sky_light"`, `"forest+no_sky_light"` (act = true),
  `"guard_exempt"` (a test fired but the starlight guard exempted the frame),
  `"none"`;
- `forest`, `no_sky_light`, `guard::Bool` — the three tests individually;
- `inputs_finite::Bool` — all four inputs were finite.

NaN / `missing` / `nothing` inputs: a comparison with NaN is false, so a missing
input means THAT TEST DOES NOT FIRE, exactly as in the prototype that produced
the approved counts. In particular a NaN `tr_dyn` (50 or fewer trace peaks, e.g. a
blank or unlit frame) does not exempt the frame: no traces is no evidence of
starlight. A NaN P(non-sky) or `oh_rel` (forest or references unavailable) can
only make the rule act LESS; check `inputs_finite` before treating `act = false`
as a positive clearance.

Wiring (AKS 2026-11-22): `exposure_predicted_bad` returns this rule's `act`
for every commanded object, so EXPFLAG_PREDICTED_BAD (0x01) carries it; the
inputs and the test that fired are recorded as `exp_class_*` fields
(`exposure_object_rule_metadata`). The frames are still reduced to 1D.
"""
function object_frame_rule(p_nonsky, oh_rel, line_z, tr_dyn;
        p_min = OBJECT_RULE_P_NONSKY_MIN, oh_max = OBJECT_RULE_OH_REL_MAX,
        z_max = OBJECT_RULE_LINE_Z_MAX, guard_min = OBJECT_RULE_TR_DYN_GUARD)
    p, oh, z, g = _rule_float(p_nonsky), _rule_float(oh_rel), _rule_float(line_z),
    _rule_float(tr_dyn)
    forest = p >= p_min
    nsl = (oh < oh_max) & (z < z_max)
    guard = g >= guard_min
    act = (forest | nsl) & !guard
    reason = act ? (forest & nsl ? "forest+no_sky_light" : forest ? "forest" : "no_sky_light") :
             ((forest | nsl) ? "guard_exempt" : "none")
    (act = act, reason = reason, forest = forest, no_sky_light = nsl, guard = guard,
        inputs_finite = all(isfinite, (p, oh, z, g)))
end

## ---------------------------------------------------------------------------
## the forest
## ---------------------------------------------------------------------------
"""
    load_exposure_classifier_2d(model_path, templates_path, light_reference_path)

Load the 2D content forest (JLD2 with `model`, `classes`, `features`) and its
two stored references. Returns `nothing`, after a clear `@warn` that the 2D
content check is DISABLED, when any of the three files is absent or unreadable.
"""
function load_exposure_classifier_2d(model_path, templates_path, light_reference_path)
    if model_path === nothing || isempty(model_path) || !isfile(model_path)
        @warn "2D content check DISABLED: forest artifact not found ($(repr(model_path))). " *
              "Rebuild with scripts/classifier/train_classifier_2d.jl."
        return nothing
    end
    refs = load_classifier2d_references(templates_path, light_reference_path)
    refs === nothing && return nothing
    d = try
        JLD2.load(model_path)
    catch e
        @warn "2D content check DISABLED: could not read the forest artifact $model_path" exception = e
        return nothing
    end
    classes = String.(d["classes"])
    "sky" in classes || error("2D forest artifact $model_path has no \"sky\" class")
    (model = d["model"], classes = classes, features = String.(d["features"]), refs = refs,
        model_path = String(model_path))
end

"""
    classify_exposure_2d(clf2d, features::AbstractDict) -> NamedTuple

Score one exposure's feature Dict (`exposure_class_features_2d`) with the 2D
forest. Returns `(pred, prob, proba, p_nonsky, p_sky)`, with
`p_nonsky = 1 - P(sky)` computed exactly as the prototype did.
"""
function classify_exposure_2d(clf2d, features::AbstractDict)
    absent = [nm for nm in clf2d.features if !haskey(features, nm)]
    isempty(absent) || error("the 2D forest artifact wants features this code does not " *
          "compute: $(absent). A v1 artifact (with the read count n_read) is not " *
          "supported since v2 (AKS 2026-11-22); use a v2 artifact.")
    x = Float64[features[nm] for nm in clf2d.features]
    x[.!isfinite.(x)] .= -9999.0
    P = vec(apply_forest_proba(clf2d.model, reshape(x, 1, :), clf2d.classes))
    maxp, mi = findmax(P)
    psky = P[findfirst(==("sky"), clf2d.classes)]
    (pred = clf2d.classes[mi], prob = maxp, proba = P, p_sky = psky, p_nonsky = 1 - psky)
end

"""
    check_exposure_2d(clf2d, fnames::AbstractDict, tele, mjd, image_type)

Read the three ar2D chip files (READ-ONLY; `fnames` = chip => path), compute
the 2D features, score them, and for a commanded `object` apply
`object_frame_rule`. `clf2d === nothing` (artifact or references absent) returns
`nothing`: the caller records "not run" rather than a verdict.

Returns `(features, score, rule)` where `rule === nothing` for non-object frames.
"""
function check_exposure_2d(clf2d, fnames::AbstractDict, tele, mjd, image_type)
    clf2d === nothing && return nothing
    prods = Dict(chip => read_classifier2d_chip_products(fnames[chip], chip)
    for chip in CHIP_LIST)
    feats = exposure_class_features_2d(prods, tele, mjd; refs = clf2d.refs)
    score = classify_exposure_2d(clf2d, feats)
    rule = lowercase(strip(String(image_type))) == "object" ?
           object_frame_rule(score.p_nonsky, feats["oh_rel"], feats["n_G_line_z"],
        feats["n_G_tr_dyn"]) : nothing
    (features = feats, score = score, rule = rule)
end
