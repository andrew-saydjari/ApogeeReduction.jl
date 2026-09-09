# Acceptance gate for the FPI-based nightly wavelength solution.
#
# Motivation (measured, 2026-09-08, 200-MJD DR21 testbed):
# On five LCO nights (59922, 59942, 59947, 60009, 60031) plus lco/apo 59894 and
# apo 60255 the FPI produced no usable light. `get_initial_fpi_peaks` seeds a
# regular peak grid from `read_fpiPeakLoc_coeffs` and *always* returns that grid,
# so the fitter came back with a full peak list whose amplitudes were 0.1-0.6
# counts (723-1039 on healthy nights) and whose widths were pinned at the 3.0
# fitting ceiling. `comb_exp_get_and_save_fpi_wavecal` then re-solved the etalon
# cavity from that noise and `fpi_medwavecal_skyline_dither` promoted the result
# to `best_wave_type = "fpi"` unconditionally and silently. The delivered 1D
# spectra on those nights carry a wavelength solution wrong by ~1.2 px RMS with
# +-3 px tails.
#
# Two guards, both required to keep the FPI solution:
#
#   G1 (primary, "is there FPI light in the FPI fibers?")
#       Runs on the `fpiPeaks` products, BEFORE the FPI wavecal is computed.
#       Per fiber, the fraction of seeded peak slots that came back as a real,
#       narrow, significant emission peak (height SNR >= FPI_PEAK_SNR_MIN and
#       fitted width clear of the FPI_WIDTH_CEILING). The night's statistic is
#       the median over fibers, minimised over exposures and chips, and further
#       minimised with the same quantity restricted to the FPI guide fibers
#       named by the almanac for THAT night.
#
#   G2 (secondary) FPI fit residual RMS, in Angstroms, over the peaks the fit
#       actually used. AR already computes this and throws it away.
#
# On failure the night is NOT failed: the nightly sky solution
# (`waveCalNightskyDither`, written by `skyline_medwavecal_skyline_dither`) is
# kept, a loud warning naming the night, the failing guard and its value is
# emitted, and the decision plus every measured statistic is written into
# `wavecalNightAve_<tele>_<mjd>.h5` under `fpi_qa/` so the choice is auditable
# after the fact.

"""Per-peak height SNR below which a fitted FPI peak is not a detection."""
const FPI_PEAK_SNR_MIN = 10.0

"""`max_widths` passed to `fit_gauss_and_bias` in `get_initial_fpi_peaks`."""
const FPI_WIDTH_CEILING = 3.0

"""A fitted width at or above this is treated as pinned at the ceiling, i.e. not a line."""
const FPI_WIDTH_MAX = 0.9 * FPI_WIDTH_CEILING

"""
G1 threshold on the lit fraction.

Measured on the 200-MJD DR21 testbed (155 telescope-nights with `fpiPeaks`):
healthy nights span 0.969-1.000 (n = 142, both telescopes, plate and FPS eras);
the eight FPI-dark telescope-nights span 0.000-0.109. The corpus is empty
between 0.109 and 0.969, so 0.5 sits in the middle of a factor-9 gap: 1.94x
below the worst healthy night and 4.6x above the best rejected one.
"""
const FPI_LIT_FRACTION_MIN = 0.5

"""
G2 threshold on the FPI fit residual RMS, Angstroms.

Measured on the same corpus (119 telescope-nights that produced a `waveCalFPI`):
healthy LCO 0.0041-0.0078 A, healthy APO 0.0201-0.0300 A (APO is intrinsically
~5x larger; one threshold still serves both). The dead-FPI nights are
0.422-0.468 A (LCO) and 0.284 A (apo 60255). 0.10 A is 3.3x above the worst
healthy night and 2.8x below the smallest rejected one, and the corpus is empty
between 0.032 and 0.284.

Note apo 59894 lands at 0.032 A and is NOT caught by G2; it is caught by G1
(lit fraction 0.070). The guards are complementary by design.
"""
const FPI_RESID_RMS_MAX = 0.10

"""
    fpi_peak_lit_fraction(fname) -> Vector{Float64}

For one `fpiPeaks_<tele>_<mjd>_<expid>_<chip>_arclamp.h5` file, the per-fiber
fraction of seeded FPI peak slots that came back as a genuine narrow emission
line. This is the peak-counting metric that separates a live FPI from a dark
one, expressed on the product the wavecal actually consumes.

Returns a length-`N_FIBERS` vector in [0, 1]; a fiber with no seeded peaks
scores 0.
"""
function fpi_peak_lit_fraction(fname)
    fpi_line_mat, fpi_line_cov_mat = h5open(fname, "r") do f
        read(f["fpi_line_mat"]), read(f["fpi_line_cov_mat"])
    end
    n_peaks, _, n_fibers = size(fpi_line_mat)
    out = zeros(Float64, n_fibers)
    for fibIndx in 1:n_fibers
        n_seed = 0
        n_lit = 0
        for peak_ind in 1:n_peaks
            # a NaN center means the slot was never populated for this fiber
            isfinite(fpi_line_mat[peak_ind, 2, fibIndx]) || continue
            n_seed += 1
            height = fpi_line_mat[peak_ind, 1, fibIndx]
            width = fpi_line_mat[peak_ind, 3, fibIndx]
            height_var = fpi_line_cov_mat[peak_ind, 1, 1, fibIndx]
            (isfinite(height_var) && (height_var > 0)) || continue
            snr = height / sqrt(height_var)
            if isfinite(snr) && (snr >= FPI_PEAK_SNR_MIN) &&
               isfinite(width) && (width < FPI_WIDTH_MAX)
                n_lit += 1
            end
        end
        out[fibIndx] = (n_seed == 0) ? 0.0 : n_lit / n_seed
    end
    return out
end

"""
    fpi_light_check(fpiPeaks_fname_list; fpi_fiberIndxs = Int[]) -> NamedTuple

G1. Evaluate every chip of every FPI exposure of the night and return

  `lit`         the gating statistic, `min(lit_all, lit_fpifib)`
  `lit_all`     worst over (exposure, chip) of the median-over-all-fibers lit fraction
  `lit_fpifib`  the same restricted to `fpi_fiberIndxs`; `NaN` when that list is empty
  `n_files`     number of `fpiPeaks` files actually read
  `pass`        `lit >= FPI_LIT_FRACTION_MIN`

`fpi_fiberIndxs` are fiber INDICES (`301 - fiber_id`), derived per night from the
almanac by `get_fpi_fiberIDs_from_almanac`; they are never hardcoded, because the
LCO FPI sat on fiber_id 142/153 during MJD 59810-59850 rather than the usual
82/213 and a gate keyed to the wrong fibers would be worse than no gate.

Both statistics are needed. Measured on the testbed corpus: `lit_fpifib` alone
misses apo 60255 (0.994 on the guide fibers while the bulk of the fibers sit at
0.085) and is undefined on lco 59894 (the almanac names no `bonus` stub for that
night); `lit_all` alone would miss a night where only the guide feed went dark.
Taking the minimum costs nothing on this corpus -- the worst healthy value of
either statistic is 0.969.
"""
function fpi_light_check(fpiPeaks_fname_list; fpi_fiberIndxs = Int[])
    lit_all = Inf
    lit_fpifib = Inf
    n_files = 0
    for fname in unique(fpiPeaks_fname_list)
        for chip in CHIP_LIST
            fname_chip = replace(fname, "_$(FIRST_CHIP)_" => "_$(chip)_")
            isfile(fname_chip) || continue
            lf = try
                fpi_peak_lit_fraction(fname_chip)
            catch e
                @warn "FPI gate: could not read $(fname_chip); treating it as unlit." exception=e
                zeros(Float64, N_FIBERS)
            end
            n_files += 1
            lit_all = min(lit_all, nanmedian(lf))
            if !isempty(fpi_fiberIndxs)
                keep = filter(i -> 1 <= i <= length(lf), fpi_fiberIndxs)
                isempty(keep) || (lit_fpifib = min(lit_fpifib, nanmedian(lf[keep])))
            end
        end
    end
    if n_files == 0
        return (lit = NaN, lit_all = NaN, lit_fpifib = NaN, n_files = 0, pass = false)
    end
    isfinite(lit_fpifib) || (lit_fpifib = NaN)
    lit = isnan(lit_fpifib) ? lit_all : min(lit_all, lit_fpifib)
    return (lit = lit, lit_all = lit_all, lit_fpifib = lit_fpifib,
        n_files = n_files, pass = isfinite(lit) && (lit >= FPI_LIT_FRACTION_MIN))
end

"""
    fpi_resid_stats(waveCalFPI_fname) -> NamedTuple

G2. RMS and MAD of the FPI wavelength-fit residuals (Angstroms) over the peaks
`comb_exp_get_and_save_fpi_wavecal` actually used, plus `n_used` and the
fraction of peaks that survived the fit. `resid_vec` is stored as
(N_FIBERS, n_peaks) while `resid_used_in_fit` is (n_peaks, N_FIBERS), hence the
transpose.

`n_used == 0` is its own failure mode and is reported separately rather than as
an RMS of NaN: the fit converged on nothing at all, so the saved solution is
whatever the seed produced. Measured example: apo 59573 in the DR21 testbed has
zero peaks used and `fpi_m0 = 0` (healthy nights: ~4409) yet was still promoted
to `best_wave_type = "fpi"`.
"""
function fpi_resid_stats(waveCalFPI_fname)
    isfile(waveCalFPI_fname) ||
        return (rms = NaN, mad = NaN, n_used = 0, frac_used = NaN, pass = false)
    resid_vec, used = h5open(waveCalFPI_fname, "r") do f
        haskey(f, "resid_vec") || return (nothing, nothing)
        (read(f["resid_vec"]),
            haskey(f, "resid_used_in_fit") ? read(f["resid_used_in_fit"]) : nothing)
    end
    isnothing(resid_vec) &&
        return (rms = NaN, mad = NaN, n_used = 0, frac_used = NaN, pass = false)
    r = vec(resid_vec)
    frac_used = NaN
    if !isnothing(used)
        u = vec(permutedims(used))
        if length(u) == length(r)
            r = r[u .!= 0]
            frac_used = count(!=(0), used) / length(used)
        else
            @warn "FPI gate: resid_used_in_fit shape does not match resid_vec in $(waveCalFPI_fname); using all residuals."
        end
    end
    r = filter(isfinite, r)
    n_used = length(r)
    n_used == 0 &&
        return (rms = NaN, mad = NaN, n_used = 0, frac_used = frac_used, pass = false)
    rms = sqrt(sum(abs2, r) / n_used)
    mad = median(abs.(r))
    return (rms = rms, mad = mad, n_used = n_used, frac_used = frac_used,
        pass = isfinite(rms) && (rms <= FPI_RESID_RMS_MAX))
end

"""
    save_fpi_qa!(wavecalNightAve_fname, qa::Dict)

Write the FPI gate verdict and every statistic behind it into the night's
`wavecalNightAve` file under `fpi_qa/`, replacing any previous copy, so the
choice of `best_wave_type` is auditable without re-deriving anything from the
peak files. Never throws: QA bookkeeping must not be able to fail a night.
"""
function save_fpi_qa!(wavecalNightAve_fname, qa::AbstractDict)
    try
        h5open(wavecalNightAve_fname, "r+") do f
            haskey(f, "fpi_qa") && delete_object(f, "fpi_qa")
            g = create_group(f, "fpi_qa")
            for k in sort(collect(keys(qa)))
                g[k] = qa[k]
            end
        end
    catch e
        @warn "FPI gate: could not record QA in $(wavecalNightAve_fname)." exception=e
    end
    return nothing
end
