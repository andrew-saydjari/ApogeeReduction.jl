using FastRunningMedian: running_median
using Distributions: cdf, Normal
using Interpolations: linear_interpolation, Line
using DataFrames

# this file contains the code needed to extract a 1D spectrum from a 2D images.
# trace_params is of size (n_x_pix, n_fibers, 3)
# the elements correspodn to flux, y-pixel, and gaussian witdth (sigma)

# Per-(tele, mjd) cache of the exposure-type check table written between the 2D
# and 1D stages, so the 1D pass reads each night's table once per worker process
# rather than once per exposure per chip. A night with no table caches the empty
# Dict, so we do not stat a missing file repeatedly.
#
# n.b. this is a per-PROCESS cache filled under pmap (Distributed), NOT a
# per-thread buffer indexed by threadid(). Do not "optimize" it into one.
# Keyed on outdir too, not just (tele, mjd): a single process can legitimately
# be pointed at two reduction directories, and (tele, mjd) alone would silently
# serve one run's verdicts for the other's exposures.
const EXP_CLASS_CHECK_CACHE = Dict{
    Tuple{String, String, String}, Dict{Int, Dict{String, Any}}}()

"""
    exposure_class_metadata_for(outdir, tele, mjd, expnum)

The `exp_class_*` metadata block for one exposure, read from the per-MJD
exposure-type check table produced between the 2D and 1D stages.

Returns the explicit-unknown block whenever the classifier did not judge this
exposure: NO `exposure_flags` key, and `status = "notrun"`. The absence of the
bitmask is what marks it unjudged — a zeroed mask would read as "fine" and must
never stand in for "we did not look". Read it with `exposure_class_verdict`.
"""
function exposure_class_metadata_for(outdir, tele, mjd, expnum)
    tbl = get!(EXP_CLASS_CHECK_CACHE, (String(outdir), String(tele), String(mjd))) do
        read_exposure_type_check(exposure_type_check_path(outdir, tele, mjd))
    end
    get(tbl, Int(expnum), exposure_class_unknown_metadata())
end

# hold off on prop ivar through until we switch to sutr_wood, also could implement a chi2 cut here
# add a condition that we should drop any x pixel where a bad bit in any of the pixels being summed is bad

# ---------------------------------------------------------------------------
# Per-FIBER relative-throughput bitmask (`bitmsk_relthrpt`)
#
# This is the authoritative, per fiber-EPOCH statement of "is this fiber
# delivering light on this exposure?". It is produced by `get_relFlux` from a
# dome flat, stored in the `relFlux_*` cal products, copied into the per-chip
# `ar1D*` products by `process_1D`, and stacked into the resampled
# `ar1Duni*` products by `reinterp_spectra`. Downstream consumers (arMADGICS
# and the QA/census scripts) are expected to read it and mask on
# `RELTHRPT_UNUSABLE_BITS`.
#
# n.b. this axis is FIBER-level and per exposure. It is deliberately separate
# from the EXPOSURE-level `exp_class_*` / engineering-carton flags that live in
# the 1D product metadata: a fiber can be dead on an otherwise perfect
# exposure, and an exposure can be junk with every fiber healthy.
# ---------------------------------------------------------------------------

"Fiber throughput is low relative to the exposure median, but still usable."
const RELTHRPT_WARN_BIT = 2^0
"Fiber throughput is below `rel_val_cut`: the fiber is dead/near-dead."
const RELTHRPT_BROKEN_BIT = 2^1
"No fluxing (dome flat) file was available; `relthrpt` was forced to 1."
const RELTHRPT_NOFILE_BIT = 2^2
"""
`relthrpt` is not finite (NaN/Inf), i.e. the fiber was all-NaN or all-zero in
the flat, or the exposure-level normalization was itself degenerate.
Always accompanied by `RELTHRPT_BROKEN_BIT`, so bit-1 consumers are correct
without changes.
"""
const RELTHRPT_NOTFINITE_BIT = 2^3
"""
Too few pixels survived `bad_pix_bits` masking for the throughput median to be a
measurement rather than noise (fewer than `RELTHRPT_MIN_GOODPIX`). `relthrpt` is
forced to NaN, so this always also sets `RELTHRPT_NOTFINITE_BIT`,
`RELTHRPT_BROKEN_BIT` and `RELTHRPT_WARN_BIT`. Bit 3 is the SYMPTOM (the value is
not finite); this bit is the CAUSE.
"""
const RELTHRPT_LOWGOODPIX_BIT = 2^4

"""
Minimum number of pixels that must contribute to a fiber's throughput median.

DERIVED FROM THE DATA, not picked round. Sub-sampling the good pixels of real
domeflat fibers (2026_09_08 `floor_study.jl`, 42 flats spanning MJD 57652-61190,
both telescopes) gives the 95th-percentile fractional error of a median over `n`
pixels:

    n      16     32     64     96    128    192    256    512   1024
    err  0.223  0.167  0.123  0.106  0.089  0.075  0.064  0.044  0.025

256 is the smallest `n` at which that error (0.064) falls below `rel_val_cut`
(0.07), the sharpest cut this function makes -- i.e. the point below which
sampling noise alone could push a healthy fiber across the "broken" threshold.

It cannot fire on healthy data: over 892,800 fiber measurements on the DR21
200-MJD testbed the MINIMUM good-pixel count was 1805 of 2048 (APO) and 1903
(LCO), 7x above this floor, and zero fibers fell below it. This is insurance,
not a filter. It is also above the proportional analogue of arMADGICS'
`INGEST_MIN_GOODPIX` (500 of 8700 -> 118 of 2048).
"""
const RELTHRPT_MIN_GOODPIX = 256

"""
Bits which mean "do not trust this fiber's flux calibration at all".

This is the aggressive mask downstream analysis should use: a fiber carrying
any of these bits is NOT flux-scaled by `process_1D`, so its flux is on an
arbitrary scale and any chi2 computed against it is meaningless.

`RELTHRPT_NOFILE_BIT` is deliberately NOT in this set: in that case `relthrpt`
is forced to exactly 1 and the spectrum is simply unfluxed-but-unscaled, which
is a known and benign state, not a broken fiber.
"""
const RELTHRPT_UNUSABLE_BITS = RELTHRPT_BROKEN_BIT | RELTHRPT_NOTFINITE_BIT |
                               RELTHRPT_LOWGOODPIX_BIT

"""
    relthrpt_fiber_unusable(bitmsk_relthrpt)

`true` wherever the per-fiber relative-throughput bitmask says the fiber's flux
calibration is unusable (dead fiber or non-finite throughput). Broadcasts over
any shape, so it works on the `(300,)` per-chip vectors in `ar1D*` and on the
`(N_CHIPS, 300)` stacks in `ar1Duni*`.
"""
relthrpt_fiber_unusable(bitmsk_relthrpt) = (bitmsk_relthrpt .& RELTHRPT_UNUSABLE_BITS) .!= 0

"""
    relthrpt_fiber_fluxable(bitmsk_relthrpt)

`true` wherever `process_1D` will actually divide the fiber by `relthrpt`.
The single source of truth for which fibers get flux-scaled: `process_1D` uses
it, and downstream code can use it to reason about what was done.
"""
function relthrpt_fiber_fluxable(bitmsk_relthrpt)
    (bitmsk_relthrpt .& (RELTHRPT_UNUSABLE_BITS | RELTHRPT_NOFILE_BIT)) .== 0
end

"""
    get_relFlux(fname; sig_cut, rel_val_cut, use_pix_mask, min_goodpix)

Per-fiber relative throughput and its quality bitmask, from a flat exposure.

`use_pix_mask` (default `true`) takes the per-fiber throughput median over
pixels carrying no `bad_pix_bits` only. `min_goodpix` is the floor below which a
fiber is flagged rather than assigned a meaningless median; see
`RELTHRPT_MIN_GOODPIX`.
"""
function get_relFlux(fname; sig_cut = 4.5, rel_val_cut = 0.07,
        use_pix_mask::Bool = true, min_goodpix::Int = RELTHRPT_MIN_GOODPIX)
    f = jldopen(fname)
    flux_1d = f["flux_1d"]
    mask_1d = f["mask_1d"]
    mask_1d_good = (mask_1d .& bad_pix_bits) .== 0
    close(f)
    metadata = read_metadata(fname)

    # `mask_1d_good` was computed here and then never used, so the throughput median
    # ran over every pixel in the fiber. `nanzeromedian` drops NaN and exact zeros,
    # but a bad-yet-finite-nonzero pixel -- cosmic ray, saturation -- still
    # contributed to the number the whole flux scale is built on. It is now honoured.
    #
    # CAVEAT, ON THE RECORD: `bad_pix_bits` (24566) does NOT include
    # `pix_not_dark_corr_bits` (2^3 = 8), which is what the coherent APO chip-G
    # defect block at columns ~512-531, rows ~1362-1377 reads. Those pixels STILL
    # contribute to throughput after this change. That is expected: landing a
    # defect-region bit is a separate approved change, and it deliberately lands
    # OUTSIDE `bad_pix_bits` first. Revisit this line when it moves inside.
    flux_for_thrpt = if use_pix_mask
        masked = copy(flux_1d)
        masked[.!mask_1d_good] .= NaN
        masked
    else
        flux_1d
    end

    # How many pixels actually enter each fiber's median. Masking can starve a fiber,
    # and the median of a handful of pixels is noise, not throughput -- worse, the
    # median of an EMPTY set is NaN, and `NaN < thresh` is false, so without the
    # non-finite guard below such a fiber would be judged GOOD and then NaN-poison
    # flux and ivar downstream. Turning the mask on without that guard would be
    # actively worse than leaving it off.
    n_contrib = dropdims(sum(.!isnanorzero.(flux_for_thrpt), dims = 1), dims = 1)

    # `absthrpt` is the un-normalized per-fiber median. No pipeline stage consumes it;
    # it is read only by make_relFlux.jl's QA plots. Kept deliberately: it costs one
    # `dropdims` (it is the value `relthrpt` is derived from anyway) and it is the only
    # record of ABSOLUTE throughput, which is what you need to tell "this fiber died"
    # from "the whole exposure was faint". Documented here so it is not mistaken for
    # dead code again.
    absthrpt = dropdims(nanzeromedian(flux_for_thrpt, 1), dims = 1)
    bitmsk_relthrpt = zeros(Int, length(absthrpt))

    # Starved fibers are NaN'd BEFORE the exposure-level normalization below, so a
    # meaningless few-pixel median cannot contaminate the median-of-fibers that every
    # other fiber is divided by.
    starved = n_contrib .< min_goodpix
    absthrpt[starved] .= NaN
    bitmsk_relthrpt[starved] .|= (RELTHRPT_LOWGOODPIX_BIT | RELTHRPT_NOTFINITE_BIT |
                                  RELTHRPT_BROKEN_BIT | RELTHRPT_WARN_BIT)

    relthrpt = copy(absthrpt)
    relthrpt ./= nanzeromedian(relthrpt)

    # `nanzeromedian` returns NaN for a fiber that is entirely NaN-or-zero, and
    # `NaN < x` is `false` in Julia, so BOTH threshold tests below silently pass a
    # NaN fiber as good. `process_1D` would then divide flux and multiply ivar by
    # NaN and poison the fiber all the way downstream with nothing flagged. Catch
    # non-finite throughput explicitly and FIRST. Same for a degenerate
    # exposure-level normalization, which makes every fiber's relthrpt non-finite
    # and must flag the whole exposure's fibers rather than none of them.
    notfinite = .!isfinite.(relthrpt)
    bitmsk_relthrpt[notfinite] .|= (RELTHRPT_NOTFINITE_BIT | RELTHRPT_BROKEN_BIT |
                                    RELTHRPT_WARN_BIT)

    # `nanzeroiqr` is itself NaN when the whole exposure is NaN-or-zero, which would
    # make `thresh` NaN and flag nothing. In that case every fiber is already caught
    # by the non-finite test above, but guard the comparison anyway so the warn bit
    # is never decided by a NaN threshold.
    iqr_relthrpt = nanzeroiqr(relthrpt)
    thresh = isfinite(iqr_relthrpt) ? (1 - sig_cut * iqr_relthrpt) : NaN
    if isfinite(thresh)
        bitmsk_relthrpt[isfinite.(relthrpt) .& (relthrpt .< thresh)] .|= RELTHRPT_WARN_BIT
    end
    bitmsk_relthrpt[isfinite.(relthrpt) .& (relthrpt .< rel_val_cut)] .|= RELTHRPT_BROKEN_BIT
    return absthrpt, relthrpt, bitmsk_relthrpt, metadata
end

"""
Regularize the trace by applying a running median filter to each param in each fiber.
Could be denoised further by fitting a low-order polynomial or similar.
"""
function regularize_trace(trace_params; window_size = 101)
    @assert isodd(window_size) # otherwise the length of the regularized array is incorrect
    n_xpix, n_fibers = size(trace_params)[1:2]
    regularized_trace = similar(trace_params)
    for fiber in 1:n_fibers, param in 1:3
        regularized_trace[:,
        fiber,
        param] = running_median(
            trace_params[:, fiber, param], window_size, :asym_trunc; nan = :ignore)
    end
    regularized_trace
end

function extract_boxcar(dimage, ivarimage, pix_bitmask, trace_params;
        boxcar_halfwidth = 2, return_resids = false)
    if return_resids
        resid_fluxes_2d = zeros(Float64, size(dimage))
        resid_ivars_2d = zeros(Float64, size(dimage))
    end
    flux_1d = extract_boxcar_core(dimage, trace_params, boxcar_halfwidth)
    var_1d = extract_boxcar_core(1 ./ ivarimage, trace_params, boxcar_halfwidth)
    ivar_1d = 1.0 ./ var_1d
    mask_1d = extract_boxcar_bitmask(pix_bitmask, trace_params, boxcar_halfwidth)
    dropped_pixel_mask_1d = zeros(Int64, N_XPIX, N_FIBERS)

    if return_resids
        return flux_1d, ivar_1d, mask_1d, dropped_pixel_mask_1d, resid_fluxes_2d, resid_ivars_2d
    else
        return flux_1d, ivar_1d, mask_1d, dropped_pixel_mask_1d
    end
end

"""
Extract a 1D spectrum using a boxcar kernel with width from trace_params. This is used twice.
Once for the flux and once for the variance.
"""
function extract_boxcar_core(dimage_in, trace_params, boxcar_halfwidth)
    out = zeros(Float64, N_XPIX, N_FIBERS)
    n_xpix, n_fibers = size(trace_params)[1:2]
    for xpix in 1:n_xpix, fib in 1:n_fibers
        _, ypixf, _ = trace_params[xpix, fib, :]
        ypix = round(Int, ypixf)
        out[xpix, fib] = sum(dimage_in[xpix, (ypix - boxcar_halfwidth):(ypix + boxcar_halfwidth)])
    end
    out
end

function extract_boxcar_bitmask(dimage_in, trace_params, boxcar_halfwidth)
    mask = zeros(Int64, N_XPIX, N_FIBERS)
    n_xpix, n_fibers = size(trace_params)[1:2]
    for xpix in 1:n_xpix, fib in 1:n_fibers
        _, ypixf, _ = trace_params[xpix, fib, :]
        ypix = round(Int, ypixf)
        mask[xpix,
        fib] = reduce(
            |, dimage_in[xpix, (ypix - boxcar_halfwidth):(ypix + boxcar_halfwidth)])
    end
    mask
end

"""
Extract a 1D spectrum using the a gaussian kernel with center and width from trace_params.

# Keyword arguments:
- `window_size` is the number of pixels to sum over.

Spacing between traces is ~ 2000 / 300 = 6.67 pixels.
"""
function extract_optimal(dimage, ivarimage, pix_bitmask, trace_params; window_half_size = 4)
    n_xpix = size(trace_params, 1)
    n_fibers = size(trace_params, 2)

    # return values to be filled
    flux_1d = Matrix{Float64}(undef, n_xpix, n_fibers)
    ivar_1d = Matrix{Float64}(undef, n_xpix, n_fibers)
    mask_1d = Matrix{Int64}(undef, n_xpix, n_fibers)

    for xpix in 1:n_xpix, fib in 1:n_fibers
        _, y_peak, y_sigma = trace_params[xpix, fib, :]

        ypixels = floor(Int, y_peak - window_half_size):ceil(
            Int, y_peak + window_half_size)
        ypix_boundaries = [ypixels .- 0.5; ypixels[end] + 0.5]
        weights = diff(cdf.(Normal(y_peak, y_sigma), ypix_boundaries))

        flux_1d[xpix,
        fib] = sum(weights .* dimage[xpix, ypixels] .* ivarimage[xpix, ypixels]) /
               sum(weights .^ 2 .* ivarimage[xpix, ypixels])
        ivar_1d[xpix, fib] = sum(weights .^ 2 .* ivarimage[xpix, ypixels])

        # bitmask
        mask_1d[xpix, fib] = reduce(|, pix_bitmask[xpix, ypixels])
    end
    flux_1d, ivar_1d, mask_1d
end

function extract_optimal_iter(dimage, ivarimage, pix_bitmask, trace_params,
        med_center_to_fiber_func, x_prof_min, x_prof_max_ind,
        n_sub, min_prof_fib, max_prof_fib, all_y_prof, all_y_prof_deriv;
        small_window_half_size = 2, fit_window_half_size = 4,
        large_window_half_size = 12, n_max_repeat = 5, flag_thresh = 0.001,
        return_resids = false, neff_thresh = 10.0)
    n_xpix = size(trace_params, 1)
    n_ypix = size(dimage, 2)
    n_fibers = size(trace_params, 2)

    # return values to be filled
    flux_1d = Matrix{Float64}(undef, n_xpix, n_fibers)
    ivar_1d = Matrix{Float64}(undef, n_xpix, n_fibers)
    mask_1d = Matrix{Int64}(undef, n_xpix, n_fibers)
    # must be zero-initialized: elements are only assigned when pixels are
    # actually dropped, and undef Int64 memory otherwise leaks recycled heap
    # contents (garbage like +/-2^62, even negative "bitmask" values) into
    # dropped_pixels_mask_1d, nondeterministically run-to-run
    dropped_pixel_mask_1d = zeros(Int64, n_xpix, n_fibers)

    good_pixels = ((pix_bitmask .& bad_pix_bits) .== 0) .& (ivarimage .> 0)

    model_flux_indv = zeros(Float64, n_fibers, n_ypix)
    comb_model_flux = zeros(Float64, n_ypix)
    new_comb_model_flux = zeros(Float64, n_ypix)

    model_var_indv = zeros(Float64, n_fibers, n_ypix)
    comb_model_var = zeros(Float64, n_ypix)
    new_comb_model_var = zeros(Float64, n_ypix)

    new_flux_1d = zeros(Float64, n_fibers)

    if return_resids
        resid_fluxes_2d = zeros(Float64, size(dimage))
        resid_ivars_2d = zeros(Float64, size(dimage))
    end

    for xpix in 1:n_xpix
        #iterate on best-fit fluxes
        #to remove crosstalk/contribution
        #from neighbours

        #reset model fluxes
        model_flux_indv .= 0
        comb_model_flux .= 0
        model_var_indv .= 0
        comb_model_var .= 0
        new_flux_1d .= 0

        for repeat_ind in 1:n_max_repeat
            new_comb_model_flux .= 0
            new_comb_model_var .= 0

            if repeat_ind == 1
                window_half_size = small_window_half_size
            else
                window_half_size = fit_window_half_size
            end

            for fib in 1:n_fibers
                _, y_peak, y_sigma = trace_params[xpix, fib, :]
		y_peak_round = ceil(Int,round(y_peak))

		if ((y_peak_round + window_half_size) < 1) | ((y_peak_round - window_half_size) > N_XPIX) 
	            #then the peak is so far off the edge that
		    #we can't measure it's flux from the wings
		    #so give 0 flux and skip
                    new_flux_1d[fib] = 0.0
                    ivar_1d[xpix, fib] = 0.0
		    #propagate flags from pixels near the edge of the chip
		    #and also add new bad flags
		    ypixels = 1:large_window_half_size
                    mask_1d[xpix, fib] = reduce(|, pix_bitmask[xpix, ypixels])
                    mask_1d[xpix, fib] |= bad_1d_no_good_pix
                    mask_1d[xpix, fib] |= bad_1d_failed_extract
		    continue
		end

                #use large window to get model fluxes
		full_ypixels = max(1,y_peak_round - large_window_half_size):min(
			              N_XPIX,y_peak_round + large_window_half_size)
		y_peak_ind = y_peak_round - full_ypixels[begin] + 1
                full_ypix_boundaries = [full_ypixels .- 0.5; full_ypixels[end] + 0.5]
                #                full_model_vals = diff(cdf.(Normal(y_peak, y_sigma), full_ypix_boundaries))
                prof_fib_ind = clamp(fib, min_prof_fib, max_prof_fib)
                full_model_vals = diff(cdf_func_indv(full_ypix_boundaries, y_peak, y_sigma,
                    prof_fib_ind, x_prof_min, x_prof_max_ind,
                    n_sub, min_prof_fib, all_y_prof, all_y_prof_deriv))

		ypixels = full_ypixels[max(1,y_peak_ind - window_half_size):min(size(full_ypixels,1),y_peak_ind + window_half_size)]
		model_vals = full_model_vals[max(1,y_peak_ind - window_half_size):min(size(full_ypixels,1),y_peak_ind + window_half_size)]
                if any(good_pixels[xpix, ypixels])
                    #then mask the bad pixels, same as setting ivar=0 there
                    model_vals[.!good_pixels[xpix, ypixels]] .= 0
                    good_flux_1d = true
                else
                    good_flux_1d = false
                end

                curr_pix_ivars = 1 ./ (1 ./ ivarimage[xpix, ypixels] .+ comb_model_var[ypixels] .-
                                  model_var_indv[fib, ypixels])
                ivar_1d[xpix, fib] = sum(model_vals .^ 2 .* curr_pix_ivars)
                flux_weights = model_vals .* curr_pix_ivars ./ ivar_1d[xpix, fib]
                new_flux_1d[fib] = sum(flux_weights .*
                                       (dimage[xpix, ypixels] .- comb_model_flux[ypixels] .+
                                        model_flux_indv[fib, ypixels]))

                # bitmask
                #		curr_good_fluxes = flux_weights .>= flag_thresh
                curr_good_fluxes = model_vals .>= flag_thresh
                if !good_flux_1d
                    curr_neff = sqrt(1 / sum(model_vals .^ 2))
                    mask_1d[xpix, fib] = reduce(|, pix_bitmask[xpix, ypixels])
                    mask_1d[xpix, fib] |= bad_1d_no_good_pix
                elseif any(curr_good_fluxes)
                    curr_neff = sqrt(1 / sum(model_vals[curr_good_fluxes] .^ 2))
                    mask_1d[xpix, fib] = reduce(|, pix_bitmask[xpix, ypixels[curr_good_fluxes]])
                    if any(.!curr_good_fluxes)
                        dropped_pixel_mask_1d[xpix, fib] = reduce(
                            |, pix_bitmask[xpix, ypixels[.!curr_good_fluxes]])
                    end
                else
                    curr_neff = sqrt(1 / sum(model_vals .^ 2))
                    mask_1d[xpix, fib] = reduce(|, pix_bitmask[xpix, ypixels])
                end

                if (!isfinite(new_flux_1d[fib])) | (ivar_1d[xpix, fib] == 0.0)
                    new_flux_1d[fib] = 0.0
                    ivar_1d[xpix, fib] = 0.0
                    mask_1d[xpix, fib] |= bad_1d_failed_extract
                end

                if curr_neff > neff_thresh
                    #                    new_flux_1d[fib] = 0.0
                    #                    ivar_1d[xpix, fib] = 0.0
                    mask_1d[xpix, fib] |= bad_1d_neff
                end

                if good_flux_1d
                    model_flux_indv[fib, full_ypixels] .= max(0, new_flux_1d[fib]) * full_model_vals
                    new_comb_model_flux[full_ypixels] .+= model_flux_indv[fib, full_ypixels]

                    if ivar_1d[xpix, fib] > 0
                        model_var_indv[fib,
                        full_ypixels] .= max(0, 1 / ivar_1d[xpix, fib]) *
                                         (full_model_vals .^ 2)
                        new_comb_model_var[full_ypixels] .+= model_var_indv[fib, full_ypixels]
                    end
                end
            end

            flux_1d[xpix, :] .= new_flux_1d
            if return_resids
                resid_fluxes_2d[xpix, :] .= dimage[xpix, :] .- new_comb_model_flux
                resid_ivars_2d[xpix, :] .= 1 ./ (1 ./ ivarimage[xpix, :] .+ new_comb_model_var)
            end

            if all(abs.(new_flux_1d .- flux_1d[xpix, :]) .< 0.01) & (repeat_ind > 1)
                break
            end

            comb_model_flux .= new_comb_model_flux
            comb_model_var .= new_comb_model_var
        end
    end

    if return_resids
        return flux_1d, ivar_1d, mask_1d, dropped_pixel_mask_1d, resid_fluxes_2d, resid_ivars_2d
    else
        return flux_1d, ivar_1d, mask_1d, dropped_pixel_mask_1d
    end
end

"""
Return the FPI guide fiber IDs for one night, derived from the configurations
themselves rather than from a hardcoded constant.

The FPI feed has no positioner and therefore no FIBERMAP row in the
confSummary; almanac synthesizes a stub with `category == "bonus"` for every
APOGEE fiber_id in 1:300 that a configuration does not account for. Those stubs
carry no fiber_type, so they must be read before any `fiber_type == "APOGEE"`
filter.

The FPI feed does not move within a night, so the first object configuration
that identifies it is used. Returns an empty vector when the night has no such
configuration (a calibration-only night, the plate era, or a read failure) --
callers should then decline to label anything rather than guess, because a
stale constant silently overriding real fiber assignments is the failure mode
this replaces.
"""
function get_fpi_fiberIDs_from_almanac(almanac_file, tele, mjd)
    fiberIDs = Int[]
    isfile(almanac_file) || return fiberIDs
    try
        f = h5open(almanac_file, "r")
        try
            df_exp = read_almanac_exp_df(f, tele, mjd)
            if "config_id" in names(df_exp)
                for row in eachrow(df_exp)
                    (row.image_type == "object") || continue
                    config_id = row.config_id
                    (config_id > 0) || continue
                    fibers_path = "raw/$(tele)/$(mjd)/fibers/$(config_id)"
                    haskey(f, fibers_path) || continue
                    df_fib = DataFrame(read(f[fibers_path]))
                    rename!(df_fib, lowercase.(names(df_fib)))
                    ("category" in names(df_fib)) || continue
                    ids = df_fib[df_fib[!, "category"].=="bonus", "fiber_id"]
                    if !isempty(ids)
                        fiberIDs = sort(unique(Int.(ids)))
                        break
                    end
                end
            end
        finally
            close(f)
        end
    catch e
        @warn "Could not derive FPI guide fibers from $(almanac_file) for $(tele)/$(mjd); " *
              "they will not be annotated."
        show(e)
    end
    return fiberIDs
end

"""
Given an open HDF.file, `f`, and the telescope, mjd, and expnum, return a dictionary
mapping fiber index (1:300 laid out on the chip) to fiber type.
"""
function get_fibTargDict(f, tele, mjd, dfindx)
    # translate confSummary/almanac terminology to AR.jl terminology
    fiber_type_names = Dict(
        # fps era
        "science" => "sci",
        "sky_boss" => "skyB",
        "standard_apogee" => "tel",
        "sky_apogee" => "sky",

        # Plate era
        "science" => "sci",
        "standard" => "tel",
        "sky" => "sky"
    )
    # TODO other fiber types:
    # "blank"s from plate era
    # FPI era "serendipitous" APOGEE fibers are those which "accidentally" point at a bright
    # star (for BOSS reasons).
    # TODO Andrew thinks the fibers with category "" might be serendipitous targets

    mjdfps2plate = get_fps_plate_divide(tele)
    configIdCol = if parse(Int, mjd) > mjdfps2plate
        "config_id"
    else
        "plate_id"
    end

    df_exp = read_almanac_exp_df(f, tele, mjd)

    if !(dfindx in df_exp.exposure)
        @warn "Exposure $(dfindx) not found in $(tele)/$(mjd)/exposures"
        return Dict(1:300 .=> "fiberTypeFail"), Dict(1:300 .=> -2)
    end
    exposure_info = df_exp[dfindx, :]
    config_id = exposure_info[configIdCol]

    # FPI guide fibers, derived from the configuration below. Empty means either
    # that no configuration was read (calibration exposure, config_id == -1, or
    # a read failure) or that the configuration accounted for all 300 APOGEE
    # fibers, leaving no "bonus" stub for the FPI feed.
    fpi_fiberIndxs = Int[]
    read_config = false

    fibtargDict, fiber_sdss_id_Dict = if exposure_info.image_type != "object"
        (Dict(1:300 .=> "cal"), Dict(1:300 .=> -2))
    else
        # Check if config_id is -1, which should not exist in the HDF5 file
        if config_id == -1
            @warn "config_id is -1 for exposure $(dfindx) in $(tele)/$(mjd). This should have been filtered out as flagged_bad=1. Returning fiberTypeFail for all fibers."
            (Dict(1:300 .=> "fiberTypeFail"), Dict(1:300 .=> -2))
        else
            try
                df_fib = DataFrame(read(f["raw/$(tele)/$(mjd)/fibers/$(config_id)"]))
                # normalizes all column names to lowercase
                rename!(df_fib, lowercase.(names(df_fib)))

                # The FPI feed is a fixed illumination source: it has no
                # positioner and therefore no FIBERMAP row in the confSummary.
                # almanac notices the gap and synthesizes a stub row with
                # category "bonus" for every APOGEE fiber_id in 1:300 that the
                # configuration does not account for. Those stubs carry no
                # fiber_type, so they must be read BEFORE the "APOGEE" filter
                # below or they are silently discarded.
                read_config = true
                fpi_fiberIndxs = if configIdCol == "config_id"
                    fiberID2fiberIndx.(df_fib[df_fib[!, "category"].=="bonus", "fiber_id"])
                else
                    Int[]
                end

                # limit to only the APOGEE fiber/hole information
                df_fib = if configIdCol == "config_id"
                    df_fib[df_fib[!, "fiber_type"].=="APOGEE",:]
                else
                    df_fib
                end

                fiber_types = map(df_fib[!, "category"]) do t
                    if t in keys(fiber_type_names)
                        fiber_type_names[t]
                    else
                        # @warn "Unknown fiber type for $(tele)/$(mjd)/fibers/$(config_id): $(repr(t))"
                        "fiberTypeFail"
                    end
                end
                fibernumvec = df_fib[!, "fiber_id"]
                fiber_sdss_id = df_fib[!, "sdss_id"]

                #this is a Hack and Andy Casey will replace this very very soon
                msknofiberdefaults = (fibernumvec .!= -1)
                fiber_types_full = repeat(["fiberTypeFail"], N_FIBERS)
                fiber_sdss_id_full = repeat([-2], N_FIBERS)
                try
                    fiber_types_full[fiberID2fiberIndx.(fibernumvec[msknofiberdefaults])] .= fiber_types[msknofiberdefaults]
                    fiber_sdss_id_full[fiberID2fiberIndx.(fibernumvec[msknofiberdefaults])] .= fiber_sdss_id[msknofiberdefaults]
                catch e
                    @warn "Problem with getting fiber type information for $(tele)/$(mjd)/fibers/$(config_id) (exposure $(dfindx)). Returning fiberTypeFail for all fibers."
                    show(e)
                    fiber_types_full .= "fiberTypeFail"
                end
                (Dict(1:N_FIBERS .=> fiber_types_full), Dict(1:N_FIBERS .=> fiber_sdss_id_full))
            catch e
                @warn "Failed to get fiber type information for $(tele)/$(mjd)/fibers/$(config_id) (exposure $(dfindx)). Returning fiberTypeFail for all fibers."
                show(e)
                (Dict(1:300 .=> "fiberTypeFail"), Dict(1:300 .=> -2))
            end
        end
    end

    # FPS era only: label the FPI guide fibers from the configuration itself.
    #
    # This used to call a hardcoded per-telescope fiber pair (removed) and
    # overwrite those two fibers unconditionally. That is wrong whenever the FPI
    # feed does not sit where the constant says. Measured on the LCO FPS
    # commissioning nights MJD 59820/59826/59827 (16 configurations): the FPI is
    # on fiber_id 142/153, while the hardcoded 82/213 carry category "science" —
    # so 133 object exposures had two real science fibers relabelled "fpiguide"
    # and discarded, and the true FPI fibers were never labelled at all.
    # Confirmed in the extracted spectra: counting narrow emission peaks per
    # fiber, 142/153 rank 1-2 at ~485 peaks against a per-fiber median of 60-100,
    # while 82/213 sit mid-pack; on later LCO nights the ranking flips to 82/213,
    # matching "bonus" on those configurations.
    #
    # There is deliberately NO fallback to the hardcoded pair: a wrong constant
    # silently destroying science fibers is worse than no FPI label. If an FPS
    # configuration yields no "bonus" stub, warn loudly and label nothing. The
    # plate era has no configurations and no FPI, so it is silent by design.
    if parse(Int, mjd) > mjdfps2plate
        if isempty(fpi_fiberIndxs)
            if read_config
                @warn "No FPI guide fibers found for $(tele)/$(mjd)/fibers/$(config_id) " *
                      "(exposure $(dfindx)): the configuration accounts for all 300 APOGEE " *
                      "fibers, so almanac synthesized no \"bonus\" row. No fiber will be " *
                      "labelled fpiguide for this exposure."
            end
        else
            for fibindx in fpi_fiberIndxs
                fibtargDict[fibindx] = "fpiguide"
            end
        end
    end
    return (fibtargDict, fiber_sdss_id_Dict)
end

# hardcoded to use chip c only for now
# must use dome flats, not quartz flats (need fiber runs to telescope)
# use full exposure_id
function get_fluxing_file(dfalmanac, parent_dir, tele, mjd, dfindx, runname; fluxing_chip = "B")
    df_mjd = sort(
        dfalmanac[(dfalmanac.mjd .== parse(Int, mjd)) .& (dfalmanac.observatory .== tele), :],
        :exposure)
    expIndex = dfindx
    cartId = df_mjd.cart_id[expIndex]
    image_type = df_mjd.image_type[expIndex]

    valid_flats4fluxing_fname = joinpath(parent_dir, "almanac/valid_domeflats4fluxing_$(runname).h5")
    if !isfile(valid_flats4fluxing_fname)
        @warn "Could not find any useful relfluxing files after looking for file $(valid_flats4fluxing_fname)"
	    return 2^2, nothing
    end
    f = h5open(valid_flats4fluxing_fname, "r")
    found_tele_mjd = false
    if tele in keys(f)
        if "$(mjd)" in keys(f[tele])
            found_tele_mjd = true
	    end
    end

    if !found_tele_mjd
        close(f)
        if ((image_type == "object") | (image_type == "domeflat"))
            @warn "Could not find any useful relfluxing files in file $(valid_flats4fluxing_fname) for tele $(tele) mjd $(mjd)"
        end
        return 2^2,nothing
    end

    cal_expid_list = read(f["$(tele)/$(mjd)"])
    close(f)

    if dfindx in cal_expid_list
        #the current files is one of the dome flats that has a relfluxing file
        return 2^0,get_fluxing_file_name(
            parent_dir, tele, mjd, df_mjd.exposure[expIndex], fluxing_chip, cartId)
    end

    expIndex_before = findlast(cal_expid_list .< dfindx)
    if !isnothing(expIndex_before)
        expIndex_before = cal_expid_list[expIndex_before]
    end
    expIndex_after = findfirst(cal_expid_list .> dfindx)
    if !isnothing(expIndex_after)
        expIndex_after = cal_expid_list[expIndex_after]
    end

    valid_before = if isnothing(expIndex_before)
        0
    elseif all(df_mjd.cart_id[expIndex_before:expIndex] .== cartId)
	    1
    elseif !isnothing(expIndex_before) & (df_mjd.cart_id[expIndex_before] .== cartId)
        2
    else
        0
    end
    valid_after = if isnothing(expIndex_after)
        0
    elseif all(df_mjd.cart_id[expIndex:expIndex_after] .== cartId)
        1
    elseif !isnothing(expIndex_after) & (df_mjd.cart_id[expIndex_after] .== cartId)
        2
    else
        0
    end

    if valid_before == 1
        return 2^0, get_fluxing_file_name(
            parent_dir, tele, mjd, df_mjd.exposure[expIndex_before], fluxing_chip, cartId)
    elseif valid_after == 1
        return 2^0, get_fluxing_file_name(
            parent_dir, tele, mjd, df_mjd.exposure[expIndex_after], fluxing_chip, cartId)
        # any of the cases below here we could consider using a global file
    elseif valid_before == 2
        return 2^1, get_fluxing_file_name(
            parent_dir, tele, mjd, df_mjd.exposure[expIndex_before], fluxing_chip, cartId)
    elseif valid_after == 2
        return 2^1, get_fluxing_file_name(
            parent_dir, tele, mjd, df_mjd.exposure[expIndex_after], fluxing_chip, cartId)
    else
        return 2^2,nothing
    end
end

# TODO: switch to meta data dict and then save wavecal flags etc.
"""
    normalize_reinterp_spectra!(outflux, outvar, cntvec) -> (outivar, outmsk)

Final normalization step of `reinterp_spectra`: divide the accumulated flux and
variance by the per-fiber frame count (the max number of contributing frames
over the uniform-grid pixels) and derive the inverse variance and the Bool
good-pixel mask. `outflux` and `outvar` are mutated in place.

Fibers with no good pixels at all (`framecnts == 0`, e.g. dead/broken fibers)
come out as flux = 0, ivar = 0, msk = false. (A4 fix: they previously came out
as flux = 0/0 = NaN, ivar = NaN, msk = (0 .== 0) = true — NaN presented as
GOOD data, which NaN-poisoned downstream consumers such as the arMADGICS
ingest.)
"""
function normalize_reinterp_spectra!(outflux, outvar, cntvec)
    framecnts = maximum(cntvec, dims = 1) #     framecnts = maximum(cntvec) # a little shocked that I throw it away if it is bad in even one frame
    outflux ./= framecnts
    outvar ./= (framecnts .^ 2)
    # need to update this to a bit mask that is all or any for the pixels contributing to the reinterpolation
    outmsk = (cntvec .== framecnts) .& (framecnts .> 0)
    outivar = 1 ./ outvar
    outivar[.!outmsk] .= 0.0
    # zero-good-pixel fibers: replace the 0/0 = NaN flux/var with zeros
    zerofibs = dropdims(framecnts .== 0, dims = 1)
    outflux[:, zerofibs] .= 0.0
    outvar[:, zerofibs] .= 0.0
    return outivar, outmsk
end

function reinterp_spectra(fname, roughwave_dict; checkpoint_mode = "commit_same", outdir = "../outdir")
    # might need to add in telluric div functionality here?
    outname = replace(replace(fname, "ar1D" => "ar1Duni"), "_$(FIRST_CHIP)_" => "_")
    if check_file(outname, mode = checkpoint_mode)
        return
    end

    sname = split(split(split(fname, "/")[end], ".h5")[1], "_")
    fnameType, tele, mjd, expnum, chip, image_type = sname[(end - 5):end]
    mjd_int = parse(Int, mjd)

    backupWave_fname = joinpath(
        outdir, "wavecal", "wavecalNightAve_$(tele)_$(mjd).h5")

    # could shift this to a preallocation step
    outflux = zeros(length(logUniWaveAPOGEE), N_FIBERS)
    outvar = zeros(length(logUniWaveAPOGEE), N_FIBERS)
    outivar = zeros(length(logUniWaveAPOGEE), N_FIBERS)
    outmsk = zeros(Int, length(logUniWaveAPOGEE), N_FIBERS)
    outDropMsk = zeros(Int, length(logUniWaveAPOGEE), N_FIBERS)
    outTraceCoords = zeros(length(logUniWaveAPOGEE), N_FIBERS, 3) #(x,y,chipInt)
    cntvec = zeros(Int, length(logUniWaveAPOGEE), N_FIBERS)

    pixvec = 1:(N_CHIPS * N_XPIX)
    xpix_stack = ((pixvec .- 1) .% N_XPIX) .+ 1
    flux_stack = zeros(N_CHIPS * N_XPIX, N_FIBERS)
    ivar_stack = zeros(N_CHIPS * N_XPIX, N_FIBERS)
    mask_stack = zeros(Int, N_CHIPS * N_XPIX, N_FIBERS)
    dropped_mask_stack = zeros(Int, N_CHIPS * N_XPIX, N_FIBERS)
    wave_stack = zeros(N_CHIPS * N_XPIX, N_FIBERS)
    trace_center_stack = zeros(N_CHIPS * N_XPIX, N_FIBERS)
    chipInt_stack = zeros(N_CHIPS * N_XPIX, N_FIBERS)
    chipBit_stack = zeros(Int, N_CHIPS * N_XPIX, N_FIBERS)
    thrpt_stack = zeros(N_CHIPS, N_FIBERS)
    bitmsk_thrpt_stack = zeros(Int, N_CHIPS, N_FIBERS)

    ingestBit = zeros(Int, N_FIBERS)

    # add a for loop over the exposures (stop thinking about "visits" for now)
    # probably just generate ap1D file names from the alamanc files

    # this was used for looping over exposures in the visit
    # outdir = "/uufs/chpc.utah.edu/common/home/u6039752/scratch1/working/2024_12_05/outdir/"
    # fname = outdir * "apred/$(mjd)/" * get_1d_name(parse(Int, last(expid,4)), df) * ".h5"

    # The night's wavelength solution is chosen by preference order, but the FPI
    # solution is only a candidate if the night's FPI acceptance gate passed (see
    # src/fpi_gate.jl). `wavecalNightAve` records that decision in
    # `best_wave_type`. Consulting it matters because this loop keys on FILE
    # PRESENCE alone: a `waveCalNightfpiDither` left behind by an earlier run (or
    # by a checkpointed resume) would otherwise silently win over the sky
    # solution the gate chose.
    wavetype_order = if isfile(backupWave_fname)
        night_wave_type = try
            h5open(backupWave_fname, "r") do f
                haskey(f, "best_wave_type") ? read(f["best_wave_type"]) : "fpi"
            end
        catch
            "fpi"
        end
        # only ever *narrows* the order, and only on an explicit "sky" verdict
        night_wave_type == "sky" ? ["sky"] : ["fpi", "sky"]
    else
        ["fpi", "sky"]
    end
    found_soln = false
    wavecal_type = ""
    for wavetype in wavetype_order
        wavecal_type = "waveCalNight$(wavetype)Dither"
        wavefname = replace(replace(fname, fnameType => wavecal_type), "_$(FIRST_CHIP)_" => "_")
        if isfile(wavefname)
            f = jldopen(wavefname)
            chipWaveSoln = f["chipWaveSoln"]
            close(f)
	    if all(isfinite.(chipWaveSoln))
                found_soln = true
                break
	    end
        end
    end

    if !found_soln
        if !isfile(backupWave_fname)
	    curr_best_wave_type = "rough"
            backupWaveSoln = nothing
	    backupWave_fname = nothing
        else
            f = h5open(backupWave_fname, "r")
            curr_best_wave_type = read(f["best_wave_type"])
	    backupWaveSoln = read(f["$(curr_best_wave_type)/nightAve_wave_soln"])
	    close(f)
        end

        #this is not a great global fallback, but it works so we get something to look at
        if isnothing(backupWave_fname) || (!all(isfinite.(backupWaveSoln)))
            chipWaveSoln = zeros(N_XPIX, N_FIBERS, N_CHIPS)
            for (chipind, chip) in enumerate(CHIP_LIST)
                chipWaveSoln[:, :, chipind] .= rough_linear_wave.(
                    1:N_XPIX, a = roughwave_dict[tele][chip][1], b = roughwave_dict[tele][chip][2])
            end
            if !(image_type in ["dark", "internalflat", "quartzflat", "domeflat"])
                println("No wavecal found for $(fname), using rough linear fallback as fallback")
                flush(stdout)
            end
            wavecal_type = "error_fixed_fallback"
        else
            chipWaveSoln = backupWaveSoln
            if !(image_type in ["dark", "internalflat", "quartzflat", "domeflat"])
                println("No wavecal found for $(fname), using $(curr_best_wave_type) nightly average as fallback")
            end
            flush(stdout)
	        wavecal_type = "error_night_ave_$(curr_best_wave_type)"
        end
    end

    metadata_lst = []
    trace_lst = []
    for (chipind, chip) in enumerate(CHIP_LIST) # This needs to be the in abc RGB order, changing that will break this section
        fnameloc = replace(fname, "_$(FIRST_CHIP)_" => "_$(chip)_")
        f = jldopen(fnameloc)
        flux_1d = f["flux_1d"]
        ivar_1d = f["ivar_1d"]
        mask_1d = f["mask_1d"]
        dropped_pixels_mask_1d = f["dropped_pixels_mask_1d"]
        extract_trace_centers = f["extract_trace_centers"]
        relthrpt = f["relthrpt"]
        bitmsk_relthrpt = f["bitmsk_relthrpt"]
        close(f)
        push!(metadata_lst, read_metadata(fnameloc))

        flux_stack[(1:N_XPIX) .+ (3 - chipind) * N_XPIX, :] .= flux_1d[end:-1:1, :]
        ivar_stack[(1:N_XPIX) .+ (3 - chipind) * N_XPIX, :] .= ivar_1d[end:-1:1, :]
        mask_stack[(1:N_XPIX) .+ (3 - chipind) * N_XPIX, :] .= mask_1d[end:-1:1, :]
        dropped_mask_stack[(1:N_XPIX) .+ (3 - chipind) * N_XPIX, :] .= dropped_pixels_mask_1d[
            end:-1:1, :]
        try
            wave_stack[(1:N_XPIX) .+ (3 - chipind) * N_XPIX, :] .= chipWaveSoln[end:-1:1, :, chipind]
        catch
            println((typeof(wave_stack), typeof(N_XPIX), typeof(chipind), typeof(chipWaveSoln),fname, mjd_int, typeof(backupWaveSoln), typeof(backupWaveSoln)))
        end
        wave_stack[(1:N_XPIX) .+ (3 - chipind) * N_XPIX, :] .= chipWaveSoln[end:-1:1, :, chipind]
        trace_center_stack[(1:N_XPIX) .+ (3 - chipind) * N_XPIX, :] .= extract_trace_centers[
            end:-1:1, :]
        chipBit_stack[(1:N_XPIX) .+ (3 - chipind) * N_XPIX, :] .+= 2^(chipind)
        chipInt_stack[(1:N_XPIX) .+ (3 - chipind) * N_XPIX, :] .= chipind
        thrpt_stack[chipind, :] .= relthrpt
        bitmsk_thrpt_stack[chipind, :] .= bitmsk_relthrpt
    end

    # should add a check all entries of metadata_lst to be equal
    metadata = metadata_lst[1]
    metadata["wavecal_type"] = wavecal_type # add wavecal type to metadata

    noBadBits = (mask_stack .& bad_pix_bits .== 0)
    chipBit_stack[.!(noBadBits)] .+= 2^4 # call pixels with bad bits thrown bad
    chipBit_stack[chipBit_stack .== 0] .+= 2^4 # call missing chips bad

    # think about adding warnings for the last two cases
    good_pix = ((noBadBits) .& (chipBit_stack .& 2^4 .== 0)) .& (.!isnan.(ivar_stack)) .&
               (ivar_stack .> (10^-20))

    ## need to propagate the bit mask
    for fiberindx in 1:N_FIBERS
        good_pix_fiber = good_pix[:, fiberindx]
        flux_fiber = flux_stack[good_pix_fiber, fiberindx]
        ivar_fiber = ivar_stack[good_pix_fiber, fiberindx]
        wave_fiber = wave_stack[good_pix_fiber, fiberindx]
        trace_center_fiber = trace_center_stack[good_pix_fiber, fiberindx]
        chipBit_fiber = chipBit_stack[good_pix_fiber, fiberindx]
        chipInt_fiber = chipInt_stack[good_pix_fiber, fiberindx]
        pixindx_fiber = pixvec[good_pix_fiber]
        xpix_fiber = xpix_stack[good_pix_fiber]

        Rinv = generateInterpMatrix_sparse_inv(
            wave_fiber, chipBit_fiber, logUniWaveAPOGEE, pixindx_fiber)
        normvec = dropdims(sum(Rinv, dims = 2), dims = 2)
        msk_inter = (normvec .!= 0)

        outflux[msk_inter, fiberindx] .+= (Rinv * flux_fiber)[msk_inter]
        outvar[msk_inter, fiberindx] .+= ((Rinv .^ 2) * (1 ./ ivar_fiber))[msk_inter]
        cntvec[:, fiberindx] .+= msk_inter

        #right now, only works for a single exposure
        if length(wave_fiber) > 1
            # Check if wave_fiber has unique, sorted values
            if length(unique(wave_fiber)) != length(wave_fiber) || !issorted(wave_fiber)
                @warn "Non-unique or unsorted wavelengths for fiber $fiberindx in $fnameType $tele $mjd $expnum $chip $image_type. Cannot interpolate, filling trace coordinates with NaN."
                outTraceCoords[:, fiberindx, 1] .= NaN
                outTraceCoords[:, fiberindx, 2] .= NaN
                outTraceCoords[:, fiberindx, 3] .= NaN
            else
                outTraceCoords[:, fiberindx, 1] .= linear_interpolation(
                    wave_fiber, xpix_fiber, extrapolation_bc = Line()).(logUniWaveAPOGEE)
                outTraceCoords[:, fiberindx, 2] .= linear_interpolation(
                    wave_fiber, trace_center_fiber, extrapolation_bc = Line()).(logUniWaveAPOGEE)
                outTraceCoords[:, fiberindx, 3] .= linear_interpolation(
                    wave_fiber, chipInt_fiber, extrapolation_bc = Line()).(logUniWaveAPOGEE)
            end
        else
            # If no good pixels or only 1 good pixel for this fiber, fill with NaN
            if length(wave_fiber) == 0
                @warn "No good pixels found for fiber $fiberindx in $fnameType $tele $mjd $expnum $chip $image_type. Filling trace coordinates with NaN."
            else
                @warn "Only 1 good pixel found for fiber $fiberindx in $fnameType $tele $mjd $expnum $chip $image_type. Cannot interpolate, filling trace coordinates with NaN."
            end
            outTraceCoords[:, fiberindx, 1] .= NaN
            outTraceCoords[:, fiberindx, 2] .= NaN
            outTraceCoords[:, fiberindx, 3] .= NaN
        end

        if all(isnanorzero.(flux_fiber)) && ((ingestBit[fiberindx] & 2^1) == 0)
            ingestBit[fiberindx] += 2^1 # ap1D exposure flux are all NaNs (for at least one of the exposures)
        elseif all(.!((chipBit_fiber .& 2^4) .== 0)) && ((ingestBit[fiberindx] & 2^2) == 0)
            ingestBit[fiberindx] += 2^2 # ap1D exposure pixels are all bad by bit mask  (for at least one of the exposures)
        elseif all(isnanorzero.(ivar_fiber)) && ((ingestBit[fiberindx] & 2^3) == 0)
            ingestBit[fiberindx] += 2^3 # ap1D exposure ivars are all NaNs or zeros (for at least one of the exposures)
        end
    end

    outivar, outmsk = normalize_reinterp_spectra!(outflux, outvar, cntvec)

    # Write reinterpolated data.
    #
    # `ingestBit` is a SECOND per-fiber quality flag (bit 1: flux all NaN/zero; bit 2:
    # every pixel bad by bitmask; bit 3: ivars all NaN/zero) that this loop computes and
    # that, until now, was thrown away at the end of the function -- exactly the same
    # class of defect as the throughput flag never reaching a consumer. It is now
    # written out, as `bitmsk_ingest`. It is a DIFFERENT axis from `bitmsk_relthrpt`:
    # that one says the fiber's dome-flat throughput is dead, this one says the
    # extracted data itself is unusable. A consumer should read both.
    #
    # n.b. named `bitmsk_ingest`, NOT `ingestBit`: arMADGICS has its own per-spectrum
    # `ingestBit` column with an entirely different bit table, and two flags with one
    # name in adjacent products is a trap.
    safe_jldsave(
        outname, metadata; flux_1d = outflux, ivar_1d = outivar, mask_1d = outmsk,
        extract_trace_coords = outTraceCoords, relthrpt = thrpt_stack,
        bitmsk_relthrpt = bitmsk_thrpt_stack, bitmsk_ingest = ingestBit)
    return
end

const logUniWaveAPOGEE = 10 .^ range((start = 4.17825), step = 6.0e-6, length = 8700);

#should add a check_file call for this one
function process_1D(fname;
        outdir::String,
        runname::String,
        extraction::String,
        relFlux::Bool,
        chip_list::Vector{String} = CHIP_LIST,
        profile_path = "./data/",
        plot_path = "../outdir/$(sjd)/plots/",
        checkpoint_mode = "commit_same",
        per_chip_relflux::Bool = false)
    sname = split(split(split(fname, "/")[end], ".h5")[1], "_")
    fnameType, tele, mjd, expnum, chip, image_type = sname[(end - 5):end]
    dfindx = parse(Int, expnum)

    outfname = if relFlux
        replace(fname, "ar2D" => "ar1D")
    else
        replace(replace(fname, "ar2D" => "ar1D"), "apred" => "apredrelflux")
    end

    if check_file(outfname, mode = checkpoint_mode)
        return true
    end

    # this seems annoying to load so often if we know we are doing a daily... need to ponder
    traceFname = outdir * "apred/$(mjd)/traceMain_$(tele)_$(mjd)_$(chip).h5"

    # how worried should I be about loading this every time?
    falm = h5open(joinpath(outdir, "almanac/$(runname).h5"))
    dfalmanac = read_almanac_exp_df(falm, tele, mjd)

    (med_center_to_fiber_func, x_prof_min, x_prof_max_ind, n_sub, min_prof_fib, max_prof_fib,
    all_y_prof, all_y_prof_deriv) = get_default_trace_hyperparams(tele, chip, profile_path = profile_path, plot_path = plot_path)

    fnamecal = if (fnameType == "ar2D")
        replace(fname, "ar2D" => "ar2Dcal")
    else
        fname
    end

    dimage = load(fname, "dimage")
    ivarimage = load(fname, "ivarimage")
    pix_bitmask = load(fnamecal, "pix_bitmask")
    metadata = read_metadata(fname)

    regularized_trace_params = try
        load(traceFname, "regularized_trace_params")
    catch
        @warn "No regularized trace params found for $(traceFname)"
        return false
    end
    trace_metadata = read_metadata(traceFname)
    metadata = merge(metadata, trace_metadata)

    # Carry the exposure-type classifier's verdict into the 1D products, so a
    # downstream consumer reading a 1D file can see WHETHER the frame was judged
    # bad and WHY without re-deriving anything from the 2D products or the
    # almanac. Written into `metadata`, so it also propagates to the
    # wavelength-reinterpolated ar1Duni*/ar1Dunical* files, which inherit their
    # metadata from the first chip's ar1D* file (see `reinterp_spectra`).
    metadata = merge(metadata, exposure_class_metadata_for(outdir, tele, mjd, dfindx))

    flux_1d, ivar_1d,
    mask_1d,
    dropped_pixels_mask_1d,
    resid_flux,
    resid_ivar = if extraction == "boxcar"
        extract_boxcar(
            dimage, ivarimage, pix_bitmask, regularized_trace_params, return_resids = true)
    elseif extraction == "optimal"
        #            extract_optimal(dimage, ivarimage, pix_bitmask, regularized_trace_params)
        extract_optimal_iter(dimage, ivarimage, pix_bitmask, regularized_trace_params,
            med_center_to_fiber_func, x_prof_min, x_prof_max_ind, n_sub,
            min_prof_fib, max_prof_fib, all_y_prof, all_y_prof_deriv, return_resids = true)
    else
        error("Extraction method $(extraction) not recognized")
    end

    resid_outfname = replace(fname, "ar2D" => "ar2Dresiduals")
    safe_jldsave(resid_outfname, metadata; resid_flux, resid_ivar, trace_used_param_fname = traceFname)
    if relFlux
        # ###################################################################
        # WARNING -- KNOWN LIMITATION, PLEASE READ BEFORE TRUSTING `relthrpt`
        #
        # `make_relFlux.jl` computes and stores a SEPARATE `relthrpt` for every
        # chip, but this call hard-selects ONE chip (`chip_list[end]`, i.e. B)
        # and the resulting throughput is applied to R, G and B alike, a few
        # lines below. The per-chip `relFlux_*_<chip>.h5` symlink created below
        # is NAMED as though it were the chip's own fluxing solution; it is not.
        # Chromatic (per-chip) throughput differences are therefore unmodelled
        # BY CONSTRUCTION, and a fiber that is dead on one chip only will not be
        # flagged unless it is also dead on B.
        #
        # Changing this is a survey-wide science change (it moves the flux scale
        # of every R and G spectrum ever reduced), so it is left as-is and made
        # explicit rather than silently "fixed": the chip actually used is now
        # recorded in the product metadata as `relflux_chip`, and
        # `per_chip_relflux` below flips the behaviour for anyone who wants to
        # measure the difference. Do not flip it in production without a
        # before/after comparison.
        # ###################################################################
        fluxing_chip = per_chip_relflux ? chip : chip_list[end]
        # this is the path to the underlying fluxing file.
        # it is symlinked below to an exposure-specific file (linkPath).
        relflux_bit,calPath = get_fluxing_file(
            dfalmanac, outdir, tele, mjd, dfindx, runname, fluxing_chip = fluxing_chip)
        fibtargDict, fiber_sdss_id_Dict = get_fibTargDict(falm, tele, mjd, dfindx)
        fiberTypeList = map(x -> fibtargDict[x], 1:300)

        if isnothing(calPath)
            # TODO uncomment this
            if (image_type == "object") | (image_type == "domeflat")
                @warn "No fluxing file available for $(tele) $(mjd) $(dfindx) $(chip)"
            end
            relthrpt = ones(size(flux_1d, 2))
            # `relthrptr` must exist on every branch: it used to be defined only in
            # the `else` below and this branch survived purely because no bit-2 fiber
            # is ever fluxable. That is an accident, not a design.
            relthrptr = reshape(relthrpt, (1, length(relthrpt)))
            bitmsk_relthrpt = RELTHRPT_NOFILE_BIT * ones(Int, size(flux_1d, 2))
        elseif !isfile(calPath)
            error("Fluxing file $(calPath) for $(tele) $(mjd) $(dfindx) $(chip) does not exist")
        else
            linkPath = abspath(joinpath(
                dirname(fname), "relFlux_$(tele)_$(mjd)_$(dfindx_fname_format(dfindx))_$(chip).h5"))
            if !islink(linkPath)
                symlink(abspath(calPath), linkPath)
            end
            relthrpt = load(linkPath, "relthrpt")
            relthrptr = reshape(relthrpt, (1, length(relthrpt)))
            bitmsk_relthrpt = load(linkPath, "bitmsk_relthrpt")
        end

        # don't flux broken fibers (don't use warn fibers for sky)
        # broken fibers (bit 1, relthrpt < rel_val_cut) always also have the low-throughput
        # warn bit 0 set, so the mask must exclude on bit 1 explicitly: their relthrpt is a
        # noise-level (possibly negative) domeflat measurement and dividing by it produces
        # arbitrarily inflated flux. Bit 2 (no fluxing file) is excluded as before, and
        # bit 3 (non-finite relthrpt) would divide by NaN.
        #
        # A fiber excluded HERE is left on an arbitrary flux scale -- it is not zeroed,
        # not NaN'd and not masked, by design, so nothing is dropped from the reduction.
        # That is exactly why `bitmsk_relthrpt` has to travel with the data: it is the
        # only record that this fiber's flux is uncalibrated, and any downstream chi2
        # computed against it is meaningless.
        msk_goodwarn = relthrpt_fiber_fluxable(bitmsk_relthrpt)
        if any(msk_goodwarn)
            flux_1d[:, msk_goodwarn] ./= relthrptr[:, msk_goodwarn]
            ivar_1d[:, msk_goodwarn] .*= relthrptr[:, msk_goodwarn] .^ 2
        end

	metadata["bitmsk_relFluxFile"] = relflux_bit
        # Which chip's throughput solution was actually applied to this chip's data.
        # Equal to `chip` only when `per_chip_relflux` is on; see the WARNING above.
        metadata["relflux_chip"] = fluxing_chip
        # we probably want to append info from the fiber dictionary from alamanac into the file name
        safe_jldsave(outfname, metadata; flux_1d, ivar_1d, mask_1d, dropped_pixels_mask_1d,
            extract_trace_centers = regularized_trace_params[:, :, 2],
            relthrpt, bitmsk_relthrpt, fiberTypeList,
            trace_used_param_fname = traceFname)
    else
        dirName = dirname(outfname)
        if !ispath(dirName)
            mkpath(dirName)
        end
        safe_jldsave(
            outfname, metadata; flux_1d, ivar_1d, mask_1d, dropped_pixels_mask_1d,
            extract_trace_centers = regularized_trace_params[:, :, 2],
            trace_used_param_fname = traceFname)
    end
    close(falm)
    return true
end
