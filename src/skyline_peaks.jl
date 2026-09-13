import FastRunningMedian: running_median
using Optim, HDF5

# wavecal_class codes stored in the skyLinePeaks product (numeric for portability):
# B = associate + fit + use in the wavelength solution;
# A = associate + fit only (association guard; NaN-masked at ingest, see
#     ingest_skyLines_file); X / unclassed rows are never selected.
const WAVECAL_CLASS_CODES = Dict("X" => 0, "A" => 1, "B" => 2)

"""
    rough_wave(pix, rw)

Rough wavelength (vacuum A) at pixel `pix` from a roughwave_dict entry `rw`.
Entries carrying the quartic model (9 elements; see
scripts/cal/make_roughwave_dict.jl) use it -- max error < 0.004 A vs the adopted
solutions, where the old linear model erred by 5-17 A with chip-dependent
curvature (the root cause of sky-line association thefts). 4-element (legacy)
entries fall back to the linear model.

This model is a SEED, never load-bearing for correctness, and its terms carry
very different trust levels (measured across 90 pass-1 nights, 2016-2025):
the ZEROPOINT (c0) moves ~0.1 A with every dither shift, 0.3-2.6 A night to
night, and up to ~10 A across the era (APO re-registrations), so it is never
trusted -- `get_sky_peaks` re-derives it per fiber/chip/exposure from the
data (global offset scan, +-15 A). The SHAPE terms are 100-1000x more stable
(c1 edge effect: ~8 mA within a night, 0.05-0.19 A night to night; c2+:
single mA within-night), which is what justifies seeding them -- and even
they are refined per fiber by the self-calibrating association correction
before anything depends on them.
"""
function rough_wave(pix, rw)
    if length(rw) >= 9
        x = (pix - (N_XPIX ÷ 2)) / N_XPIX
        return rw[5] + rw[6] * x + rw[7] * x^2 + rw[8] * x^3 + rw[9] * x^4
    else
        return rough_linear_wave(pix; a = rw[1], b = rw[2])
    end
end

"""
    rough_dispersion(pix, rw)

Local dispersion d(lambda)/d(pix) (A per pixel, negative) of the rough wavelength
model at pixel `pix` -- the quartic's analytic derivative, or the global linear
slope for legacy 4-element entries. Used to convert catalog doublet separations
to pixels at the line's own location (the global slope is up to ~5% off locally).
"""
function rough_dispersion(pix, rw)
    if length(rw) >= 9
        x = (pix - (N_XPIX ÷ 2)) / N_XPIX
        return (rw[6] + 2 * rw[7] * x + 3 * rw[8] * x^2 + 4 * rw[9] * x^3) / N_XPIX
    else
        return rw[2]
    end
end

"""
    assign_segments(pred_wavs, bright_wavs; accept_radius = 10.0)

Nearest-catalog-line assignment of detected peak segments given their
predicted wavelengths: each segment goes to its nearest catalog line, each
line keeps only its closest segment (`is_closest` deduplication), and
assignments farther than `accept_radius` are rejected. Returns
`(nearest_idx, d2th_wavs, msk2use)`.
"""
function assign_segments(pred_wavs, bright_wavs; accept_radius = 10.0)
    n = length(pred_wavs)
    nearest_idx = zeros(Int, n)
    d2th_wavs = zeros(Float64, n)
    for (i, w) in enumerate(pred_wavs)
        nearest_idx[i] = argmin(abs.(w .- bright_wavs))
        d2th_wavs[i] = abs(w - bright_wavs[nearest_idx[i]])
    end
    is_closest = zeros(Bool, n)
    for nind in unique(nearest_idx)
        matches = findall(x -> x == nind, nearest_idx)
        is_closest[matches] .= false
        is_closest[matches[argmin(d2th_wavs[matches])]] = true
    end
    return nearest_idx, d2th_wavs, is_closest .& (d2th_wavs .<= accept_radius)
end

function get_sky_peaks(flux_vec, tele, chip, roughwave_dict, df_sky_lines;
				 med_flux_window = 31, flux_thresh = 97,
				 max_pix_sep = 5, n_pad = 5, min_seg_length = 2,
				 self_calibrate = true, max_assoc_iter = 5, assoc_loo_nsigma = 7.0,
				 assoc_only = false)

    #use running median to help identify skyline peaks in data
    #(especially for bright stars)
    med_flux_vec = running_median(
        flux_vec, med_flux_window, :asym_trunc, nan = :ignore)
    med_flux_vec[med_flux_vec .== 0.0] .= 1.0
    scaled_flux_vec = flux_vec ./ med_flux_vec

    # Find indices where flux is above flux_thresh percentile
    thresh = nanzeropercentile(scaled_flux_vec, percent_vec = [flux_thresh])[1]
    above_thresh = findall(x -> x > thresh, scaled_flux_vec)

    if size(above_thresh,1) < 1
        return nothing, nothing, nothing, (0, true, 0)
    end

    # Group indices into segments, combining those less than max_pix_sep pixels apart
    segments = []
    current_segment = [above_thresh[1]]
    for i in 2:length(above_thresh)
        if above_thresh[i] - above_thresh[i - 1] <= max_pix_sep
            push!(current_segment, above_thresh[i])
        else
            push!(segments, current_segment)
            current_segment = [above_thresh[i]]
        end
    end
    push!(segments, current_segment)
    mean_x = mean.(segments)
    length_segs = length.(segments)
    segments = segments[(mean_x .> 64) .& (mean_x .< 1984) .& (length_segs .>= min_seg_length)]

    # Preallocate array for segment fluxes
    segment_fluxes = zeros(length(segments))

    # For each segment, compute flux in padded range (using n_pad) after subtracting median
    for (i, segment) in enumerate(segments)
        # Get range with n_pad pixel padding on each side
        start_idx = maximum([1, minimum(segment) - n_pad])
        end_idx = minimum([length(flux_vec), maximum(segment) + n_pad])

        # Get flux in range
        flux_range = flux_vec[start_idx:end_idx]

        # Subtract median of range
        med = median(flux_range)
        flux_range .-= med

        # Sum total flux in range
        segment_fluxes[i] = sum(flux_range)
    end

    rw = roughwave_dict[tele][chip]
    seg_rough_wavs = map(x -> rough_wave(mean(x), rw), segments)

    ######

    wav_range = (rw[4], rw[3])
    # Select the association set: classed rows of data/APOGEE_sky_linelist.csv
    # (B = fit + use in solution, A = association guard, fit only) inside this
    # chip's wavelength range. The strength/cleanliness cuts are frozen in the
    # generated file (see scripts/cal/make_sky_linelist.jl), not applied here.
    bright_lines = sort(
        filter(
            row -> (wav_range[1] <= row.wave_cen_ang <= wav_range[2]) &&
                   (coalesce(row.wavecal_class, "") in ("A", "B")),
            df_sky_lines),
        :wave_cen_ang
    )

    bright_wavs = collect(bright_lines.wave_cen_ang)

    # Find optimal offset by trying different values
    function calc_total_diff(offset)
        # Add offset to rough wavelengths
        adjusted_wavs = seg_rough_wavs .+ offset

        # For each adjusted wavelength, find distance to nearest bright line
        diffs = zeros(Float64, length(adjusted_wavs))
        ninds = zeros(Int, length(adjusted_wavs))
        for (i, w) in enumerate(adjusted_wavs)
            ninds[i] = argmin(abs.(w .- bright_wavs))
            diffs[i] = abs(w - bright_wavs[ninds[i]])
        end

        is_closest = zeros(Bool, length(diffs))
        for nind in unique(ninds)
            matches = findall(x -> x == nind, ninds)
            is_closest[matches] .= false
            closest_match = matches[argmin(diffs[matches])]
            is_closest[closest_match] = true
        end

        # Calculate sigma using nanzeroiqr
        sigma = jack_std(diffs[(abs.(diffs) .<= 10) .& is_closest])

        # Filter out differences more than 10 sigma and more than 15 Å
        good_diffs = is_closest .& (abs.(diffs) .<= 10 * sigma) .& (abs.(diffs) .<= 10)

        # Return number of good differences and sum of squares
        # return count(good_diffs), mean(diffs[good_diffs].^2)
        return count(good_diffs), mean(diffs[good_diffs] .^ 2), sigma
    end

    # Try range of offsets around 0
    offsets = -10:0.1:10
    offsets = -15:0.1:15
    n_good_diffs = zeros(Int, length(offsets))
    total_diffs = zeros(Float64, length(offsets))

    # Calculate metrics for each offset
    for (i, offset) in enumerate(offsets)
        n_good_diffs[i], total_diffs[i] = calc_total_diff(offset)
    end

    # Find offsets with maximum number of good differences
    # zero_off_good = n_good_diffs[Int(ceil(length(offsets)/2))]
    zero_off_good = maximum(n_good_diffs)
    max_good_indices = findall(x -> x == zero_off_good, n_good_diffs)

    # Among those, find the one with minimum total difference
    best_idx = max_good_indices[argmin(total_diffs[max_good_indices])]
    best_offset = offsets[best_idx]
#    println(best_offset,chip,tele)

    # println("Best offset: $best_offset Å with $(max_good) good differences")

    seg_pix = mean.(segments)
    pred_wavs = seg_rough_wavs .+ best_offset
    nearest_idx, d2th_wavs, msk2use = assign_segments(pred_wavs, bright_wavs)

    # ---- self-calibrating association -------------------------------------
    # The rough model (quartic or legacy linear) is only a SEED: after the
    # initial nearest-line assignment, a low-order correction to the rough
    # model is fit to THIS fiber's own (segment pixel, assigned wavelength)
    # pairs -- with leave-one-out rejection of inconsistent pairs (thresholds
    # relative to the pairs' own scatter, no absolute wavelength numbers) --
    # and all segments are re-associated under the corrected, data-driven
    # model. Iterating to a fixed point makes the association depend on the
    # night's own data, not on the hardcoded seed: a wrong seed (curvature
    # error, instrument re-registration) is corrected here instead of
    # silently mislabeling lines. Convergence = the segment->line assignment
    # repeats; failure to converge within `max_assoc_iter` is flagged loudly
    # (assoc_converged = false in the skyLinePeaks product).
    assoc_niter = 0
    assoc_converged = true
    assoc_ndrop = 0
    if self_calibrate
        xnorm = (seg_pix .- (N_XPIX ÷ 2)) ./ N_XPIX
        # scatter floor: segment centroids are means of integer pixel indices
        # (quantized at the ~0.3 px level), expressed in wavelength through the
        # model's own local dispersion -- a resolution scale, not a position
        # lookup
        floor_ang = 0.3 * abs(rough_dispersion(N_XPIX ÷ 2, rw))
        prev_assign = nearest_idx .* msk2use
        assoc_converged = false
        for iter in 1:max_assoc_iter
            assoc_niter = iter
            use = findall(msk2use)
            if length(use) < 5
                # too few pairs to self-calibrate; keep the seed association
                assoc_converged = true
                break
            end
            # smooth wavelength-space correction to the seed model, robust to
            # misassigned pairs via sequential LOO rejection
            dlam = bright_wavs[nearest_idx[use]] .- seg_rough_wavs[use]
            keep, cq, _ = loo_poly_reject(xnorm[use], dlam; porder = 2,
                nsigma = assoc_loo_nsigma, max_drop_frac = 1 / 3,
                scatter_floor = floor_ang)
            assoc_ndrop = length(use) - count(keep)
            pred_wavs = seg_rough_wavs .+ positional_poly_mat(xnorm, porder = 2) * cq
            new_nearest, new_d2, new_msk = assign_segments(pred_wavs, bright_wavs)
            new_assign = new_nearest .* new_msk
            converged_now = (new_assign == prev_assign)
            nearest_idx, d2th_wavs, msk2use = new_nearest, new_d2, new_msk
            if converged_now
                assoc_converged = true
                break
            end
            prev_assign = new_assign
        end
        # final consistency cut: pairs inconsistent with the converged
        # data-driven model are removed from the emitted association entirely.
        # (A stolen label whose own feature is absent re-associates to the
        # thief peak whenever the thief sits inside the acceptance radius --
        # e.g. the 5.7 A 15546->15540 mode; it fails this cut instead, so the
        # label goes unmeasured rather than silently wrong.)
        use = findall(msk2use)
        if length(use) >= 5
            dlam = bright_wavs[nearest_idx[use]] .- seg_rough_wavs[use]
            keep, _, _ = loo_poly_reject(xnorm[use], dlam; porder = 2,
                nsigma = assoc_loo_nsigma, max_drop_frac = 1 / 3,
                scatter_floor = floor_ang)
            assoc_ndrop = length(use) - count(keep)
            msk2use[use[.!keep]] .= false
        end
    end
    assoc_qa = (assoc_niter, assoc_converged, assoc_ndrop)

    th_norm_flux = [bright_lines.model_rel_int[nearest_idx[i]]
                    for i in eachindex(nearest_idx)]
    if count(msk2use) < 1
        return nothing, nothing, nothing, assoc_qa
    end

    if assoc_only
        # association-stage result without the (expensive) profile fits:
        # rows are [segment centroid pixel; assigned catalog centroid; linindx]
        sel = findall(msk2use)
        amat = zeros(Float64, 3, length(sel))
        for (k, i) in enumerate(sel)
            ref = nearest_idx[i]
            amat[:, k] .= [seg_pix[i], bright_lines.wave_cen_ang[ref],
                bright_lines.linindx[ref]]
        end
        return amat, best_offset, length(sel), assoc_qa
    end
    max_obs_flux = maximum(segment_fluxes[msk2use])
    obs_norm_flux = segment_fluxes ./ max_obs_flux
    th_norm_flux = th_norm_flux ./ maximum(th_norm_flux[msk2use])
    # for (i,w) in enumerate(adjusted_wavs)
    #     if msk2use[i]
    #         println("Adjusted: $(round(w,digits=1)) Å -> Nearest line: $(round(bright_wavs[nearest_idx[i]],digits=1)) Å (th_norm_flux: $(round(th_norm_flux[i],digits=3)), obs_norm_flux: $(round(obs_norm_flux[i],digits=3)))")
    #     end
    # end

    # every selected (class A/B) row is a clean fit-eligible Lambda-doublet by
    # construction of the generated list, so no per-row cleanliness mask is needed
    # msk2use .&= (abs.(obs_norm_flux .- th_norm_flux) .<= 0.3)
    cen_pixs = seg_pix[msk2use]

    ######

    # Fixed-splitting, fixed-ratio Lambda-doublet fit ("V1"): the component
    # separation is a spectroscopic constant known to ~0.01 A (Brooke et al. 2016)
    # -- far better than any per-fiber fit -- and the e/f component ratio is
    # 1.000 +/- 0.002 and temperature-independent. Freeing the separation railed
    # at its bound in 50-91% of unresolved-line fits and carried a degeneracy
    # jump mode of 150-235 mA under LSF-width mismatch; fixing both leaves a
    # 5-parameter fit (amplitude, centroid, width, background slope + offset).
    function get_subline_params(lindx)
        cpix = cen_pixs[lindx]
        ref_indx = nearest_idx[msk2use][lindx]

        # sublines are ascending in wavelength in the generated list
        subline_wavs = [bright_lines.subwave_1_ang[ref_indx],
            bright_lines.subwave_2_ang[ref_indx]]
        subLine_weight = [bright_lines.weight_1[ref_indx],
            bright_lines.weight_2[ref_indx]]
        subLine_weight = subLine_weight ./ sum(subLine_weight)

        outwave = subline_wavs' * subLine_weight

        # fixed splitting in pixels, converted with the LOCAL dispersion of the
        # rough model at this peak (the global linear slope is up to ~5% off
        # locally, which matters for the resolved doublets)
        sepd = (subline_wavs[2] - subline_wavs[1]) / abs(rough_dispersion(cpix, rw))

        # dispersion is NEGATIVE (wavelength falls with pixel), so the BLUE
        # component (subline_wavs[1]) sits at the HIGHER pixel: pair the weights
        # accordingly. (The previous free-separation fit paired weight 1 with the
        # lower pixel, leaving outpix and outwave mutually inconsistent by
        # (w1 - w2) * sep -- harmless at 1:1 weights, first-order otherwise.)
        w_lo_pix = subLine_weight[2]   # red component -> lower pixel
        w_hi_pix = subLine_weight[1]   # blue component -> higher pixel

        peak_range = collect(Int(floor(cpix - 10)):Int(floor(cpix + 10)))
        x = peak_range
        x_edges = [peak_range .- 0.5; peak_range[end] + 0.5]
        y = flux_vec[peak_range]
        yamp = max(abs(maximum(y)) / sqrt(2 * pi),0.001)

        # p = [amplitude, centroid pix, width pix, bkg slope, bkg offset]
        function gfit(p; exclude_idx = nothing)
            lam1 = p[2] - sepd / 2
            lam2 = p[2] + sepd / 2

            model = p[1] * w_lo_pix * diff(normal_cdf.(x_edges .- lam1, p[3])) .+
                    p[1] * w_hi_pix * diff(normal_cdf.(x_edges .- lam2, p[3])) .+
                    p[4] .* (x .- cpix) .+ p[5]
            residuals = (y .- model) .^ 2
            if !isnothing(exclude_idx)
                return sum(residuals[1:(exclude_idx - 1)]) + sum(residuals[(exclude_idx + 1):end])
            else
                return sum(residuals)
            end
        end

        function gfit2(p)
            lam1 = p[2] - sepd / 2
            lam2 = p[2] + sepd / 2
            model = p[1] * w_lo_pix * diff(normal_cdf.(x_edges .- lam1, p[3])) .+
                    p[1] * w_hi_pix * diff(normal_cdf.(x_edges .- lam2, p[3])) .+
                    p[4] .* (x .- cpix) .+ p[5]
            return model
        end

        yoffset = abs(nanzeromedian(y))
        p0 = [yamp, cpix, 1, 0.0, yoffset]
        lb = [0.0, cpix - 10, 0.3, -5, 0]
        ub = [20 * yamp, cpix + 10, 3.0, 5, 2 * yoffset]

        opt_prob = Optim.optimize(gfit, lb, ub, p0, Fminbox(NelderMead()))
        fitparams = Optim.minimizer(opt_prob)

        model_pred = gfit2(fitparams)
        mod_res = y .- model_pred
        exclude_idx = argmax(abs.(mod_res))

        gfit_excl = (p) -> gfit(p; exclude_idx = exclude_idx)
        opt_prob = Optim.optimize(gfit_excl, lb, ub, p0, Fminbox(NelderMead()))
        fitparams = Optim.minimizer(opt_prob)

        # weight-averaged component position with the corrected pairing
        # (equals fitparams[2] for 1:1 weights)
        outpix = fitparams[2] + (w_hi_pix - w_lo_pix) * sepd / 2
        # keep the saved 9-row layout: the (now fixed) separation is written into
        # the old free-separation slot so downstream consumers are untouched
        return [outpix, outwave, fitparams[1], fitparams[2], fitparams[3], sepd,
            fitparams[4], fitparams[5], bright_lines.linindx[ref_indx]]
    end

    # make sure to save the pamt and the tparams (BOTH!)
    tout = get_subline_params.(1:length(cen_pixs))
    pmat = hcat(tout...)
    th_norm_flux = map(x -> x[3], tout)
    th_norm_flux ./= maximum(th_norm_flux)
    mskFlux = ones(Bool, length(cen_pixs)) #.& (abs.(obs_norm_flux[msk2use] .- th_norm_flux) .<= 0.15) # chip b needs this, but chip a/c hates it
    return pmat[:, mskFlux], best_offset, count(mskFlux), assoc_qa
end

function get_and_save_sky_peaks(fname, roughwave_dict, df_sky_lines; checkpoint_mode = "commit_same")
    outname = replace(replace(fname, "ar1Dcal" => "skyLinePeaks"), "ar1D" => "skyLinePeaks")
    if check_file(outname, mode = checkpoint_mode)
        return
    end

    sname = split(split(fname, "/")[end], "_")
    fnameType, tele, mjd, expnum, chip, image_type = sname[(end - 5):end]

    f = jldopen(fname, "r+")
    flux_1d = f["flux_1d"]
    extract_trace_centers = f["extract_trace_centers"]
    close(f)
    function get_sky_peaks_partial(flux_1d)
        get_sky_peaks(flux_1d, tele, chip, roughwave_dict, df_sky_lines)
    end
    pout = map(get_sky_peaks_partial, eachcol(flux_1d))
    boff = map(x -> replace_data_1(x[2]), pout)
    # self-calibrating association QA (see get_sky_peaks): iterations used,
    # convergence, and number of LOO-rejected pairs, per fiber
    assoc_niter = map(x -> x[4][1], pout)
    assoc_converged = map(x -> x[4][2], pout)
    assoc_ndrop = map(x -> x[4][3], pout)
    n_unconverged = count(.!assoc_converged)
    if n_unconverged > 0
        println("ASSOC NOT CONVERGED for $(n_unconverged) fibers in $(fname): " *
                "self-calibrating sky-line association hit max iterations")
    end

    unique_skyline_inds = sort(unique(vcat(map(x -> get_last_ind(x[1]), pout)...)))
    # wavecal class per line (aligned with the sky_line_mat rows), recorded in the
    # product so wavelength-solution ingest can mask class-A guards (code 1)
    # without re-reading the line list
    line_class = Int[get(WAVECAL_CLASS_CODES,
                         String(coalesce(df_sky_lines.wavecal_class[i], "")), 0)
                     for i in unique_skyline_inds]
    # println("$(tele) $(chip): $(length(unique_skyline_inds)) unique sky lines")
    # println(unique_skyline_inds)
    first_valid_idx = findfirst(x -> !isnothing(x[1]), pout)
    if isnothing(first_valid_idx)
        println("No valid sky lines found for ANY fiber in$(fname)")
        return
    end
    sky_line_mat = zeros(
        Float64, length(unique_skyline_inds), size(pout[first_valid_idx][1], 1), N_FIBERS)
    sky_trace_centers = zeros(
        Float64, length(unique_skyline_inds), N_FIBERS)
    fill!(sky_line_mat, NaN)
    x_pixels = collect(1:N_XPIX)
    for i in 1:N_FIBERS
        for (eindx, skyindx) in enumerate(unique_skyline_inds)
            if !isnothing(pout[i][1])
                locIndx = findfirst(pout[i][1][end, :] .== skyindx)
                if !isnothing(locIndx)
                    sky_line_mat[eindx, :, i] .= pout[i][1][:, locIndx]
                end
            end
        end
	sky_trace_centers[:, i] .= linear_interpolation(x_pixels, extract_trace_centers[:, i], extrapolation_bc = Line()).(sky_line_mat[:,1,i])
    end

    medx_detect = running_median.(eachrow(sky_line_mat[:, 1, :]), 31, :asym_trunc, nan = :ignore)
    medx_detect_mat = hcat(medx_detect...)
    sigma_detect = nanzeroiqr(sky_line_mat[:, 1, :] .- medx_detect_mat', 2)

    sky_line_mat_clean = copy(sky_line_mat)
    for i in 1:size(sky_line_mat, 1)
        msk = abs.(sky_line_mat[i, 1, :] .- medx_detect[i]) .<= 3 * sigma_detect[i]
        sky_line_mat_clean[i, :, .!msk] .= NaN
    end
    
    # f = h5open(outname, "w")
    # # Write cleaned data
    # write(f, "sky_line_mat_clean", sky_line_mat_clean)
    # attrs(f["sky_line_mat_clean"])["axis_1"] = "skylineID"
    # attrs(f["sky_line_mat_clean"])["axis_2"] = "fit_info"
    # attrs(f["sky_line_mat_clean"])["axis_3"] = "fibers"

    # # Write original data
    # write(f, "sky_line_mat", sky_line_mat)
    # attrs(f["sky_line_mat"])["axis_1"] = "skylineID"
    # attrs(f["sky_line_mat"])["axis_2"] = "fit_info"
    # attrs(f["sky_line_mat"])["axis_3"] = "fibers"

    # write(f, "sky_line_trace_centers", sky_trace_centers)
    # attrs(f["sky_line_trace_centers"])["axis_1"] = "skylineID"
    # attrs(f["sky_line_trace_centers"])["axis_2"] = "fibers"

    # # Write boff data
    # write(f, "boff", boff)
    # attrs(f["boff"])["axis_1"] = "fibers"
    # close(f)

    safe_jldsave(outname, sky_line_mat_clean = sky_line_mat_clean, sky_line_mat = sky_line_mat, sky_line_trace_centers = sky_trace_centers, boff = boff, line_class = line_class,
        assoc_niter = collect(assoc_niter), assoc_converged = collect(assoc_converged),
        assoc_ndrop = collect(assoc_ndrop), no_metadata = true)

    return
end

function rough_linear_wave(pix; a = 16156.8, b = -0.282)
    return a + (pix - (N_XPIX ÷ 2)) * b
end

function replace_data(x)
    if isnothing(x)
        return [NaN, NaN, NaN]
    else
        return x
    end
end

function replace_data_1(x)
    if isnothing(x)
        return NaN
    else
        return x
    end
end
