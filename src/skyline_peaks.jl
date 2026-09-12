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

function get_sky_peaks(flux_vec, tele, chip, roughwave_dict, df_sky_lines;
				 med_flux_window = 31, flux_thresh = 97,
				 max_pix_sep = 5, n_pad = 5, min_seg_length = 2)

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
        return nothing, nothing, nothing
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

    adjusted_wavs = seg_rough_wavs .+ best_offset
    d2th_wavs = zeros(Float64, length(adjusted_wavs))
    nearest_idx = zeros(Int, length(adjusted_wavs))
    th_norm_flux = zeros(Float64, length(adjusted_wavs))
    for (i, w) in enumerate(adjusted_wavs)
        nearest_idx[i] = argmin(abs.(w .- bright_wavs))
        d2th_wavs[i] = abs(w - bright_wavs[nearest_idx[i]])
        th_norm_flux[i] = bright_lines.model_rel_int[nearest_idx[i]]
    end

    is_closest = zeros(Bool, length(d2th_wavs))
    for n in unique(nearest_idx)
        matches = findall(x -> x == n, nearest_idx)
        is_closest[matches] .= false
        closest_match = matches[argmin(d2th_wavs[matches])]
        is_closest[closest_match] = true
    end

    msk2use = is_closest .& (d2th_wavs .<= 10)
    if count(msk2use) < 1
        return nothing, nothing, nothing
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
    cen_pixs = mean.(segments[msk2use])

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
    return pmat[:, mskFlux], best_offset, count(mskFlux)
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

    safe_jldsave(outname, sky_line_mat_clean = sky_line_mat_clean, sky_line_mat = sky_line_mat, sky_line_trace_centers = sky_trace_centers, boff = boff, line_class = line_class, no_metadata = true)

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
