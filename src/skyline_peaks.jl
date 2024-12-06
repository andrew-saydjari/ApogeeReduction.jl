import FastRunningMedian: running_median
using Optim

function get_sky_peaks(flux_vec, tele, chip, roughwave_dict, df_sky_lines)
    # Find indices where flux is above 99th percentile
    thresh = percentile(flux_vec, 97.5)
    above_thresh = findall(x -> x > thresh, flux_vec)

    # Group indices into segments, combining those less than 10 pixels apart
    segments = []
    current_segment = [above_thresh[1]]
    for i in 2:length(above_thresh)
        if above_thresh[i] - above_thresh[i - 1] <= 10
            push!(current_segment, above_thresh[i])
        else
            push!(segments, current_segment)
            current_segment = [above_thresh[i]]
        end
    end
    push!(segments, current_segment)
    mean_x = mean.(segments)
    length_segs = length.(segments)
    segments = segments[(mean_x .> 64) .& (mean_x .< 1984) .& (length_segs .> 2)]

    # Preallocate array for segment fluxes
    segment_fluxes = zeros(length(segments))

    # For each segment, compute flux in padded range after subtracting median
    for (i, segment) in enumerate(segments)
        # Get range with 10 pixel padding on each side
        start_idx = maximum([1, minimum(segment) - 10])
        end_idx = minimum([length(flux_vec), maximum(segment) + 10])

        # Get flux in range
        flux_range = flux_vec[start_idx:end_idx]

        # Subtract median of range
        med = median(flux_range)
        flux_range .-= med

        # Sum total flux in range
        segment_fluxes[i] = sum(flux_range)
    end

    seg_rough_wavs = map(
        x -> rough_linear_wave(
            mean(x); a = roughwave_dict[tele][chip][1], b = roughwave_dict[tele][chip][2]),
        segments)

    ######

    wav_range = (roughwave_dict[tele][chip][4], roughwave_dict[tele][chip][3])
    # Filter lines within wavelength range and sort by emission strength
    bright_lines0 = sort(
        filter(
            row -> (wav_range[1] <= row.wav * 10 <= wav_range[2]) .& (row.number_of_lines .== 2),
            df_sky_lines),
        :intensity,
        rev = true
    )[1:8, :]
    rescale_int = maximum(bright_lines0.intensity)
    bright_lines0.norm_intensity .= (bright_lines0.intensity ./ rescale_int)
    bright_lines0.norm_sigma_intensity .= (bright_lines0.sigma_intensity ./ rescale_int)
    bright_lines = bright_lines0[bright_lines0.norm_intensity .> 0.12, :] #prev 12

    # Convert bright_lines wavelengths to Angstroms and get as vector
    bright_wavs = bright_lines.wav .* 10

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

    # println("Best offset: $best_offset Å with $(max_good) good differences")

    adjusted_wavs = seg_rough_wavs .+ best_offset
    d2th_wavs = zeros(Float64, length(adjusted_wavs))
    nearest_idx = zeros(Int, length(adjusted_wavs))
    th_norm_flux = zeros(Float64, length(adjusted_wavs))
    for (i, w) in enumerate(adjusted_wavs)
        nearest_idx[i] = argmin(abs.(w .- bright_wavs))
        d2th_wavs[i] = abs(w - bright_wavs[nearest_idx[i]])
        th_norm_flux[i] = bright_lines.intensity[nearest_idx[i]]
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

    # msk2use .&= (bright_lines.clean[nearest_idx] .== true) .& (abs.(obs_norm_flux .- th_norm_flux) .<= 0.3)
    msk2use .&= (bright_lines.clean[nearest_idx] .== true)
    cen_pixs = mean.(segments[msk2use])

    ######

    function get_subline_params(lindx)
        cpix = cen_pixs[lindx]
        ref_indx = nearest_idx[msk2use][lindx]

        subline_wav_diffs = (parse.(Float64,
            split(
                replace(replace(bright_lines[ref_indx, :subline_wav][2:(end - 1)], "  " => ","),
                    " " => ","),
                ",")[1:2]) .- bright_lines[ref_indx, :wav]) .* 10
        subline_wav_pix_diffs = subline_wav_diffs ./ (roughwave_dict[tele][chip][2])
        subLine_weight = parse.(Float64,
            split(replace(bright_lines[ref_indx, :subline_I][2:(end - 1)], " " => ","), ","))
        subLine_ratio = subLine_weight[2] / subLine_weight[1]

        peak_range = collect(Int(floor(cpix - 10)):Int(floor(cpix + 10)))
        x = peak_range
        y = flux_vec[peak_range]
        yamp = abs(maximum(y)) / sqrt(2 * pi)

        function gfit(p; exclude_idx = nothing)
            lam1 = p[2] - p[4] / 2
            lam2 = p[2] + p[4] / 2
            model = p[1] * subLine_weight[1] * normal_pdf.(x .- lam1, p[3]) .+
                    p[1] * subLine_weight[2] * normal_pdf.(x .- lam2, p[3]) .+
                    p[5] .* (x .- cpix) .+ p[6]
            residuals = (y .- model) .^ 2
            if !isnothing(exclude_idx)
                return sum(residuals[1:(exclude_idx - 1)]) + sum(residuals[(exclude_idx + 1):end])
            else
                return sum(residuals)
            end
        end

        function gfit2(p)
            lam1 = p[2] - p[4] / 2
            lam2 = p[2] + p[4] / 2
            model = p[1] * subLine_weight[1] * normal_pdf.(x .- lam1, p[3]) .+
                    p[1] * subLine_weight[2] * normal_pdf.(x .- lam2, p[3]) .+
                    p[5] .* (x .- cpix) .+ p[6]
            return model
        end

        sepd = subline_wav_pix_diffs[2] - subline_wav_pix_diffs[1]
        yoffset = abs(nanzeromedian(y))
        p0 = [yamp, cpix, 1, sepd, 0.0, yoffset]
        lb = [0.0, cpix - 10, 0.3, 0.8 * sepd, -5, 0]
        ub = [20 * yamp, cpix + 10, 3.0, 1.5 * sepd, 5, 2 * yoffset]

        opt_prob = Optim.optimize(gfit, lb, ub, p0, Fminbox(NelderMead()))
        fitparams = Optim.minimizer(opt_prob)

        model_pred = gfit2(fitparams)
        mod_res = y .- model_pred
        exclude_idx = argmax(abs.(mod_res))

        gfit_excl = (p) -> gfit(p; exclude_idx = exclude_idx)
        opt_prob = Optim.optimize(gfit_excl, lb, ub, p0, Fminbox(NelderMead()))
        fitparams = Optim.minimizer(opt_prob)

        outpix = [fitparams[2] - fitparams[5] / 2, fitparams[2] + fitparams[5] / 2]' *
                 subLine_weight
        outwave = bright_lines[ref_indx, :wav] * 10
        return [outpix, outwave, fitparams..., bright_lines.linindx[ref_indx]]
    end

    # make sure to save the pamt and the tparams (BOTH!)
    tout = get_subline_params.(1:length(cen_pixs))
    pmat = hcat(tout...)
    th_norm_flux = map(x -> x[3], tout)
    th_norm_flux ./= maximum(th_norm_flux)
    mskFlux = ones(Bool, length(cen_pixs)) #.& (abs.(obs_norm_flux[msk2use] .- th_norm_flux) .<= 0.15) # chip b needs this, but chip a/c hates it
    return pmat[:, mskFlux], best_offset, count(mskFlux)
end

function get_and_save_sky_peaks(fname, roughwave_dict, df_sky_lines)
    sname = split(fname, "_")
    tele, mjd, chip, expid = sname[(end - 4):(end - 1)]
    f = jldopen(fname, "r+")
    flux_1d = f["flux_1d"]
    close(f)
    function get_sky_peaks_partial(flux_1d)
        get_sky_peaks(flux_1d, tele, chip, roughwave_dict, df_sky_lines)
    end
    pout = map(get_sky_peaks_partial, eachcol(flux_1d))
    boff = map(x -> replace_data_1(x[2]), pout)

    unique_skyline_inds = sort(unique(vcat(map(x -> get_last_ind(x[1]), pout)...)))
    # println("$(tele) $(chip): $(length(unique_skyline_inds)) unique sky lines")
    # println(unique_skyline_inds)
    sky_line_mat = zeros(Float64, length(unique_skyline_inds), size(pout[1][1], 1), 300)
    fill!(sky_line_mat, NaN)
    for i in 1:300
        for (eindx, skyindx) in enumerate(unique_skyline_inds)
            if !isnothing(pout[i][1])
                locIndx = findfirst(pout[i][1][end, :] .== skyindx)
                if !isnothing(locIndx)
                    sky_line_mat[eindx, :, i] .= pout[i][1][:, locIndx]
                end
            end
        end
    end

    medx_detect = running_median.(eachrow(sky_line_mat[:, 1, :]), 31, :asym_trunc, nan = :ignore)
    medx_detect_mat = hcat(medx_detect...)
    sigma_detect = nanzeroiqr(sky_line_mat[:, 1, :] .- medx_detect_mat', 2)

    sky_line_mat_clean = copy(sky_line_mat)
    for i in 1:size(sky_line_mat, 1)
        msk = abs.(sky_line_mat[i, 1, :] .- medx_detect[i]) .<= 3 * sigma_detect[i]
        sky_line_mat_clean[i, :, .!msk] .= NaN
    end
    outname = replace(fname, "ap1D" => "skyLine_peaks")
    f = h5open(outname, "w")
    # Write cleaned data
    write(f, "sky_line_mat_clean", sky_line_mat_clean)
    attrs(f["sky_line_mat_clean"])["axis_1"] = "skylineID"
    attrs(f["sky_line_mat_clean"])["axis_2"] = "fit_info"
    attrs(f["sky_line_mat_clean"])["axis_3"] = "fibers"

    # Write original data
    write(f, "sky_line_mat", sky_line_mat)
    attrs(f["sky_line_mat"])["axis_1"] = "skylineID"
    attrs(f["sky_line_mat"])["axis_2"] = "fit_info"
    attrs(f["sky_line_mat"])["axis_3"] = "fibers"

    # Write boff data
    write(f, "boff", boff)
    attrs(f["boff"])["axis_1"] = "fibers"
    close(f)
end

function rough_linear_wave(pix; a = 16156.8, b = -0.282)
    return a + (pix - 1024) * b
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
