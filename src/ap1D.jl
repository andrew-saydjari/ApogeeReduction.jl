using FastRunningMedian: running_median
using Distrubtions: cdf, Normal

# this file contains the code needed to extract a 1D spectrum from a 2D images.
# trace_params is of size (n_x_pix, n_fibers, 3)
# the elements correspodn to flux, y-pixel, and gaussian witdth (sigma)

# hold off on prop ivar through until we switch to sutr_wood, also could implement a chi2 cut here
# add a condition that we should drop any x pixel where a bad bit in any of the pixels being summed is bad
function extract_boxcar!(extract_out, dimage_in, trace_params; boxcar_halfwidth = 2)
    fill!(extract_out, 0)
    for xpix in 1:N_XPIX
        for fib in 1:N_FIBERS
            _, ypixf, _ = trace_params[xpix, fib, :]
            ypix = round(Int, ypixf)
            extract_out[xpix, fib] = sum(dimage_in[
                xpix, (ypix - boxcar_halfwidth):(ypix + boxcar_halfwidth)])
        end
    end
end

function extract_boxcar_bitmask!(extract_out, dimage_in, trace_params; boxcar_halfwidth = 2)
    fill!(extract_out, 0)
    for xpix in 1:N_XPIX
        for fib in 1:N_FIBERS
            _, ypixf, _ = trace_params[xpix, fib, :]
            ypix = round(Int, ypixf)
            extract_out[xpix, fib] = reduce(
                |, dimage_in[xpix, (ypix - boxcar_halfwidth):(ypix + boxcar_halfwidth)])
        end
    end
end

"""
Regularize the trace by applying a running median filter to each param in each fiber.
Could be denoised further by fitting a low-order polynomial or similar.
"""
function regularize_trace(trace_params; window_size = 101)
    @assert isodd(window_size) # otherwise the length of the regularized array is incorrect

    regularized_trace = similar(trace_params)
    for fiber in 1:300, param in 1:3
        regularized_trace[:, fiber, param] = running_median(
            trace_params[:, fiber, param], window_size, :asym_trunc; nan = :ignore)
    end
    regularized_trace
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
    @assert size(dimage, 1) == n_xpix
    @assert size(ivarimage, 1) == n_xpix

    # return values to be filled
    flux_1d = Vector{Float64}(undef, n_xpix, n_fibers)
    var_1d = Vector{Float64}(undef, n_xpix, n_fibers)
    mask_1d = Vector{Int64}(undef, n_xpix, n_fibers)

    for xpix in 1:n_xpix, fib in 1:n_fibers
        _, y_peak, y_sigma = trace_params[xpix, fib, :]

        ypixels = floor(Int, y_peak - window_half_size):ceil(
            Int, y_peak + window_half_size)
        ypix_boundaries = [ypixels .- 0.5; ypixels[end] + 0.5]
        unnormalized_weights = diff(cdf.(Normal(y_peak, y_sigma), ypix_boundaries))
        weights = unnormalized_weights ./ sum(unnormalized_weights)

        flux_1d[xpix, fib] = sum(weights .* dimage[xpix, ypixels] .* ivarimage[xpix, ypixels]) /
                             sum(weights .^ 2 * ivarimage[xpix, ypixels])
        var_1d[xpix, fib] = 1 / sum(weights .^ 2 * ivarimage[xpix, ypixels])

        # bitmask
        mask_1d[xpix, fib] = reduce(|, pix_bitmask[xpix, ypixels])
    end
    flux_1d, var_1d, mask_1d
end
