using FastRunningMedian

# this file contains the code needed to extract a 1D spectrum from a 2D images.
# trace_params is of size (n_x_pix, n_fibers, 3)
# the elements correspodn to flux, y-pixel, and gaussian witdth (sigma)

# hold off on prop ivar through until we switch to sutr_wood, also could implement a chi2 cut here
# add a condition that we should drop any x pixel where a bad bit in any of the pixels being summed is bad
function extract_boxcar!(extract_out, dimage_in, trace_params; boxcar_halfwidth = 2)
    fill!(extract_out, 0)
    for xpix in 1:2048
        for fib in 1:300
            _, ypixf, _ = trace_params[xpix, fib, :]
            ypix = round(Int, ypixf)
            extract_out[xpix, fib] = sum(dimage_in[
                xpix, (ypix - boxcar_halfwidth):(ypix + boxcar_halfwidth)])
        end
    end
end

function extract_boxcar_bitmask!(extract_out, dimage_in, trace_params; boxcar_halfwidth = 2)
    fill!(extract_out, 0)
    for xpix in 1:2048
        for fib in 1:300
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
"""
function extract_gaussian!(extract_out, dimage_in, trace_params; window_size = 19)
    for xpix in 1:2048, fib in 1:300
        _, ypixf, sigf = trace_params[xpix, fib, :]
        #TODO
    end
end
