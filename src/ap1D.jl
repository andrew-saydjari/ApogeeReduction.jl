using FastRunningMedian: running_median
using Distributions: cdf, Normal

# this file contains the code needed to extract a 1D spectrum from a 2D images.
# trace_params is of size (n_x_pix, n_fibers, 3)
# the elements correspodn to flux, y-pixel, and gaussian witdth (sigma)

# hold off on prop ivar through until we switch to sutr_wood, also could implement a chi2 cut here
# add a condition that we should drop any x pixel where a bad bit in any of the pixels being summed is bad

"""
Regularize the trace by applying a running median filter to each param in each fiber.
Could be denoised further by fitting a low-order polynomial or similar.
"""
function regularize_trace(trace_params; window_size = 101)
    @assert isodd(window_size) # otherwise the length of the regularized array is incorrect
    n_xpix, n_fibers = size(trace_params)[1:2]
    regularized_trace = similar(trace_params)
    for fiber in 1:n_fibers, param in 1:3
        regularized_trace[:, fiber, param] = running_median(
            trace_params[:, fiber, param], window_size, :asym_trunc; nan = :ignore)
    end
    regularized_trace
end

function extract_boxcar(dimage, ivarimage, pix_bitmask, trace_params; boxcar_halfwidth = 2)
    flux_1d = extract_boxcar_core(dimage, trace_params, boxcar_halfwidth)
    var_1d = extract_boxcar_core(1 ./ ivarimage, trace_params, boxcar_halfwidth)
    ivar_1d = 1.0 ./ var_1d
    mask_1d = extract_boxcar_bitmask(pix_bitmask, trace_params, boxcar_halfwidth)

    flux_1d, ivar_1d, mask_1d
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
        mask[xpix, fib] = reduce(
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

        flux_1d[xpix, fib] = sum(weights .* dimage[xpix, ypixels] .* ivarimage[xpix, ypixels]) /
                             sum(weights .^ 2 .* ivarimage[xpix, ypixels])
        ivar_1d[xpix, fib] = sum(weights .^ 2 .* ivarimage[xpix, ypixels])

        # bitmask
        mask_1d[xpix, fib] = reduce(|, pix_bitmask[xpix, ypixels])
    end
    flux_1d, ivar_1d, mask_1d
end

function get_fibTargDict(f, tele, mjd, expid)
    # worry about the read-in overhead per expid?
    df_exp = DataFrame(read(f["$(tele)/$(mjd)/exposures"]))
    mjdfps2plate = get_fps_plate_divide(tele)
    configName = (mjd > mjdfps2plate) ? "fps" : "plates"
    configIDname = (mjd > mjdfps2plate) ? "configid" : "plateid"
    configNumStr = df_exp[expid, Symbol(configIDname)]

    fibtargDict = if (df_exp.exptype[expid] == "OBJECT")
        try
            df_fib = DataFrame(read(f["$(tele)/$(mjd)/fibers/$(configName)/$(configNumStr)"]))
            Dict(parse.(Int, df_fib[!, "fiber_id"]) .=> df_fib[!, "target_type"])
        catch
            Dict(1:300 .=> "fiberTypeFail")
        end
    else
        Dict(1:300 .=> "cal")
    end

    if mjd > mjdfps2plate
        fpifib1, fpifib2 = get_fpi_guide_fiberID(tele)
        fibtargDict[fpifib1] = "fpiguide"
        fibtargDict[fpifib2] = "fpiguide"
    end
    return fibtargDict
end