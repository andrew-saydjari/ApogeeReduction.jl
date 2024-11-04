
# hold off on prop ivar through until we switch to sutr_wood, also could implement a chi2 cut here
# add a condition that we should drop any x pixel where a bad bit in any of the pixels being summed is bad
function extract_boxcar!(extract_out, dimage_in, trace_params; widy = 2)
    fill!(extract_out, 0)
    for xpix in 1:2048
        for fib in 1:300
            flux, ypixf, sig = trace_params[xpix, fib, :]
            ypix = round(Int, ypixf)
            extract_out[xpix, fib] = sum(dimage_in[xpix, (ypix - widy):(ypix + widy)])
        end
    end
end

function extract_boxcar_bitmask!(extract_out, dimage_in, trace_params; widy = 2)
    fill!(extract_out, 0)
    for xpix in 1:2048
        for fib in 1:300
            flux, ypixf, sig = trace_params[xpix, fib, :]
            ypix = round(Int, ypixf)
            extract_out[xpix, fib] = reduce(|,dimage_in[xpix, (ypix - widy):(ypix + widy)])
        end
    end
end