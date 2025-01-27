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

function reinterp_spectra(fname ;wavecal_type = "wavecal_skyline")
    # might need to add in telluric div functionality here?

    sname = split(split(fname, "/")[end], "_")
    fnameType, tele, mjd, chip, expid = sname[(end - 5):(end - 1)]

    # could shift this to a preallocation step
    outflux = zeros(length(logUniWaveAPOGEE),300)
    outvar = zeros(length(logUniWaveAPOGEE),300)
    outmsk = zeros(Int,length(logUniWaveAPOGEE),300);
    cntvec = zeros(Int,length(logUniWaveAPOGEE),300);

    pixvec = 1:(3*2048)
    flux_stack = zeros(3*2048,300)
    ivar_stack = zeros(3*2048,300)
    mask_stack = zeros(Int,3*2048,300)
    wave_stack = zeros(3*2048,300)
    chipBit_stack = zeros(Int,3*2048,300);

    ingestBit = zeros(Int,300)

    # add a for loop over the exposures (stop thinking about "visits" for now)
    # probably just generate ap1D file names from the alamanc files

    # this was used for looping over exposures in the visit
    # outdir = "/uufs/chpc.utah.edu/common/home/u6039752/scratch1/working/2024_12_05/outdir/"
    # fname = outdir * "apred/$(mjd)/" * get_1d_name(parse(Int, last(expid,4)), df) * ".jld2"

    f = jldopen(replace(replace(fname, fnameType => wavecal_type), "_a_" => "_"))
    chipWaveSoln = f["chipWaveSoln"]
    close(f)

    for (chipind,chip) in enumerate(["a","b","c"])
        fnameloc = replace(fname,"_a_"=>"_$(chip)_")
        f = jldopen(fnameloc)
        flux_1d = f["flux_1d"]
        ivar_1d = f["ivar_1d"]
        mask_1d = f["mask_1d"]
        close(f)

        flux_stack[(1:2048).+(3-chipind)*2048,:] .= flux_1d[end:-1:1,:]
        ivar_stack[(1:2048).+(3-chipind)*2048,:] .= ivar_1d[end:-1:1,:]
        mask_stack[(1:2048).+(3-chipind)*2048,:] .= mask_1d[end:-1:1,:]
        wave_stack[(1:2048).+(3-chipind)*2048,:] .= chipWaveSoln[end:-1:1,:,chipind]
        chipBit_stack[(1:2048).+(3-chipind)*2048,:] .+= 2^(chipind)
    end

    noBadBits = (mask_stack .& bad_pix_bits .== 0)
    chipBit_stack[.!(noBadBits)] .+= 2^4 # call pixels with bad bits thrown bad
    chipBit_stack[chipBit_stack.==0] .+= 2^4 # call missing chips bad

    # think about adding warnings for the last two cases
    good_pix = ((noBadBits) .& (chipBit_stack .& 2^4 .== 0)) .& (.!isnan.(ivar_stack)) .& (ivar_stack .> (10^-20))

    ## need to propagate the bit mask
    for fiberindx in 1:300
        good_pix_fiber = good_pix[:,fiberindx]
        flux_fiber = flux_stack[good_pix_fiber,fiberindx]
        ivar_fiber = ivar_stack[good_pix_fiber,fiberindx]
        wave_fiber = wave_stack[good_pix_fiber,fiberindx]
        chipBit_fiber = chipBit_stack[good_pix_fiber,fiberindx]
        pixindx_fiber = pixvec[good_pix_fiber]

        Rinv = generateInterpMatrix_sparse_inv(wave_fiber,chipBit_fiber,logUniWaveAPOGEE,pixindx_fiber)
        normvec = dropdims(sum(Rinv,dims=2),dims=2)
        msk_inter = (normvec.!=0)

        outflux[msk_inter,fiberindx] .+= (Rinv*flux_fiber)[msk_inter]
        outvar[msk_inter,fiberindx] .+= ((Rinv.^2)*(1 ./ ivar_fiber))[msk_inter]
        cntvec[:,fiberindx] .+= msk_inter

        if all(isnanorzero.(flux_fiber)) && ((ingestBit[fiberindx] & 2^1)==0)
            ingestBit[fiberindx] += 2^1 # ap1D exposure flux are all NaNs (for at least one of the exposures)
        elseif all(.!((chipBit_fiber .& 2^4).==0)) && ((ingestBit[fiberindx] & 2^2)==0)
            ingestBit[fiberindx] += 2^2 # ap1D exposure pixels are all bad by bit mask  (for at least one of the exposures)
        elseif all(isnanorzero.(ivar_fiber)) && ((ingestBit[fiberindx] & 2^3)==0)
            ingestBit[fiberindx] += 2^3 # ap1D exposure ivars are all NaNs or zeros (for at least one of the exposures)
        end
    end

    framecnts = maximum(cntvec,dims=1) #     framecnts = maximum(cntvec) # a little shocked that I throw it away if it is bad in even one frame
    outflux./=framecnts
    outvar./=(framecnts.^2)
    # need to update this to a bit mask that is all or any for the pixels contributing to the reinterpolation
    outmsk = (cntvec .== framecnts)

    outname = replace(replace(replace(fname, "ap1D" => "ar1Duni"),"_a_" => "_"),".jld2" => ".h5")
    f = h5open(outname, "w")
    # Write reinterpolated data
    write(f, "flux_1d", outflux)
    write(f, "ivar_1d", 1 ./outvar)
    write(f, "mask_1d", outmsk)
    close(f)
    return
end

logUniWaveAPOGEE = 10 .^range((start=4.17825),step=6.0e-6,length=8700);