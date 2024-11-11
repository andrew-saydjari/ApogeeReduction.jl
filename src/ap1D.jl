
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
            extract_out[xpix, fib] = reduce(|, dimage_in[xpix, (ypix - widy):(ypix + widy)])
        end
    end
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