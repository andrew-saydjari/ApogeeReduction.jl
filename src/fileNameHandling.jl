function getUtahBase(release_dir, redux_ver)
    if occursin("dr", release_dir)
        dr_number = parse(Int, match(r"dr(\d+)", release_dir).captures[1])
        if (10 <= dr_number <= 17)
            return "/uufs/chpc.utah.edu/common/home/sdss/$(release_dir)/apogee/spectro/redux/$(redux_ver)/"
        elseif (18 <= dr_number)
            return "/uufs/chpc.utah.edu/common/home/sdss/$(release_dir)/spectro/apogee/redux/$(redux_ver)/"
        end
    end
    if release_dir[1:3] == "ipl"
        return "/uufs/chpc.utah.edu/common/home/sdss/$(release_dir)/spectro/apogee/redux/$(redux_ver)/"
    end
    # the last case catches the dev versions under daily/trial versions
    return "/uufs/chpc.utah.edu/common/home/sdss/$(release_dir)/apogee/spectro/redux/$(redux_ver)/"
end

function build_raw_path(tele, chip, mjd, exposure_id)
    base = "/uufs/chpc.utah.edu/common/home/sdss/sdsswork/data/apogee" #the raw data is NOT version dependent
    fname = if tele == "apo"
        "apR-$chip-$exposure_id.apz"
    elseif tele == "lco"
        "asR-$chip-$exposure_id.apz"
    else
        error("Unknown telescope")
    end
    return join([base, tele, mjd, fname], "/")
end

function short_expid_to_long(mjd, expid)
    # this should probably be done with math
    return lpad(parse(Int, string(mjd - 55562), 4, "0") * lpad(expid, 4, "0"))
end

function long_expid_to_short(mjd, expid)
    return expid - (mjd - 55562) * 10000
end

function get_cal_file(parent_dir, tele, mjd, expnum, chip, exptype; use_cal = false)
    if use_cal
        fname_type = "ar2Dcal"
    else
        fname_type = "ar2D"
    end
    return parent_dir *
           "apred/$(mjd)/$(fname_type)_$(tele)_$(mjd)_$(lpad(expnum, 4, "0"))_$(chip)_$(exptype).h5"
end

function get_fluxing_file_name(parent_dir, tele, chip, mjd, exposure_id, cartid)
    return parent_dir *
           "dome_flats/$(mjd)/domeFlux_$(tele)_$(mjd)_$(chip)_$(lpad(exposure_id, 8, "0"))_DOMEFLAT_$(cartid).h5"
end

function get_1d_name(expid, df; cal = false)
    fnameType = if cal
        "ar1Dcal"
    else
        "ar1D"
    end
    return join([fnameType, df.observatory[expid], df.mjd[expid], last(df.exposure_str[expid],4), df.chip[expid], df.exptype[expid]], "_")
end

function fiberIndx2fiberID(fibindx)
    return 301 - fibindx
end

function fiberID2fiberIndx(fiberID)
    return 301 - fiberID
end

function adjFiberIndx(teleind, fibindx)
    return (teleind - 1) * 300 + fibindx
end

function adjFiberIndx2FiberIndx(adjfibindx)
    return mod1(adjfibindx, 300)
end

function get_fpi_guide_fiberID(tele)
    if (tele == "apo")
        return 75, 225
    elseif (tele == "lco")
        return 82, 213
    end
end

function get_fps_plate_divide(tele)
    if (tele == "apo")
        return 59423
    elseif (tele == "lco")
        return 59810 # still need to confirm this
    end
end
