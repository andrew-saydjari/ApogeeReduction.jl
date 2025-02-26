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

function build_raw_path(obs, mjd, chip, expid)
    base = "/uufs/chpc.utah.edu/common/home/sdss/sdsswork/data/apogee" #the raw data is NOT version dependent
    fname = if obs == "apo"
        "apR-$chip-$expid.apz"
    elseif obs == "lco"
        "asR-$chip-$expid.apz"
    else
        error("Unknown telescope")
    end
    return join([base, obs, mjd, fname], "/")
end

function get_cal_file(parent_dir, tele, mjd, expid, chip, imtype; use_cal = false)
    expid_adj = string(mjd - 55562) * lpad(expid, 4, "0") # what a fun forced hard code!
    if use_cal
        fname_type = "ar2Dcal"
    else
        fname_type = "ar2D"
    end
    return parent_dir *
           "apred/$(mjd)/$(fname_type)_$(tele)_$(mjd)_$(chip)_$(expid_adj)_$(imtype).jld2"
end

function get_1d_name(expid, df)
    return join(
        ["ar1D", df.observatory[expid], df.mjd[expid],
            df.chip[expid], df.exposure[expid], df.exptype[expid]],
        "_")
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
