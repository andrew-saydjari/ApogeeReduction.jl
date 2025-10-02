
chipRaw2Redux = Dict(
    "a" => "R",
    "b" => "G",
    "c" => "B",
)

chipRedux2Raw = Dict(
    "R" => "a",
    "G" => "b",
    "B" => "c",
)

function build_raw_path(tele, chip, mjd, exposure_id; cluster = "sdss", suppress_warning = false)
    base = if cluster == "sdss"
        "/uufs/chpc.utah.edu/common/home/sdss/sdsswork/data/apogee" #the raw data is NOT version dependent
    elseif cluster == "cca"
        "/mnt/ceph/users/sdssv/raw/APOGEE"
    else
        !suppress_warning && warn("Unknown cluster, interpreting as a path: $cluster")
        cluster
    end
    fname = if tele == "apo"
        "apR-$(chipRedux2Raw[chip])-$exposure_id.apz"
    elseif tele == "lco"
        "asR-$(chipRedux2Raw[chip])-$exposure_id.apz"
    else
        error("Unknown telescope")
    end
    return join([base, tele, mjd, fname], "/")
end

function dfindx_fname_format(dfindx)
    return lpad(dfindx, 4, "0")
end

function get_cal_file(parent_dir, tele, mjd, dfindx, chip, image_type; use_cal = false)
    if use_cal
        fname_type = "ar2Dcal"
    else
        fname_type = "ar2D"
    end
    return parent_dir *
           "apred/$(mjd)/$(fname_type)_$(tele)_$(mjd)_$(dfindx_fname_format(dfindx))_$(chip)_$(lowercase(image_type)).h5"
end

function get_fluxing_file_name(parent_dir, tele, mjd, dfindx, chip, cartid)
    return parent_dir *
           "dome_flats/$(mjd)/domeFlux_$(tele)_$(mjd)_$(dfindx_fname_format(dfindx))_$(chip)_domeflat_$(cartid).h5"
end

function get_1d_name(dfindx, df; cal = false)
    fnameType = if cal
        "ar1Dcal"
    else
        "ar1D"
    end
    return join([fnameType, df.observatory[dfindx], df.mjd[dfindx], dfindx_fname_format(df.exposure[dfindx]), FIRST_CHIP, df.image_type[dfindx]], "_")
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
