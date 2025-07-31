
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
        "/mnt/ceph/users/asaydjari/apogee/raw"
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

# both inputs are strings, as is the output
function short_expid_to_long(mjd, expnum)
    return lpad(string(parse(Int,mjd) - 55562), 4, "0") * lpad(expnum, 4, "0")
end

# both inputs are integers, as is the output
function long_expid_to_short(mjd, expnum)
    return expnum - (mjd - 55562) * 10000
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

function get_fluxing_file_name(parent_dir, tele, mjd, exposure, chip, cartid)
    return parent_dir *
           "dome_flats/$(mjd)/domeFlux_$(tele)_$(mjd)_$(exposure)_$(chip)_DOMEFLAT_$(cartid).h5"
end

function get_1d_name(expid, df; cal = false)
    fnameType = if cal
        "ar1Dcal"
    else
        "ar1D"
    end
    return join([fnameType, df.observatory[expid], df.mjd[expid], last(df.exposure_str[expid],4), chipRaw2Redux[df.chip[expid]], df.exptype[expid]], "_")
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
