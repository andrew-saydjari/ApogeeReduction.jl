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

function build_raw_path(release_dir,obs,mjd,chip,expid)
    base = "/uufs/chpc.utah.edu/common/home/sdss/$(release_dir)/data/apogee"
    fname = "apR-$chip-$expid.apz"
    return join([base,obs,mjd,fname],"/")
end