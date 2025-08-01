using DataFrames

function calStrip(fname)
    sname = split(split(split(fname, "/")[end],".h5")[1], "_")
    calType, tele, chip, mjdstart, mjdend = sname[(end - 4):end]
    return vcat(abspath(fname), calType, tele, chip, parse(Int, mjdstart), parse(Int, mjdend))
end

function cal2df(flist)
    df = DataFrame(path = String[], calType = String[], telescope = String[],
        chip = String[], mjdstart = Int[], mjdend = Int[])
    for row in calStrip.(flist)
        push!(df, row)
    end
    sort!(df, :mjdstart)
    return df
end

function get_cal_path(df, tele, mjd, chip)
    mskTeleChip = (df.telescope .== tele) .& (df.chip .== chip)
    msk0 = mskTeleChip .& (df.mjdstart .<= mjd .<= df.mjdend)
    msk1 = mskTeleChip .& (mjd .<= df.mjdstart)
    
    if count(msk0) > 0
        return df.path[findlast(msk0)], 0
    elseif count(msk1) > 0
        return df.path[findlast(msk1)], 2^1 # flag for using cal from after exposure
    elseif any(mskTeleChip)
        return df.path[findfirst(mskTeleChip)], 2^2 # flag for using cal from before exposure
    else
        error("No calibration files found for telescope $tele, chip $chip. Please ensure calibration files exist in the calibration directory.")
    end
end
