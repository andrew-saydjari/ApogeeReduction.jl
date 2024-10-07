using DataFrames

function calStrip(fname)
    calType, telescope, chip, mjdstart, mjdend = split(split(split(fname,"/")[end],".")[1],"_")
    return vcat(abspath(fname),calType, telescope, chip, parse(Int,mjdstart), parse(Int,mjdend))
end

function cal2df(flist)
    df = DataFrame(path = String[], calType = String[], telescope = String[], chip = String[], mjdstart = Int[], mjdend = Int[])
    for row in calStrip.(flist)
        push!(df, row)
    end
    sort!(df, :mjdstart)
    return df
end

function get_cal_path(df,tele,mjd,chip)
    mskTeleChip = (df.telescope .== tele) .& (df.chip .== chip) 
    msk0 = mskTeleChip .& (df.mjdstart .<= mjd .<= df.mjdend)
    msk1 = mskTeleChip .& (mjd .<= df.mjdstart)
    msk2 = mskTeleChip .& (df.mjdend .<= mjd)
    if count(msk0) > 0
        return df.path[findlast(msk0)], 0
    elseif count(msk1) > 0
        return df.path[findlast(msk1)], 2^1 # flag for using cal from after exposure
    else count(msk2) > 0
        return df.path[findfirst(msk2)], 2^2 # flag for using cal from before exposure
    end
end