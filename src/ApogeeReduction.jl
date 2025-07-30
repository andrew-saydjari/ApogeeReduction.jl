module ApogeeReduction

function __init__()
    # ENV["SLACK_CHANNEL"] = "C08B7FKMP16" #apogee-reduction-jl
    if !haskey(ENV, "SLACK_CHANNEL")
        ENV["SLACK_CHANNEL"] = "C07KQ7BJY5P" #apogee-reduction-jl-dev
    end
end

export N_FIBERS, N_XPIX, N_CHIPS, CHIP_LIST, FIRST_CHIP, LAST_CHIP
const N_FIBERS = 300
const N_XPIX = 2048
const N_CHIPS = 3
const CHIP_LIST = ["R", "G", "B"]
const FIRST_CHIP = CHIP_LIST[1]
const LAST_CHIP = CHIP_LIST[end]

include("ar3D.jl")
include("ar2Dcal.jl")
include("ar1D.jl")
include("fileNameHandling.jl")
include("error_bits.jl")
include("utils.jl")

include("wavecal.jl")
include("skyline_peaks.jl")
include("arclamp_peaks.jl")
include("spectraInterpolation.jl")
include("traceExtract_GH.jl")

#export get_fibTargDict, fiberID2fiberIndx
end
