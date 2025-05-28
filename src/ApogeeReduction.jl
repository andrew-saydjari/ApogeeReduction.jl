module ApogeeReduction

export N_FIBERS, N_XPIX, N_CHIPS, CHIP_LIST, FIRST_CHIP
const N_FIBERS = 300
const N_XPIX = 2048
const N_CHIPS = 3
const CHIP_LIST = ["a","b","c"]
const FIRST_CHIP = CHIP_LIST[1]

include("utils.jl")
include("ar3D.jl")
include("ar1D.jl")
end
