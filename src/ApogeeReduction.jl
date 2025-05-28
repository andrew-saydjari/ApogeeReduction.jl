module ApogeeReduction

export N_FIBERS, N_XPIX, N_CHIPS
const N_FIBERS = 300
const N_XPIX = 2048
const N_CHIPS = 3
const CHIP_LIST = ["a","b","c"]
const FIRST_CHIP = "a"

include("utils.jl")
include("ar3D.jl")
include("ar1D.jl")
end
