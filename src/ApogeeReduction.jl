module ApogeeReduction

export N_FIBERS, N_XPIX, N_CHIPS, CHIP_LST
const N_FIBERS = 300
const N_XPIX = 2048
const N_CHIPS = 3
const CHIP_LST = ["R", "G", "B"]

include("utils.jl")
include("ar3D.jl")
include("ar1D.jl")
end
