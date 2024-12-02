module ApogeeReduction

const N_FIBERS = 300
const N_XPIX = 2048

ENV["SLACK_CHANNEL"] = "C07KQ7BJY5P"

include("fileNameHandling.jl") # where files are located on utah
include("utils.jl")            # general utilities
include("ap3D.jl")             # 3D -> 2D
include("ap2Dcal.jl")          # calibration file stuff (???)
include("ap1D.jl")             # 2D -> 1D
include("plotutils.jl")        # setup for automatic plots on slack
end
