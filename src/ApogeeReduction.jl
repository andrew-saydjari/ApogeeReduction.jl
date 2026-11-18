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
include("exposureClassifier.jl")
include("ar1D.jl")
include("fileNameHandling.jl")
include("utils.jl")

# acceptance gate for the FPI nightly wavelength solution
# (uses nanmedian from utils.jl and fiberID2fiberIndx from fileNameHandling.jl)
include("fpi_gate.jl")
# corpus-level FPI fiber-index homogeneity survey
# (uses get_fps_plate_divide from fileNameHandling.jl and
#  read_almanac_exp_df from utils.jl, included above)
include("fpi_homogeneity.jl")
# comb-measured inter-chip transform calibration (read by the sky wavelength fit)
include("chip_transform.jl")
# per-epoch trace templates: calibration file, validation, lookup by MJD, identity
# (read by make_traces_from_flats.jl and make_traceMain in trace_qa.jl)
include("trace_template_cal.jl")
include("wavecal.jl")
include("skyline_peaks.jl")
include("arclamp_peaks.jl")
include("spectraInterpolation.jl")
# per-flat trace QA, gate, night trace choice, whole-slit identification, trace bits
# (uses utils.jl helpers; called by traceExtract_GH.jl and ar1D.jl)
include("trace_qa.jl")
include("traceExtract_GH.jl")
# epoch trace templates: registration, template-seeded per-column refinement, fallback levels
# (uses trace_qa.jl and traceExtract_GH.jl)
include("slit_model.jl")
include("trace_levels.jl")
# lsf.jl depends on traceExtract_GH.jl (int_gauss_hermite_term) and
# fileNameHandling.jl (adjFiberIndx2FiberIndx), included above
include("lsf.jl")

#export get_fibTargDict, fiberID2fiberIndx
end
