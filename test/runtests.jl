using ApogeeReduction, Test, Random, Statistics
using FITSIO: FITS

# The file structure in test roughly mirrors that of src.  Each file is included below.
@testset verbose=true "ApogeeReduction.jl" begin
    include("ar3D.jl")
    include("exposureClassifier.jl")
    include("ar1D.jl")
    include("fiber_types.jl")
    include("wavecal.jl")
    include("chip_transform.jl")
    include("chip_transform_bootstrap.jl")
    include("trace_template_cal.jl")
    include("fpi_gate.jl")
    include("fpi_homogeneity.jl")
    include("trace_qa.jl")
    include("safe_jldsave.jl")
    include("lsf.jl")
    include("readonly_opens.jl")
    include("scripts_parse.jl")
end
