using ApogeeReduction, Test, Random, Statistics

src_dir = "../"
include(src_dir * "src/ap3D.jl")

# The file structure in test roughly mirrors that of src.  Each file is included below.
@testset "ApogeeReduction.jl" begin
    include("ap3D.jl")
end
