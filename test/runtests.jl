using ApogeeReduction, Test, Random, Statistics

# The file structure in test roughly mirrors that of src.  Each file is included below.
@testset verbose=true "ApogeeReduction.jl" begin
    include("ap3D.jl")
end
