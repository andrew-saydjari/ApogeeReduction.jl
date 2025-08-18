using HDF5
using ApogeeReduction: safe_jldsave

@testset "safe_jldsave" begin
    path = tempdir() * "/test.h5"
    metadata = Dict{String, Any}()

    data = ones(Int, 10, 10)

    safe_jldsave(path, metadata; data)

    git_branch = h5read(path, "metadata/git_branch")
    git_commit = h5read(path, "metadata/git_commit")
    git_clean = h5read(path, "metadata/git_clean")

    @test git_branch == ApogeeReduction.git_branch
    @test git_commit == ApogeeReduction.git_commit
    @test git_clean == ApogeeReduction.git_clean

    # would be good to test here that the h5 file is "clean"
end
