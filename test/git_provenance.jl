const _LG = ApogeeReduction.LibGit2   # stdlib via the package (not a test dependency)
@testset "git provenance read at runtime" begin
    # must not be a precompile-time constant
    @test !isconst(ApogeeReduction, :git_commit) && !isconst(ApogeeReduction, :git_clean)
    root = pkgdir(ApogeeReduction)
    if isdir(joinpath(root, ".git")) || isfile(joinpath(root, ".git"))
        @test ApogeeReduction.git_commit == string(_LG.head(root * "/"))
        @test ApogeeReduction.git_clean == !_LG.isdirty(_LG.GitRepo(root))
        # refreshing reads the current state again
        b, c, k = ApogeeReduction.refresh_git_provenance!()
        @test c == string(_LG.head(root * "/"))
    end
end
