# Every Julia script outside src/ must at least parse. The package tests never load the pipeline
# entry points (pipeline.jl, pipeline_2d_1d.jl, scripts/...), so a syntax error in one of them
# (2026-11-18: unescaped quotes in a help string of pipeline_2d_1d.jl) otherwise passes CI and
# only shows up when a night is run.
@testset "scripts parse" begin
    root = dirname(@__DIR__)
    files = String[]
    for f in readdir(root; join = true)
        endswith(f, ".jl") && push!(files, f)
    end
    for (dir, _, fs) in walkdir(joinpath(root, "scripts"))
        for f in fs
            endswith(f, ".jl") && push!(files, joinpath(dir, f))
        end
    end
    @test !isempty(files)
    for f in files
        ex = Meta.parseall(read(f, String); filename = f)
        bad = [e for e in ex.args if e isa Expr && e.head in (:error, :incomplete)]
        isempty(bad) || @info "does not parse" file = relpath(f, root)
        @test isempty(bad)
    end
end
