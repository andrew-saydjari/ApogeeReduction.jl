using HDF5
using SHA: sha256

# Functions that only READ a product file must open it read-only. Opening with
# "r+" leaves the bytes identical but bumps the file's mtime, cannot run against
# a read-only filesystem or snapshot, and takes a write lock. Each reader below
# is called on a temp copy that has been chmod'ed 0444; the call must succeed
# and leave both the bytes and the mtime untouched.

function _freeze(path)
    chmod(path, 0o444)
    return (sha = sha256(read(path)), mtime = stat(path).mtime)
end

function _unchanged(path, before)
    return sha256(read(path)) == before.sha && stat(path).mtime == before.mtime
end

@testset "read-only opens" begin
    dir = mktempdir()
    # the directory stays writable (some readers write a separate output file
    # next to their input); only the input product itself is read-only
    n_fib = ApogeeReduction.N_FIBERS
    n_pix = ApogeeReduction.N_XPIX

    @testset "ingest_fpiLines_file" begin
        p = joinpath(dir, "fpiPeaks_apo_60000_0001_R_arclamp.h5")
        mat = fill(1.0, 5, 4, n_fib)
        cov = fill(0.01, 5, 4, 4, n_fib)
        h5open(p, "w") do f
            write(f, "fpi_line_mat", mat)
            write(f, "fpi_line_cov_mat", cov)
        end
        before = _freeze(p)
        xlst, errs = ApogeeReduction.ingest_fpiLines_file(p)
        @test size(xlst) == (5, n_fib)
        @test all(errs .≈ 0.1)
        @test _unchanged(p, before)
    end

    @testset "ingest_skyLines_file" begin
        p = joinpath(dir, "skyLinePeaks_apo_60000_0002_R_object.h5")
        h5open(p, "w") do f
            write(f, "sky_line_mat_clean", fill(1024.0, 4, 3, n_fib))
        end
        before = _freeze(p)
        xlst, wlst = ApogeeReduction.ingest_skyLines_file(p)
        @test all(xlst .== 0.0)
        @test all(wlst .== 1024.0)
        @test _unchanged(p, before)
    end

    @testset "get_ave_night_wave_soln" begin
        paths = String[]
        for expid in ("0003", "0004")
            p = joinpath(dir, "waveCalNightskyDither_apo_60000_$(expid)_object.h5")
            linParams = repeat([16000.0 300.0 1.0], n_fib)
            nlParams = zeros(n_fib, 4)
            h5open(p, "w") do f
                write(f, "linParams", linParams)
                write(f, "nlParams", nlParams)
            end
            push!(paths, p)
        end
        befores = map(_freeze, paths)
        med_lin, med_nl, soln = ApogeeReduction.get_ave_night_wave_soln(copy(paths))
        @test med_lin ≈ repeat([16000.0 300.0 1.0], n_fib)
        @test size(soln) == (n_pix, n_fib, ApogeeReduction.N_CHIPS)
        @test all(_unchanged.(paths, befores))
    end

    # synthetic ar1Dcal product (all-zero flux): exercises the read at the top
    # of each peak finder; with no signal they find no peaks and return early
    function _write_ar1d(p; n_fibers = n_fib)
        ApogeeReduction.safe_jldsave(p;
            flux_1d = zeros(n_pix, n_fibers),
            ivar_1d = ones(n_pix, n_fibers),
            mask_1d = zeros(Int, n_pix, n_fibers),
            extract_trace_centers = repeat(collect(1.0:n_fibers)', n_pix),
            no_metadata = true)
        return p
    end

    @testset "get_and_save_sky_peaks" begin
        p = _write_ar1d(joinpath(dir, "ar1Dcal_apo_60000_0005_R_object.h5"))
        before = _freeze(p)
        @test isnothing(ApogeeReduction.get_and_save_sky_peaks(
            p, nothing, nothing; checkpoint_mode = "clobber"))
        @test _unchanged(p, before)
    end

    @testset "get_and_save_arclamp_peaks" begin
        p = _write_ar1d(joinpath(dir, "ar1Dcal_apo_60000_0006_R_arclamp.h5"); n_fibers = 2)
        before = _freeze(p)
        ApogeeReduction.get_and_save_arclamp_peaks(p; checkpoint_mode = "clobber")
        @test _unchanged(p, before)
    end

    @testset "get_and_save_fpi_peaks" begin
        p = _write_ar1d(joinpath(dir, "ar1Dcal_apo_60000_0007_R_arclamp.h5"); n_fibers = 2)
        before = _freeze(p)
        ApogeeReduction.get_and_save_fpi_peaks(p;
            data_path = joinpath(pkgdir(ApogeeReduction), "data"), checkpoint_mode = "clobber")
        @test _unchanged(p, before)
    end

    # Tripwire: every remaining read-write open of an existing file must be a
    # deliberate writer. If you add one, add it here with the reason.
    @testset "no stray read-write opens" begin
        root = pkgdir(ApogeeReduction)
        allowed = Dict(
            "src/fpi_gate.jl" => 1,                       # save_fpi_qa!: rewrites fpi_qa/
            "src/utils.jl" => 2,                          # safe_jldsave metadata; modify_saved_dict
            "src/wavecal.jl" => 1,                        # FPI gate promotes best_wave_type
            "src/trace_qa.jl" => 1,                       # make_traceMain: bitmsk_trace into its own temp traceMain
            "scripts/cal/decorate_almanac_exptype.jl" => 1 # writes exposure_class/
        )
        files = String[]
        for d in ("src", "scripts")
            for (r, _, fs) in walkdir(joinpath(root, d)), fn in fs
                endswith(fn, ".jl") && push!(files, joinpath(r, fn))
            end
        end
        append!(files, filter(f -> endswith(f, ".jl"), readdir(root; join = true)))
        rw = r"(h5open|jldopen)\([^\n]*\"(r\+|a|a\+)\""
        found = Dict{String, Int}()
        for f in files
            n = count(l -> !startswith(lstrip(l), "#") && occursin(rw, l), eachline(f))
            n > 0 && (found[relpath(f, root)] = n)
        end
        @test found == allowed
    end
end
