using HDF5

# Build a synthetic fpiPeaks file. `lit` selects a live FPI (bright, narrow,
# high-SNR peaks) or a dark one (the seed grid returned unchanged: amplitudes of
# a fraction of a count, widths pinned at the FPI_WIDTH_CEILING).
function _write_fpipeaks(path; n_peaks = 130, n_fibers = 300, lit = true,
        lit_fibers = nothing, height = 900.0, width = 1.02, height_err = 4.0)
    mat = fill(NaN, n_peaks, 4, n_fibers)
    cov = fill(NaN, n_peaks, 4, 4, n_fibers)
    for i in 1:n_fibers
        fiber_lit = isnothing(lit_fibers) ? lit : (i in lit_fibers)
        for p in 1:n_peaks
            mat[p, 2, i] = 5.6 + 15.7 * (p - 1)          # the seed grid: always populated
            if fiber_lit
                mat[p, 1, i] = height
                mat[p, 3, i] = width
                mat[p, 4, i] = 11.0
                cov[p, 1, 1, i] = height_err^2
            else
                mat[p, 1, i] = 0.3
                mat[p, 3, i] = ApogeeReduction.FPI_WIDTH_CEILING
                mat[p, 4, i] = 0.0
                cov[p, 1, 1, i] = 1.4^2
            end
        end
    end
    h5open(path, "w") do f
        write(f, "fpi_line_mat", mat)
        write(f, "fpi_line_cov_mat", cov)
    end
    return path
end

function _write_wavecalfpi(path; n_fibers = 8, n_peaks = 20, resid = 0.004, n_used = nothing)
    r = fill(resid, n_fibers, n_peaks)
    used = ones(Int, n_peaks, n_fibers)
    if !isnothing(n_used)
        used .= 0
        idx = 1
        for j in 1:n_fibers, p in 1:n_peaks
            idx > n_used && break
            used[p, j] = 1
            idx += 1
        end
    end
    h5open(path, "w") do f
        write(f, "resid_vec", r)
        write(f, "resid_used_in_fit", used)
    end
    return path
end

@testset "fpi_gate" begin
    dir = mktempdir()

    @testset "fpi_peak_lit_fraction" begin
        live = _write_fpipeaks(joinpath(dir, "fpiPeaks_lco_60077_0009_R_arclamp.h5"), lit = true)
        dark = _write_fpipeaks(joinpath(dir, "fpiPeaks_lco_60031_0008_R_arclamp.h5"), lit = false)

        lf_live = ApogeeReduction.fpi_peak_lit_fraction(live)
        lf_dark = ApogeeReduction.fpi_peak_lit_fraction(dark)

        @test length(lf_live) == 300
        @test all(lf_live .== 1.0)
        # the dark case is exactly the observed failure: the peak list is full,
        # but nothing in it is a detection
        @test all(lf_dark .== 0.0)
    end

    @testset "fpi_light_check: live FPI passes" begin
        base = joinpath(dir, "live", "fpiPeaks_lco_60077_0009_R_arclamp.h5")
        mkpath(dirname(base))
        for chip in ApogeeReduction.CHIP_LIST
            _write_fpipeaks(replace(base, "_R_" => "_$(chip)_"), lit = true)
        end
        res = ApogeeReduction.fpi_light_check([base]; fpi_fiberIndxs = [88, 219])
        @test res.n_files == 3
        @test res.lit_all == 1.0
        @test res.lit_fpifib == 1.0
        @test res.pass
    end

    @testset "fpi_light_check: dark FPI fails" begin
        base = joinpath(dir, "dark", "fpiPeaks_lco_60031_0008_R_arclamp.h5")
        mkpath(dirname(base))
        for chip in ApogeeReduction.CHIP_LIST
            _write_fpipeaks(replace(base, "_R_" => "_$(chip)_"), lit = false)
        end
        res = ApogeeReduction.fpi_light_check([base]; fpi_fiberIndxs = [88, 219])
        @test res.lit == 0.0
        @test !res.pass
    end

    @testset "fpi_light_check: no almanac fibers falls back to all-fiber" begin
        base = joinpath(dir, "dark", "fpiPeaks_lco_60031_0008_R_arclamp.h5")
        res = ApogeeReduction.fpi_light_check([base]; fpi_fiberIndxs = Int[])
        @test isnan(res.lit_fpifib)
        @test res.lit == res.lit_all == 0.0
        @test !res.pass       # lco 59894: almanac names no bonus stub, gate must still fire
    end

    @testset "fpi_light_check: guide fibers lit but bulk dark still fails" begin
        # measured shape of apo 60255: lit_fpifib 0.994, lit_all 0.085
        base = joinpath(dir, "mixed", "fpiPeaks_apo_60255_0009_R_arclamp.h5")
        mkpath(dirname(base))
        for chip in ApogeeReduction.CHIP_LIST
            _write_fpipeaks(replace(base, "_R_" => "_$(chip)_"), lit_fibers = Set([88, 219]))
        end
        res = ApogeeReduction.fpi_light_check([base]; fpi_fiberIndxs = [88, 219])
        @test res.lit_fpifib == 1.0
        @test res.lit_all == 0.0
        @test res.lit == 0.0
        @test !res.pass
    end

    @testset "fpi_light_check: bulk lit but guide fibers dark still fails" begin
        base = joinpath(dir, "mixed2", "fpiPeaks_lco_60077_0009_R_arclamp.h5")
        mkpath(dirname(base))
        lit_fibers = Set(setdiff(1:300, [88, 219]))
        for chip in ApogeeReduction.CHIP_LIST
            _write_fpipeaks(replace(base, "_R_" => "_$(chip)_"), lit_fibers = lit_fibers)
        end
        res = ApogeeReduction.fpi_light_check([base]; fpi_fiberIndxs = [88, 219])
        @test res.lit_all == 1.0
        @test res.lit_fpifib == 0.0
        @test !res.pass
    end

    @testset "fpi_light_check: no files at all" begin
        res = ApogeeReduction.fpi_light_check([joinpath(dir, "nope", "fpiPeaks_lco_1_0001_R_arclamp.h5")])
        @test res.n_files == 0
        @test !res.pass
    end

    @testset "fpi_resid_stats" begin
        good = _write_wavecalfpi(joinpath(dir, "waveCalFPI_lco_60077_arclamp.h5"), resid = 0.0043)
        bad = _write_wavecalfpi(joinpath(dir, "waveCalFPI_lco_60031_arclamp.h5"), resid = 0.45)
        none = _write_wavecalfpi(joinpath(dir, "waveCalFPI_apo_59573_arclamp.h5"), resid = 0.01, n_used = 0)

        g = ApogeeReduction.fpi_resid_stats(good)
        @test isapprox(g.rms, 0.0043; atol = 1e-9)
        @test g.n_used == 160
        @test g.pass

        b = ApogeeReduction.fpi_resid_stats(bad)
        @test isapprox(b.rms, 0.45; atol = 1e-9)
        @test !b.pass

        # apo 59573: the fit used zero peaks and the seed solution was promoted
        n = ApogeeReduction.fpi_resid_stats(none)
        @test n.n_used == 0
        @test isnan(n.rms)
        @test !n.pass

        @test !ApogeeReduction.fpi_resid_stats(joinpath(dir, "does_not_exist.h5")).pass
    end

    @testset "thresholds bracket the measured corpus" begin
        # DR21 200-MJD testbed, 155 telescope-nights with fpiPeaks (2026-09-08)
        worst_healthy_lit = 0.96899
        best_rejected_lit = 0.08462
        worst_healthy_rms = 0.03001
        best_rejected_rms = 0.28364
        @test best_rejected_lit < ApogeeReduction.FPI_LIT_FRACTION_MIN < worst_healthy_lit
        @test worst_healthy_rms < ApogeeReduction.FPI_RESID_RMS_MAX < best_rejected_rms
    end

    @testset "save_fpi_qa! round trip" begin
        p = joinpath(dir, "wavecalNightAve_lco_60031.h5")
        h5open(p, "w") do f
            f["best_wave_type"] = "sky"
        end
        qa = Dict{String, Any}("gate_pass" => false, "failed_guard" => "fpi_light",
            "lit_fraction" => 0.0, "mjd" => 60031, "best_wave_type" => "sky")
        ApogeeReduction.save_fpi_qa!(p, qa)
        got = h5open(p, "r") do f
            (read(f["best_wave_type"]), read(f["fpi_qa/failed_guard"]),
                read(f["fpi_qa/lit_fraction"]), read(f["fpi_qa/gate_pass"]))
        end
        @test got[1] == "sky"          # the gate must not disturb the sky verdict
        @test got[2] == "fpi_light"
        @test got[3] == 0.0
        @test got[4] == false

        # a second call must overwrite rather than throw
        qa["failed_guard"] = "fpi_resid_rms"
        ApogeeReduction.save_fpi_qa!(p, qa)
        @test h5open(f -> read(f["fpi_qa/failed_guard"]), p, "r") == "fpi_resid_rms"

        # QA bookkeeping must never be able to fail a night
        @test ApogeeReduction.save_fpi_qa!(joinpath(dir, "missing.h5"), qa) === nothing
    end
end
