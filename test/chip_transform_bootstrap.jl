using HDF5

# The chip-transform calibration must be rebuildable FROM SCRATCH: the rebuild driver's derivation
# stages (collect, derive, write) may never need an existing chip-transform file. This runs them on a
# synthetic reduction root (comb fits only) with no calibration file anywhere, the pipeline default
# pointing at a path that does not exist, and the SLURM-script variable pointing at another one.
@testset "chip_transform bootstrap (no calibration file)" begin
    mktempdir() do dir
        root = joinpath(dir, "reduction")
        fib = 1:N_FIBERS
        truth_R = 1 - 2.6e-4 .+ 3e-5 .* sin.(fib ./ 7)          # per-fiber structure, as the comb shows
        truth_B = 1 + 1.1e-4 .+ 2e-5 .* cos.(fib ./ 11)
        rng = MersenneTwister(11)
        mjds = 60000:60029
        npk = 330                                                # 110 comb peaks per chip
        for mjd in mjds
            d = joinpath(root, "apred", string(mjd)); mkpath(d)
            nl = hcat(fill(-1.0707, N_FIBERS), truth_R .+ 2e-6 .* randn(rng, N_FIBERS),
                fill(1.0752, N_FIBERS), truth_B .+ 2e-6 .* randn(rng, N_FIBERS))
            h5open(joinpath(d, "waveCalFPI_apo_$(mjd)_arclamp.h5"), "w") do f
                f["nlParams"] = nl
                f["resid_vec"] = 0.002 .* randn(rng, N_FIBERS, npk)
                f["resid_used_in_fit"] = ones(UInt8, npk, N_FIBERS)
                f["resid_chipInts"] = repeat(vcat(fill(1, 110), fill(2, 110), fill(3, 110)), 1, N_FIBERS)
            end
            mkpath(joinpath(root, "wavecal"))
            h5open(joinpath(root, "wavecal", "wavecalNightAve_apo_$(mjd).h5"), "w") do f
                f["best_wave_type"] = "fpi"
                f["fpi_qa/gate_pass"] = true
                f["fpi_qa/failed_guard"] = "none"
            end
        end
        wd = joinpath(dir, "build")
        driver = joinpath(pkgdir(ApogeeReduction), "scripts", "cal", "chip_transform", "build_chip_transform.jl")
        # no chip-transform file exists anywhere; both configuration routes point at nothing
        @test !isfile(joinpath(dir, "nonexistent", "chiptransform.h5"))
        cmd = addenv(`$(Base.julia_cmd()) --project=$(pkgdir(ApogeeReduction)) $(driver)
                      --stage collect,derive,write --root $(root) --tele apo --workdir $(wd)
                      --version bootstrap_test --cover_from 59000 --cover_to 99999 --qa_plots false`,
            "AR_CHIP_TRANSFORM_CAL" => joinpath(dir, "nonexistent", "chiptransform.h5"),
            "JULIA_NUM_THREADS" => "1", "OPENBLAS_NUM_THREADS" => "1")
        @test success(pipeline(cmd; stdout = devnull, stderr = devnull))
        calp = joinpath(wd, "bootstrap_test.h5")
        @test isfile(calp)
        cal = ApogeeReduction.load_chip_transform_cal(calp)
        @test length(cal.blocks["apo"]) == 1                      # one all-time block by default
        b = cal.blocks["apo"][1]
        @test (b.mjd_start, b.mjd_end, b.n_nights) == (59000, 99999, 30)
        @test maximum(abs.(b.scale[:, 1] .- truth_R)) < 3e-6      # per-fiber scales recovered
        @test maximum(abs.(b.scale[:, 2] .- truth_B)) < 3e-6
        @test startswith(cal.provenance["needs_chip_transform_file"], "no")
        @test isfile(joinpath(wd, "stage2_apo_detected_change_points.csv"))
        @test ApogeeReduction.chip_scales(cal, "apo", 60010) == b.scale
    end
end
