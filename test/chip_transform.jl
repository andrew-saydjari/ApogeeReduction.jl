using HDF5

# Synthetic calibration with one block per telescope over [mjd_start, mjd_end].
function synthetic_chip_cal(path; blocks_apo = [(59000, 59999, 0.99974, 1.00011), (60000, 61000, 0.99973, 1.00010)],
        version = "test_v1")
    mk(ms, me, bR, bB) = ApogeeReduction.ChipTransformBlock(ms, me, ms, me, 10,
        hcat(fill(bR, N_FIBERS), fill(bB, N_FIBERS)), fill(1e-6, N_FIBERS, 2),
        hcat(fill(-1.0707, N_FIBERS), fill(1.0752, N_FIBERS)), falses(N_FIBERS))
    blocks = Dict("apo" => [mk(b...) for b in blocks_apo])
    prov = Dict("code_commit" => "abc123", "input_root" => "/some/root", "created_utc" => "2026-11-17T00:00:00")
    ApogeeReduction.write_chip_transform_cal(path, blocks, prov; calibration_version = version)
end

# per-chip RMS (m/s) of a wavelength error after removing a per-chip zero point and stretch
function shape_rms(G, T)
    u = ((1:N_XPIX) .- 1024.5) ./ 1024
    X = hcat(ones(N_XPIX), u)
    [begin
         v = 299792458.0 .* (G[:, c] .- T[:, c]) ./ T[:, c]
         r = v .- X * (X \ v)
         sqrt(sum(abs2, r[65:1984]) / 1920)
     end for c in 1:3]
end

@testset "chip_transform" begin
    @testset "file round trip and provenance" begin
        mktempdir() do dir
            p = joinpath(dir, "cal.h5")
            synthetic_chip_cal(p)
            cal = ApogeeReduction.load_chip_transform_cal(p)
            @test cal.calibration_version == "test_v1"
            @test cal.path == p
            @test cal.provenance["code_commit"] == "abc123"
            @test cal.provenance["input_root"] == "/some/root"
            @test length(cal.blocks["apo"]) == 2
            b = cal.blocks["apo"][2]
            @test (b.mjd_start, b.mjd_end, b.n_nights) == (60000, 61000, 10)
            @test all(b.scale[:, 1] .== 0.99973) && all(b.scale[:, 2] .== 1.00010)
            @test all(b.offset[:, 1] .== -1.0707)
            @test !any(b.fiber_filled)
            # a configured but missing file is a hard error, never a silent b = 1
            @test_throws ArgumentError ApogeeReduction.load_chip_transform_cal(joinpath(dir, "missing.h5"))
            # wrong format version
            h5open(p, "r+") do f
                delete_attribute(f, "format_version")
                attributes(f)["format_version"] = 99
            end
            @test_throws ArgumentError ApogeeReduction.load_chip_transform_cal(p)
            # overlapping blocks are rejected at write time
            @test_throws ArgumentError synthetic_chip_cal(joinpath(dir, "bad.h5");
                blocks_apo = [(59000, 60000, 1.0, 1.0), (60000, 61000, 1.0, 1.0)])
        end
    end

    @testset "block lookup by MJD" begin
        mktempdir() do dir
            p = joinpath(dir, "cal.h5")
            synthetic_chip_cal(p; blocks_apo = [(59000, 59499, 0.99974, 1.00011), (59600, 61000, 0.99973, 1.00010)])
            cal = ApogeeReduction.load_chip_transform_cal(p)
            cs = ApogeeReduction.chip_scales
            @test cs(cal, "apo", 59000)[1, 1] == 0.99974      # first day inclusive
            @test cs(cal, "apo", 59499)[1, 1] == 0.99974      # last day inclusive
            @test cs(cal, "apo", 59600)[1, 2] == 1.00010
            @test cs(cal, "apo", 61000)[1, 1] == 0.99973
            @test isnothing(cs(cal, "apo", 59550))             # gap between blocks
            @test isnothing(cs(cal, "apo", 58000))             # before the first block
            @test isnothing(cs(cal, "apo", 61001))             # after the last block
            @test isnothing(cs(cal, "lco", 60000))             # telescope not in the file
            @test isnothing(cs(nothing, "apo", 60000))
        end
    end

    @testset "data-driven block boundaries" begin
        rng = MersenneTwister(3)
        mjds = collect(1:200)
        y = 1.0 .+ 2e-6 .* randn(rng, 200, 2)
        @test isempty(ApogeeReduction.detect_blocks(mjds, y; min_nights = 20, nsigma = 5, min_step = [1e-5, 1e-5]))
        y2 = copy(y); y2[121:end, 1] .+= 3e-5                 # a step in b_R at night 121
        @test ApogeeReduction.detect_blocks(mjds, y2; min_nights = 20, nsigma = 5, min_step = [1e-5, 1e-5]) == [121]
        # a significant step smaller than min_step does not open a block
        y3 = copy(y); y3[121:end, 2] .+= 5e-6
        @test isempty(ApogeeReduction.detect_blocks(mjds, y3; min_nights = 20, nsigma = 5, min_step = [1e-5, 1e-5]))
        # two steps
        y4 = copy(y2); y4[61:end, 2] .-= 4e-5
        @test ApogeeReduction.detect_blocks(mjds, y4; min_nights = 20, nsigma = 5, min_step = [1e-5, 1e-5]) == [61, 121]
    end

    @testset "sky fit with a known chip scale (end to end)" begin
        # truth: one quartic in the green-chip coordinate shared by the chips, with the
        # APO comb-measured chip transform (b_R - 1 = -2.6e-4, b_B - 1 = +1.1e-4)
        aR, bR, aB, bB = -1.0707, 1 - 2.6e-4, 1.0752, 1 + 1.1e-4
        q = [16000.0, -590.0, -25.0, 3.0, 0.5]
        lam(xp) = sum(q[k + 1] * xp^k for k in 0:4)
        xprime(x, c) = c == 1 ? aR + bR * x : c == 3 ? aB + bB * x : x
        rng = MersenneTwister(7)
        nl = 40
        pix = sort(30 .+ 1990 .* rand(rng, nl))       # same pixels on every chip and fiber
        truth = [lam(xprime((p - 1024) / 2048, c)) for p in 1:N_XPIX, c in 1:3]
        mktempdir() do dir
            fname = joinpath(dir, "skyLinePeaks_apo_60500_0001_R_object.h5")
            for (c, chip) in enumerate(CHIP_LIST)
                M = fill(NaN, nl, 9, N_FIBERS)
                for f in 1:N_FIBERS, l in 1:nl
                    # tiny fiber-to-fiber jitter (1e-4 px) so the across-slit smoothing is not degenerate
                    M[l, 1, f] = pix[l] + 1e-4 * randn(rng)
                    M[l, 2, f] = lam(xprime((M[l, 1, f] - 1024) / 2048, c))
                end
                h5open(replace(fname, "_R_" => "_$(chip)_"), "w") do fo
                    fo["sky_line_mat_clean"] = M
                end
            end
            # old model: chip scale forced to 1
            out = ApogeeReduction.get_and_save_sky_wavecal(fname; checkpoint_mode = "clobber")
            G0, nl0 = h5open(f -> (read(f["chipWaveSoln"]), read(f["nlParams"])), out)
            @test size(nl0, 2) == 2                       # unchanged cporder-0 layout
            r0 = shape_rms(G0[:, 1, :], truth)
            @test maximum(r0) > 10.0                      # the forced scale bends the chips (m/s)

            # new model: chip scale fixed to the calibration file's value
            calp = joinpath(dir, "cal.h5")
            synthetic_chip_cal(calp; blocks_apo = [(60000, 61000, bR, bB)])
            cal = ApogeeReduction.load_chip_transform_cal(calp)
            out1 = ApogeeReduction.get_and_save_sky_wavecal(fname; checkpoint_mode = "clobber", chip_transform_cal = cal)
            G1, nl1, meta = h5open(f -> (read(f["chipWaveSoln"]), read(f["nlParams"]), read(f["metadata"])), out1)
            @test size(nl1, 2) == 4                       # cporder-1 layout a_R, b_R, a_B, b_B
            @test all(nl1[:, 2] .== bR) && all(nl1[:, 4] .== bB)   # fixed, not smoothed or refitted
            @test isapprox(median(nl1[:, 1]), aR; atol = 1e-6)
            @test isapprox(median(nl1[:, 3]), aB; atol = 1e-6)
            r1 = shape_rms(G1[:, 1, :], truth)
            @test maximum(r1) < 0.5                       # no curvature left (m/s)
            v1 = 299792458.0 .* (G1[65:1984, 1, :] .- truth[65:1984, :]) ./ truth[65:1984, :]
            @test maximum(abs.(v1)) < 1.0                 # and the absolute solution is right
            @test meta["chip_transform_cal"] == calp
            @test meta["chip_transform_version"] == "test_v1"

            @test meta["chip_transform_id"] == ApogeeReduction.chip_transform_id(cal)
            @test meta["chip_transform_fingerprint"] == ApogeeReduction.file_fingerprint(calp)

            # checkpointing ("commit_exists": the stored product passes check_file whatever the commit):
            # same calibration -> reused; different calibration, or none -> recomputed
            t1 = mtime(out1)
            sleep(1.1)
            ApogeeReduction.get_and_save_sky_wavecal(fname; checkpoint_mode = "commit_exists", chip_transform_cal = cal)
            @test mtime(out1) == t1                      # reused
            out0 = ApogeeReduction.get_and_save_sky_wavecal(fname; checkpoint_mode = "commit_exists")
            @test mtime(out0) > t1                       # off: recomputed with b = 1
            @test size(h5open(f -> read(f["nlParams"]), out0), 2) == 2
            @test h5open(f -> read(f["metadata"]), out0)["chip_transform_id"] == "none"
            calp3 = joinpath(dir, "cal3.h5")
            synthetic_chip_cal(calp3; blocks_apo = [(60000, 61000, bR + 1e-6, bB)], version = "test_v1")
            cal3 = ApogeeReduction.load_chip_transform_cal(calp3)
            @test cal3.fingerprint != cal.fingerprint    # same version string, different content
            ApogeeReduction.get_and_save_sky_wavecal(fname; checkpoint_mode = "commit_exists", chip_transform_cal = cal)
            t2 = mtime(out1); sleep(1.1)
            ApogeeReduction.get_and_save_sky_wavecal(fname; checkpoint_mode = "commit_exists", chip_transform_cal = cal3)
            @test mtime(out1) > t2                       # changed file -> recomputed
            @test h5open(f -> read(f["nlParams"]), out1)[1, 2] == bR + 1e-6

            # an MJD the file does not cover is an error, never a silent b = 1
            synthetic_chip_cal(calp; blocks_apo = [(59000, 59999, bR, bB)])
            cal2 = ApogeeReduction.load_chip_transform_cal(calp)
            @test_throws ArgumentError ApogeeReduction.get_and_save_sky_wavecal(fname; checkpoint_mode = "clobber",
                chip_transform_cal = cal2)
        end
    end
end
