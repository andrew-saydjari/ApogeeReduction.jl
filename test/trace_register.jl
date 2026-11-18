# Per-exposure registration of the night trace to the science light (src/trace_register.jl).
# Uses the synthetic-image helpers of test/slit_model.jl (sm_template, sm_truth, sm_image), which
# runtests.jl includes first.
using Random, Statistics
using ApogeeReduction: TemplateRegistration, slit_profile, trace_reference_positions,
                       register_exposure_trace, fit_exposure_registration, exposure_trace_offsets,
                       apply_exposure_registration, exposure_registration_metadata,
                       exposure_registration_bits, TRACE_EXPOSURE_REGISTERED,
                       TRACE_EXPOSURE_REG_FAILED, TRACE_BIT_NAMES, trace_bit_names, TraceBits, bad_pix_bits

@testset "per-exposure trace registration" begin
    R = trace_reference_positions("apo", "G"; data_path = SM_DATA)
    P = slit_profile("apo", "G"; data_path = SM_DATA)
    SM_P[] = P
    rng = MersenneTwister(11)
    href = 0.6 .+ 0.8 .* rand(rng, 300)
    t = sm_template(R; href = href)
    night = sm_truth(t, TemplateRegistration([0.4, 0.05, -0.1], [0.0002, 0.0]); heights = 1000.0 .* href)
    # the science frame: the same slit moved by a known rigid amount; a third of the fibers empty
    # (sky-faint) and a few very bright (stars)
    moved(a, b, c) = (T = copy(night);
        for x in axes(T, 1), f in axes(T, 2)
            y = night[x, f, 2]
            T[x, f, 2] = y + a + b * (x - 1024) / 1024 + c * (y - 1024) / 1024
        end; T)
    h = 300.0 .* href
    h[rand(MersenneTwister(2), 1:300, 100)] .= 0.0
    h[1:25:300] .*= 20

    @testset "recovers shift and rotation (and scale when asked)" begin
        truth = moved(0.07, -0.03, 0.02)
        truth[:, :, 1] .= reshape(h, 1, :)
        img, iv, m = sm_image(truth; noise = 3.0)
        reg = register_exposure_trace(img, iv, m, night, P)              # default :shift_rot
        @test reg.outcome == "ok" && reg.model == :shift_rot
        @test abs(reg.a - 0.07) < 0.006 && abs(reg.b + 0.03) < 0.006 && reg.c == 0
        @test reg.nfib >= 150 && reg.n >= 1000 && reg.half_a < 0.01
        reg3 = register_exposure_trace(img, iv, m, night, P; model = :shift_rot_scale)
        @test reg3.outcome == "ok"
        @test abs(reg3.a - 0.07) < 0.004 && abs(reg3.b + 0.03) < 0.004 && abs(reg3.c - 0.02) < 0.004
        # the moved trace sits on the light
        T3 = apply_exposure_registration(night, reg3)
        @test maximum(abs.(T3[17:2032, h .> 0, 2] .- truth[17:2032, h .> 0, 2])) < 0.01
        @test T3[:, :, [1, 3]] == night[:, :, [1, 3]]                    # heights and widths kept
        md = exposure_registration_metadata(reg3)
        @test md["trace_expreg_outcome"] == "ok" && md["trace_expreg_model"] == "shift_rot_scale"
        @test md["trace_expreg_a"] == reg3.a && haskey(md, "trace_expreg_half_a")
        @test exposure_registration_bits(reg3) == TRACE_EXPOSURE_REGISTERED
        @test exposure_registration_bits(nothing) == zero(TraceBits)
        # an unmoved frame registers to ~0
        night_img = copy(night); night_img[:, :, 1] .= reshape(h, 1, :)
        img0, iv0, m0 = sm_image(night_img; noise = 3.0, rng = MersenneTwister(9))
        r0 = register_exposure_trace(img0, iv0, m0, night, P)
        @test r0.outcome == "ok" && abs(r0.a) < 0.004 && abs(r0.b) < 0.004
        # masked pixels are ignored (bad_pix_bits): a masked half chip still registers
        mm = copy(m); mm[1:1024, :] .= bad_pix_bits
        rm_ = register_exposure_trace(img, iv, mm, night, P)
        @test rm_.outcome == "ok" && abs(rm_.a - 0.07) < 0.01 && abs(rm_.b + 0.03) < 0.01
    end

    @testset "refusals fall back to the night trace" begin
        # no light: too few fiber-bands
        dark = copy(night); dark[:, :, 1] .= 0
        img, iv, m = sm_image(dark; noise = 3.0)
        reg = register_exposure_trace(img, iv, m, night, P)
        @test reg.outcome == "too_few"
        @test apply_exposure_registration(night, reg) == night
        @test exposure_registration_bits(reg) == TRACE_EXPOSURE_REG_FAILED
        @test exposure_registration_metadata(reg)["trace_expreg_outcome"] == "too_few"
        # a shift beyond the accepted range is refused, not applied
        truth = moved(0.07, 0.0, 0.0); truth[:, :, 1] .= reshape(h, 1, :)
        img, iv, m = sm_image(truth; noise = 3.0)
        reg = register_exposure_trace(img, iv, m, night, P; amax = 0.05)
        @test reg.outcome == "too_large" && abs(reg.a - 0.07) < 0.006
        @test apply_exposure_registration(night, reg) == night
        # inconsistent halves: even and odd fibers disagree (no rigid motion explains the light)
        n = 4000
        rng2 = MersenneTwister(4)
        f = rand(rng2, 1:300, n)
        mfake = (b = rand(rng2, 1:32, n), f = f, u = 2 .* rand(rng2, n) .- 1, v = 2 .* rand(rng2, n) .- 1,
            r = ifelse.(iseven.(f), 0.05, -0.05) .+ 0.002 .* randn(rng2, n), sig = fill(0.003, n),
            snr = fill(100.0, n))
        reg = fit_exposure_registration(mfake)
        @test reg.outcome == "inconsistent_halves" && reg.half_a > 0.04
        @test_throws ErrorException fit_exposure_registration(mfake; model = :bow)
        # fibers that do not move together (per-fiber scatter 0.15 px): not a rigid motion, refused
        off = 0.15 .* randn(MersenneTwister(8), 300)
        mscat = (b = mfake.b, f = f, u = mfake.u, v = mfake.v, r = 0.05 .+ off[f] .+ 0.002 .* randn(rng2, n),
            sig = fill(0.003, n), snr = fill(100.0, n))
        reg = fit_exposure_registration(mscat)
        @test reg.outcome == "not_rigid" && reg.rms > 0.08
        @test fit_exposure_registration(mscat; rms_max = 1.0, half_max = 1.0).outcome == "ok"
    end

    @testset "trace bits" begin
        @test TRACE_BIT_NAMES[TRACE_EXPOSURE_REGISTERED] == "TRACE_EXPOSURE_REGISTERED"
        @test trace_bit_names(TRACE_EXPOSURE_REG_FAILED) == ["TRACE_EXPOSURE_REG_FAILED"]
    end
end
