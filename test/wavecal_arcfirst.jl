# In-order ("arc-first") wavelength chain: ordering and fallback logic (src/wavecal_arcfirst.jl).
using HDF5

@testset "wavecal arc_first" begin
    AR = ApogeeReduction

    @testset "order switch" begin
        @test AR.arcfirst_check_order("legacy") == "legacy"
        @test AR.arcfirst_check_order("arc_first") == "arc_first"
        @test_throws ArgumentError AR.arcfirst_check_order("sky_first")
    end

    @testset "stage table is in order (no stage reads its own or a later stage's product)" begin
        S = AR.ARCFIRST_STAGES
        written_by = Dict{String, Int}()
        for (i, st) in enumerate(S), w in st.writes
            @test !haskey(written_by, w)          # each product has one writer
            written_by[w] = i
        end
        for (i, st) in enumerate(S), r in st.reads
            haskey(written_by, r) || continue     # static input (almanac, line list, calibration files)
            @test written_by[r] < i
        end
        # the previous run's products the archive refit reads are not inputs of any in-order stage
        for st in S, r in st.reads
            @test !(r in ("wavecalNightAve", "waveCalNightfpiDither", "waveCalNightskyDither", "waveCalFPI",
                          "ar1Duni", "ar1Dunical"))
        end
        @test [s.stage for s in S] == ["peaks", "arc", "comb", "fpifib", "assemble", "resample"]
        # no stage reads another night's product (no cross-night step, AKS 2026-11-18)
        @test !any(st -> any(r -> r in ("anchorCheck", "crossnight_flags"), st.reads), S)
        @test !isdefined(AR, :refit_night_anchorcheck) && !isdefined(AR, :arcfirst_crossnight)
        @test !(:history_roots in Base.kwarg_decl(first(methods(AR.arcfirst_night))))
    end

    @testset "night route and comb decision" begin
        @test AR.arcfirst_night_route(true) == :arcfirst
        @test AR.arcfirst_night_route(false) == :legacy
        @test_throws ArgumentError AR.arcfirst_night_route(false, "arcwin")
        @test AR.arcfirst_comb_decision(0, false, false, false, false) == "nocomb"
        @test AR.arcfirst_comb_decision(4, false, true, true, true) == "nocomb_gated"   # G1 first
        @test AR.arcfirst_comb_decision(4, true, false, false, false) == "nocomb_gated" # provisional fit failed
        @test AR.arcfirst_comb_decision(4, true, true, false, true) == "nocomb_gated"   # G2
        @test AR.arcfirst_comb_decision(4, true, true, true, false) == "nocomb_unmatched"
        @test AR.arcfirst_comb_decision(4, true, true, true, true) == "ok_arc"
    end

    @testset "bits" begin
        @test length(AR.WAVE_BIT_NAMES) == 39
        for (i, n) in enumerate(AR.WAVE_BIT_NAMES)
            @test getfield(AR, Symbol(n)) == 1 << (i - 1)
        end
        nm, b = AR.arcfirst_legacy_type("waveCalNightskyDither")
        @test nm == "sky_skydith" && (b & AR.WAVE_NOT_ARCFIRST) != 0 && (b & AR.WAVE_SKY) != 0
        nm, b = AR.arcfirst_legacy_type("error_fixed_fallback")
        @test nm == "none_rough" && (b & AR.WAVE_ROUGH) != 0 && (b & AR.WAVE_NOT_ARCFIRST) != 0
    end

    @testset "exposure list from file names" begin
        fn = ["/x/ar1Dcal_apo_60110_0012_R_object.h5", "/x/ar1Dcal_apo_60110_0012_G_object.h5",
              "/x/ar1Dcal_apo_60110_0003_R_arclamp.h5", "/x/ar1Dcal_lco_60110_0004_R_object.h5",
              "/x/ar1Dcal_apo_60111_0005_R_object.h5"]
        @test AR.arcfirst_exps_from_names(fn, "apo", 60110) == [(3, "arclamp"), (12, "object")]
    end

    @testset "waveCalExp writer fills non-finite fibers and flags them" begin
        d = mktempdir()
        G = zeros(N_XPIX, N_FIBERS, N_CHIPS)
        for c in 1:N_CHIPS, f in 1:N_FIBERS
            G[:, f, c] .= 15000 .+ 1000c .+ 0.3 .* (1:N_XPIX) .+ 0.001f
        end
        fn = joinpath(d, "waveCalExp_apo_60110_0012_object.h5")
        @test AR.arcfirst_write_exp(fn, G, "arc_fpifib", 0, 1.5f0; tele = "apo", mjd = 60110, expnum = 12,
            kind = "object") == 0
        G[:, 7, 2] .= NaN
        add = AR.arcfirst_write_exp(fn, G, "arc_fpifib", AR.WAVE_FPIFIB_ROTATION, 1.5f0; tele = "apo",
            mjd = 60110, expnum = 12, kind = "object")
        @test add == AR.WAVE_FIBER_FILLED
        W, wt, b = h5open(f -> (read(f["chipWaveSoln"]), read(f["wavecal_type"]), read(f["bitmsk_wavecal"])), fn, "r")
        @test all(isfinite, W) && wt == "arc_fpifib"
        @test b == AR.WAVE_FPIFIB_ROTATION | AR.WAVE_FIBER_FILLED
        @test W[:, 7, 2] == G[:, 6, 2]            # nearest finite fiber (6 before 8 on a tie)
    end

    @testset "anchor exclusion falls to the next candidate" begin
        # pure selection kernel: excluding the first-ranked sequence (no comb for it) picks the next one
        tseq = [0.1, 0.5, 0.9]; q = [0, 0, 0]; sci = [0.45]
        cdt = [0.2, 0.3, 0.4]; ck = [1, 2, 3]
        S1 = AR.refit_anchor_select("nearest_science", tseq, q, cdt, ck, sci, falses(3), trues(3), falses(3))
        @test S1.ok && S1.j == 2
        ck2 = copy(ck); ck2[2] = 0; cdt2 = copy(cdt); cdt2[2] = Inf   # what exclude_seqs = [2] does
        S2 = AR.refit_anchor_select("nearest_science", tseq, q, cdt2, ck2, sci, falses(3), trues(3), falses(3))
        @test S2.ok && S2.j != 2 && (S2.bits & AR.WAVE_ANCHOR_FALLBACK) != 0
    end

    @testset "end-of-night check" begin
        d = mktempdir()
        mkpath(joinpath(d, "exp")); mkpath(joinpath(d, "comb"))
        asm = AR.refit_path(d, "assemble", "apo", 60110)
        combf = AR.refit_path(d, "comb", "apo", 60110)
        v = fill(80.0, 3, 4); v[:, 4] .= NaN
        h5open(asm, "w") do f
            f["wavecal_type"] = ["arc_fpifib", "arc_fpifib", "arc_combtime", "arc_fpifib"]
            f["sky_minus_model_ms"] = v
        end
        h5open(f -> (f["n_mode_slip"] = 0), combf, "w")
        c = AR.arcfirst_check_night(d, "apo", 60110)
        @test c.pass && c.verifiable && c.n_exp == 2 && c.metric_ms ≈ 80.0
        v[2, 1:2] .= 2600.0                       # a dither step on chip G: wrong anchor block
        h5open(f -> (delete_object(f, "sky_minus_model_ms"); f["sky_minus_model_ms"] = v), asm, "r+")
        c = AR.arcfirst_check_night(d, "apo", 60110)
        @test !c.pass && c.metric_ms ≈ 2600.0
        h5open(f -> (delete_object(f, "sky_minus_model_ms"); f["sky_minus_model_ms"] = fill(NaN, 3, 4)), asm, "r+")
        c = AR.arcfirst_check_night(d, "apo", 60110)
        @test c.pass && !c.verifiable             # no sky on FPI exposures: unverifiable, not failed
        h5open(f -> (delete_object(f, "n_mode_slip"); f["n_mode_slip"] = 3), combf, "r+")
        @test !AR.arcfirst_check_night(d, "apo", 60110).pass
    end

    @testset "arc fiber fill keeps fitted fibers" begin
        lay = AR.ArcGlobalLayout(1, [1], [0.0])
        TH = zeros(N_FIBERS, 9)
        for f in 1:N_FIBERS
            TH[f, :] .= [-1.07, 1.0, 1.075, 1.0, 16000.0 + 0.01f, 600.0, -5.0, 0.1, 0.0]
        end
        okf = trues(N_FIBERS); okf[[10, 150]] .= false; TH[150, :] .= NaN
        L = (seq = Int[], fiber = Int[], chip = Int[], id = Int[], pix = Float64[], pixvar = Float64[], keep = Bool[])
        A = AR.RefitArc("x", TH, okf, zeros(lay.n), lay, [1], [1], [0.0], [0.0], [1], 0.0, 0, [0], "fixed", L)
        B, n = AR.refit_arc_fill_fibers(A)
        @test n == 2 && all(B.okf) && all(isfinite, B.TH)
        @test B.TH[setdiff(1:N_FIBERS, [10, 150]), :] == TH[setdiff(1:N_FIBERS, [10, 150]), :]
        @test isapprox(B.TH[150, 5], 16000.0 + 1.5; atol = 0.05)   # the regulariser value (quadratic + block step)
        C, n0 = AR.refit_arc_fill_fibers(AR.RefitArc("x", ones(N_FIBERS, 9), trues(N_FIBERS), zeros(lay.n), lay,
            [1], [1], [0.0], [0.0], [1], 0.0, 0, [0], "fixed", L))
        @test n0 == 0
    end
end

module CrossnightQA
include(joinpath(@__DIR__, "..", "scripts", "arcfirst", "crossnight_comb_qa.jl"))
end

@testset "post-hoc cross-night comb QA kernel (standalone script, trailing only)" begin
    K = CrossnightQA
    mj = collect(60000:60010)
    L = fill(3.7363e7, 11); lam = fill(16000.0, 11)
    @test isempty(K.crossnight_comb_flags(mj, lam, L))
    lam2 = copy(lam); lam2[8:end] .+= lam[8]^2 / (2 * L[8])          # one FSR from night 60007 on
    f = K.crossnight_comb_flags(mj, lam2, L)
    @test length(f) == 1 && f[1][1] == 60007 && f[1][2] == "fsr_continuity"
    L2 = copy(L); L2[11] += 100.0                                     # cavity jump on the last night
    f = K.crossnight_comb_flags(mj, lam, L2)
    @test length(f) == 1 && f[1][1] == 60010 && f[1][2] == "cavity_jump_8sigma"
    @test K.crossnight_comb_flags(reverse(mj), reverse(lam), reverse(L2)) == f
end
