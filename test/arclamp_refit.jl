using ApogeeReduction: RefitConfig, refit_check_config, refit_config_dict, RefitArc, ArcGlobalLayout,
                       refit_arc_state, refit_arc_lin_nl, refit_prod_wave, arcfit_eval, arcfit_seq_state,
                       refit_mode_slips, refit_order_change, refit_fit_ab, refit_comb_delta,
                       refit_fpifib_transform, refit_pooled_sky, refit_sky_terms, refit_apply_sky_terms!,
                       refit_coeff_basis, refit_grid_to_coeffs, refit_coeffs_to_grid, refit_pass1_type,
                       refit_arc_night_bits, refit_interp_dt, RefitTemplate, refit_template_lambda,
                       refit_ref_wave_disp, refit_ref_invert, refit_vs_pass1, refit_comb_vs_arc,
                       refit_fiber_sky_terms, WAVE_SKY, WAVE_ROUGH, WAVE_NIGHT_AVE, WAVE_ARC_SEQ_GATED,
                       WAVE_ARC_MULTIPOS, WAVE_CHIPSCALE_FREE, WAVE_FPIFIB_ROTATION, WAVE_FPIFIB_UNAVAILABLE,
                       WAVE_FPIFIB_ONE_FIBER, ARCSEQ_GATE_FAIL, ARCNIGHT_MULTIPOS, ARCNIGHT_CHIPSCALE_FREE,
                       N_FIBERS, N_XPIX

# a synthetic arc night: 3 used sequences, two position groups, time terms on
function _rf_arc(; rng = MersenneTwister(3))
    S = 3
    grp = [1, 2, 1]
    tau = [-0.6, 0.0, 0.6]
    lay = ArcGlobalLayout(S, grp, tau)
    g = 1e-4 .* randn(rng, lay.n)
    g[lay.idx_d[2]] = 0.5 / 2048          # half-pixel step for sequence 2
    TH = zeros(N_FIBERS, 9)
    for f in 1:N_FIBERS
        TH[f, :] .= (-1.07 + 1e-4 * randn(rng), 1.0 - 2.6e-4, 1.0755, 1.0 + 1.1e-4,
            16200.0 + 0.01f, 1300.0, -40.0, 3.0, -0.5)
    end
    lines = (seq = Int[], fiber = Int[], chip = Int[], id = Int[], pix = Float64[], pixvar = Float64[], keep = Bool[])
    RefitArc("synthetic", TH, trues(N_FIBERS), g, lay, [1, 2, 3], grp, tau, 60000.0 .+ [0.1, 0.2, 0.4],
        [1, 2, 1], 60000.25, 0, [0, ARCSEQ_GATE_FAIL, 0], "fixed", lines)
end

@testset "arclamp_refit" begin
    @testset "config" begin
        cfg = RefitConfig()
        @test refit_check_config(cfg)
        @test cfg.exp_terms == "chipzp_stretch_tilt" && cfg.prefps_time == "regularised" && cfg.neon == "none"
        @test cfg.sky_exp_terms == "chipzp_stretch_tilt" && cfg.reg_weight_apo == 0.825 && cfg.reg_weight_lco == 0.65
        @test_throws ArgumentError refit_check_config(RefitConfig(exp_terms = "bogus"))
        @test_throws ArgumentError refit_check_config(RefitConfig(prefps_time = "both"))
        d = refit_config_dict(cfg)
        @test d["sky_fiber_term"] == "full" && haskey(d, "chip_transform_file")
    end

    @testset "bits and naming" begin
        @test refit_pass1_type("waveCalNightfpiDither") == ("sky_skydith", WAVE_SKY)
        nm, b = refit_pass1_type("error_night_ave_fpi")
        @test nm == "sky_nightpos" && (b & WAVE_NIGHT_AVE) != 0
        @test refit_pass1_type("error_fixed_fallback")[2] == WAVE_ROUGH
        b = refit_arc_night_bits(ARCNIGHT_MULTIPOS | ARCNIGHT_CHIPSCALE_FREE, [0, ARCSEQ_GATE_FAIL])
        @test (b & WAVE_ARC_SEQ_GATED) != 0 && (b & WAVE_ARC_MULTIPOS) != 0 && (b & WAVE_CHIPSCALE_FREE) != 0
    end

    @testset "arc state -> production comb form is exact" begin
        A = _rf_arc()
        @test refit_interp_dt([1.0, 3.0], [0.0, 2.0], [1.0, 1.0], 2.0) == (1.0, 1.0)
        @test refit_interp_dt([1.0, 3.0], [0.0, 2.0], [1.0, 1.0], 5.0) == (2.0, 1.0)
        px = collect(1.0:37.0:2048.0)
        for cl in (1, 2), t in (60000.05, 60000.3)
            st = refit_arc_state(A, t, cl)
            lin, nl, nfill = refit_arc_lin_nl(A, st)
            @test nfill == 0
            for f in (1, 150, 300), c in 1:3
                a = arcfit_eval(A.TH[f, :], A.g, A.lay, f, c, px; grpv = st[1], d = st[2], t = st[3], tau = st[4])
                b = refit_prod_wave(lin, nl, f, c, px)
                @test maximum(abs.(a .- b)) < 1e-8
            end
        end
        # a fiber without an arc fit takes its neighbour's seed
        A.okf[10] = false
        lin, nl, nfill = refit_arc_lin_nl(A, refit_arc_state(A, 60000.2, 1))
        @test nfill == 1 && lin[10, :] == lin[9, :] || lin[10, :] == lin[11, :]
    end

    @testset "mode slips and order changes" begin
        d = fill(0.02, N_FIBERS, 3)
        d[5, 2] = 1.05
        d[6, 1] = -0.97
        d[7, 3] = 0.5                       # far from any integer: not a slip
        n, m = refit_mode_slips(d)
        @test n == 2 && m[5, 2] && m[6, 1] && !m[7, 3]
        pk_old = repeat(collect(1:20), 1, N_FIBERS)
        chip = repeat(vcat(fill(1, 7), fill(2, 7), fill(3, 6)), 1, N_FIBERS)
        pk_new = copy(pk_old)
        pk_new[8:14, 42] .+= 1              # fiber 42, chip G: one FSR
        oc = refit_order_change(pk_new, pk_old, chip, trues(20, N_FIBERS))
        @test oc[42, 2] == 1 && count(!=(0), oc) == 1
    end

    @testset "comb displacement and FPI-fiber transform" begin
        rng = MersenneTwister(5)
        L = 3.7363e7
        eps = 0.22
        lin = repeat([16200.0 1300.0 -40.0 3.0 -0.5], N_FIBERS)
        nl = repeat([-1.07 1.0 1.0755 1.0], N_FIBERS)
        # true comb peaks (m integers) displaced by delta0(p) = 0.31 + 0.05 t px on fiber 76, chip G
        f, c = 76, 2
        px = collect(30.0:1.0:2020.0)
        lam = refit_prod_wave(lin, nl, f, c, px)
        ms = unique(round.(2L ./ lam .- eps))
        xs = Float64[]
        for m in ms
            lm = 2L / (m + eps)
            i = argmin(abs.(lam .- lm))
            (5 < i < length(px) - 5) || continue
            x = px[i]
            for _ in 1:20
                l = refit_prod_wave(lin, nl, f, c, [x])[1]
                dl = refit_prod_wave(lin, nl, f, c, [x + 0.5])[1] - refit_prod_wave(lin, nl, f, c, [x - 0.5])[1]
                x -= (l - lm) / dl
            end
            push!(xs, x + 0.31 + 0.05 * (x - 1024) / 1024 + 0.002 * randn(rng))
        end
        dl, mf = refit_comb_delta(xs, f, c, lin, nl, L, eps)
        r = refit_fit_ab(xs, dl, fill(0.002, length(xs)))
        @test abs(r[1] - 0.31) < 0.003 && abs(r[2] - 0.05) < 0.005
        # transform: two fibers with different zero points -> tilt
        P = (; a = [0.30 0.31 0.32; 0.34 0.36 0.38], b = zeros(2, 3), ea = fill(0.001, 2, 3),
            eb = fill(0.002, 2, 3), n = fill(100, 2, 3), rms = fill(0.02, 2, 3), fibers = [76, 226])
        fn, nf, tl, fl = refit_fpifib_transform(P, "chipzp_stretch_tilt")
        @test nf == 2 && isapprox(fn(76, 2, 1000.0), 0.31) && isapprox(fn(226, 2, 1000.0), 0.36)
        @test isapprox(fn(151, 2, 1000.0), 0.335) && (fl & WAVE_FPIFIB_ROTATION) != 0
        fz, _, _, _ = refit_fpifib_transform(P, "chipzp")
        @test isapprox(fz(1, 3, 10.0), 0.35)
        fs, _, _, _ = refit_fpifib_transform(P, "shift")
        @test isapprox(fs(1, 1, 10.0), mean(P.a))
        P1 = merge(P, (; n = [100 100 100; 3 3 3]))
        f1, nf1, _, fl1 = refit_fpifib_transform(P1, "chipzp_stretch_tilt")
        @test nf1 == 1 && (fl1 & WAVE_FPIFIB_ONE_FIBER) != 0 && isapprox(f1(200, 1, 1024.0), 0.30)
        P0 = merge(P, (; n = zeros(Int, 2, 3)))
        @test refit_fpifib_transform(P0, "chipzp")[4] == WAVE_FPIFIB_UNAVAILABLE
    end

    @testset "pooled sky robust to misassigned OH lines" begin
        rng = MersenneTwister(11)
        rows = (fiber = Int[], chip = Int[], line = Int[], pix = Float64[], lam = Float64[])
        v = Float64[]
        truth = [40.0, -25.0, 10.0]
        for c in 1:3, l in 1:8, f in 1:N_FIBERS
            p = 150.0 + 220.0 * l
            u = (p - 1024.5) / 1024
            ph = (f - 150.5) / 149.5
            x = truth[c] + 30.0 * u + 60.0 * ph + 15.0 * randn(rng)
            (l == 3 && f % 5 == 0) && (x += 1800.0)    # wrong OH component on 20 % of fibers
            (c == 2 && l == 6) && (x -= 2500.0)        # a whole misassigned line
            push!(rows.fiber, f); push!(rows.chip, c); push!(rows.line, 100c + l)
            push!(rows.pix, p); push!(rows.lam, 16000.0); push!(v, x)
        end
        zp, nl, keep = refit_pooled_sky(rows, v)
        @test nl == [8, 7, 8]          # the wholly misassigned chip-G line is rejected
        @test maximum(abs.(zp .- truth .- 30.0 .* mean((150.0 .+ 220.0 .* (1:8) .- 1024.5) ./ 1024))) < 8
        th = refit_sky_terms(rows, v, keep, zp, nl, "chipzp_stretch_tilt")
        @test maximum(abs.(th[1, :] .- truth)) < 8 && all(abs.(th[3, :] .- 60.0) .< 8)
        @test all(abs.(th[2, :] .- 30.0) .< 10)
        thz = refit_sky_terms(rows, v, keep, zp, nl, "chipzp")
        @test thz[1, :] == zp && all(thz[2:3, :] .== 0)
        z, s = refit_fiber_sky_terms(rows, v, keep)
        @test abs(median(filter(isfinite, z[:, 1]))) < 1e-6        # fiber-relative only
        G = fill(16000.0, N_XPIX, N_FIBERS, 3)
        refit_apply_sky_terms!(G, [100.0 0 0; 0 0 0; 0 0 0])
        @test isapprox(G[1, 1, 1], 16000.0 * (1 + 100 / 299792458.0))
        @test G[1, 1, 2] == 16000.0
    end

    @testset "stored coefficients round-trip" begin
        B = refit_coeff_basis(9)
        x = ((1:N_XPIX) .- 1024) ./ 2048
        G = zeros(N_XPIX, N_FIBERS, 3)
        for c in 1:3, f in 1:N_FIBERS
            xt = (c - 2) * 1.07 .+ (1 + 1e-4 * c) .* (x .+ 1e-5 .* x .^ 2)
            G[:, f, c] .= 16200 .+ 1300 .* xt .- 40 .* xt .^ 2 .+ 3 .* xt .^ 3 .- 0.5 .* xt .^ 4 .+ 0.01f
        end
        G[:, 7, 2] .= NaN
        C, mx = refit_grid_to_coeffs(G, B)
        @test mx < 1e-8
        G2 = refit_coeffs_to_grid(C)
        @test all(isnan, G2[:, 7, 2]) && maximum(abs.(filter(isfinite, G2 .- G))) < 1e-8
        dv, zp, sh, rs = refit_vs_pass1(G[65:20:1985, :, :] .* (1 + 50 / 299792458.0), G[65:20:1985, :, :])
        @test all(abs.(zp .- 50) .< 1e-3) && all(sh .< 1e-3)
    end

    @testset "comb shape template" begin
        lin = repeat([16200.0 1300.0 -40.0 3.0 -0.5], N_FIBERS)
        nl = repeat([-1.07 1.0 1.0755 1.0], N_FIBERS)
        Tp = RefitTemplate("apo", "all", zeros(7, 3, N_FIBERS), lin, nl, "x")
        @test refit_template_lambda(Tp, zeros(4), 5, 1, 777.0) == refit_ref_wave_disp(Tp, 5, 1, 777.0)[1]
        l = refit_ref_wave_disp(Tp, 5, 3, 1234.5)[1]
        @test abs(refit_ref_invert(Tp, 5, 3, l, 1200.0) - 1234.5) < 1e-6
        # a chip zero point of +0.5 px moves the template like the reference evaluated 0.5 px later
        @test isapprox(refit_template_lambda(Tp, [0.0, 0.5, 0.0, 0.0], 5, 2, 900.0),
            refit_ref_wave_disp(Tp, 5, 2, 900.5)[1])
        # comb vs arc QA: identical solutions -> no slip, zero velocity difference
        arcG = zeros(97, N_FIBERS, 3)
        for f in 1:N_FIBERS, c in 1:3
            arcG[:, f, c] .= refit_prod_wave(lin, nl, f, c, collect(65.0:20.0:1985.0))
        end
        dfsr, dv = refit_comb_vs_arc(lin, nl, arcG, 3.7363e7)
        @test maximum(abs.(dfsr)) < 1e-9 && maximum(abs.(dv)) < 1e-6
    end
end

using ApogeeReduction: refit_is_morning, refit_order_change_lambda, refit_anchor_rank, refit_anchor_quality,
                       refit_era, WAVE_ANCHOR_FALLBACK, WAVE_ANCHOR_NOT_SAME_BLOCK, WAVE_ANCHOR_NOT_PREFERRED,
                       WAVE_ANCHOR_UNVERIFIED, WAVE_FPS_NO_COMB, WAVE_BIT_NAMES, ARCSEQ_UNGATED, ARCSEQ_POS_OUTLIER,
                       ARCSEQ_FEW_LINES
@testset "arclamp_refit anchors" begin
    @test refit_is_morning("apo", 60000.40) && !refit_is_morning("apo", 60000.08)
    @test refit_is_morning("lco", 60000.30) && !refit_is_morning("lco", 60000.05)
    @test !refit_is_morning("lco", 60000.96) && !refit_is_morning("apo", 60000.99)   # evening before UT midnight
    lin = repeat([16200.0 1300.0 -40.0 3.0 -0.5], N_FIBERS)
    nl = repeat([-1.07 1.0 1.0755 1.0], N_FIBERS)
    L = 3.7363e7
    oc = refit_order_change_lambda(lin, nl, lin, nl, L)
    @test all(oc .== 0)
    # shift the pass-1 solution by one FSR on fiber 9, chip R; plus a half-pixel dither everywhere
    p1 = copy(lin)
    nl2 = copy(nl); nl2[:, 1] .+= 0.5 / 2048; nl2[:, 3] .+= 0.5 / 2048
    lam = refit_prod_wave(lin, nl, 9, 1, [1024.0])[1]
    p1lin = copy(lin)
    oc2 = refit_order_change_lambda(lin, nl, p1lin, nl2, L)
    @test all(oc2 .== 0)                         # a dither step is not an order change
    fsr = lam^2 / (2L)
    p1lin[9, 1] -= fsr
    oc3 = refit_order_change_lambda(lin, nl, p1lin, nl, L)
    @test oc3[9, 1] == 1 && oc3[9, 2] == 1 && count(!=(0), oc3) == 3   # fiber 9 shifted on all chips
    @test cfg_ok = refit_check_config(RefitConfig(neon = "lab_line_offset", plate_shape = "template",
        anchor_block = "any", anchor_pairing = "position"))
    @test refit_check_config(RefitConfig(anchor_block = "morning", anchor_prefer = "none"))
    @test_throws ArgumentError refit_check_config(RefitConfig(plate_shape = "comb"))
    @test_throws ArgumentError refit_check_config(RefitConfig(anchor_block = "evening"))
    c = RefitConfig()
    @test c.plate_shape == "arc" && c.anchor_block == "nearest_science" && c.anchor_prefer == "strict" &&
          c.anchor_pairing == "block" && c.neon == "none"
    @test !hasfield(RefitConfig, :plate_shape_apo) && !hasfield(RefitConfig, :plate_shape_lco)
    # era: canonical plate/FPS boundary (FPS = mjd > get_fps_plate_divide), the same for any FPI light
    @test refit_era("apo", 59393) == "plate" && refit_era("apo", 59558) == "fps"
    @test refit_era("lco", 59808) == "plate" && refit_era("lco", 59809) == "fps"
    @test length(WAVE_BIT_NAMES) == 30 && WAVE_BIT_NAMES[26] == "WAVE_ANCHOR_FALLBACK" &&
          WAVE_BIT_NAMES[30] == "WAVE_FPS_NO_COMB" && WAVE_FPS_NO_COMB == 1 << 29
    # anchor quality tiers
    @test refit_anchor_quality(0) == 0 && refit_anchor_quality(ARCSEQ_UNGATED) == 1 &&
          refit_anchor_quality(ARCSEQ_FEW_LINES) == 1 && refit_anchor_quality(ARCSEQ_GATE_FAIL) == 2 &&
          refit_anchor_quality(ARCSEQ_POS_OUTLIER | ARCSEQ_UNGATED) == 2
    # ranking (times in days): science 0.10-0.30; seq 1 evening (0.05, leading), seq 2 mid-night after
    # an object (0.205), seq 3 morning (0.34, closes the object sequence), seq 4 mid-night before an object
    tseq = [0.05, 0.205, 0.34, 0.149]
    sci = collect(0.10:0.05:0.30)
    q0 = [0, 0, 0, 0]
    r = refit_anchor_rank(tseq, q0, sci, falses(4))
    @test r.trailing == [false, true, true, false]
    @test r.dt_sci_h ≈ [1.2, 0.12, 0.96, 0.024]
    @test r.order == [2, 3, 4, 1]                # preferred (closing) first by distance, then the others
    @test refit_anchor_rank(tseq, q0, sci, falses(4); prefer = "none").order == [4, 2, 3, 1]
    @test refit_anchor_rank(tseq, q0, sci, [false, false, false, true]).order == [4, 2, 3, 1]  # same position
    # a bad preferred block drops behind the verified ones, in the same sensible order
    @test refit_anchor_rank(tseq, [0, 2, 0, 1], sci, falses(4)).order == [3, 1, 4, 2]
    # no science: nothing preferred, every dt Inf, almanac order
    r0 = refit_anchor_rank(tseq, q0, Float64[], falses(4))
    @test all(isinf, r0.dt_sci_h) && !any(r0.pref) && r0.order == [1, 2, 3, 4]
end

using ApogeeReduction: refit_choose_anchor, ArcObjRef, ARCGATE_FIBERS
@testset "arclamp_refit choose anchor" begin
    A = _rf_arc()
    A.g .*= 0.01                          # small time / chip terms: same-group states agree to << 0.25 px
    A.g[A.lay.idx_d[2]] = 0.5 / 2048      # sequence 2 half a pixel off (the other dither position)
    function ref_at(j, e, t)
        d, tt = arcfit_seq_state(A.g, A.lay, j)
        px = collect(880.0:1170.0)
        W = zeros(length(px), length(ARCGATE_FIBERS))
        for (jf, f) in enumerate(ARCGATE_FIBERS)
            W[:, jf] .= arcfit_eval(view(A.TH, f, :), A.g, A.lay, f, 2, px; grpv = A.grp_u[j], d = d, t = tt,
                tau = A.tau_u[j])
        end
        return ArcObjRef(e, t, "fpi", 880, ARCGATE_FIBERS, W)
    end
    # sequences at 60000.1 (evening, leading), .2 (other position), .4 (closes the object sequence)
    combs = [(; ref = ref_at(1, 5, 60000.105), t = 60000.105), (; ref = ref_at(3, 40, 60000.41), t = 60000.41)]
    scis = [(; ref = ref_at(1, 20 + i, t), t) for (i, t) in enumerate((60000.25, 60000.30, 60000.35))]
    ok0 = zeros(Int, 3)
    P = refit_choose_anchor(A, ok0, combs, "apo"; scis)
    @test P.ok && P.j_seq == 3 && P.k_comb == 2 && P.bits == 0 && P.rank == 1 &&
          P.kind == "nearest_science:same_position" && P.dt_sci_h ≈ 1.2
    # the preferred block unverified (ungated): the next verified block in rank anchors (same position, no flag)
    P = refit_choose_anchor(A, [0, 0, ARCSEQ_UNGATED], combs, "apo"; scis)
    @test P.ok && P.j_seq == 1 && P.k_comb == 1 && P.bits == 0 && P.quality == 0
    # first-ranked block without a same-block comb exposure: fallback to the next with one
    P = refit_choose_anchor(A, [0, 0, ARCSEQ_UNGATED], combs[2:2], "apo"; scis)
    @test P.ok && P.j_seq == 3 && P.rank == 3 && P.bits == WAVE_ANCHOR_FALLBACK | WAVE_ANCHOR_UNVERIFIED
    # no science exposures: "any" behaviour, flagged
    P = refit_choose_anchor(A, ok0, combs, "apo")
    @test P.ok && P.kind == "nearest_science:no_science" && P.j_seq == 1 &&
          P.bits == WAVE_ANCHOR_FALLBACK | WAVE_ANCHOR_NOT_PREFERRED
    # morning switch (APO local midnight = 0.29 UT day): the 0.4 sequence; without it, an evening one
    P = refit_choose_anchor(A, ok0, combs, "apo"; scis, anchor_block = "morning")
    @test P.ok && P.j_seq == 3 && P.kind == "morning" && P.bits == 0
    P = refit_choose_anchor(A, [0, 0, ARCSEQ_UNGATED], combs, "apo"; scis, anchor_block = "morning")
    @test P.ok && P.j_seq == 1 && P.kind == "evening" && P.bits == WAVE_ANCHOR_FALLBACK | WAVE_ANCHOR_NOT_PREFERRED
    P = refit_choose_anchor(A, ok0, combs, "apo"; scis, anchor_block = "any")
    @test P.ok && P.j_seq == 1 && P.kind == "any"
    # the only comb exposure sits at the other dither position: that block anchors, flagged
    P = refit_choose_anchor(A, ok0, [(; ref = ref_at(2, 7, 60000.21), t = 60000.21)], "apo"; scis,
        pos_match_px = 0.1)
    @test P.ok && P.j_seq == 2 && (P.bits & WAVE_ANCHOR_NOT_PREFERRED) != 0
end

@testset "arclamp_refit anchor samepos switch" begin
    @test refit_check_config(RefitConfig(anchor_samepos = "none")) && RefitConfig().anchor_samepos == "pass1"
end
