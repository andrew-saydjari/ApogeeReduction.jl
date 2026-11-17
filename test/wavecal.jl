using HDF5

@testset "wavecal" begin
    @testset "get_fpi_exp_m0s empty in_exp (A5)" begin
        # 3 fibers, 2 exposures. Fiber 1 has peaks in both exposures, fiber 2
        # only in exposure 1, fiber 3 has NO ingested FPI peaks at all.
        # (Pre-fix, exp_m0s was an Int array and the peakless entries were
        # assigned NaN -> InexactError, killing the night's FPI wavecal.)
        n_fnames = 2
        # peak-major layout: rows = peaks, cols = fibers (0 = unused row)
        fpi_line_expInt = [1 1 0
                           1 0 0
                           2 0 0
                           2 0 0]
        peak_ints = Float64[1204 1210 0
                            1203 0 0
                            1205 0 0
                            1202 0 0]

        exp_m0s = ApogeeReduction.get_fpi_exp_m0s(fpi_line_expInt, peak_ints, n_fnames)

        @test size(exp_m0s) == (2, 3)
        @test eltype(exp_m0s) == Float64
        @test exp_m0s[1, 1] == 1203.0
        @test exp_m0s[2, 1] == 1202.0
        @test exp_m0s[1, 2] == 1210.0
        @test isnan(exp_m0s[2, 2])   # fiber with zero peaks in exposure 2
        @test all(isnan.(exp_m0s[:, 3]))  # fiber with zero peaks everywhere

        # downstream median-m0 computation must skip the NaNs and not throw
        med_m0 = round(ApogeeReduction.nanmedian(
            ApogeeReduction.nanmedian(exp_m0s, 2), 1)[1, 1])
        @test med_m0 == round(median([1203.0 + 1210.0, 2 * 1202.0]) / 2)

        # degenerate case: no fiber has any peaks -> all NaN, still no throw
        exp_m0s_empty = ApogeeReduction.get_fpi_exp_m0s(zeros(Int, 4, 3), peak_ints, n_fnames)
        @test all(isnan.(exp_m0s_empty))
        med_empty = ApogeeReduction.nanmedian(
            ApogeeReduction.nanmedian(exp_m0s_empty, 2), 1)[1, 1]
        @test isnan(med_empty)

        # NaN peak_ints (bad peaks) propagate to NaN m0 instead of throwing
        peak_ints_nan = copy(peak_ints)
        peak_ints_nan[1, 2] = NaN
        exp_m0s_nan = ApogeeReduction.get_fpi_exp_m0s(fpi_line_expInt, peak_ints_nan, n_fnames)
        @test isnan(exp_m0s_nan[1, 2])
        @test exp_m0s_nan[1, 1] == 1203.0
    end

    @testset "deterministic night-frame anchor (Method A')" begin
        AR = ApogeeReduction
        rng = Random.MersenneTwister(2026)

        # ABBA synthetic night, labels HONEST per the verified corpus-wide
        # sign convention at LCO (direction +1: smaller DITHPIX <-> smaller
        # measured offset): group at ~0 px labeled 12.995, group at
        # ~+0.478 px labeled 13.496. The anchor default direction is +1.
        n_pairs = 6
        offs_A = 0.0 .+ 0.02 .* randn(rng, n_pairs)
        offs_B = 0.478 .+ 0.02 .* randn(rng, n_pairs)
        offsets = collect(Iterators.flatten(zip(offs_A, offs_B)))
        labels = collect(Iterators.flatten(zip(fill(12.995, n_pairs), fill(13.496, n_pairs))))

        @testset "(a) normal ABBA: smallest-DITHPIX cluster wins, deterministic" begin
            res = AR.night_frame_anchor(offsets, labels)
            @test res.mode == "dithpix_group"
            @test res.canonical_dithpix ≈ 12.995
            @test isapprox(res.shift_pix, 0.0, atol = 0.05)
            @test length(res.cluster_n) == 2
            @test !res.qa_dithpix_frozen_but_motion
            @test !res.qa_dithpix_changed_but_no_motion
            @test !res.qa_dithpix_order_inconsistent
            # every exposure labeled 12.995 is in the canonical cluster
            canon_members = res.exp_cluster .== res.canonical_cluster
            @test all(labels[canon_members] .== 12.995)
            @test all(labels[.!canon_members] .== 13.496)

            # COINCIDENCE property: with honest labels the label rule and the
            # measured-position tie-break select the SAME cluster -- rerunning
            # with the labels withheld must pick the identical membership
            res_nolabel = AR.night_frame_anchor(offsets, fill(NaN, length(offsets)))
            @test res_nolabel.mode == "measured_position"
            @test res_nolabel.canonical_cluster == res.canonical_cluster
            @test res_nolabel.exp_cluster == res.exp_cluster
            @test res_nolabel.shift_pix == res.shift_pix

            # deterministic under 0.001-px perturbations of the measured offsets
            for trial in 1:20
                pert = offsets .+ 0.001 .* randn(rng, length(offsets))
                res_p = AR.night_frame_anchor(pert, labels)
                @test res_p.mode == "dithpix_group"
                @test res_p.canonical_dithpix ≈ 12.995
                @test isapprox(res_p.shift_pix, res.shift_pix, atol = 0.01)
            end

            # deterministic under exposure reordering
            for trial in 1:5
                perm = Random.randperm(rng, length(offsets))
                res_r = AR.night_frame_anchor(offsets[perm], labels[perm])
                @test res_r.mode == "dithpix_group"
                @test res_r.canonical_dithpix ≈ 12.995
                @test res_r.shift_pix == res.shift_pix
            end

            # a NaN offset (failed dither fit) is left unclustered, no crash
            offsets_nan = copy(offsets)
            offsets_nan[3] = NaN
            res_n = AR.night_frame_anchor(offsets_nan, labels)
            @test res_n.exp_cluster[3] == 0
            @test res_n.mode == "dithpix_group"
            @test res_n.canonical_dithpix ≈ 12.995
        end

        @testset "(b) AAAA: single cluster is canonical" begin
            offs1 = 0.13 .+ 0.02 .* randn(rng, 8)
            labs1 = fill(12.995, 8)
            res = AR.night_frame_anchor(offs1, labs1)
            @test res.mode == "single_cluster"
            @test length(res.cluster_n) == 1
            @test isapprox(res.shift_pix, 0.13, atol = 0.05)
            @test !res.qa_dithpix_frozen_but_motion
            @test !res.qa_dithpix_changed_but_no_motion
            @test !res.qa_dithpix_order_inconsistent
        end

        @testset "(c) frozen labels + real motion: measured-position ordering + QA flag" begin
            frozen = fill(13.496, length(offsets))
            res = AR.night_frame_anchor(offsets, frozen)
            @test res.mode == "measured_position"
            @test res.qa_dithpix_frozen_but_motion
            @test !res.qa_dithpix_changed_but_no_motion
            @test !res.qa_dithpix_order_inconsistent
            # direction +1 (LCO): smallest measured position wins (cluster at ~0 px)
            @test res.canonical_cluster == 1
            @test isapprox(res.shift_pix, 0.0, atol = 0.05)
            # direction -1 (APO): the HIGHEST measured position maps to the
            # smallest DITHPIX, so the other cluster wins
            res_apo = AR.night_frame_anchor(offsets, frozen; dithpix_direction = -1)
            @test res_apo.mode == "measured_position"
            @test res_apo.canonical_cluster == 2
            @test isapprox(res_apo.shift_pix, 0.478, atol = 0.05)
            # same behavior when labels are entirely missing, but then the
            # "frozen" QA flag must NOT fire (nothing was frozen -- just absent)
            res_nan = AR.night_frame_anchor(offsets, fill(NaN, length(offsets)))
            @test res_nan.mode == "measured_position"
            @test !res_nan.qa_dithpix_frozen_but_motion
            @test res_nan.canonical_cluster == 1
            # DITHPIX = 0.0 sentinel labels count as missing
            res_zero = AR.night_frame_anchor(offsets, fill(0.0, length(offsets)))
            @test res_zero.mode == "measured_position"
            @test res_zero.canonical_cluster == 1
        end

        @testset "(d) labels ABBA but no motion: single cluster + QA flag" begin
            offs1 = 0.0 .+ 0.02 .* randn(rng, length(labels))
            res = AR.night_frame_anchor(offs1, labels)
            @test res.mode == "single_cluster"
            @test length(res.cluster_n) == 1
            @test res.qa_dithpix_changed_but_no_motion
            @test !res.qa_dithpix_frozen_but_motion
            @test !res.qa_dithpix_order_inconsistent
        end

        @testset "APO direction (-1): honest APO labels coincide too" begin
            # honest APO night per the verified sign convention (direction
            # -1): the SMALLER DITHPIX (12.994) sits on the HIGHER-position
            # cluster (as measured on 133/133 label-discriminating APO nights)
            labels_apo = collect(Iterators.flatten(
                zip(fill(13.499, n_pairs), fill(12.994, n_pairs))))
            res = AR.night_frame_anchor(offsets, labels_apo; dithpix_direction = -1)
            @test res.mode == "dithpix_group"
            @test res.canonical_dithpix ≈ 12.994
            @test res.canonical_cluster == 2
            @test isapprox(res.shift_pix, 0.478, atol = 0.05)
            @test !res.qa_dithpix_order_inconsistent
            # coincidence at APO: withholding the labels picks the same cluster
            res_nl = AR.night_frame_anchor(offsets, fill(NaN, length(offsets));
                dithpix_direction = -1)
            @test res_nl.mode == "measured_position"
            @test res_nl.canonical_cluster == res.canonical_cluster
            @test res_nl.shift_pix == res.shift_pix
        end

        @testset "label ordering inconsistent with measured direction: override + QA flag" begin
            # labels violate the telescope's verified sign convention (here
            # direction +1 with the smaller DITHPIX on the HIGHER-position
            # cluster). Measurement beats label: measured ordering wins and
            # the inconsistency flags.
            labels_swapped = collect(Iterators.flatten(
                zip(fill(13.496, n_pairs), fill(12.995, n_pairs))))
            res = AR.night_frame_anchor(offsets, labels_swapped)
            @test res.mode == "measured_position_override"
            @test res.qa_dithpix_order_inconsistent
            @test !res.qa_dithpix_frozen_but_motion
            @test !res.qa_dithpix_changed_but_no_motion
            # canonical = measured pick under direction +1 (the ~0 px
            # cluster), which here carries the LARGER label
            @test res.canonical_cluster == 1
            @test isapprox(res.shift_pix, 0.0, atol = 0.05)
            @test res.canonical_dithpix ≈ 13.496
            # the same disagreement in the APO direction: honest-LCO-looking
            # labels under direction -1 are an inconsistency there
            res_apo = AR.night_frame_anchor(offsets, labels; dithpix_direction = -1)
            @test res_apo.mode == "measured_position_override"
            @test res_apo.qa_dithpix_order_inconsistent
            @test res_apo.canonical_cluster == 2
        end

        @testset "singleton stray cluster is not a canonical candidate" begin
            offs_stray = vcat(offsets, -0.9)          # one bad fit far away
            labs_stray = vcat(labels, 12.995)
            res = AR.night_frame_anchor(offs_stray, labs_stray)
            @test length(res.cluster_n) == 3
            @test res.cluster_n[res.canonical_cluster] >= 2
            @test res.mode == "dithpix_group"
            @test res.canonical_dithpix ≈ 12.995
            @test isapprox(res.shift_pix, 0.0, atol = 0.05)
            @test !res.qa_dithpix_order_inconsistent
        end

        @testset "(e) pin_night_frame! leaves composed per-exposure solutions invariant" begin
            NXP = ApogeeReduction.N_XPIX
            n_lin_coeffs = 5
            n_fib = ApogeeReduction.N_FIBERS
            med_linParams = zeros(n_fib, n_lin_coeffs)
            # APOGEE-like: ~16000 A zero point, ~580 A/unit-xt slope, small curvature
            for fib in 1:n_fib
                med_linParams[fib, :] .= [16150.0 + 0.1 * fib, 580.0, -20.0, 3.0, -0.5]
            end

            shift_pix = 0.478
            pinned = ApogeeReduction.pin_night_frame!(copy(med_linParams), shift_pix)

            # per-exposure composed solution: wave(xt) = P(off + scale * xt)
            # against the pinned frame the refit dither is (off - s, scale);
            # the composition must be numerically unchanged
            xt = ((1:NXP) .- (NXP ÷ 2)) ./ NXP
            s = shift_pix / NXP
            for (off, scale) in ((2.1e-4, 1.0), (-3.7e-4, 0.999985), (0.0, 1.00002))
                for fib in (1, 150, 300)
                    p_old = evalpoly.(off .+ scale .* xt, Ref(med_linParams[fib, :]))
                    p_new = evalpoly.((off - s) .+ scale .* xt, Ref(pinned[fib, :]))
                    @test maximum(abs.(p_new .- p_old)) < 1e-9
                end
            end

            # offset + SCALE pin: against a frame pinned with (s, c) the refit
            # dither is ((off - s) / c, scale / c); composition unchanged
            c_pin = 1.0 + 1.5e-5
            pinned_sc = ApogeeReduction.pin_night_frame!(copy(med_linParams), shift_pix, c_pin)
            for (off, scale) in ((2.1e-4, 1.0), (-3.7e-4, 0.999985), (0.0, 1.00002))
                for fib in (1, 150, 300)
                    p_old = evalpoly.(off .+ scale .* xt, Ref(med_linParams[fib, :]))
                    p_new = evalpoly.(((off - s) / c_pin) .+ (scale / c_pin) .* xt,
                        Ref(pinned_sc[fib, :]))
                    @test maximum(abs.(p_new .- p_old)) < 1e-9
                end
            end
            # a non-finite scale falls back to the offset-only pin
            @test ApogeeReduction.pin_night_frame!(copy(med_linParams), shift_pix, NaN) == pinned

            # the anchor's scale_pin is the canonical cluster's median scale
            # (1.0 when no scales are given)
            exp_scale = [isodd(i) ? 1.0 + 2e-7 * i : 1.0 + 1.5e-5 for i in eachindex(offsets)]
            res_sc = AR.night_frame_anchor(offsets, labels; exp_scale = exp_scale)
            canon_members = findall(res_sc.exp_cluster .== res_sc.canonical_cluster)
            @test res_sc.scale_pin == median(exp_scale[canon_members])
            @test res_sc.shift_pix == AR.night_frame_anchor(offsets, labels).shift_pix
            @test AR.night_frame_anchor(offsets, labels).scale_pin == 1.0

            # NaN shift (no clusterable exposures) leaves the frame untouched
            unpinned = ApogeeReduction.pin_night_frame!(copy(med_linParams), NaN)
            @test unpinned == med_linParams

            # empty input: mode "none", NaN shift, unit scale, no crash
            res_none = AR.night_frame_anchor(Float64[], Float64[])
            @test res_none.mode == "none"
            @test isnan(res_none.shift_pix)
            @test res_none.scale_pin == 1.0
        end
    end

    @testset "night frame does not inherit the reference exposure's dither scale" begin
        # Two runs of one synthetic ABBA night that differ only by a tiny noise
        # realization in the per-exposure sky solutions -- enough to
        # flip the argmin reference exposure into the OTHER dither cluster. The
        # sky-line peaks are identical. The clusters differ in dither SCALE by
        # 1.5e-5 (the size measured on apo 58500), so a night frame that pinned
        # only the offset would still differ by ~1e-2 A at the chip edges.
        AR = ApogeeReduction
        nfib, nx = AR.N_FIBERS, AR.N_XPIX
        # per-fiber wavelength polynomial in the registered coordinate
        # xt = (pix - 1024) / 2048 (chip G untransformed), ~ -0.28 A/pix
        base = hcat(16000.0 .+ 0.01 .* (1:nfib), fill(-580.0, nfib),
            fill(-3.0, nfib), fill(0.5, nfib))
        nl = repeat([-1.07 1.0 1.0755 1.0], nfib)          # [R0 R1 B0 B1], cporder 1
        chip_tr = (R = (-1.07, 1.0), G = (0.0, 1.0), B = (1.0755, 1.0))
        # true per-exposure dither xt -> d1 + d2 xt: cluster A at 0 px, cluster
        # B at +0.478 px; small fixed within-cluster scale spread
        isB = [false, true, true, false, false, true, true, false]
        n_exp = length(isB)
        d1 = [(b ? 0.478 : 0.0) / nx for b in isB]
        d2 = 1.0 .+ 1.5e-5 .* isB .+ 1e-7 .* [0, 1, -1, 2, -2, 0, 2, 1]
        compose(lin, a, b) = AR.Polynomial(lin)(AR.Polynomial([a, b])).coeffs
        true_lin = [permutedims(hcat([compose(base[f, :], d1[k], d2[k]) for f in 1:nfib]...))
                    for k in 1:n_exp]
        # a run's per-exposure sky solutions: the truth with a tiny per-fiber
        # position jitter (sigma 5e-4 px, below the ~1e-3 px seen between real
        # runs); the jitter realization is all that differs between two runs
        function noisy_lin(seed)
            r = MersenneTwister(seed)
            return [permutedims(hcat([compose(base[f, :], d1[k] + 5e-4 / nx * randn(r), d2[k])
                                      for f in 1:nfib]...)) for k in 1:n_exp]
        end
        # sky-line peaks: the same pixel positions in every fiber (passes the
        # cross-fiber scatter cut), wavelengths from each exposure's TRUE solution
        xpix = collect(range(30.0, 2018.0, length = 40))
        function write_night(dir, lins)
            fnames = String[]
            for k in 1:n_exp
                tag = "apo_60000_$(lpad(k, 8, '0'))"
                wf = joinpath(dir, "waveCalSkyLine_$(tag)_object.h5")
                h5open(wf, "w") do f
                    write(f, "linParams", lins[k])
                    write(f, "nlParams", nl)
                end
                for (chip, (a, b)) in pairs(chip_tr)
                    xt = a .+ b .* ((xpix .- nx ÷ 2) ./ nx)
                    mat = fill(NaN, length(xpix), 2, nfib)
                    for f in 1:nfib
                        mat[:, 1, f] .= xpix
                        mat[:, 2, f] .= AR.Polynomial(base[f, :]).(d1[k] .+ d2[k] .* xt)
                    end
                    h5open(joinpath(dir, "skyLinePeaks_$(tag)_$(chip)_object.h5"), "w") do f
                        write(f, "sky_line_mat_clean", mat)
                    end
                end
                push!(fnames, wf)
            end
            return fnames
        end
        # the reference (argmin) exposure get_ave_night_wave_soln seeds from
        function ref_exposure(lins)
            L = cat(lins..., dims = 3)
            m0 = AR.nanmedian(L, 3)
            off = AR.nanmedian(L .- m0, 2)
            m = AR.nanmedian(L .- off, 3)
            return argmin(vec(sum((L .- m) .^ 2, dims = (1, 2))))
        end
        lins1 = noisy_lin(1)
        ref1 = ref_exposure(lins1)
        seed2 = findfirst(s -> isB[ref_exposure(noisy_lin(s))] != isB[ref1], 2:200) + 1
        lins2 = noisy_lin(seed2)
        @test isB[ref_exposure(lins2)] != isB[ref1]     # the reference really flips

        runs = map((lins1, lins2)) do lins
            AR.get_ave_night_wave_soln(write_night(mktempdir(), lins);
                fit_dither = true, dithpix_direction = 1)
        end
        (lin1, _, wave1, anc1), (lin2, _, wave2, anc2) = runs
        canon = findall(.!isB)       # direction +1: lowest measured position = A
        for anc in (anc1, anc2)
            @test anc.mode == "measured_position"
            @test findall(anc.exp_cluster .== anc.canonical_cluster) == canon
        end
        # the same night frame in both runs: compare the fiber-median wavelength
        # at every pixel of every chip (the per-fiber scatter is the runs'
        # genuine jitter difference; an offset-only pin leaves ~1e-2 A here)
        dwave_fibmed = maximum(abs, AR.nanmedian(wave1 .- wave2, 2))
        println("scale-pin test: ref ", ref1, " -> ", ref_exposure(lins2), " (seed ", seed2,
            "), max |fiber-median dwave| = ", dwave_fibmed, " A, max |dwave| = ",
            maximum(abs, wave1 .- wave2), " A")
        @test dwave_fibmed < 1e-4                   # ~2 m/s
        @test maximum(abs, wave1 .- wave2) < 2e-3   # per-fiber jitter only
        # ... and that frame is the canonical cluster's own (median offset,
        # median scale) solution, whichever exposure seeded the refinement
        target = permutedims(hcat([compose(base[f, :], median(d1[canon]), median(d2[canon]))
                                   for f in 1:nfib]...))
        for L in (lin1, lin2)
            @test abs(median(L[:, 1] .- target[:, 1])) < 1e-5
            @test abs(median(L[:, 2] .- target[:, 2])) < 1e-6
        end
        @test anc1.scale_pin ≈ median(d2[canon]) / d2[ref_exposure(lins1)] rtol = 1e-6
    end

    @testset "interpolate_wave_params: no slit-edge bias, outliers rejected" begin
        # Per-fiber wavelength parameters whose true fiber dependence is CURVED
        # (quadratic in fiber), plus noise, a few gross outlier fibers and one
        # NaN fiber. The pre-fix regulariser fitted linParams[:, 3:end] and all
        # nlParams LINEARLY in fiber, which leaves a systematic bias at the slit
        # edges (fibers 1-30 and 271-300) of ~ curvature * (u^2 - <u^2>).
        rng = MersenneTwister(20261117)
        nfib = ApogeeReduction.N_FIBERS
        fibInds = collect(1:nfib)
        u = (fibInds .- (nfib + 1) / 2) ./ ((nfib - 1) / 2)   # -1 .. 1
        edge = vcat(1:30, (nfib - 29):nfib)
        centre = 101:200
        # (offset, slope, curvature, noise) per column; curvature = 25 noise
        lin_spec = [(16000.0, 0.5, 0.25, 0.01), (130.0, 0.02, 0.05, 0.002),
            (-3.0, 0.002, 0.005, 2e-4), (0.5, 1e-3, 2.5e-3, 1e-4), (-0.2, 5e-4, 1e-3, 4e-5)]
        nl_spec = [(-1.070, 2e-4, 5e-4, 2e-5), (1.0755, -1e-4, 5e-4, 2e-5)]
        truth(s) = s[1] .+ s[2] .* u .+ s[3] .* u .^ 2
        lin_true = reduce(hcat, truth.(lin_spec))
        nl_true = reduce(hcat, truth.(nl_spec))
        lin_noise = [s[4] for s in lin_spec]'
        nl_noise = [s[4] for s in nl_spec]'
        linParams = lin_true .+ lin_noise .* randn(rng, nfib, length(lin_spec))
        nlParams = nl_true .+ nl_noise .* randn(rng, nfib, length(nl_spec))
        outliers = [3, 17, 150, 222, 296]
        for (k, f) in enumerate(outliers)
            sgn = isodd(k) ? 1 : -1
            linParams[f, :] .+= sgn .* 200 .* vec(lin_noise)
            nlParams[f, :] .-= sgn .* 200 .* vec(nl_noise)
        end
        linParams[60, :] .= NaN
        nlParams[60, :] .= NaN

        interp_nl, interp_lin = ApogeeReduction.interpolate_wave_params(
            fibInds, nlParams, linParams)

        @test size(interp_lin) == size(linParams)
        @test size(interp_nl) == size(nlParams)
        @test all(isfinite, interp_lin)
        @test all(isfinite, interp_nl)

        for (P, T, N) in ((interp_lin, lin_true, vec(lin_noise)),
            (interp_nl, nl_true, vec(nl_noise)))
            for j in axes(P, 2)
                r = (P[:, j] .- T[:, j]) ./ N[j]   # in units of the per-fiber noise
                # no slit-edge bias relative to the centre (the linear fit gave
                # ~ 0.55 * curvature = ~14 noise units here)
                @test abs(mean(r[edge]) - mean(r[centre])) < 1.0
                # outlier and NaN fibers are replaced by the smooth solution
                @test maximum(abs.(r[vcat(outliers, 60)])) < 3.0
                # the regularised solution is smooth: no fiber far off the truth
                @test maximum(abs.(r)) < 5.0
            end
        end
    end

    @testset "interpolate_wave_params: per-MTP-block steps, joint across parameters" begin
        # Curved fiber dependence as above, plus a constant step per 30-fiber MTP
        # block that is shared (with fixed per-parameter signs, i.e. correlated)
        # by every non-constant parameter, plus noise, gross outliers, a NaN fiber
        # and one block with only 4 finite fibers. A smooth-only regulariser leaves
        # each block off by ~ its step; the block steps must recover them, and the
        # nearly empty block must fall back to the smooth solution.
        rng = MersenneTwister(20261118)
        nfib = ApogeeReduction.N_FIBERS
        fibInds = collect(1:nfib)
        u = (fibInds .- (nfib + 1) / 2) ./ ((nfib - 1) / 2)
        blk = (fibInds .- 1) .÷ 30 .+ 1
        lin_spec = [(16000.0, 0.5, 0.25, 0.01), (130.0, 0.02, 0.05, 0.002),
            (-3.0, 0.002, 0.005, 2e-4), (0.5, 1e-3, 2.5e-3, 1e-4), (-0.2, 5e-4, 1e-3, 4e-5)]
        nl_spec = [(-1.070, 2e-4, 5e-4, 2e-5), (1.0755, -1e-4, 5e-4, 2e-5)]
        truth(s) = s[1] .+ s[2] .* u .+ s[3] .* u .^ 2
        lin_noise = [s[4] for s in lin_spec]
        nl_noise = [s[4] for s in nl_spec]
        # block step amplitudes in units of the per-fiber noise; one pattern for all
        # parameters, with per-parameter signs (correlated, as in the real fits)
        step_units = [6.0, -5.0, 4.0, -6.0, 5.0, -4.0, 6.0, -5.0, 3.0, -6.0]
        lin_sign = [0.0, 1.0, -1.0, 1.0, -1.0]   # the per-MTP constant term has no step here
        nl_sign = [1.0, -1.0]
        lin_true = reduce(hcat, truth.(lin_spec)) .+
                   step_units[blk] .* (lin_sign .* lin_noise)'
        nl_true = reduce(hcat, truth.(nl_spec)) .+ step_units[blk] .* (nl_sign .* nl_noise)'
        linParams = lin_true .+ lin_noise' .* randn(rng, nfib, length(lin_spec))
        nlParams = nl_true .+ nl_noise' .* randn(rng, nfib, length(nl_spec))
        outliers = [3, 17, 150, 222, 296]
        for (k, f) in enumerate(outliers)
            sgn = isodd(k) ? 1 : -1
            linParams[f, :] .+= sgn .* 200 .* lin_noise
            nlParams[f, :] .-= sgn .* 200 .* nl_noise
        end
        linParams[60, :] .= NaN
        nlParams[60, :] .= NaN
        sparse_block = 7                     # fibers 181-210: only 4 finite fibers
        dead = setdiff(findall(blk .== sparse_block), [183, 191, 199, 207])
        linParams[dead, :] .= NaN
        nlParams[dead, :] .= NaN

        interp_nl, interp_lin = ApogeeReduction.interpolate_wave_params(
            fibInds, nlParams, linParams)
        smooth_nl, smooth_lin = ApogeeReduction.interpolate_wave_params(
            fibInds, nlParams, linParams; block_steps = false)
        @test all(isfinite, interp_lin) && all(isfinite, interp_nl)

        # block-mean error in noise units, per non-constant parameter and block
        function block_errors(P_nl, P_lin)
            E = Float64[]
            for (P, T, N, cols) in ((P_nl, nl_true, nl_noise, 1:2), (P_lin, lin_true, lin_noise, 2:5))
                for j in cols, b in 1:10
                    b == sparse_block && continue
                    push!(E, mean((P[blk .== b, j] .- T[blk .== b, j]) ./ N[j]))
                end
            end
            E
        end
        E_new = block_errors(interp_nl, interp_lin)
        E_smooth = block_errors(smooth_nl, smooth_lin)
        # the steps are recovered: every block mean within 1 noise unit of the truth
        @test maximum(abs.(E_new)) < 1.0
        # and the smooth-only (pre-block-step) regulariser fails that same check
        @test maximum(abs.(E_smooth)) > 3.0
        # the nearly empty block keeps the smooth solution (no step from 4 fibers)
        in_sparse = blk .== sparse_block
        @test interp_nl[in_sparse, :] == smooth_nl[in_sparse, :]
        @test interp_lin[in_sparse, 2:end] == smooth_lin[in_sparse, 2:end]
        # fibers with too few sky lines do not contribute to the steps: a block whose
        # fibers all have 6 lines keeps the smooth solution, the others are unchanged
        n_lines = fill(18, nfib)
        n_lines[blk .== 2] .= 6
        lines_nl, lines_lin = ApogeeReduction.interpolate_wave_params(
            fibInds, nlParams, linParams; n_lines = n_lines)
        @test lines_nl[blk .== 2, :] == smooth_nl[blk .== 2, :]
        @test lines_lin[blk .== 2, 2:end] == smooth_lin[blk .== 2, 2:end]
        @test lines_nl[blk .!= 2, :] == interp_nl[blk .!= 2, :]
        # outlier and NaN fibers get their block's solution, not their own values
        for (P, T, N, cols) in ((interp_nl, nl_true, nl_noise, 1:2),
            (interp_lin, lin_true, lin_noise, 2:5))
            for j in cols
                r = (P[:, j] .- T[:, j]) ./ N[j]
                @test maximum(abs.(r[vcat(outliers, 60)])) < 3.0
                @test maximum(abs.(r[.!in_sparse])) < 5.0
            end
        end
        # the block step itself: joint fiber selection, nothing below min_fibers
        R = randn(rng, 30, 3) .+ [1.0 -2.0 3.0]
        R[5, 2] = 50.0                       # outlier in ONE parameter drops the whole fiber
        R[9, :] .= NaN
        st = ApogeeReduction.mtp_block_steps(R)
        kept = setdiff(1:30, [5, 9])
        @test st ≈ vec(mean(R[kept, :], dims = 1)) atol = 1e-12
        @test isnothing(ApogeeReduction.mtp_block_steps(R[1:9, :]; min_fibers = 10))
    end
end

# ---------------------------------------------------------------------------
# FPI night-average grid (comb_exp_get_and_save_fpi_wavecal).
#
# The FPI night grid -- saved to waveCalFPI, returned, and stored by
# fpi_medwavecal_skyline_dither as `fpi/nightAve_wave_soln` -- must be the
# fitted night solution at the reference (first-exposure) dither. It used to be
# the LAST FPI exposure's dithered grid, because the per-exposure loop that
# writes the waveCalNightfpiDither files overwrote the same buffer.

# Independent evaluation of the wavecal model (not using the package helpers):
# x = (pix - 1024)/2048; chip G untransformed, R/B registered by nlParams
# [R0, R1, B0, B1] (cporder 1) or [R0, B0] (cporder 0); optional dither poly
# on the registered x; lambda = sum lin[k] xt^(k-1).
function _fpi_test_grid(lin, nl; dith = nothing)
    nf = size(lin, 1)
    cp = size(nl, 2) ÷ 2 - 1
    x = ((1:2048) .- 1024) ./ 2048
    W = zeros(2048, nf, 3)
    for f in 1:nf, c in 1:3
        xt = if c == 2
            collect(x)
        else
            p = c == 1 ? nl[f, 1:(cp + 1)] : nl[f, (cp + 2):end]
            cp == 0 ? p[1] .+ x : p[1] .+ p[2] .* x
        end
        if !isnothing(dith)
            d = dith[f, :]
            xt = sum(d[k] .* xt .^ (k - 1) for k in eachindex(d))
        end
        W[:, f, c] .= evalpoly.(xt, Ref(lin[f, :]))
    end
    return W
end

# Write fpiPeaks_{R,G,B} files for one exposure whose FPI lines (lambda_m =
# 2 cav / (m + m_off)) fall where `lin`/`nl` put them, shifted by `shift_pix`.
function _write_fpipeaks_model(fname_R, lin, nl; shift_pix = 0.0,
        cav = 3.7363125e7, m_off = 0.225)
    W = _fpi_test_grid(lin, nl)
    for (c, chip) in enumerate(("R", "G", "B"))
        cols = Vector{Vector{Float64}}(undef, size(lin, 1))
        for f in axes(lin, 1)
            lam = W[:, f, c]               # decreasing with pixel
            mlo = ceil(Int, 2cav / lam[1] - m_off)
            mhi = floor(Int, 2cav / lam[end] - m_off)
            xs = Float64[]
            for m in mlo:mhi
                lm = 2cav / (m + m_off)
                i = findlast(>=(lm), lam)
                (isnothing(i) || i >= 2048) && continue
                push!(xs, i + (lam[i] - lm) / (lam[i] - lam[i + 1]) + shift_pix)
            end
            cols[f] = filter(x -> 3 <= x <= 2046, xs)
        end
        n_peaks = maximum(length, cols)
        mat = fill(NaN, n_peaks, 4, size(lin, 1))
        cov = fill(NaN, n_peaks, 4, 4, size(lin, 1))
        for f in axes(lin, 1), (p, xp) in enumerate(cols[f])
            mat[p, :, f] .= (800.0, xp, 1.25, 9.0)
            cov[p, 2, 2, f] = 1e-4
        end
        h5open(replace(fname_R, "_R_" => "_$(chip)_"), "w") do fh
            write(fh, "fpi_line_mat", mat)
            write(fh, "fpi_line_cov_mat", cov)
        end
    end
    return fname_R
end

@testset "FPI night-average grid is the reference-dither solution" begin
    dir = mktempdir()
    nf = ApogeeReduction.N_FIBERS
    # seed (sky-like, cporder 0) and truth (FPI, cporder 1), typical APO values
    lin = repeat([16150.35 -579.42 -43.24 1.05 -0.015], nf)
    nl_sky = repeat([-1.0702 1.0754], nf)
    nl_true = repeat([-1.0702 1.0 1.0754 1.0], nf)
    # exposure 1 = reference; exposure 2 is dithered by 0.5 pix
    f1 = _write_fpipeaks_model(joinpath(dir, "fpiPeaks_apo_60000_0009_R_arclamp.h5"),
        lin, nl_true; shift_pix = 0.0)
    f2 = _write_fpipeaks_model(joinpath(dir, "fpiPeaks_apo_60000_0010_R_arclamp.h5"),
        lin, nl_true; shift_pix = 0.5)

    out = ApogeeReduction.comb_exp_get_and_save_fpi_wavecal(
        [f1, f2], lin, nl_sky; cporder = 1, wporder = 4, dporder = 2,
        n_sigma = 4, max_ang_sigma = 0.2, max_iter = 1, verbose = false,
        checkpoint_mode = "clobber")
    @test !isnothing(out)
    outname, night_lin, night_nl, night_grid = out

    saved_grid = h5open(f -> read(f["chipWaveSoln"]), outname)
    dith1 = joinpath(dir, "waveCalNightfpiDither_apo_60000_0009_arclamp.h5")
    dith2 = joinpath(dir, "waveCalNightfpiDither_apo_60000_0010_arclamp.h5")
    exp1_grid, exp2_grid, exp2_dith = h5open(f -> read(f["chipWaveSoln"]), dith1),
    h5open(f -> read(f["chipWaveSoln"]), dith2), h5open(f -> read(f["ditherParams"]), dith2)

    ref_grid = _fpi_test_grid(night_lin, night_nl)
    # what is returned (and so stored as nightAve_wave_soln) is what is saved ...
    @test saved_grid == night_grid
    # ... and it is the fitted night solution at the reference dither
    @test maximum(abs, night_grid .- ref_grid) < 1e-8
    @test maximum(abs, night_grid .- exp1_grid) < 1e-8
    # the dithered exposure's own grid is still its parameters + its dither ...
    @test maximum(abs, exp2_grid .- _fpi_test_grid(night_lin, night_nl; dith = exp2_dith)) < 1e-8
    # ... it really is dithered (~0.5 pix x ~0.28 A/pix) ...
    @test 0.05 < ApogeeReduction.nanmedian(abs.(exp2_grid .- ref_grid)) < 0.3
    # ... and the night grid is NOT the last exposure's dithered grid (pre-fix bug)
    @test ApogeeReduction.nanmedian(abs.(night_grid .- exp2_grid)) > 0.05
end
