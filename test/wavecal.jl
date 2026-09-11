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

        # a fixed synthetic fiducial: wavelength of the frame shifted by s px,
        # with an APOGEE-like positive dispersion (0.284 A/px at chip G center)
        fid_pos(s) = 16151.8 + 0.284 * s
        fid_neg(s) = 16151.8 - 0.284 * s

        # ABBA synthetic night: group at ~0 px labeled 13.496, group at
        # ~+0.478 px labeled 12.995 (smallest DITHPIX -> canonical)
        n_pairs = 6
        offs_A = 0.0 .+ 0.02 .* randn(rng, n_pairs)
        offs_B = 0.478 .+ 0.02 .* randn(rng, n_pairs)
        offsets = collect(Iterators.flatten(zip(offs_A, offs_B)))
        labels = collect(Iterators.flatten(zip(fill(13.496, n_pairs), fill(12.995, n_pairs))))

        @testset "(a) normal ABBA: smallest-DITHPIX cluster wins, deterministic" begin
            res = AR.night_frame_anchor(offsets, labels, fid_pos)
            @test res.mode == "dithpix_group"
            @test res.canonical_dithpix ≈ 12.995
            @test isapprox(res.shift_pix, 0.478, atol = 0.05)
            @test length(res.cluster_n) == 2
            @test !res.qa_dithpix_frozen_but_motion
            @test !res.qa_dithpix_changed_but_no_motion
            # every exposure labeled 12.995 is in the canonical cluster
            canon_members = res.exp_cluster .== res.canonical_cluster
            @test all(labels[canon_members] .== 12.995)
            @test all(labels[.!canon_members] .== 13.496)

            # deterministic under 0.001-px perturbations of the measured offsets
            for trial in 1:20
                pert = offsets .+ 0.001 .* randn(rng, length(offsets))
                res_p = AR.night_frame_anchor(pert, labels, fid_pos)
                @test res_p.mode == "dithpix_group"
                @test res_p.canonical_dithpix ≈ 12.995
                @test isapprox(res_p.shift_pix, res.shift_pix, atol = 0.01)
            end

            # deterministic under exposure reordering
            for trial in 1:5
                perm = Random.randperm(rng, length(offsets))
                res_r = AR.night_frame_anchor(offsets[perm], labels[perm], fid_pos)
                @test res_r.mode == "dithpix_group"
                @test res_r.canonical_dithpix ≈ 12.995
                @test res_r.shift_pix == res.shift_pix
            end

            # a NaN offset (failed dither fit) is left unclustered, no crash
            offsets_nan = copy(offsets)
            offsets_nan[3] = NaN
            res_n = AR.night_frame_anchor(offsets_nan, labels, fid_pos)
            @test res_n.exp_cluster[3] == 0
            @test res_n.mode == "dithpix_group"
            @test res_n.canonical_dithpix ≈ 12.995
        end

        @testset "(b) AAAA: single cluster is canonical" begin
            offs1 = 0.13 .+ 0.02 .* randn(rng, 8)
            labs1 = fill(12.995, 8)
            res = AR.night_frame_anchor(offs1, labs1, fid_pos)
            @test res.mode == "single_cluster"
            @test length(res.cluster_n) == 1
            @test isapprox(res.shift_pix, 0.13, atol = 0.05)
            @test !res.qa_dithpix_frozen_but_motion
            @test !res.qa_dithpix_changed_but_no_motion
        end

        @testset "(c) frozen labels + real motion: absolute tie-break + QA flag" begin
            frozen = fill(13.496, length(offsets))
            res = AR.night_frame_anchor(offsets, frozen, fid_pos)
            @test res.mode == "fiducial_tiebreak"
            @test res.qa_dithpix_frozen_but_motion
            @test !res.qa_dithpix_changed_but_no_motion
            # smallest fiducial wavelength wins: with positive dispersion the
            # cluster at smaller offset; with negative dispersion the other one
            @test isapprox(res.shift_pix, 0.0, atol = 0.05)
            res_neg = AR.night_frame_anchor(offsets, frozen, fid_neg)
            @test res_neg.mode == "fiducial_tiebreak"
            @test isapprox(res_neg.shift_pix, 0.478, atol = 0.05)
            # same behavior when labels are entirely missing, but then the
            # "frozen" QA flag must NOT fire (nothing was frozen -- just absent)
            res_nan = AR.night_frame_anchor(offsets, fill(NaN, length(offsets)), fid_pos)
            @test res_nan.mode == "fiducial_tiebreak"
            @test !res_nan.qa_dithpix_frozen_but_motion
            # DITHPIX = 0.0 sentinel labels count as missing
            res_zero = AR.night_frame_anchor(offsets, fill(0.0, length(offsets)), fid_pos)
            @test res_zero.mode == "fiducial_tiebreak"
        end

        @testset "(d) labels ABBA but no motion: single cluster + QA flag" begin
            offs1 = 0.0 .+ 0.02 .* randn(rng, length(labels))
            res = AR.night_frame_anchor(offs1, labels, fid_pos)
            @test res.mode == "single_cluster"
            @test length(res.cluster_n) == 1
            @test res.qa_dithpix_changed_but_no_motion
            @test !res.qa_dithpix_frozen_but_motion
        end

        @testset "singleton stray cluster is not a canonical candidate" begin
            offs_stray = vcat(offsets, -0.9)          # one bad fit far away
            labs_stray = vcat(labels, 12.995)
            res = AR.night_frame_anchor(offs_stray, labs_stray, fid_pos)
            @test length(res.cluster_n) == 3
            @test res.cluster_n[res.canonical_cluster] >= 2
            @test res.mode == "dithpix_group"
            @test res.canonical_dithpix ≈ 12.995
            @test isapprox(res.shift_pix, 0.478, atol = 0.05)
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
            med_nlParams = zeros(n_fib, 2)
            med_nlParams[:, 1] .= -1.07   # chip R offset
            med_nlParams[:, 2] .= 1.07    # chip B offset

            shift_pix = 0.478
            pinned = ApogeeReduction.pin_night_frame!(copy(med_linParams), shift_pix)

            # the fiducial wavelength of the pinned frame equals the
            # shifted-fiducial wavelength of the unpinned frame
            w_pred = ApogeeReduction.night_frame_fiducial_wave(
                med_linParams, med_nlParams, shift_pix)
            w_pinned = ApogeeReduction.night_frame_fiducial_wave(
                pinned, med_nlParams, 0.0)
            @test isapprox(w_pred, w_pinned, atol = 1e-10)

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

            # NaN shift (no clusterable exposures) leaves the frame untouched
            unpinned = ApogeeReduction.pin_night_frame!(copy(med_linParams), NaN)
            @test unpinned == med_linParams

            # empty input: mode "none", NaN shift, no crash
            res_none = AR.night_frame_anchor(Float64[], Float64[], fid_pos)
            @test res_none.mode == "none"
            @test isnan(res_none.shift_pix)
        end
    end
end
