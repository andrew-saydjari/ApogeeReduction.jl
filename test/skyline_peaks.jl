using CSV, DataFrames, JLD2

@testset "skyline_peaks / sky line list" begin
    AR = ApogeeReduction
    repo = dirname(@__DIR__)
    list_file = joinpath(repo, "data", "APOGEE_sky_linelist.csv")
    df = CSV.read(list_file, DataFrame, comment = "#")

    @testset "generated line-list schema" begin
        for col in [:label, :chip, :wavecal_class, :wave_cen_ang, :subwave_1_ang,
            :subwave_2_ang, :sep_ang, :weight_1, :weight_2, :model_rel_int,
            :meas_peak_apo, :meas_peak_lco, :class_note]
            @test col in propertynames(df)
        end
        cls = coalesce.(df.wavecal_class, "")
        @test count(cls .== "B") == 18
        @test count(cls .== "A") == 4
        @test count(cls .== "X") == 1
        # sublines ascending, separation and centroid self-consistent
        @test all(df.subwave_1_ang .< df.subwave_2_ang)
        @test all(abs.(df.sep_ang .- (df.subwave_2_ang .- df.subwave_1_ang)) .< 2e-4)
        @test all(abs.(df.wave_cen_ang .-
                       (df.weight_1 .* df.subwave_1_ang .+ df.weight_2 .* df.subwave_2_ang)) .<
                  2e-4)
        # Lambda-doublet component ratio fixed 1:1
        @test all(df.weight_1 .== 0.5) && all(df.weight_2 .== 0.5)
        @test all(15125.0 .<= df.wave_cen_ang .<= 16960.0)
        @test all(in(["B", "G", "R"]).(df.chip))
        # every classed row carries a provenance note
        @test all(.!ismissing.(df.class_note[cls .!= ""]))

        # the approved class-A guards and class-X exclusion, at their Brooke positions
        a = df[cls .== "A", :]
        @test any(abs.(a.subwave_1_ang .- 16702.0357) .< 1e-3) # OH(5-3)Q2(1.5) f
        @test any(abs.(a.subwave_2_ang .- 16703.2436) .< 1e-3) # OH(5-3)Q2(1.5) e
        @test any(abs.(a.wave_cen_ang .- 15546.14) .< 0.05)    # OH(4-2)R1(3.5)
        @test any(abs.(a.wave_cen_ang .- 16414.74) .< 0.05)    # OH(5-3)R2(2.5)
        @test any(abs.(a.wave_cen_ang .- 16689.20) .< 0.05)    # OH(5-3)Q2(0.5)
        x = df[cls .== "X", :]
        @test nrow(x) == 1 && abs(x.wave_cen_ang[1] - 16442.16) < 0.05

        # class-B chip counts: 7 B / 6 G / 5 R (old selection minus 16442.15)
        b = df[cls .== "B", :]
        @test sum(b.chip .== "B") == 7 && sum(b.chip .== "G") == 6 && sum(b.chip .== "R") == 5
    end

    @testset "class-B centroids match previous production list at mA level" begin
        old = CSV.read(joinpath(repo, "data", "APOGEE_lines.csv"), DataFrame)
        parsevec(s) = parse.(Float64,
            filter(!isempty, split(replace(s[2:(end - 1)], r"\s+" => ","), ",")))
        cls = coalesce.(df.wavecal_class, "")
        for r in eachrow(df[cls .== "B", :])
            d = abs.(old.wav .* 10 .- r.wave_cen_ang)
            i = argmin(d)
            @test d[i] < 0.3
            # compare against the centroid the old pipeline actually consumed
            # (top-2 subline weighted mean), not the raw wav column
            w = parsevec(old.subline_I[i])
            lam = parsevec(old.subline_wav[i])
            si = sortperm(w, rev = true)[1:2]
            w2 = w[si] ./ sum(w[si])
            old_consumed = sum(lam[si] .* w2) * 10
            @test abs(r.wave_cen_ang - old_consumed) < 3e-3
        end
    end

    # a quartic roughwave entry mimicking apo chip R (see data/rough_poly_coeffs.csv)
    rw_quartic = (16725.0346, -0.23544, 16960, 16454,
        16719.943094, -483.275754, -46.582182, 1.197764, 0.022371)
    rw_linear = (16725.0346, -0.23544, 16960, 16454)

    @testset "rough_wave / rough_dispersion" begin
        # quartic evaluation and analytic derivative vs finite differences
        for pix in [100.0, 700.0, 1024.0, 1500.0, 1900.0]
            fd = (AR.rough_wave(pix + 0.5, rw_quartic) - AR.rough_wave(pix - 0.5, rw_quartic))
            @test abs(AR.rough_dispersion(pix, rw_quartic) - fd) < 1e-6
            @test AR.rough_dispersion(pix, rw_quartic) < 0  # wavelength falls with pixel
        end
        @test abs(AR.rough_wave(1024, rw_quartic) - 16719.943094) < 1e-8
        # legacy 4-element entries fall back to the linear model
        @test AR.rough_wave(1200, rw_linear) ==
              AR.rough_linear_wave(1200, a = rw_linear[1], b = rw_linear[2])
        @test AR.rough_dispersion(1200, rw_linear) == rw_linear[2]
    end

    # --- synthetic end-to-end helpers ------------------------------------------
    rwd = Dict("apo" => Dict("R" => rw_quartic))
    grid = [AR.rough_wave(p, rw_quartic) for p in 1:2048]
    function wave2pix(w)
        i = argmin(abs.(grid .- w))
        return i + (w - grid[i]) / AR.rough_dispersion(Float64(i), rw_quartic)
    end
    function synth_flux(lines; sigma = 1.0, amp = 4000.0)
        # lines: vector of (wave, weight) component tuples. The background carries
        # small deterministic pixel-to-pixel structure so the percentile threshold
        # of the peak finder behaves as it does on real data (an exactly constant
        # background makes the threshold a floating-point tie).
        flux = [50.0 + 0.3 * sin(7.7 * p) for p in 1:2048]
        for (w, wt) in lines
            pc = wave2pix(w)
            for p in max(1, Int(floor(pc - 8))):min(2048, Int(ceil(pc + 8)))
                flux[p] += amp * wt * exp(-0.5 * ((p - pc) / sigma)^2)
            end
        end
        return flux
    end
    function mini_df(rows)
        d = DataFrame(rows)
        d.linindx = 1:nrow(d)
        return d
    end
    line_row(label, w1, w2, wt1, wt2, int, cls) = (label = label,
        wave_cen_ang = wt1 * w1 + wt2 * w2, subwave_1_ang = w1, subwave_2_ang = w2,
        weight_1 = wt1, weight_2 = wt2, model_rel_int = int, wavecal_class = cls)

    @testset "class selection + association (synthetic apo R)" begin
        d = mini_df([
            line_row("bright_B", 16690.0, 16690.2, 0.5, 0.5, 1.0, "B"),
            line_row("resolved_B", 16550.0, 16551.2, 0.5, 0.5, 0.5, "B"),
            line_row("guard_A", 16700.0, 16701.2, 0.5, 0.5, 0.05, "A"),
            line_row("excluded_X", 16600.0, 16600.4, 0.5, 0.5, 0.5, "X"),
            line_row("inert", 16630.0, 16630.4, 0.5, 0.5, 0.5, missing)])
        # synthesize ONLY the A and B features (X/inert features absent too)
        comps = vcat([[(r.subwave_1_ang, r.weight_1), (r.subwave_2_ang, r.weight_2)]
                      for r in eachrow(d[1:3, :])]...)
        flux = synth_flux(comps)
        pmat, boff, cnt = AR.get_sky_peaks(flux, "apo", "R", rwd, d)
        fit_inds = sort(Int.(pmat[end, :]))
        @test fit_inds == [1, 2, 3]           # A and B fit; X and inert never selected
        @test abs(boff) < 0.3                 # quartic model needs no offset
        for j in 1:size(pmat, 2)
            li = Int(pmat[end, j])
            # recovered position consistent with the injected wavelength
            @test abs(AR.rough_wave(pmat[1, j], rw_quartic) - pmat[2, j]) < 0.02
            # fixed separation written into the old free-separation slot (row 6)
            expect_sep = (d.subwave_2_ang[li] - d.subwave_1_ang[li]) /
                         abs(AR.rough_dispersion(pmat[4, j], rw_quartic))
            @test abs(pmat[6, j] - expect_sep) < 0.05
        end
    end

    @testset "fixed splitting is robust to LSF-width mismatch (jump mode)" begin
        # Under the old free-separation fit, a +-15% LSF-width error could move
        # unresolved-line centroids by 150-235 mA (sep-width degeneracy ridge).
        # With the separation fixed, the free width absorbs the mismatch and the
        # centroid must not move.
        d = mini_df([line_row("unres_B", 16690.0, 16690.3, 0.5, 0.5, 1.0, "B")])
        for sig in [0.85, 1.0, 1.15]
            flux = synth_flux(
                [(16690.0, 0.5), (16690.3, 0.5)]; sigma = sig)
            pmat, _, cnt = AR.get_sky_peaks(flux, "apo", "R", rwd, d)
            @test cnt == 1
            # < 10 mA recovered-centroid error at every true width
            @test abs(AR.rough_wave(pmat[1, 1], rw_quartic) - pmat[2, 1]) < 0.010
        end
    end

    @testset "weight-pixel pairing fix (asymmetric doublet)" begin
        # blue component (subwave_1) carries weight 0.7 and sits at the HIGHER
        # pixel (dispersion is negative). outpix must be the weighted component
        # position such that rough_wave(outpix) matches outwave; the pre-fix
        # pairing would miss by (w1 - w2) * sep = 0.48 A here.
        w1, w2 = 16640.0, 16641.2
        d = mini_df([line_row("asym", w1, w2, 0.7, 0.3, 1.0, "B")])
        flux = synth_flux([(w1, 0.7), (w2, 0.3)])
        pmat, _, cnt = AR.get_sky_peaks(flux, "apo", "R", rwd, d)
        @test cnt == 1
        @test abs(pmat[2, 1] - (0.7 * w1 + 0.3 * w2)) < 1e-6   # outwave
        @test abs(AR.rough_wave(pmat[1, 1], rw_quartic) - pmat[2, 1]) < 0.05
    end

    # ---- self-calibrating association (stage 1 of the LOO design) ---------
    # 7 class-B lines with irregular spacing like the real chip-R list
    sc_waves = [16502.4, 16553.8, 16620.0, 16692.4, 16708.9, 16800.0, 16903.7]
    sc_df = mini_df([line_row("L$i", w - 0.15, w + 0.15, 0.5, 0.5, 1.0, "B")
                     for (i, w) in enumerate(sc_waves)])
    all_comps(d; skip = "") = vcat([[(r.subwave_1_ang, 0.5), (r.subwave_2_ang, 0.5)]
                                    for r in eachrow(d) if r.label != skip]...)
    function assigned_ok(pmat, d)
        # is every fitted line's recovered pixel consistent with its assigned
        # wavelength under the TRUE (quartic) model?
        nbad = 0
        for j in 1:size(pmat, 2)
            truepix_wav = AR.rough_wave(pmat[1, j], rw_quartic)
            nbad += abs(truepix_wav - pmat[2, j]) > 0.5
        end
        return nbad
    end
    rwd_q = Dict("apo" => Dict("R" => rw_quartic))
    rwd_l = Dict("apo" => Dict("R" => rw_linear))

    @testset "theft reproduced without self-calibration, killed with it" begin
        # victim (L4) feature absent; unlisted thief pair 10.3 A away -- the
        # 16702 -> 16692 mode. Without self-calibration the thief peak takes
        # L4's label; with it, L4 goes unmeasured (correct: NaN, not wrong).
        comps = all_comps(sc_df; skip = "L4")
        push!(comps, (16702.036, 0.15))
        push!(comps, (16703.244, 0.15))
        flux = synth_flux(comps)
        for rwd in (rwd_l, rwd_q)
            pmat0, _, _, qa0 = AR.get_sky_peaks(flux, "apo", "R", rwd, sc_df;
                self_calibrate = false)
            @test assigned_ok(pmat0, sc_df) == 1        # the theft
            pmat1, _, _, qa1 = AR.get_sky_peaks(flux, "apo", "R", rwd, sc_df;
                self_calibrate = true)
            @test assigned_ok(pmat1, sc_df) == 0        # no wrong assignment
            @test !(4 in Int.(pmat1[end, :]))           # L4 unmeasured, not stolen
            @test qa1[2]                                # converged
        end
    end

    @testset "close thief (5.4 A, inside acceptance radius): emission cut" begin
        # the 15546 -> 15540 geometry: the thief is INSIDE the 10 A radius,
        # so re-association alone cannot reject it -- the final LOO
        # consistency cut must
        comps = all_comps(sc_df; skip = "L4")
        push!(comps, (16697.65, 0.3))
        push!(comps, (16697.95, 0.3))
        flux = synth_flux(comps)
        pmat0, _, _, _ = AR.get_sky_peaks(flux, "apo", "R", rwd_q, sc_df;
            self_calibrate = false)
        @test assigned_ok(pmat0, sc_df) == 1
        pmat1, _, _, qa1 = AR.get_sky_peaks(flux, "apo", "R", rwd_q, sc_df;
            self_calibrate = true)
        @test assigned_ok(pmat1, sc_df) == 0
        @test !(4 in Int.(pmat1[end, :]))
        @test qa1[2] && qa1[3] >= 1                     # converged, >=1 pair cut
    end

    @testset "healthy fiber: converges immediately, no drops" begin
        flux = synth_flux(all_comps(sc_df))
        pmat, _, cnt, qa = AR.get_sky_peaks(flux, "apo", "R", rwd_q, sc_df;
            self_calibrate = true)
        @test cnt == 7 && assigned_ok(pmat, sc_df) == 0
        @test qa == (1, true, 0)                        # one pass, converged, no drops
    end

    @testset "iteration cap respected and flagged" begin
        comps = all_comps(sc_df; skip = "L4")
        push!(comps, (16702.036, 0.15))
        push!(comps, (16703.244, 0.15))
        flux = synth_flux(comps)
        # this scenario needs 2 iterations to converge; capping at 1 must
        # return the loud not-converged flag rather than pretending
        _, _, _, qa = AR.get_sky_peaks(flux, "apo", "R", rwd_l, sc_df;
            self_calibrate = true, max_assoc_iter = 1)
        @test qa[1] == 1
        @test !qa[2]
    end

    @testset "seed decomposition: zeropoint is never load-bearing" begin
        # dithers (and re-registrations) move the ZEROPOINT of the true
        # solution while the shape terms are stable (measured: c0 moves
        # ~0.1 A per dither and up to ~10 A across the era; c1/c2 edge
        # effects are 100-1000x smaller). The association must therefore be
        # invariant to zeropoint errors of the SEED, and tolerant of
        # shape-term errors at many times their measured drift.
        flux = synth_flux(all_comps(sc_df))
        ref, _, _, _ = AR.get_sky_peaks(flux, "apo", "R", rwd_q, sc_df)
        for dc0 in [-8.0, -3.0, 3.0, 8.0]
            rw_shift = (rw_quartic[1], rw_quartic[2], rw_quartic[3], rw_quartic[4],
                rw_quartic[5] + dc0, rw_quartic[6], rw_quartic[7],
                rw_quartic[8], rw_quartic[9])
            pmat, boff, cnt, qa = AR.get_sky_peaks(flux, "apo", "R",
                Dict("apo" => Dict("R" => rw_shift)), sc_df)
            @test cnt == size(ref, 2)
            @test sort(Int.(pmat[end, :])) == sort(Int.(ref[end, :]))
            @test maximum(abs.(sort(pmat[1, :]) .- sort(ref[1, :]))) < 0.01
        end
        # shape terms perturbed at ~10x the measured night-to-night drift
        rw_shape = (rw_quartic[1], rw_quartic[2], rw_quartic[3], rw_quartic[4],
            rw_quartic[5], rw_quartic[6] + 2.0, rw_quartic[7] + 0.3,
            rw_quartic[8], rw_quartic[9])
        pmat, _, cnt, qa = AR.get_sky_peaks(flux, "apo", "R",
            Dict("apo" => Dict("R" => rw_shape)), sc_df)
        @test cnt == size(ref, 2)
        @test sort(Int.(pmat[end, :])) == sort(Int.(ref[end, :]))
        @test qa[2]
    end

    @testset "assoc_only fast path matches the full path's associations" begin
        flux = synth_flux(all_comps(sc_df))
        amat, boff, cnt, qa = AR.get_sky_peaks(flux, "apo", "R", rwd_q, sc_df;
            assoc_only = true)
        pmat, _, _, _ = AR.get_sky_peaks(flux, "apo", "R", rwd_q, sc_df)
        @test size(amat, 1) == 3 && cnt == size(pmat, 2)
        @test sort(Int.(amat[3, :])) == sort(Int.(pmat[end, :]))
    end

    @testset "ingest_skyLines_file masks non-class-B lines" begin
        mktempdir() do dir
            nfib = 8
            mat = fill(NaN, 3, 9, nfib)
            for i in 1:3, f in 1:nfib
                mat[i, 1, f] = 500.0 * i    # pixel
                mat[i, 2, f] = 16500.0 + i  # wavelength
            end
            # with line_class: A (code 1) rows masked, B (code 2) kept
            fn = joinpath(dir, "skyLinePeaks_test.h5")
            AR.safe_jldsave(fn; sky_line_mat_clean = mat, line_class = [2, 1, 2],
                no_metadata = true)
            xlst, wlst = AR.ingest_skyLines_file(fn)
            @test all(.!isnan.(xlst[1, :])) && all(.!isnan.(xlst[3, :]))
            @test all(isnan.(xlst[2, :])) && all(isnan.(wlst[2, :]))
            @test all(.!isnan.(wlst[1, :])) && all(.!isnan.(wlst[3, :]))

            # legacy file without line_class: nothing masked
            fn2 = joinpath(dir, "skyLinePeaks_legacy.h5")
            AR.safe_jldsave(fn2; sky_line_mat_clean = mat, no_metadata = true)
            xlst2, wlst2 = AR.ingest_skyLines_file(fn2)
            @test all(.!isnan.(xlst2)) && all(.!isnan.(wlst2))
        end
    end
end
