using ApogeeReduction: load_arclamp_linelist, ARCLAMP_LINELIST_DEFAULT, get_arclamp_peaks_local,
                       arclamp_associate, arclamp_rough_lambda, arclamp_seed_dispersion,
                       ArcGlobalLayout, ArcJointData, arcfit_joint, arcfit_eval, arcfit_seq_state,
                       arclamp_position_clusters, fit_arclamp_night, arclamp_night_grid,
                       arclamp_sequence_reject, ARCSEQ_REJECTED, ARCSEQ_QUALITY, ARCSEQ_SHAPE,
                       ARCSEQ_NOT_COMPARED, ARCNIGHT_CHIPSCALE_FREE, ARCNIGHT_SEQ_REJECTED,
                       ARCNIGHT_FEW_SEQ, ARCLAMP_SEQFLAG_BITS, ARCLAMP_NIGHTFLAG_BITS,
                       ARCREJ_DEFAULTS, ARCSEQ_CHIPOFF, ARCNIGHT_ABS_FAIL,
                       bad_fully_saturated

# --- a synthetic instrument built from the frozen line list's own seed --------------------
# true per-fiber model: chip transform (a_R, b_R, a_B, b_B) + one quartic in z; a sequence s
# adds a shift d_s (x units) to z.
const _AW_CT = [-1.07, 1.0, 1.0755, 1.0]
_aw_xt(ct, c, x) = c == 2 ? x : (c == 1 ? ct[1] + ct[2] * x : ct[3] + ct[4] * x)
_aw_poly(p, z) = p[1] + z * (p[2] + z * (p[3] + z * (p[4] + z * p[5])))
_aw_dpoly(p, z) = p[2] + z * (2p[3] + z * (3p[4] + z * 4p[5]))
function _aw_truth_poly(ll)
    zs = Float64[]
    ls = Float64[]
    for c in 1:3, pix in 40.0:40.0:2008.0
        x = (pix - 1024) / 2048
        push!(zs, _aw_xt(_AW_CT, c, x))
        push!(ls, arclamp_rough_lambda(ll, c, [pix])[1])
    end
    A = hcat([zs .^ j for j in 0:4]...)
    return A \ ls
end
"pixel at which fiber model (ct, p) with total shift `sh` reaches wavelength w on chip c"
function _aw_pix(ct, p, c, w, sh, ll)
    l0 = arclamp_rough_lambda(ll, c, collect(1.0:2048.0))
    pix = argmin(abs.(l0 .- w)) * 1.0
    x = (pix - 1024) / 2048
    for _ in 1:30
        z = _aw_xt(ct, c, x) + sh
        dzdx = c == 2 ? 1.0 : (c == 1 ? ct[2] : ct[4])
        x -= (_aw_poly(p, z) - w) / (_aw_dpoly(p, z) * dzdx)
    end
    return x * 2048 + 1024
end
"line table of one sequence (ids of the fit selection) for fibers `fibs` at shift `sh`"
function _aw_lines(ll, p, fibs, sh; ct = _AW_CT, phi_tilt = 0.0)
    rows = (fiber = Int16[], chip = Int8[], lamp = Int8[], id = Int16[], pix = Float32[],
        pixvar = Float32[], height = Float32[])
    df = ll.df
    for f in fibs, i in 1:size(df, 1)
        id = df.lab_id[i]
        (id > 0 && id in ll.fit_ids) || continue
        c = df.chip_idx[i]
        s = sh + phi_tilt * (f - 150.5) / 150.0
        px = _aw_pix(ct, p, c, ll.labwave[id], s, ll)
        (20 <= px <= 2028) || continue
        push!(rows.fiber, f)
        push!(rows.chip, c)
        push!(rows.lamp, df.lamp[i] == "ThAr" ? 1 : 2)
        push!(rows.id, id)
        push!(rows.pix, px)
        push!(rows.pixvar, 1e-4)
        push!(rows.height, 100.0)
    end
    return rows
end

@testset "arclamp_wavecal" begin
    @testset "frozen line list" begin
        ll = load_arclamp_linelist()
        @test isfile(ARCLAMP_LINELIST_DEFAULT)
        @test length(ll.fit_ids) == 93
        @test length(ll.fingerprint) == 16
        df = ll.df
        # the default selection excludes neon and the DRP re-derived (blended) lines
        @test !any(df.use_default .& df.neon)
        @test !any(df.use_default .& df.drp_rederived)
        @test all(df.store[df.use_default])
        @test count(df.drp_rederived) == 3
        # switching on neon / blended lines is a one-line selection change
        lln = load_arclamp_linelist(; selection = "use_with_neon")
        @test length(lln.fit_ids) > length(ll.fit_ids)
        llb = load_arclamp_linelist(; selection = "use_with_blended")
        @test length(llb.fit_ids) == length(ll.fit_ids) + 3
        # a per-lamp neon offset moves only the neon wavelengths
        llo = load_arclamp_linelist(; selection = "use_with_neon",
            neon_offset_A = Dict("ThAr" => 0.021, "UNe" => 0.006))
        nid = df.lab_id[df.neon .& (df.lamp .== "ThAr") .& df.store]
        tid = df.lab_id[df.use_default]
        @test all(llo.labwave[i] - lln.labwave[i] ≈ 0.021 for i in nid)
        @test all(llo.labwave[i] == lln.labwave[i] for i in tid)
        @test occursin("neon_offset", llo.selection)
        @test_throws ArgumentError load_arclamp_linelist(; selection = "no_such_column")
        # seed dispersion is data-derived and close to the known APOGEE values
        d = arclamp_seed_dispersion(ll)
        @test all(abs.(collect(d) .- [0.2359, 0.2828, 0.3265]) .< 0.02)
    end

    @testset "local-noise peak finder" begin
        rng = MersenneTwister(11)
        n = 2048
        centers = sort(100 .+ 1800 .* rand(rng, 40))
        centers = centers[[true; diff(centers) .> 8]]
        flux = fill(50.0, n)
        x = collect(1.0:n)
        for c in centers
            flux .+= 2000 .* exp.(-0.5 .* ((x .- c) ./ 1.1) .^ 2)
        end
        sig = sqrt.(flux)
        flux .+= sig .* randn(rng, n)
        ivar = 1 ./ sig .^ 2
        mask = zeros(Int, n)
        _, prm, cov, info = get_arclamp_peaks_local(flux, ivar, mask; nsig = 8.0)
        @test size(prm, 1) == length(centers)
        @test maximum(abs.(prm[:, 2] .- centers)) < 0.05
        @test all(prm[:, 3] .> 0.9) && all(prm[:, 3] .< 1.3)
        # a fully saturated line is flagged, not fitted
        k = round(Int, centers[5])
        mask2 = copy(mask)
        mask2[(k - 1):(k + 1)] .|= bad_fully_saturated
        _, prm2, _, info2 = get_arclamp_peaks_local(flux, ivar, mask2; nsig = 8.0)
        @test size(prm2, 1) == length(centers) - 1
        @test length(info2.sat) == 1
        # pure noise gives (almost) nothing at 8 sigma
        noise = 50 .+ sqrt(50) .* randn(rng, n)
        _, prm3, _, _ = get_arclamp_peaks_local(noise, fill(1 / 50, n), mask; nsig = 8.0)
        @test size(prm3, 1) <= 1
    end

    @testset "association" begin
        ll = load_arclamp_linelist()
        cs = ll.slices[(2, 3)]                     # UNe, chip B
        # pixels of the catalog lines under the seed shifted by 4 A (an instrument move)
        lam_all = arclamp_rough_lambda(ll, 3, collect(1.0:2048.0))
        sel = findall(cs.store)[1:2:end]
        pix = [argmin(abs.(lam_all .+ 4.0 .- cs.wave[j])) + 0.0 for j in sel]
        lam0 = arclamp_rough_lambda(ll, 3, pix)
        a = arclamp_associate(lam0, (pix .- 1024) ./ 2048, cs.wave; tol = 1.0)
        @test count(a .== sel) >= 0.9 * length(sel)
        @test count((a .!= 0) .& (a .!= sel)) == 0
    end

    @testset "positions and arc-only rejection rule" begin
        posout, cl, ncl = arclamp_position_clusters([0.0, 0.51, 0.02, 0.50, 29.5])
        @test posout == [false, false, false, false, true]
        @test ncl == 2
        @test cl[1] == cl[3] && cl[2] == cl[4] && cl[1] != cl[2]
        st(r, k, sh, co = 0.0; n = 20000) = (n_rows = n, n_keep = isfinite(k) ? round(Int, k * n) : 0, keep_frac = k,
            rsig_px = r, shape_px = sh, chipoff_px = co, chip_zp_px = zeros(3),
            chip_stretch_px = zeros(3))
        good = [st(0.050, 0.997, 0.002), st(0.051, 0.996, 0.003), st(0.052, 0.997, 0.002)]
        q, sh, co, ab, sc, cmp = arclamp_sequence_reject(good)
        @test cmp && !any(q) && !any(sh) && !any(co) && !any(ab) && all(sc .< 1)
        # a sequence with a much larger line scatter, one with most lines clipped, one whose
        # residual pattern differs from the night's, one with a chip-differential offset, one
        # with no kept rows
        bad = vcat(good, [st(0.090, 0.997, 0.002), st(0.050, 0.60, 0.002), st(0.050, 0.997, 0.060),
            st(0.050, 0.997, 0.002, 0.08), st(NaN, NaN, NaN; n = 0)])
        q, sh, co, ab, sc, cmp = arclamp_sequence_reject(bad)
        @test q == [false, false, false, true, true, false, false, true]
        @test sh == [false, false, false, false, false, true, false, false]
        @test co == [false, false, false, false, false, false, true, false]
        @test ab[5] && ab[8] && !any(ab[[1, 2, 3, 6, 7]])
        @test all(sc[4:8] .> 1) && all(sc[1:3] .< 1)
        # differences below the absolute excess floors are never rejected, however tight the night
        tight = [st(0.0500, 0.997, 0.002), st(0.0500, 0.997, 0.002),
            st(0.0500 + 0.5 * ARCREJ_DEFAULTS.rsig_floor_px, 0.997, 0.002 + 0.5 * ARCREJ_DEFAULTS.shape_floor_px)]
        q, sh, co, _, _, _ = arclamp_sequence_reject(tight)
        @test !any(q) && !any(sh) && !any(co)
        # nights with fewer than min_seq_compare sequences: no relative test, absolute test only
        q, sh, co, ab, sc, cmp = arclamp_sequence_reject(bad[[1, 5]])
        @test !cmp && !any(q) && !any(sh) && !any(co) && ab == [false, true]
        @test length(ARCLAMP_SEQFLAG_BITS) == 10 && length(ARCLAMP_NIGHTFLAG_BITS) == 10
        @test ARCSEQ_REJECTED == ARCSEQ_QUALITY | ARCSEQ_SHAPE | ARCSEQ_CHIPOFF
        # retired gate-B bit positions are not reused
        @test ARCSEQ_REJECTED & 0x0006 == 0 && ARCNIGHT_SEQ_REJECTED & 0x0018 == 0
    end

    @testset "joint night fit: synthetic recovery, rejection and fixed chip scales" begin
        ll = load_arclamp_linelist()
        p = _aw_truth_poly(ll)
        fibs = 120:149
        # sequences: s1 position 1 (evening), s2 position 2 (+0.5 px), s3 position 1 but
        # displaced by 0.05 px (an unrecorded settling error, like the LCO evening block),
        # s4 position 1 (morning). Every sequence has its own measured shift, so the displaced
        # sequence is absorbed: it is NOT rejected and does not bias the solution.
        shifts_px = [0.0, 0.5, 0.05, 0.0]
        lines = [_aw_lines(ll, p, fibs, s / 2048) for s in shifts_px]
        tseq = [60000.10, 60000.11, 60000.12, 60000.40]
        R = fit_arclamp_night(lines, tseq, "apo", ll)
        @test R.n_clusters == 2
        @test R.dpx1[2] ≈ 0.5 atol = 1e-4
        @test R.dpx1[3] ≈ 0.05 atol = 1e-4
        @test R.compared
        @test all((R.seqflags .& ARCSEQ_REJECTED) .== 0)
        @test all((R.seqflags .& ARCSEQ_NOT_COMPARED) .== 0)
        @test R.used_seq == [1, 2, 3, 4]
        @test eltype(R.seqflags) == UInt16 && R.night_flags isa UInt16
        @test (R.night_flags & ARCNIGHT_CHIPSCALE_FREE) != 0
        @test (R.night_flags & ARCNIGHT_SEQ_REJECTED) == 0
        @test count(R.okf) == length(fibs)
        # the fitted model reproduces the truth at the reference state
        G = arclamp_night_grid(R)
        px = collect(100.0:100.0:1900.0)
        err = Float64[]
        for f in fibs, c in 1:3
            tr = [_aw_poly(p, _aw_xt(_AW_CT, c, (q - 1024) / 2048)) for q in px]
            append!(err, abs.(G[Int.(px), f, c] .- tr))
        end
        @test maximum(err) < 1e-4                 # Å (≈ 2 m/s)

        # a sequence whose lines carry a shape the others do not (a curvature of 0.1 px across
        # each chip) is rejected and refitted out; the solution still reproduces the truth
        linesb = deepcopy(lines)
        linesb[3] = merge(linesb[3], (pix = Float32[q + 0.1 * (6 * ((q - 1024) / 2048)^2 - 0.5)
                                                    for q in linesb[3].pix],))
        Rb = fit_arclamp_night(linesb, tseq, "apo", ll)
        @test (Rb.seqflags[3] & ARCSEQ_REJECTED) != 0
        @test all((Rb.seqflags[[1, 2, 4]] .& ARCSEQ_REJECTED) .== 0)
        @test Rb.used_seq == [1, 2, 4]
        @test (Rb.night_flags & ARCNIGHT_SEQ_REJECTED) != 0
        Gb = arclamp_night_grid(Rb)
        errb = [abs(Gb[Int(q), f, c] - _aw_poly(p, _aw_xt(_AW_CT, c, (q - 1024) / 2048)))
                for f in fibs, c in 1:3, q in px]
        @test maximum(errb) < 1e-4

        # fixed chip scales (the comb-measured b_R, b_B) are held exactly
        ct2 = [-1.07, 1.0 - 2.6e-4, 1.0755, 1.0 + 1.1e-4]
        lines2 = [_aw_lines(ll, p, fibs, s / 2048; ct = ct2) for s in (0.0, 0.5)]
        scales = repeat([ct2[2] ct2[4]], 300, 1)
        R2 = fit_arclamp_night(lines2, [60000.1, 60000.11], "apo", ll; chip_scales = scales)
        @test all(R2.TH[fibs, 2] .== ct2[2]) && all(R2.TH[fibs, 4] .== ct2[4])
        @test (R2.night_flags & ARCNIGHT_CHIPSCALE_FREE) == 0
        # two sequences: nothing to compare against, flagged and never rejected
        @test !R2.compared
        @test all((R2.seqflags .& ARCSEQ_NOT_COMPARED) .!= 0)
        @test all((R2.seqflags .& ARCSEQ_REJECTED) .== 0)
        @test (R2.night_flags & ARCNIGHT_FEW_SEQ) != 0
        G2 = arclamp_night_grid(R2)
        err2 = [abs(G2[Int(q), f, 1] - _aw_poly(p, _aw_xt(ct2, 1, (q - 1024) / 2048)))
                for f in fibs, q in px]
        @test maximum(err2) < 1e-4
    end
end
