using ApogeeReduction: load_arclamp_linelist, ARCLAMP_LINELIST_DEFAULT, get_arclamp_peaks_local,
                       arclamp_associate, arclamp_rough_lambda, arclamp_seed_dispersion,
                       ArcGlobalLayout, ArcJointData, arcfit_joint, arcfit_eval, arcfit_seq_state,
                       arclamp_position_clusters, arclamp_gate, ArcObjRef, fit_arclamp_night,
                       arclamp_night_grid, ARCSEQ_GATE_FAIL, ARCSEQ_UNGATED, ARCNIGHT_CHIPSCALE_FREE,
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

    @testset "positions and gate" begin
        posout, cl, ncl = arclamp_position_clusters([0.0, 0.51, 0.02, 0.50, 29.5])
        @test posout == [false, false, false, false, true]
        @test ncl == 2
        @test cl[1] == cl[3] && cl[2] == cl[4] && cl[1] != cl[2]
        dev, gated, flagged = arclamp_gate([0.001, -0.002, NaN, 0.045]; gate_px = 0.02)
        @test gated == [true, true, false, true]
        @test flagged == [false, false, false, true]
    end

    @testset "joint night fit: synthetic recovery, gate and fixed chip scales" begin
        ll = load_arclamp_linelist()
        p = _aw_truth_poly(ll)
        fibs = 120:149
        # sequences: s1 position 1 (evening), s2 position 2 (+0.5 px), s3 position 1 but
        # displaced by 0.05 px (an unrecorded settling error, like the LCO evening block),
        # s4 position 1 (morning)
        shifts_px = [0.0, 0.5, 0.05, 0.0]
        lines = [_aw_lines(ll, p, fibs, s / 2048) for s in shifts_px]
        tseq = [60000.10, 60000.11, 60000.12, 60000.40]
        # object references: the true solution at position 1 and position 2
        function ref_at(sh_px, t, e)
            W = zeros(291, 101)
            for (jf, f) in enumerate(100:200), (i, pix) in enumerate(880:1170)
                W[i, jf] = _aw_poly(p, _aw_xt(_AW_CT, 2, (pix - 1024) / 2048) + sh_px / 2048)
            end
            ArcObjRef(e, t, "fpi", 880, 100:200, W)
        end
        refs = [ref_at(0.0, 60000.13, 21), ref_at(0.5, 60000.14, 22), ref_at(0.0, 60000.39, 40)]
        R = fit_arclamp_night(lines, tseq, "apo", ll; objrefs = refs)
        @test R.n_clusters == 2
        @test R.dpx1[2] ≈ 0.5 atol = 1e-4
        @test R.dpx1[3] ≈ 0.05 atol = 1e-4
        @test (R.seqflags[3] & ARCSEQ_GATE_FAIL) != 0
        @test all((R.seqflags[[1, 2, 4]] .& ARCSEQ_GATE_FAIL) .== 0)
        @test all((R.seqflags .& ARCSEQ_UNGATED) .== 0)
        @test R.used_seq == [1, 2, 4]
        @test (R.night_flags & ARCNIGHT_CHIPSCALE_FREE) != 0
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

        # fixed chip scales (the comb-measured b_R, b_B) are held exactly
        ct2 = [-1.07, 1.0 - 2.6e-4, 1.0755, 1.0 + 1.1e-4]
        lines2 = [_aw_lines(ll, p, fibs, s / 2048; ct = ct2) for s in (0.0, 0.5)]
        scales = repeat([ct2[2] ct2[4]], 300, 1)
        R2 = fit_arclamp_night(lines2, [60000.1, 60000.11], "apo", ll; chip_scales = scales)
        @test all(R2.TH[fibs, 2] .== ct2[2]) && all(R2.TH[fibs, 4] .== ct2[4])
        @test (R2.night_flags & ARCNIGHT_CHIPSCALE_FREE) == 0
        @test all((R2.seqflags .& ARCSEQ_UNGATED) .!= 0)    # no object references given
        G2 = arclamp_night_grid(R2)
        err2 = [abs(G2[Int(q), f, 1] - _aw_poly(p, _aw_xt(ct2, 1, (q - 1024) / 2048)))
                for f in fibs, q in px]
        @test maximum(err2) < 1e-4
    end
end
