@testset "ar1D" begin
    @testset "normalize_reinterp_spectra! zero-good-pixel fiber (A4)" begin
        npix, nfib = 12, 4
        cntvec = zeros(Int, npix, nfib)
        cntvec[:, 1] .= 2           # healthy fiber, 2 frames everywhere
        cntvec[:, 2] .= 0           # dead fiber: no good pixels in any frame
        cntvec[:, 3] .= 2           # fiber with partial coverage:
        cntvec[3, 3] = 1            #   one pixel seen in only 1 of 2 frames
        cntvec[:, 4] .= 1           # single-frame fiber

        # accumulated (pre-normalization) flux and variance, as built by
        # reinterp_spectra's accumulation loop: zero where nothing contributed
        outflux = zeros(npix, nfib)
        outvar = zeros(npix, nfib)
        outflux[:, 1] .= 20.0
        outvar[:, 1] .= 8.0
        outflux[:, 3] .= 20.0
        outvar[:, 3] .= 8.0
        outflux[3, 3] = 7.0
        outvar[3, 3] = 3.0
        outflux[:, 4] .= 5.0
        outvar[:, 4] .= 2.0

        outivar, outmsk = ApogeeReduction.normalize_reinterp_spectra!(outflux, outvar, cntvec)

        # A4: dead fiber must come out flux = 0, ivar = 0, msk = false
        # (previously NaN / NaN / true)
        @test all(outflux[:, 2] .== 0.0)
        @test all(outivar[:, 2] .== 0.0)
        @test all(.!outmsk[:, 2])
        @test !any(isnan.(outflux)) && !any(isnan.(outivar))

        # healthy fibers unchanged by the fix: flux/2, ivar = 4/var, msk true
        @test all(outflux[:, 1] .== 10.0)
        @test all(outivar[:, 1] .== 0.5)
        @test all(outmsk[:, 1])
        @test all(outflux[:, 4] .== 5.0)
        @test all(outivar[:, 4] .== 0.5)
        @test all(outmsk[:, 4])

        # partial-coverage pixel: masked, ivar zeroed, flux still normalized
        @test outmsk[3, 3] == false
        @test outivar[3, 3] == 0.0
        @test outflux[3, 3] == 3.5
        @test all(outmsk[[1:2; 4:npix], 3])
        @test all(outivar[[1:2; 4:npix], 3] .== 0.5)
    end
end

@testset "get_relFlux throughput bitmask" begin
    using ApogeeReduction: RELTHRPT_WARN_BIT, RELTHRPT_BROKEN_BIT, RELTHRPT_NOFILE_BIT,
                           RELTHRPT_NOTFINITE_BIT, RELTHRPT_LOWGOODPIX_BIT,
                           RELTHRPT_MIN_GOODPIX, RELTHRPT_UNUSABLE_BITS,
                           relthrpt_fiber_unusable, relthrpt_fiber_fluxable, get_relFlux,
                           safe_jldsave

    # bit table is disjoint and the aggressive set is exactly broken|notfinite
    @test RELTHRPT_UNUSABLE_BITS ==
          (RELTHRPT_BROKEN_BIT | RELTHRPT_NOTFINITE_BIT | RELTHRPT_LOWGOODPIX_BIT)
    @test (RELTHRPT_UNUSABLE_BITS & RELTHRPT_NOFILE_BIT) == 0
    @test relthrpt_fiber_unusable([RELTHRPT_LOWGOODPIX_BIT]) == [true]

    # a warn-only fiber is still fluxed; broken / non-finite are not
    @test relthrpt_fiber_fluxable([0, RELTHRPT_WARN_BIT]) == [true, true]
    @test relthrpt_fiber_fluxable([RELTHRPT_WARN_BIT | RELTHRPT_BROKEN_BIT]) == [false]
    @test relthrpt_fiber_fluxable([RELTHRPT_NOTFINITE_BIT]) == [false]
    @test relthrpt_fiber_fluxable([RELTHRPT_NOFILE_BIT]) == [false]
    @test relthrpt_fiber_unusable([0, RELTHRPT_WARN_BIT, RELTHRPT_NOFILE_BIT]) ==
          [false, false, false]
    @test relthrpt_fiber_unusable([RELTHRPT_BROKEN_BIT]) == [true]
    # broadcasts over the (N_CHIPS, N_FIBERS) stack shape used in ar1Duni*
    @test relthrpt_fiber_unusable([0 RELTHRPT_BROKEN_BIT; 0 0]) == [false true; false false]

    # --- get_relFlux against a synthetic 1D product ---
    # npix is the real 2048: anything below RELTHRPT_MIN_GOODPIX would (correctly)
    # trip the good-pixel floor on every fiber and mask what these cases are testing.
    mktempdir() do dir
        npix, nfib = 2048, 8
        function write_1d(fn, flux, mask = zeros(Int, npix, nfib))
            safe_jldsave(fn, Dict{String, Any}(); flux_1d = flux, mask_1d = mask)
            fn
        end

        flux = fill(100.0, npix, nfib)
        flux[:, 2] .= 0.5      # relthrpt 0.005 -> broken
        flux[:, 3] .= 70.0     # 0.7 -> low but well above rel_val_cut
        flux[:, 4] .= -100.0   # negative throughput -> broken
        fn = write_1d(joinpath(dir, "ar1Dcal_apo_57652_0010_B_domeflat.h5"), flux)
        absthrpt, relthrpt, bits, _ = get_relFlux(fn)

        @test length(relthrpt) == nfib
        @test bits[1] == 0
        @test (bits[2] & RELTHRPT_BROKEN_BIT) != 0
        @test (bits[4] & RELTHRPT_BROKEN_BIT) != 0
        @test (bits[3] & RELTHRPT_BROKEN_BIT) == 0   # 0.7 is not broken
        @test !any(relthrpt_fiber_unusable(bits)[[1, 3, 5, 6, 7, 8]])

        # THE NaN BUG: `NaN < x` is false in Julia, so before the fix an all-NaN
        # (or all-zero) fiber passed BOTH threshold tests as good and was then
        # divided into flux/ivar, silently NaN-poisoning the fiber downstream.
        flux_nan = fill(100.0, npix, nfib)
        flux_nan[:, 5] .= NaN
        flux_nan[:, 6] .= 0.0   # nanzeromedian of an all-zero fiber is also NaN
        fn2 = write_1d(joinpath(dir, "ar1Dcal_apo_57652_0011_B_domeflat.h5"), flux_nan)
        _, relthrpt2, bits2, _ = get_relFlux(fn2)

        @test isnan(relthrpt2[5]) && isnan(relthrpt2[6])
        for f in (5, 6)
            @test (bits2[f] & RELTHRPT_NOTFINITE_BIT) != 0
            @test (bits2[f] & RELTHRPT_BROKEN_BIT) != 0  # legacy bit-1 consumers stay correct
            @test relthrpt_fiber_unusable(bits2)[f]
            @test !relthrpt_fiber_fluxable(bits2)[f]
        end
        @test !any(relthrpt_fiber_unusable(bits2)[[1, 2, 3, 4, 7, 8]])

        # a wholly dead exposure must flag every fiber, not none: the exposure-level
        # normalization is NaN, so the warn threshold is NaN and the old code flagged
        # nothing at all
        fn3 = write_1d(joinpath(dir, "ar1Dcal_apo_57652_0012_B_domeflat.h5"),
            fill(NaN, npix, nfib))
        _, _, bits3, _ = get_relFlux(fn3)
        @test all(relthrpt_fiber_unusable(bits3))

        # use_pix_mask is ON by default and honours bad_pix_bits; the old behaviour
        # is still reachable for before/after comparisons
        mask = zeros(Int, npix, nfib)
        mask[1:(npix ÷ 2), 7] .= ApogeeReduction.bad_pix_bits
        flux_m = fill(100.0, npix, nfib)
        flux_m[1:(npix ÷ 2), 7] .= 1.0  # the bad half drags the unmasked median down
        fn4 = write_1d(joinpath(dir, "ar1Dcal_apo_57652_0013_B_domeflat.h5"), flux_m, mask)
        _, rel_off, _, _ = get_relFlux(fn4, use_pix_mask = false)
        _, rel_on, _, _ = get_relFlux(fn4)   # default is ON
        @test rel_off[7] < rel_on[7]
        @test rel_on[7] ≈ 1.0
        @test rel_off[1] ≈ rel_on[1]   # unaffected fibers unchanged
    end
end

@testset "get_relFlux pixel mask + good-pixel floor" begin
    using ApogeeReduction: RELTHRPT_WARN_BIT, RELTHRPT_BROKEN_BIT, RELTHRPT_NOTFINITE_BIT,
                           RELTHRPT_LOWGOODPIX_BIT, RELTHRPT_MIN_GOODPIX,
                           relthrpt_fiber_unusable, relthrpt_fiber_fluxable, get_relFlux,
                           safe_jldsave, bad_pix_bits

    # the floor is data-derived, not round-number: 256 is the smallest sub-sample size
    # whose p95 median error (0.064) is below rel_val_cut (0.07)
    @test RELTHRPT_MIN_GOODPIX == 256
    # ... and it is far below anything the real data ever shows (min 1805 of 2048)
    @test RELTHRPT_MIN_GOODPIX < 1805 / 7

    mktempdir() do dir
        npix, nfib = 2048, 8
        write_1d = function (fn, flux, mask)
            safe_jldsave(fn, Dict{String, Any}(); flux_1d = flux, mask_1d = mask)
            fn
        end

        # ---- (a) THE ALL-MASKED FIBER MUST FLAG, NOT PASS -----------------------
        # Masking can leave a fiber with ZERO good pixels. median(empty) is NaN and
        # `NaN < thresh` is false, so before the non-finite fix such a fiber would be
        # judged GOOD and then divided into flux/ivar as NaN. This is the case that
        # makes turning the mask on safe rather than actively worse.
        flux = fill(100.0, npix, nfib)
        mask = zeros(Int, npix, nfib)
        mask[:, 2] .= bad_pix_bits            # fiber 2: EVERY pixel bad
        mask[1:(npix - 10), 3] .= bad_pix_bits # fiber 3: only 10 good pixels left
        fn = write_1d(joinpath(dir, "ar1Dcal_apo_57652_0020_B_domeflat.h5"), flux, mask)
        absthrpt, relthrpt, bits, _ = get_relFlux(fn)

        @test isnan(relthrpt[2])                             # not silently "good"
        @test (bits[2] & RELTHRPT_NOTFINITE_BIT) != 0
        @test (bits[2] & RELTHRPT_BROKEN_BIT) != 0           # legacy bit-1 consumers
        @test (bits[2] & RELTHRPT_LOWGOODPIX_BIT) != 0       # and the CAUSE is recorded
        @test relthrpt_fiber_unusable(bits)[2]
        @test !relthrpt_fiber_fluxable(bits)[2]              # never divided into flux

        # ---- (b) THE FLOOR: a starved-but-nonempty fiber also flags -------------
        @test isnan(relthrpt[3])
        @test (bits[3] & RELTHRPT_LOWGOODPIX_BIT) != 0
        @test !relthrpt_fiber_fluxable(bits)[3]
        # everyone else is untouched and unflagged
        @test all(bits[[1, 4, 5, 6, 7, 8]] .== 0)
        @test all(relthrpt[[1, 4, 5, 6, 7, 8]] .≈ 1.0)

        # a starved fiber must NOT contaminate the exposure-level normalization:
        # its junk median is NaN'd before the median-of-fibers is taken
        flux2 = fill(100.0, npix, nfib)
        flux2[:, 3] .= 1.0e6                  # would wreck the normalizer if it counted
        fn2 = write_1d(joinpath(dir, "ar1Dcal_apo_57652_0021_B_domeflat.h5"), flux2, mask)
        _, rel2, bits2, _ = get_relFlux(fn2)
        @test all(rel2[[1, 4, 5, 6, 7, 8]] .≈ 1.0)
        @test (bits2[3] & RELTHRPT_LOWGOODPIX_BIT) != 0

        # exactly at the floor is acceptable; one below is not
        for (ngood, want_flag) in ((RELTHRPT_MIN_GOODPIX, false),
            (RELTHRPT_MIN_GOODPIX - 1, true))
            mk = zeros(Int, npix, nfib)
            mk[1:(npix - ngood), 5] .= bad_pix_bits
            fnx = write_1d(joinpath(dir, "ar1Dcal_apo_57652_002$(ngood % 7)_B_domeflat.h5"),
                fill(100.0, npix, nfib), mk)
            _, _, bx, _ = get_relFlux(fnx)
            @test ((bx[5] & RELTHRPT_LOWGOODPIX_BIT) != 0) == want_flag
        end

        # ---- the floor never fires on healthy data ------------------------------
        clean = write_1d(joinpath(dir, "ar1Dcal_apo_57652_0030_B_domeflat.h5"),
            fill(100.0, npix, nfib), zeros(Int, npix, nfib))
        _, _, bclean, _ = get_relFlux(clean)
        @test all(bclean .== 0)

        # ---- bad-yet-finite-nonzero pixels are what this change removes ---------
        # nanzeromedian already dropped NaN and exact zeros; a cosmic-ray pixel that
        # is finite and non-zero used to contribute to the throughput estimate.
        fluxc = fill(100.0, npix, nfib)
        maskc = zeros(Int, npix, nfib)
        fluxc[1:1200, 6] .= 5.0               # finite, non-zero, and flagged bad
        maskc[1:1200, 6] .= bad_pix_bits      # a majority of the fiber, so the median moves
        fnc = write_1d(joinpath(dir, "ar1Dcal_apo_57652_0031_B_domeflat.h5"), fluxc, maskc)
        _, rel_off, boff, _ = get_relFlux(fnc, use_pix_mask = false)
        _, rel_on, bon, _ = get_relFlux(fnc)
        @test rel_off[6] < 0.5 * rel_on[6]    # the old median was dragged to ~5/100
        @test rel_on[6] ≈ 1.0
        # and the consequence: unmasked, a perfectly healthy fiber is called BROKEN
        @test (boff[6] & RELTHRPT_BROKEN_BIT) != 0
        @test bon[6] == 0
        # 848 good pixels is still comfortably above the floor
        @test (bon[6] & RELTHRPT_LOWGOODPIX_BIT) == 0
    end
end
