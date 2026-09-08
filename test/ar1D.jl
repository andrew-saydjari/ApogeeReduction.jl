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
                           RELTHRPT_NOTFINITE_BIT, RELTHRPT_UNUSABLE_BITS,
                           relthrpt_fiber_unusable, relthrpt_fiber_fluxable, get_relFlux,
                           safe_jldsave

    # bit table is disjoint and the aggressive set is exactly broken|notfinite
    @test RELTHRPT_UNUSABLE_BITS == (RELTHRPT_BROKEN_BIT | RELTHRPT_NOTFINITE_BIT)
    @test (RELTHRPT_UNUSABLE_BITS & RELTHRPT_NOFILE_BIT) == 0

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
    mktempdir() do dir
        npix, nfib = 64, 8
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

        # use_pix_mask is off by default (historical numbers preserved) and honours
        # bad_pix_bits when on
        mask = zeros(Int, npix, nfib)
        mask[1:(npix ÷ 2), 7] .= ApogeeReduction.bad_pix_bits
        flux_m = fill(100.0, npix, nfib)
        flux_m[1:(npix ÷ 2), 7] .= 1.0  # the bad half drags the unmasked median down
        fn4 = write_1d(joinpath(dir, "ar1Dcal_apo_57652_0013_B_domeflat.h5"), flux_m, mask)
        _, rel_off, _, _ = get_relFlux(fn4)
        _, rel_on, _, _ = get_relFlux(fn4, use_pix_mask = true)
        @test rel_off[7] < rel_on[7]
        @test rel_on[7] ≈ 1.0
        @test rel_off[1] ≈ rel_on[1]   # unaffected fibers unchanged
    end
end
