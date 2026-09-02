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
