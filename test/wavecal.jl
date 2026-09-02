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
end
