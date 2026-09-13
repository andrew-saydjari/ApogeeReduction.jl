using Statistics: median, std
using Random

# Leave-one-out validation of sky-line -> wavelength assignments (stage 2 of
# the self-validating association design; see loo_validate_fiber docstring).
@testset "wavecal LOO assignment validation" begin
    AR = ApogeeReduction

    @testset "loo_linear_residuals is exact" begin
        Random.seed!(1)
        x = randn(20)
        A = AR.positional_poly_mat(x, porder = 3)
        y = A * [1.0, 0.5, -0.2, 0.05] .+ 0.01 .* randn(20)
        r, h = AR.loo_linear_residuals(A, y)
        for i in [1, 7, 20]
            idx = setdiff(1:20, i)
            ci = A[idx, :] \ y[idx]
            @test abs((y[i] - (A[i:i, :] * ci)[1]) - r[i]) < 1e-10
        end
        @test all(0 .<= h .<= 1)
    end

    @testset "lms_poly_classify survives heavy contamination" begin
        Random.seed!(2)
        x = collect(range(-0.5, 0.5, length = 8))
        y = 2.0 .+ 3.0 .* x .- 1.5 .* x .^ 2 .+ 0.005 .* randn(8)
        y[2] += 7.0
        y[5] -= 9.0
        keep, _, s = AR.lms_poly_classify(x, y; porder = 2, nsigma = 7.0)
        @test .!keep == [false, true, false, false, true, false, false, false]
        @test s < 0.1
    end

    @testset "loo_poly_reject drops outliers, respects scatter floor" begin
        Random.seed!(3)
        x = collect(range(-0.5, 0.5, length = 10))
        y = 1.0 .- 2.0 .* x .+ 0.01 .* randn(10)
        y[4] += 5.0
        keep, coeffs, capped = AR.loo_poly_reject(x, y; porder = 2, nsigma = 7.0)
        @test !keep[4] && count(keep) == 9 && !capped
        # a deviation below nsigma x floor must NOT be dropped even when the
        # rest of the points are unnaturally consistent
        y2 = 1.0 .- 2.0 .* x
        y2[4] += 0.3
        keep2, _, _ = AR.loo_poly_reject(x, y2; porder = 2, nsigma = 7.0,
            scatter_floor = 0.1)
        @test all(keep2)
    end

    # ---- synthetic per-fiber joint 3-chip solution ------------------------
    truth_lin = [16150.0, -1370.0, -30.0, 5.0, 1.0]
    offs = [-1.070, 0.0, 1.0755]
    function mkfiber(; nper = 6, noise = 0.010)
        xv = Float64[]
        yv = Float64[]
        ci = Int[]
        for c in 1:3, k in 1:nper
            x = -0.45 + 0.9 * (k - 1) / (nper - 1) + 0.01 * randn()
            xt = x + offs[c]
            push!(xv, x)
            push!(yv,
                (AR.positional_poly_mat([xt], porder = 4) * truth_lin)[1] + noise * randn())
            push!(ci, c)
        end
        return xv, yv, ci
    end
    cpp0 = zeros(3, 1)
    cpp0[1, 1] = -1.070
    cpp0[3, 1] = 1.0755

    @testset "healthy fiber: no drops, status OK" begin
        Random.seed!(4)
        for t in 1:10
            xv, yv, ci = mkfiber()
            d, rl, st, nl, lin, rv = AR.loo_validate_fiber(xv, yv, ci, cpp0)
            @test count(d) == 0
            @test st == AR.LOO_OK
            @test maximum(abs.(filter(!isnan, rl))) < 0.2
        end
    end

    @testset "single misassignment: detected, dropped, solution restored" begin
        Random.seed!(5)
        xv, yv, ci = mkfiber()
        yv[8] += 10.0   # ~the 16702->16692 theft distance
        d, rl, st, nl, lin, rv = AR.loo_validate_fiber(xv, yv, ci, cpp0)
        @test d[8] && count(d) == 1
        @test st == AR.LOO_DROPPED
        @test abs(rl[8] - 10.0) < 0.5          # LOO residual ~ the injected error
        @test isnan(rv[8])                     # excluded from the returned fit
        @test maximum(abs.(filter(!isnan, rv))) < 0.1  # remaining fit clean
    end

    @testset "historical pattern (one theft per chip): all detected" begin
        Random.seed!(6)
        for t in 1:5
            xv, yv, ci = mkfiber()
            bads = [2, 8, 14]
            for (b, s) in zip(bads, [-9.65, +27.4, -5.73])  # R, G, B theft modes
                yv[b] += s
            end
            d, rl, st, = AR.loo_validate_fiber(xv, yv, ci, cpp0)
            @test all(d[bads])
            @test count(d) == 3
        end
    end

    @testset "too few lines: untested flag, no silent verdict" begin
        Random.seed!(7)
        xv, yv, ci = mkfiber(nper = 3)   # 9 lines < nparam + 3
        d, rl, st, = AR.loo_validate_fiber(xv, yv, ci, cpp0)
        @test st == AR.LOO_UNTESTED
        @test count(d) == 0
    end

    @testset "get_sky_wavecal: drops applied, junk fiber flagged wholesale" begin
        Random.seed!(8)
        nline = 18
        nf = AR.N_FIBERS
        ux = fill(NaN, nline, nf)
        fw = fill(NaN, nline, nf)
        cI = ones(Int, nline, nf)
        for f in 1:nf
            xv, yv, ci = mkfiber()
            ux[:, f] .= xv
            fw[:, f] .= yv
            cI[:, f] .= ci
        end
        fw[5, 3] += 8.0                      # one misassignment on fiber 3
        fw[:, 7] .+= 3.0 .* randn(nline)     # fiber 7 is junk everywhere
        cpp0all = zeros(nf, 3, 1)
        cpp0all[:, 1, 1] .= -1.070
        cpp0all[:, 3, 1] .= 1.0755
        linP, nlP, rv, lr, ld, ls = AR.get_sky_wavecal(
            ux, fw, cI, cpp0all; cporder = 0, wporder = 4)
        @test ld[3, 5] && ls[3] == AR.LOO_DROPPED
        @test isnan(rv[3, 5])
        # the junk fiber cannot be caught from inside itself (its lines are
        # consistently wrong) -- the cross-fiber tripwire must flag it
        @test ls[7] == AR.LOO_WHOLESALE
        @test count(ls .== AR.LOO_OK) >= nf - 5
        # backward-compatible destructuring still works
        a, b, c = AR.get_sky_wavecal(ux, fw, cI, cpp0all; cporder = 0, wporder = 4)
        @test size(a) == size(linP)
    end
end
