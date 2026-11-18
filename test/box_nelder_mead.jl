using ApogeeReduction: BoxAffineSimplexer, box_nelder_mead_fit, normal_cdf
using Optim

# Regression for the Optim pin: the sky-line profile fits use Fminbox(NelderMead()). With
# Optim's default simplex the centroid vertex starts ~500 px outside its +-10 px box, and
# since Optim 1.12 that makes Nelder-Mead stop at iteration 0 (1.12-1.13: out-of-box
# minimizer; 2.x: the centroid never leaves the seed). The in-box simplex must give a real
# fit on every Optim version.
@testset "box Nelder-Mead (sky-line fit optimizer)" begin
    lb = [0.0, 990.4, 0.3, -5.0, 0.0]
    ub = [900.0, 1010.4, 3.0, 5.0, 40.0]
    x0 = [45.0, 1000.4, 1.0, 0.0, 20.0]
    S = Optim.simplexer(BoxAffineSimplexer(lb, ub), x0)
    @test length(S) == length(x0) + 1
    @test S[1] == x0
    for (j, v) in enumerate(S[2:end])
        @test all(lb .< v .< ub)          # every vertex strictly inside the box
        @test count(v .!= x0) == 1        # one coordinate moved per vertex
        @test v[j] != x0[j]
    end
    # vertices that Optim's default simplex already puts inside the box are unchanged
    @test S[2][1] ≈ 1.5 * x0[1] + 0.025
    @test S[4][3] ≈ 1.5 * x0[3] + 0.025

    # synthetic fixed-splitting doublet, same model/bounds/start as get_sky_peaks
    cpix_true, amp, wid, bkg, sepd = 1000.37, 300.0, 1.1, 20.0, 0.2
    cpix = 1000.4
    x = collect(Int(floor(cpix - 10)):Int(floor(cpix + 10)))
    x_edges = [x .- 0.5; x[end] + 0.5]
    mdl(p) = p[1] * 0.5 * diff(normal_cdf.(x_edges .- (p[2] - sepd / 2), p[3])) .+
             p[1] * 0.5 * diff(normal_cdf.(x_edges .- (p[2] + sepd / 2), p[3])) .+
             p[4] .* (x .- cpix) .+ p[5]
    y = mdl([amp, cpix_true, wid, 0.0, bkg])     # noiseless
    f(p) = sum((y .- mdl(p)) .^ 2)
    yamp = maximum(y) / sqrt(2pi)
    yoff = sort(y)[11]
    p0 = [yamp, cpix, 1, 0.0, yoff]
    lbf = [0.0, cpix - 10, 0.3, -5, 0]
    ubf = [20 * yamp, cpix + 10, 3.0, 5, 2 * yoff]
    p = box_nelder_mead_fit(f, lbf, ubf, p0)
    @test all(lbf .<= p .<= ubf)
    @test abs(p[2] - cpix_true) < 1e-3           # the centroid actually moved from the 1000.4 seed
    @test abs(p[1] / amp - 1) < 1e-3
    @test abs(p[3] - wid) < 1e-3
end
