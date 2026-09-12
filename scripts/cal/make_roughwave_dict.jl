# make_roughwave_dict.jl -- regenerate data/roughwave_dict.jld2 with the quartic
# rough wavelength model appended to each (telescope, chip) entry.
#
# Entry layout (tuple):
#   old: (a, b, lam_max, lam_min)                    -- linear model lambda = a + (pix - 1024) * b
#   new: (a, b, lam_max, lam_min, c0, c1, c2, c3, c4) -- quartic lambda = sum_k c_k x^k,
#        x = (pix - 1024)/2048 (coefficients from data/rough_poly_coeffs.csv; fit from
#        the adopted FPI nightAve solutions of reference nights apo/lco 60105, max
#        residual < 0.004 A vs 9-17 A for the linear model)
#
# Elements 1-4 are UNCHANGED, so every existing consumer that indexes [1]..[4]
# (chip wavelength ranges, linear fallback) behaves identically. The quartic is
# consumed by ApogeeReduction.rough_wave / rough_dispersion (src/skyline_peaks.jl).
#
# Usage (from the repo root):  julia --project=. scripts/cal/make_roughwave_dict.jl
# Idempotent: entries that already carry 9 elements have elements 5-9 replaced.

using CSV, DataFrames, JLD2, Printf

const REPO = dirname(dirname(@__DIR__))
const DICTFILE = joinpath(REPO, "data", "roughwave_dict.jld2")
const COEFFILE = joinpath(REPO, "data", "rough_poly_coeffs.csv")

coeffs = CSV.read(COEFFILE, DataFrame; comment = "#")
old_dict = load(DICTFILE, "roughwave_dict")
# rebuild with a widened value type (the stored Dict is typed to the old 4-tuples)
roughwave_dict = Dict{String, Dict{String, Tuple}}(
    tele => Dict{String, Tuple}(chip => v for (chip, v) in d) for (tele, d) in old_dict)

for r in eachrow(coeffs)
    old = roughwave_dict[r.tele][r.chip]
    lin = old[1:4]
    new = (lin..., r.c0, r.c1, r.c2, r.c3, r.c4)
    # consistency checks: the quartic must agree with the linear model it extends
    # (same reference geometry) at the few-A level of the linear model's curvature error
    quart(pix) = (x = (pix - 1024) / 2048; r.c0 + r.c1 * x + r.c2 * x^2 + r.c3 * x^3 +
                                           r.c4 * x^4)
    dc = abs(quart(1024) - old[1])
    db = abs(r.c1 / 2048 - old[2])
    dc < 8.0 || error("$(r.tele) $(r.chip): quartic center $(quart(1024)) vs linear a $(old[1])")
    db < 0.001 || error("$(r.tele) $(r.chip): quartic mean dispersion inconsistent with linear b")
    roughwave_dict[r.tele][r.chip] = new
    @printf("%s %s: center %.3f A (linear a %.3f), dispersion %.5f A/pix (linear b %.5f)\n",
        r.tele, r.chip, quart(1024), old[1], r.c1 / 2048, old[2])
end

jldsave(DICTFILE; roughwave_dict = roughwave_dict)
println("wrote $(DICTFILE)")
