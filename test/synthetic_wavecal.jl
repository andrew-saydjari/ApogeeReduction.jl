# Synthetic sky-line and FPI peak products for the wavelength-fit tests: a known per-fiber
# wavelength solution (wporder 4 polynomial in the transformed coordinate xt, chip offsets
# a_R / a_B with chip scale 1, the cporder-0 layout of the per-exposure sky fit), an optional
# per-exposure dither xt' = d(xt), and the files the package readers expect.
module SyntheticWavecal

using HDF5
using ApogeeReduction: N_FIBERS, N_XPIX, CHIP_LIST

const LIN0 = [16000.0, -580.0, 3.0, 0.5, 0.1]
const A_R = -1.070
const A_B = 1.0755
const CAV_APO = 3.7363125e7     # first-guess cavity size of the comb fit (APO), Angstrom
const MOFF_APO = 0.225

lin_true(k) = LIN0 .+ [0.01 * k, 1e-4 * k, 0.0, 0.0, 0.0]
xnorm(x) = (x - (N_XPIX ÷ 2)) / N_XPIX
chip_offset(c) = c == 1 ? A_R : c == 3 ? A_B : 0.0

"""wavelength at pixel x of chip c (1 R, 2 G, 3 B) of fiber k, dither polynomial d (coefficients)."""
lam(x, c, k; d = [0.0, 1.0]) = evalpoly(evalpoly(chip_offset(c) + xnorm(x), d), lin_true(k))

"""pixel of wavelength w on chip c of fiber k (bisection; wavelength decreases with pixel)."""
function pix(w, c, k; d = [0.0, 1.0])
    lo, hi = 1.0, Float64(N_XPIX)
    for _ in 1:80
        mid = (lo + hi) / 2
        lam(mid, c, k; d = d) > w ? (lo = mid) : (hi = mid)
    end
    return (lo + hi) / 2
end

chip_range(c, k; d = [0.0, 1.0]) = (lam(Float64(N_XPIX), c, k; d = d), lam(1.0, c, k; d = d))

"""
    write_sky_peaks(dir, expid, nlines_per_fiber; d, tele, mjd) -> path of the R-chip file

skyLinePeaks files for the 3 chips: fiber k gets `nlines_per_fiber[k]` lines spread over the
chips (round robin R, G, B), at the pixels of the true solution with dither `d`.
"""
function write_sky_peaks(dir, expid, nlines_per_fiber; d = [0.0, 1.0], tele = "apo", mjd = 60000)
    nmax = 12   # line slots per chip
    paths = String[]
    for (c, chip) in enumerate(CHIP_LIST)
        m = fill(NaN, nmax, 2, N_FIBERS)
        for k in 1:N_FIBERS
            n = nlines_per_fiber[k]
            nc = count(j -> mod1(j, 3) == c, 1:n)
            nc == 0 && continue
            lo, hi = chip_range(c, 150)
            for j in 1:nc
                # same catalog wavelengths for every fiber (slot j), well inside the chip
                w = lo + (hi - lo) * (j - 0.5) / nmax
                m[j, 1, k] = pix(w, c, k; d = d)
                m[j, 2, k] = w
            end
        end
        p = joinpath(dir, "skyLinePeaks_$(tele)_$(mjd)_$(expid)_$(chip)_object.h5")
        h5open(p, "w") do f
            f["sky_line_mat_clean"] = m
        end
        push!(paths, p)
    end
    return paths[1]
end

"""
    write_fpi_peaks(dir, expid, npeaks_max_per_fiber; d, lit, tele, mjd) -> path of the R-chip file

fpiPeaks files for the 3 chips: every FPI order m with 2 L / (m + eps) on the chip, at the pixels
of the true solution with dither `d`. Fibers with `lit[k] == false` have no peaks; fiber k keeps
at most `npeaks_max_per_fiber[k]` peaks over all chips (the first ones of each chip).
"""
function write_fpi_peaks(dir, expid, npeaks_max_per_fiber; d = [0.0, 1.0], lit = trues(N_FIBERS),
        tele = "apo", mjd = 60000, cav = CAV_APO, moff = MOFF_APO)
    paths = String[]
    for (c, chip) in enumerate(CHIP_LIST)
        lo, hi = chip_range(c, 150; d = d)
        ms = ceil(Int, 2cav / (hi - 3.0) - moff):floor(Int, 2cav / (lo + 3.0) - moff)
        npk = length(ms)
        mat = fill(NaN, npk, 3, N_FIBERS)
        cov = zeros(npk, 3, 3, N_FIBERS)
        for k in 1:N_FIBERS
            lit[k] || continue
            nkeep = min(npk, cld(npeaks_max_per_fiber[k], 3))
            for (j, mm) in enumerate(ms)
                j > nkeep && break
                w = 2cav / (mm + moff)
                x = pix(w, c, k; d = d)
                (1.5 <= x <= N_XPIX - 0.5) || continue
                mat[j, 2, k] = x
                cov[j, 2, 2, k] = 1e-4
            end
        end
        p = joinpath(dir, "fpiPeaks_$(tele)_$(mjd)_$(expid)_$(chip)_arclamp.h5")
        h5open(p, "w") do f
            f["fpi_line_mat"] = mat
            f["fpi_line_cov_mat"] = cov
        end
        push!(paths, p)
    end
    return paths[1]
end

"""Night solution in the cporder-0 layout: linParams (N_FIBERS x 5), nlParams (N_FIBERS x 2)."""
night_params() = (permutedims(reduce(hcat, [lin_true(k) for k in 1:N_FIBERS])),
    repeat([A_R A_B], N_FIBERS))

end
