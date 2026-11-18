# Comb-measured inter-chip transform: calibration file, lookup and derivation helpers.
#
# Per fiber, a wavelength solution is ONE polynomial in the green-chip coordinate shared by
# the three chips; chips R and B are mapped onto it by x' = a_c + b_c * x (`transform_x_chips`).
# The sky-line fit (`get_and_save_sky_wavecal`, cporder 0) forces the chip scale b_c = 1, but
# the FPI comb (cporder 1, ~500 uniform peaks per chip) measures b_R - 1 ~ -2.6e-4 and
# b_B - 1 ~ +1.1e-4 at APO (LCO -0.9e-4 / -0.4e-4), stable to <~1e-5. Forcing b = 1 makes the
# shared polynomial absorb the scale error as 65-95 m/s of per-chip curvature.
#
# This file holds everything the pipeline needs to use a comb-derived chip scale WITHOUT any
# number baked into the source: the values live in a versioned HDF5 calibration file written by
# `scripts/cal/chip_transform/build_chip_transform.jl` (which re-derives them from the comb fits
# of a whole reduction), and the pipeline is pointed at that file by path
# (`pipeline_2d_1d.jl --chip_transform_cal`), exactly like the other calibration inputs.
# With no file configured the pipeline behaves as before (b = 1).
#
# File layout (format version 1), one group per telescope:
#   /<tele>/mjd_start, /<tele>/mjd_end   Int64 [n_blocks]   inclusive MJD range each block serves
#   /<tele>/data_mjd_min, data_mjd_max   Int64 [n_blocks]   first/last comb night actually used
#   /<tele>/n_nights                     Int64 [n_blocks]
#   /<tele>/scale                        Float64 [N_FIBERS, 2, n_blocks]  (b_R, b_B)
#   /<tele>/scale_err                    Float64 [N_FIBERS, 2, n_blocks]  robust s.e. of the block median
#   /<tele>/offset                       Float64 [N_FIBERS, 2, n_blocks]  (a_R, a_B), ximport units; QA only
#   /<tele>/fiber_filled                 Bool    [N_FIBERS, n_blocks]     fiber had no usable comb data;
#                                                                         value from the across-slit fit
#   /provenance/<key>                    String attributes of the build (code commit, input root,
#                                        block rule and its parameters, counts, date, ...)
#   attribute "format_version", "calibration_version"

using HDF5
using Statistics: median, mean

const CHIP_TRANSFORM_FORMAT_VERSION = 1

struct ChipTransformBlock
    mjd_start::Int
    mjd_end::Int
    data_mjd_min::Int
    data_mjd_max::Int
    n_nights::Int
    scale::Matrix{Float64}      # N_FIBERS x 2 (b_R, b_B)
    scale_err::Matrix{Float64}  # N_FIBERS x 2
    offset::Matrix{Float64}     # N_FIBERS x 2 (a_R, a_B)
    fiber_filled::Vector{Bool}
end

struct ChipTransformCal
    path::String
    calibration_version::String
    fingerprint::String     # FNV-1a 64-bit hash of the file's bytes (see `file_fingerprint`)
    provenance::Dict{String, String}
    blocks::Dict{String, Vector{ChipTransformBlock}}
end

"""
    file_fingerprint(path) -> String

FNV-1a 64-bit hash of a file's bytes, as 16 hex digits. Identifies the exact calibration file
a product was made with (two files with the same version string but different contents get
different fingerprints). Dependency-free; the calibration files are ~100 KB.
"""
function file_fingerprint(path::AbstractString)
    h = 0xcbf29ce484222325
    for b in read(path)
        h = (h ⊻ UInt64(b)) * 0x00000100000001b3
    end
    return string(h, base = 16, pad = 16)
end

"""
    chip_transform_id(cal) -> String

The identity of the chip-transform input recorded in (and checked against) the metadata of
every sky wavelength product: "none" without a calibration, else "<version>:<fingerprint>".
"""
chip_transform_id(::Nothing) = "none"
chip_transform_id(cal::ChipTransformCal) = "$(cal.calibration_version):$(cal.fingerprint)"

function Base.show(io::IO, cal::ChipTransformCal)
    nb = join(["$(t): $(length(b)) block(s)" for (t, b) in sort(collect(cal.blocks))], ", ")
    print(io, "ChipTransformCal($(cal.calibration_version); $(nb); $(cal.path))")
end

"""
    write_chip_transform_cal(path, blocks::Dict{String, Vector{ChipTransformBlock}},
        provenance::Dict{String, String}; calibration_version)

Write a chip-transform calibration file (format above). Blocks of one telescope must be
sorted, non-overlapping MJD ranges. Written to a temporary name and moved into place.
"""
function write_chip_transform_cal(path::AbstractString, blocks::Dict{String, Vector{ChipTransformBlock}},
        provenance::Dict{String, String}; calibration_version::AbstractString)
    for (tele, bl) in blocks
        check_chip_transform_blocks(tele, bl)
    end
    tmp = path * ".tmp$(getpid())"
    h5open(tmp, "w") do f
        HDF5.attributes(f)["format_version"] = CHIP_TRANSFORM_FORMAT_VERSION
        HDF5.attributes(f)["calibration_version"] = String(calibration_version)
        for (tele, bl) in blocks
            g = create_group(f, tele)
            g["mjd_start"] = [b.mjd_start for b in bl]
            g["mjd_end"] = [b.mjd_end for b in bl]
            g["data_mjd_min"] = [b.data_mjd_min for b in bl]
            g["data_mjd_max"] = [b.data_mjd_max for b in bl]
            g["n_nights"] = [b.n_nights for b in bl]
            g["scale"] = cat([b.scale for b in bl]...; dims = 3)
            g["scale_err"] = cat([b.scale_err for b in bl]...; dims = 3)
            g["offset"] = cat([b.offset for b in bl]...; dims = 3)
            g["fiber_filled"] = UInt8.(hcat([b.fiber_filled for b in bl]...))
        end
        p = create_group(f, "provenance")
        for (k, v) in provenance
            HDF5.attributes(p)[k] = String(v)
        end
    end
    mv(tmp, path; force = true)
    return path
end

function check_chip_transform_blocks(tele, bl::Vector{ChipTransformBlock})
    isempty(bl) && throw(ArgumentError("telescope $(tele): no blocks"))
    for (i, b) in enumerate(bl)
        b.mjd_start <= b.mjd_end ||
            throw(ArgumentError("telescope $(tele) block $(i): mjd_start > mjd_end"))
        size(b.scale) == (N_FIBERS, 2) ||
            throw(ArgumentError("telescope $(tele) block $(i): scale must be $(N_FIBERS) x 2"))
        all(isfinite, b.scale) ||
            throw(ArgumentError("telescope $(tele) block $(i): non-finite chip scale"))
        i > 1 && bl[i - 1].mjd_end >= b.mjd_start &&
            throw(ArgumentError("telescope $(tele): blocks $(i-1) and $(i) overlap or are unsorted"))
    end
    return nothing
end

"""
    load_chip_transform_cal(path) -> ChipTransformCal

Read a calibration file written by `write_chip_transform_cal`. Throws if the file does not
exist or has an unknown format version (a configured but unreadable calibration input is a
hard error, never a silent fallback to b = 1).
"""
function load_chip_transform_cal(path::AbstractString)
    isfile(path) || throw(ArgumentError("chip-transform calibration file not found: $(path)"))
    h5open(path, "r") do f
        fv = read(HDF5.attributes(f)["format_version"])
        fv == CHIP_TRANSFORM_FORMAT_VERSION ||
            throw(ArgumentError("$(path): chip-transform format version $(fv), expected $(CHIP_TRANSFORM_FORMAT_VERSION)"))
        version = read(HDF5.attributes(f)["calibration_version"])
        prov = Dict{String, String}()
        if haskey(f, "provenance")
            a = HDF5.attributes(f["provenance"])
            for k in keys(a)
                prov[k] = string(read(a[k]))
            end
        end
        blocks = Dict{String, Vector{ChipTransformBlock}}()
        for tele in keys(f)
            tele == "provenance" && continue
            g = f[tele]
            ms = read(g["mjd_start"]); me = read(g["mjd_end"])
            dmin = read(g["data_mjd_min"]); dmax = read(g["data_mjd_max"])
            nn = read(g["n_nights"])
            sc = read(g["scale"]); se = read(g["scale_err"]); of = read(g["offset"])
            ff = read(g["fiber_filled"]) .!= 0
            blocks[tele] = [ChipTransformBlock(ms[i], me[i], dmin[i], dmax[i], nn[i],
                                sc[:, :, i], se[:, :, i], of[:, :, i], ff[:, i])
                            for i in eachindex(ms)]
            check_chip_transform_blocks(tele, blocks[tele])
        end
        return ChipTransformCal(String(path), version, file_fingerprint(path), prov, blocks)
    end
end

"""
    chip_transform_block(cal, tele, mjd) -> ChipTransformBlock or nothing

The block whose inclusive [mjd_start, mjd_end] contains `mjd`, or `nothing` if the file has no
block for this telescope and MJD (the caller then keeps the b = 1 behaviour and says so).
Which MJDs a block serves (e.g. whether the first FPS block is extended back over the plate
era) is decided when the file is built and recorded in its provenance, not here. The sky fit
treats an MJD without a block as an error (see `get_and_save_sky_wavecal`).
"""
function chip_transform_block(cal::ChipTransformCal, tele::AbstractString, mjd::Integer)
    bl = get(cal.blocks, String(tele), nothing)
    isnothing(bl) && return nothing
    for b in bl
        b.mjd_start <= mjd <= b.mjd_end && return b
    end
    return nothing
end

"""
    chip_scales(cal, tele, mjd) -> Matrix{Float64} (N_FIBERS x 2: b_R, b_B) or nothing
"""
function chip_scales(cal::ChipTransformCal, tele::AbstractString, mjd::Integer)
    b = chip_transform_block(cal, tele, mjd)
    return isnothing(b) ? nothing : copy(b.scale)
end
chip_scales(::Nothing, tele, mjd) = nothing

"""
    scale_chip_x(x, chipInt, scales_fiber) -> Vector

Apply a FIXED chip scale to green-chip-relative pixel coordinates: x -> b_c * x for chips R
(index 1) and B (index 3), unchanged for G. A cporder-0 fit (x' = a_c + x) on the scaled
coordinates is then exactly the cporder-1 model x' = a_c + b_c * x with b_c held fixed.
`scales_fiber` = (b_R, b_B) for this fiber.
"""
function scale_chip_x(x::AbstractVector, chipInt::AbstractVector, scales_fiber)
    out = copy(x)
    for i in eachindex(out)
        if chipInt[i] == 1
            out[i] *= scales_fiber[1]
        elseif chipInt[i] == 3
            out[i] *= scales_fiber[2]
        end
    end
    return out
end

"""
    expand_nlParams_fixed_scale(nlParams0, scales) -> nlParams1

Convert cporder-0 chip parameters (a_R, a_B) fitted on scaled coordinates into the cporder-1
layout (a_R, b_R, a_B, b_B) used everywhere downstream (`params2ChipPolyParams!`), so that
night averaging, dither fits, the FPI initial guess and the stored chipWaveSoln all apply the
fixed scale.
"""
function expand_nlParams_fixed_scale(nlParams0::AbstractMatrix, scales::AbstractMatrix)
    size(nlParams0, 2) == 2 || throw(ArgumentError("expected cporder-0 nlParams (a_R, a_B)"))
    return hcat(nlParams0[:, 1], scales[:, 1], nlParams0[:, 2], scales[:, 2])
end

# ------------------------------------------------------------------------------------------
# derivation helpers (used by scripts/cal/chip_transform/build_chip_transform.jl)
# ------------------------------------------------------------------------------------------

"Robust scale: 1.4826 * median absolute deviation over the finite entries (NaN if < 3)."
function robust_sigma(x)
    y = filter(isfinite, x)
    length(y) < 3 && return NaN
    m = median(y)
    return 1.4826 * median(abs.(y .- m))
end

"""
    detect_blocks(mjds, series; min_nights, nsigma, min_step) -> Vector{Int}

Data-driven block boundaries for a night-ordered set of measurements. `series` is
n_nights x n_quantities (e.g. the fiber-median b_R and b_B of each night).

1. Noise: per quantity, sigma = robust_sigma(successive differences) / sqrt 2, estimated once
   over the whole series (insensitive to a few steps).
2. Outlier nights are screened out of the statistic: a night more than `nsigma` sigma from the
   running median of the `min_nights` nights around it.
3. Binary segmentation: within a segment, every split leaving >= `min_nights` screened nights
   on each side is scored, per quantity, by
       z = |mean(right) - mean(left)| / (sigma * sqrt(1/n_left + 1/n_right)).
   The best split is accepted if z > `nsigma` AND |mean(right) - mean(left)| > `min_step[q]`
   for that quantity; both halves are then searched again.

Returns the indices (into the sorted nights) at which a new block STARTS, excluding 1.
"""
function detect_blocks(mjds::AbstractVector, series::AbstractMatrix; min_nights::Integer,
        nsigma::Real, min_step::AbstractVector)
    issorted(mjds) || throw(ArgumentError("mjds must be sorted"))
    n, nq = size(series)
    length(min_step) == nq || throw(ArgumentError("min_step needs one entry per quantity"))
    sig = [robust_sigma(diff(series[isfinite.(series[:, q]), q])) / sqrt(2) for q in 1:nq]
    Y = copy(Float64.(series))
    h = max(1, min_nights ÷ 2)
    for q in 1:nq, i in 1:n
        w = filter(isfinite, series[max(1, i - h):min(n, i + h), q])
        (isfinite(Y[i, q]) && length(w) >= 3 && isfinite(sig[q])) || continue
        abs(Y[i, q] - median(w)) > nsigma * sig[q] && (Y[i, q] = NaN)
    end
    starts = Int[]
    function search(lo, hi)
        hi - lo + 1 < 2 * min_nights && return
        best = (0.0, 0)
        for k in (lo + min_nights):(hi - min_nights + 1)
            for q in 1:nq
                (isfinite(sig[q]) && sig[q] > 0) || continue
                L = filter(isfinite, view(Y, lo:(k - 1), q)); R = filter(isfinite, view(Y, k:hi, q))
                (length(L) >= min_nights && length(R) >= min_nights) || continue
                d = mean(R) - mean(L)
                abs(d) > min_step[q] || continue
                z = abs(d) / (sig[q] * sqrt(1 / length(L) + 1 / length(R)))
                z > best[1] && (best = (z, k))
            end
        end
        if best[1] > nsigma
            push!(starts, best[2])
            search(lo, best[2] - 1)
            search(best[2], hi)
        end
    end
    search(1, n)
    return sort(starts)
end

"""
    robust_block_median(V; nsigma) -> (median, stderr, n_used)

Median over nights (columns of a fibers x nights matrix, NaN = missing) per fiber, after
rejecting entries more than `nsigma` robust sigma from that fiber's median. stderr is the
large-sample standard error of a median, 1.2533 * sigma / sqrt(n).
"""
function robust_block_median(V::AbstractMatrix; nsigma::Real)
    nf = size(V, 1)
    med = fill(NaN, nf); se = fill(NaN, nf); nu = zeros(Int, nf)
    for f in 1:nf
        y = filter(isfinite, V[f, :])
        length(y) == 0 && continue
        m = median(y); s = robust_sigma(y)
        if isfinite(s) && s > 0
            y = y[abs.(y .- m) .<= nsigma * s]
            m = median(y); s = robust_sigma(y)
        end
        med[f] = m; nu[f] = length(y)
        se[f] = length(y) >= 3 ? 1.2533 * s / sqrt(length(y)) : NaN
    end
    return med, se, nu
end

"""
    fiber_poly_fit(fibers, values, deg; nsigma, max_iter = 5) -> (coefficients, keep)

Robust least-squares polynomial of degree `deg` in u = (fiber - 150.5) / 150 over the finite
entries, iteratively clipping at `nsigma` robust sigma. Evaluate with `fiber_poly_eval`.
"""
function fiber_poly_fit(fibers::AbstractVector, values::AbstractVector, deg::Integer;
        nsigma::Real, max_iter::Integer = 5)
    u = (fibers .- 150.5) ./ 150
    keep = isfinite.(values)
    c = zeros(deg + 1)
    for _ in 1:max_iter
        A = hcat([u[keep] .^ p for p in 0:deg]...)
        c = A \ values[keep]
        r = values .- fiber_poly_eval(c, fibers)
        s = robust_sigma(r[keep])
        new_keep = isfinite.(values) .& (abs.(r) .<= nsigma * s)
        (isfinite(s) && s > 0) || break
        new_keep == keep && break
        keep = new_keep
    end
    return c, keep
end
fiber_poly_eval(c, fibers) = [sum(c[p + 1] * ((f - 150.5) / 150)^p for p in 0:(length(c) - 1)) for f in fibers]
