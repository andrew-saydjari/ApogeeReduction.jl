"""
    DefectFinder

Coherent-detector-defect region finder for ApogeeReduction 2D calibration products.

The design premise is that the defects AR currently misses are *low contrast but
spatially coherent*: a per-pixel sigma clip never fires because no single pixel is a
large outlier, but a contiguous block of ~10^2 pixels all deviate in the same sense.
The detector here is therefore a **matched filter on a robustly detrended residual**,
followed by connected-component labelling, run independently on several products and
several epochs, with agreement across products/epochs used as the confirmation.

Coordinate convention (verified against AR source):
  * `dimage`, `ivarimage`, `pix_bitmask`, `dark_rate` are (2560, 2048); the first
    2048 columns are the detector, columns 2049:2560 are the reference array.
  * `flat_im` is (2040, 2040) and corresponds to `dimage[5:2044, 5:2044]`
    (`AR/src/ar3D.jl:714`).
  * 1D extraction indexes `dimage[xpix, ypix]` directly with xpix in 1:2048 and
    `extract_trace_centers` in detector row units (`AR/src/ar1D.jl:797`), so 1D column
    index == detector column and trace centre == detector row.
  * `ar1Dunical/extract_trace_coords[:,:,1]` is the *stacked* x, related to the
    detector column by `col = 2049 - x_stack` (`AR/src/ar1D.jl:545-546, 633-644`);
    `[:,:,3]` is chip index 1=R, 2=G, 3=B.

Internally everything is carried in **science coordinates** `s = 1:2040`, which map to
detector coordinates as `det = s + 4`.  All reported extents are DETECTOR coordinates.

No AR code is imported and nothing here writes into a reduction.
"""
module DefectFinder

using Statistics

export detrend_x, zmap, matched_filter_detect, standardize_rows!, highpass, label_components, dilate1,
       regions_from_mask, RegionBox, SCI_OFF, SCI_N, sci_view, DEFAULT_SCALES

const SCI_OFF = 4          # detector index = science index + SCI_OFF
const SCI_N = 2040

"""Rectangular matched-filter scales `(wx, wy)`.

`(3,3)`…`(33,33)` catch blobs; the thin ones catch bad columns/rows and thin
streaks.  Scale `(1,1)` is deliberately absent: single-pixel outliers are what AR's
existing per-pixel cuts already handle, and including them would swamp the region
list with isolated hot pixels.
"""
const DEFAULT_SCALES = [(3, 3), (5, 5), (9, 9), (17, 17), (33, 33),
    (1, 9), (1, 33), (9, 1), (33, 1),
    (3, 17), (17, 3), (5, 65), (65, 5)]

"""Science-area view of a raw (2560,2048) or (2040,2040) product."""
function sci_view(A::AbstractMatrix)
    if size(A, 1) >= 2048 && size(A, 2) >= 2048
        return Array{Float64}(A[(SCI_OFF + 1):(SCI_OFF + SCI_N), (SCI_OFF + 1):(SCI_OFF + SCI_N)])
    elseif size(A) == (SCI_N, SCI_N)
        return Array{Float64}(A)
    else
        error("unexpected product size $(size(A))")
    end
end

# ---------------------------------------------------------------------------
# robust detrending
# ---------------------------------------------------------------------------

@inline function _median!(buf, v)
    n = length(v)
    resize!(buf, n)
    copyto!(buf, v)
    sort!(buf)
    isodd(n) ? buf[(n + 1) ÷ 2] : 0.5 * (buf[n ÷ 2] + buf[n ÷ 2 + 1])
end

"""
    detrend_x(A; medblk=64, sigblk=256) -> (M, S)

Model each row of `A` (first index = the direction being detrended) as a slowly
varying level `M` plus noise of scale `S`.

`M` is a piecewise-linear interpolation of block medians (block width `medblk`), so it
follows the smooth spectral/illumination shape along the dispersion direction while
being insensitive to a defect occupying a minority of a block.

`S` is estimated from *first differences* within blocks of width `sigblk`,
`1.4826*median(|ΔA|)/√2`.  The difference estimator is immune to the smooth trend that
`M` is modelling, and — importantly — a flat-valued defect block produces small
differences, so a defect does not inflate its own noise estimate.
"""
function detrend_x(A::Matrix{Float64}; medblk::Int = 64, sigblk::Int = 256)
    nx, ny = size(A)
    nbm = cld(nx, medblk)
    nbs = cld(nx, sigblk)
    M = Array{Float64}(undef, nx, ny)
    S = Array{Float64}(undef, nx, ny)
    cm = Vector{Float64}(undef, nbm)   # block centres (median blocks)
    vm = Vector{Float64}(undef, nbm)
    cs = Vector{Float64}(undef, nbs)
    vs = Vector{Float64}(undef, nbs)
    buf = Float64[]
    dbuf = Float64[]
    for b in 1:nbm
        lo = (b - 1) * medblk + 1
        hi = min(nx, b * medblk)
        cm[b] = 0.5 * (lo + hi)
    end
    for b in 1:nbs
        lo = (b - 1) * sigblk + 1
        hi = min(nx, b * sigblk)
        cs[b] = 0.5 * (lo + hi)
    end
    @inbounds for j in 1:ny
        col = view(A, :, j)
        for b in 1:nbm
            lo = (b - 1) * medblk + 1
            hi = min(nx, b * medblk)
            vm[b] = _median!(buf, view(col, lo:hi))
        end
        for b in 1:nbs
            lo = (b - 1) * sigblk + 1
            hi = min(nx, b * sigblk)
            n = hi - lo
            if n < 8
                vs[b] = NaN
            else
                resize!(dbuf, n)
                for k in 1:n
                    dbuf[k] = abs(col[lo + k] - col[lo + k - 1])
                end
                sort!(dbuf)
                m = isodd(n) ? dbuf[(n + 1) ÷ 2] : 0.5 * (dbuf[n ÷ 2] + dbuf[n ÷ 2 + 1])
                vs[b] = 1.4826 * m / sqrt(2)
            end
        end
        _interp_into!(view(M, :, j), cm, vm)
        _interp_into!(view(S, :, j), cs, vs)
    end
    M, S
end

"""Piecewise-linear interpolation of `v` sampled at `c` onto 1:length(out), flat
outside the sampled range."""
function _interp_into!(out, c::Vector{Float64}, v::Vector{Float64})
    n = length(out)
    nb = length(c)
    if nb == 1
        out .= v[1]
        return out
    end
    b = 1
    @inbounds for i in 1:n
        x = Float64(i)
        while b < nb - 1 && x > c[b + 1]
            b += 1
        end
        if x <= c[1]
            out[i] = v[1]
        elseif x >= c[nb]
            out[i] = v[nb]
        else
            t = (x - c[b]) / (c[b + 1] - c[b])
            out[i] = (1 - t) * v[b] + t * v[b + 1]
        end
    end
    out
end

"""
    zmap(A; medblk, sigblk, sigfloor_q=0.05) -> Z

Robust per-pixel deviation of `A` from its own smooth row-wise model, in units of the
local noise.  `Z` is NaN where the noise scale is undefined.  The noise scale is
floored at `sigfloor_q` of its own chip-level median so that a locally dead (perfectly
constant) row cannot generate infinite significance.
"""
function zmap(A::Matrix{Float64}; medblk::Int = 64, sigblk::Int = 256,
        sigfloor_q::Float64 = 0.05)
    M, S = detrend_x(A; medblk = medblk, sigblk = sigblk)
    sref = median(filter(x -> isfinite(x) && x > 0, vec(S)))
    floorv = sigfloor_q * sref
    Z = Array{Float64}(undef, size(A))
    @inbounds for i in eachindex(A)
        s = S[i]
        s = (isfinite(s) && s > floorv) ? s : floorv
        z = (A[i] - M[i]) / s
        Z[i] = isfinite(z) ? z : NaN
    end
    Z
end

# ---------------------------------------------------------------------------
# matched filter
# ---------------------------------------------------------------------------

struct Integral
    I::Matrix{Float64}
    C::Matrix{Float64}
end

"""
    Integral(Z; zclip=4.0)

Integral images of `Z` and of its finite-pixel count, with `Z` **winsorised** to
`±zclip` first.

The clip is essential, not cosmetic.  Measured on an APO chip-G `dark_rate` map the
raw residual reaches z ≈ 2.7e4 at isolated hot pixels and the 99.9th percentile is
z ≈ 3.9e3; unclipped, one such pixel dominates every box that contains it and smears a
spurious detection over its whole 33×33 neighbourhood — 30 % of the array lit up.
Winsorising caps any single pixel's contribution at 4, so a 33×33 box can gain at most
0.12 from it, while a coherent 20×14 block in which *every* pixel sits at z = −4 still
returns Σz/√N = −67.  That is precisely the discrimination the method exists to make:
isolated extreme pixels are AR's existing per-pixel cuts' business, coherent blocks
are this tool's.
"""
function Integral(Z::Matrix{Float64}; zclip::Float64 = 4.0)
    nx, ny = size(Z)
    I = zeros(Float64, nx + 1, ny + 1)
    C = zeros(Float64, nx + 1, ny + 1)
    @inbounds for j in 1:ny, i in 1:nx
        z = Z[i, j]
        ok = isfinite(z)
        zc = ok ? clamp(z, -zclip, zclip) : 0.0
        I[i + 1, j + 1] = zc + I[i, j + 1] + I[i + 1, j] - I[i, j]
        C[i + 1, j + 1] = (ok ? 1.0 : 0.0) + C[i, j + 1] + C[i + 1, j] - C[i, j]
    end
    Integral(I, C)
end

"""Normalised box sum, `Σz / √N`, over a `wx × wy` window centred on each pixel.

For white unit-variance noise this has unit variance at every scale, so one threshold
serves all scales; for a coherent block of mean offset `μ` it grows as `μ√N`.  That
growth is exactly why this finds what a per-pixel cut misses.
"""
function box_significance(itg::Integral, nx::Int, ny::Int, wx::Int, wy::Int)
    hx = (wx - 1) ÷ 2
    hy = (wy - 1) ÷ 2
    I = itg.I
    C = itg.C
    S = Array{Float64}(undef, nx, ny)
    @inbounds for j in 1:ny
        y0 = max(1, j - hy)
        y1 = min(ny, j + hy)
        for i in 1:nx
            x0 = max(1, i - hx)
            x1 = min(nx, i + hx)
            s = I[x1 + 1, y1 + 1] - I[x0, y1 + 1] - I[x1 + 1, y0] + I[x0, y0]
            c = C[x1 + 1, y1 + 1] - C[x0, y1 + 1] - C[x1 + 1, y0] + C[x0, y0]
            S[i, j] = c > 0 ? s / sqrt(c) : 0.0
        end
    end
    S
end

"""
    standardize_rows!(S; stride=4)

Replace each row of `S` by `(S - median)/(1.4826*MAD)` of that row.

This is the step that makes the method work on real detector data.  The white-noise
normalisation `Σz/√N` assumes the detrended residual is independent pixel to pixel; it
is not (inter-pixel capacitance correlates neighbours, and the block-median trend model
leaves low-level curvature that is coherent across rows).  Measured on an APO chip-G
dome flat, the raw `Σz/√N` exceeds 6 over 59 % of the array — i.e. the analytic
threshold is meaningless.  Standardising each scale's filter output against *the rest
of the same row* asks the only question that is actually well posed: is this box sum an
outlier with respect to the same statistic everywhere else at the same scale and the
same illumination level?  A defect occupying a few percent of a row cannot bias the
row's own median/MAD, so this does not suppress the signal.

Statistics are estimated from every `stride`-th pixel (≈510 samples per row), which is
ample and keeps the sweep cheap.
"""
function standardize_rows!(S::Matrix{Float64}; stride::Int = 4)
    nx, ny = size(S)
    samp = Vector{Float64}(undef, length(1:stride:nx))
    dev = similar(samp)
    @inbounds for j in 1:ny
        k = 0
        for i in 1:stride:nx
            k += 1
            samp[k] = S[i, j]
        end
        m = _median!(dev, view(samp, 1:k))
        for t in 1:k
            dev[t] = abs(samp[t] - m)
        end
        sd = 1.4826 * _median!(samp, view(dev, 1:k))
        sd = (isfinite(sd) && sd > 0) ? sd : 1.0
        for i in 1:nx
            S[i, j] = (S[i, j] - m) / sd
        end
    end
    S
end

"""
    highpass(Z; zclip=4.0, hp=129) -> Zhp

Winsorise `Z` and subtract its own `hp × hp` local mean.

Both detectors carry genuine coherent structure on scales of hundreds of pixels —
illumination roll-off near the array edges, broad dark-current gradients.  It is real,
but it is not a *defect region*, and left in place it dominates the ranking: without
this step the twelve largest regions across the six telescope-chips are all
array-spanning strips of 40000-150000 px, ten of them the roll-off in the first or last
~250 columns.

`hp` must be comfortably larger than the largest matched-filter scale (65 here), so a
defect the filter bank can actually resolve loses only a few percent of its amplitude:
the 460-pixel block recovered on APO chip G occupies 2.8 % of a 129 × 129 window, so
its matched-filter output drops by 2.8 %.  Structure much broader than `hp` is removed by
construction — that is the point, and it is the reason the scale ladder stops at 65.
"""
function highpass(Z::Matrix{Float64}; zclip::Float64 = 4.0, hp::Int = 129)
    nx, ny = size(Z)
    itg = Integral(Z; zclip = zclip)
    I = itg.I
    C = itg.C
    h = (hp - 1) ÷ 2
    O = Array{Float64}(undef, nx, ny)
    @inbounds for j in 1:ny
        y0 = max(1, j - h); y1 = min(ny, j + h)
        for i in 1:nx
            x0 = max(1, i - h); x1 = min(nx, i + h)
            s = I[x1 + 1, y1 + 1] - I[x0, y1 + 1] - I[x1 + 1, y0] + I[x0, y0]
            c = C[x1 + 1, y1 + 1] - C[x0, y1 + 1] - C[x1 + 1, y0] + C[x0, y0]
            z = Z[i, j]
            O[i, j] = isfinite(z) ? clamp(z, -zclip, zclip) - (c > 0 ? s / c : 0.0) : NaN
        end
    end
    O
end

"""
    matched_filter_detect(Z; thresh=6.0, scales=DEFAULT_SCALES, standardize=true,
                          zclip=4.0, hp=129) -> (det, smax)

`det[i,j]` is true when at least one scale's row-standardised filter output reaches
|·| ≥ `thresh` at that pixel; `smax` carries the signed most-extreme value over scales,
for reporting.  Pass `hp = 0` to skip the high-pass.
"""
function matched_filter_detect(Z::Matrix{Float64}; thresh::Float64 = 6.0,
        scales = DEFAULT_SCALES, standardize::Bool = true, zclip::Float64 = 4.0,
        hp::Int = 129)
    nx, ny = size(Z)
    Zw = hp > 0 ? highpass(Z; zclip = zclip, hp = hp) : Z
    itg = Integral(Zw; zclip = zclip)
    det = falses(nx, ny)
    smax = zeros(Float64, nx, ny)
    for (wx, wy) in scales
        S = box_significance(itg, nx, ny, wx, wy)
        standardize && standardize_rows!(S)
        @inbounds for k in eachindex(S)
            s = S[k]
            if abs(s) >= thresh
                det[k] = true
            end
            if abs(s) > abs(smax[k])
                smax[k] = s
            end
        end
    end
    det, smax
end

# ---------------------------------------------------------------------------
# morphology / connected components
# ---------------------------------------------------------------------------

"""3x3 binary dilation."""
function dilate1(M::BitMatrix)
    nx, ny = size(M)
    O = falses(nx, ny)
    @inbounds for j in 1:ny, i in 1:nx
        M[i, j] || continue
        for dj in -1:1, di in -1:1
            p = i + di
            q = j + dj
            (1 <= p <= nx && 1 <= q <= ny) && (O[p, q] = true)
        end
    end
    O
end

"""8-connectivity labelling. Returns `(labels, nlabels)`."""
function label_components(M::BitMatrix)
    nx, ny = size(M)
    lab = zeros(Int32, nx, ny)
    n = 0
    stack = Tuple{Int, Int}[]
    @inbounds for j in 1:ny, i in 1:nx
        (M[i, j] && lab[i, j] == 0) || continue
        n += 1
        lab[i, j] = n
        push!(stack, (i, j))
        while !isempty(stack)
            (a, b) = pop!(stack)
            for dj in -1:1, di in -1:1
                (di == 0 && dj == 0) && continue
                p = a + di
                q = b + dj
                (1 <= p <= nx && 1 <= q <= ny) || continue
                if M[p, q] && lab[p, q] == 0
                    lab[p, q] = n
                    push!(stack, (p, q))
                end
            end
        end
    end
    lab, n
end

struct RegionBox
    id::Int
    x0::Int   # detector coordinates, inclusive
    x1::Int
    y0::Int
    y1::Int
    area::Int          # detected pixels
    pixels::Vector{Tuple{Int, Int}}   # science coords
end

"""
    regions_from_mask(M; min_area=4, close=false) -> Vector{RegionBox}

Connected regions of `M` (science coords) returned in detector coords.

`close=true` labels the 1-pixel dilation instead, so that a defect broken into
touching-but-not-adjacent fragments is reported as one region (the recorded pixels are
still the undilated ones).  It defaults to **off**: at the few-percent detection
densities these products actually produce, a 1-pixel closure chains unrelated
detections into array-spanning blobs — on an APO chip-G flat it merged a 115-column
edge strip and everything near it into a single 92 000-pixel "region".
"""
function regions_from_mask(M::BitMatrix; min_area::Int = 4, close::Bool = false)
    L, n = label_components(close ? dilate1(M) : M)
    px = [Tuple{Int, Int}[] for _ in 1:n]
    nx, ny = size(M)
    @inbounds for j in 1:ny, i in 1:nx
        M[i, j] || continue
        l = L[i, j]
        l > 0 && push!(px[l], (i, j))
    end
    out = RegionBox[]
    for l in 1:n
        p = px[l]
        length(p) >= min_area || continue
        xs = first.(p)
        ys = last.(p)
        push!(out,
            RegionBox(length(out) + 1,
                minimum(xs) + SCI_OFF, maximum(xs) + SCI_OFF,
                minimum(ys) + SCI_OFF, maximum(ys) + SCI_OFF,
                length(p), p))
    end
    sort!(out, by = r -> -r.area)
    [RegionBox(i, r.x0, r.x1, r.y0, r.y1, r.area, r.pixels) for (i, r) in enumerate(out)]
end

end # module
