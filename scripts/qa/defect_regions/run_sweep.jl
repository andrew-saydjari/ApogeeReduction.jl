# Coherent-defect-region sweep over all chips and both telescopes.
#
#   AR_APRED=<outdir>/apred AR_QA_OUT=<workdir> julia --project=. run_sweep.jl            # full sweep
#   NMJD=4 AR_APRED=<outdir>/apred AR_QA_OUT=<workdir> julia --project=. run_sweep.jl     # quick smoke sweep
#
# Reads only; writes regions.jls and masks_<tele>_<chip>.h5 into AR_QA_OUT.
# Nothing in the reduction is touched.

include("DefectFinder.jl")
include("types.jl")
using .DefectFinder
using HDF5, Statistics, Printf, Serialization

const APRED = get(ENV, "AR_APRED") do
    error("set AR_APRED to a reduction's apred directory, e.g. <outdir>/apred")
end
const OUT = get(ENV, "AR_QA_OUT", pwd())
const NMJD = parse(Int, get(ENV, "NMJD", "14"))
const THRESH = 6.0
const MARGINAL = 3.0
const MIN_AREA = 6
const BAD_PIX_BITS = 24566        # AR/src/utils.jl, verified 2026-09-08

log_io = open(joinpath(OUT, "sweep.log"), "a")
P(a...) = (println(log_io, a...); flush(log_io); println(a...); flush(stdout))

mjds() = sort(parse.(Int, filter(x -> occursin(r"^\d+$", x), readdir(APRED))))

"""First matching exposure file of the given kind, or nothing."""
function findexp(mjd, tele, chip, kind)
    d = joinpath(APRED, string(mjd))
    isdir(d) || return nothing
    pre = "ar2Dcal_$(tele)_$(mjd)_"
    suf = "_$(chip)_$(kind).h5"
    c = sort(filter(x -> startswith(x, pre) && endswith(x, suf), readdir(d)))
    isempty(c) ? nothing : joinpath(d, c[1])
end

function find_static(tele, chip, pre, key)
    for m in reverse(mjds())
        d = joinpath(APRED, string(m))
        isdir(d) || continue
        c = sort(filter(x -> startswith(x, "$(pre)_$(tele)_$(chip)_"), readdir(d)))
        isempty(c) || return joinpath(d, c[1])
    end
    nothing
end

readsci(path, key) = h5open(path) do h
    d = h[key]
    size(d, 1) >= 2048 ? Float64.(d[5:2044, 5:2044]) : Float64.(d[:, :])
end
readsci_int(path, key) = h5open(path) do h
    d = h[key]
    size(d, 1) >= 2048 ? Int64.(d[5:2044, 5:2044]) : Int64.(d[:, :])
end

"""Evenly spaced subset of `v` of length ≤ n, always keeping the endpoints."""
function spread(v, n)
    length(v) <= n && return collect(v)
    idx = unique(round.(Int, range(1, length(v), length = n)))
    collect(v)[idx]
end

# --------------------------------------------------------------------------
# per-(telescope, chip) sweep
# --------------------------------------------------------------------------

function sweep(tele, chip, mjdlist)
    P("--- $tele $chip : $(length(mjdlist)) MJDs ---")
    nx = ny = DefectFinder.SCI_N

    # ---- static products -------------------------------------------------
    static = Dict{String, NamedTuple}()
    dkf = find_static(tele, chip, "darkRate", "dark_rate")
    if dkf !== nothing
        A = readsci(dkf, "dark_rate")
        Z = zmap(A)
        d, s = matched_filter_detect(Z; thresh = THRESH)
        static["dark_rate"] = (Z = Z, det = d, smax = s, file = basename(dkf))
        P("    dark_rate  $(basename(dkf))  det=$(count(d)) px")
    end
    ffn = find_static(tele, chip, "flatFraction", "flat_im")
    if ffn !== nothing
        A = readsci(ffn, "flat_im")
        Z = zmap(A)
        d, s = matched_filter_detect(Z; thresh = THRESH)
        static["flat_im"] = (Z = Z, det = d, smax = s, file = basename(ffn))
        P("    flat_im    $(basename(ffn))  det=$(count(d)) px")
    end

    # ---- 2D flat-frame channel, epoch by epoch ---------------------------
    epochs = Tuple{Int, String}[]         # (mjd, kind)
    epdet = BitMatrix[]
    sumz = zeros(Float32, nx, ny)
    maxS = zeros(Float32, nx, ny)
    nz = 0
    for m in mjdlist, kind in ("domeflat", "quartzflat")
        f = findexp(m, tele, chip, kind)
        f === nothing && continue
        A = try
            readsci(f, "dimage")
        catch e
            P("    !! read failed $(basename(f)): $e"); continue
        end
        Z = zmap(A)
        d, s = matched_filter_detect(Z; thresh = THRESH)
        push!(epochs, (m, kind))
        push!(epdet, d)
        @inbounds for k in eachindex(Z)
            v = Z[k]
            isfinite(v) && (sumz[k] += Float32(clamp(v, -20, 20)))
            abs(s[k]) > abs(maxS[k]) && (maxS[k] = Float32(s[k]))
        end
        nz += 1
    end
    P("    2D flat channel: $(length(epochs)) epoch-frames")
    length(epochs) == 0 && nz == 0 && isempty(static) && return Region[], nothing
    meanz2d = nz > 0 ? sumz ./ nz : sumz
    cnt2d = zeros(Int16, nx, ny)
    for d in epdet
        @inbounds for k in eachindex(d)
            d[k] && (cnt2d[k] += Int16(1))
        end
    end

    # ---- union mask ------------------------------------------------------
    # A 2D-channel pixel must fire in >= 2 independent frames: a cosmic ray or a
    # one-off artefact cannot survive that, a detector defect always does.
    U = falses(nx, ny)
    haskey(static, "dark_rate") && (U .|= static["dark_rate"].det)
    haskey(static, "flat_im") && (U .|= static["flat_im"].det)
    U .|= (cnt2d .>= 2)
    regs = regions_from_mask(U; min_area = MIN_AREA)
    P("    union mask: $(count(U)) px, $(length(regs)) regions (min_area=$MIN_AREA)")

    # ---- AR's own view ---------------------------------------------------
    arbad = falses(nx, ny)
    repfile = nothing
    for m in reverse(mjdlist)
        f = findexp(m, tele, chip, "domeflat")
        f === nothing && (f = findexp(m, tele, chip, "quartzflat"))
        f === nothing && continue
        repfile = f
        B = readsci_int(f, "pix_bitmask")
        arbad .|= (B .& BAD_PIX_BITS) .!= 0
        break
    end
    P("    AR bad_pix_bits coverage of the science area: $(round(100*count(arbad)/length(arbad),digits=3)) %")

    # ---- trace centres ---------------------------------------------------
    tc = nothing
    for m in reverse(mjdlist)
        d = joinpath(APRED, string(m))
        c = sort(filter(
            x -> startswith(x, "ar1Dcal_$(tele)_$(m)_") && endswith(x, "_$(chip)_domeflat.h5"),
            readdir(d)))
        isempty(c) && (c = sort(filter(
            x -> startswith(x, "ar1Dcal_$(tele)_$(m)_") && endswith(x, "_$(chip)_quartzflat.h5"),
            readdir(d))))
        isempty(c) && continue
        tc = h5open(joinpath(d, c[1])) do h
            h["extract_trace_centers"][:, :]
        end
        break
    end

    out = Region[]
    for r in regs
        ev = ProductEvidence[]
        for nm in ("dark_rate", "flat_im")
            haskey(static, nm) || continue
            st = static[nm]
            zs = [st.Z[p[1], p[2]] for p in r.pixels]
            ss = [st.smax[p[1], p[2]] for p in r.pixels]
            df = mean(Bool[st.det[p[1], p[2]] for p in r.pixels])
            push!(ev,
                ProductEvidence(nm, median(filter(isfinite, zs)),
                    ss[argmax(abs.(ss))], df))
        end
        if nz > 0
            zs = [Float64(meanz2d[p[1], p[2]]) for p in r.pixels]
            ss = [Float64(maxS[p[1], p[2]]) for p in r.pixels]
            df = mean([cnt2d[p[1], p[2]] >= 2 for p in r.pixels])
            push!(ev, ProductEvidence("flat2d", median(filter(isfinite, zs)),
                ss[argmax(abs.(ss))], df))
        end
        # epoch membership: an epoch "has" the region if it detects >= 25 % of it
        eflags = [mean(Bool[d[p[1], p[2]] for p in r.pixels]) >= 0.25 for d in epdet]
        got = findall(eflags)
        mlo = isempty(got) ? -1 : minimum(epochs[i][1] for i in got)
        mhi = isempty(got) ? -1 : maximum(epochs[i][1] for i in got)
        af = mean(Bool[arbad[p[1], p[2]] for p in r.pixels])
        fibs = Int[]
        if tc !== nothing
            xs = clamp(r.x0, 1, 2048):clamp(r.x1, 1, 2048)
            for f in 1:size(tc, 2)
                any(x -> (r.y0 - 1) <= tc[x, f] <= (r.y1 + 1), xs) && push!(fibs, f)
            end
        end
        push!(out,
            Region(tele, chip, r.id, r.x0, r.x1, r.y0, r.y1, r.area, ev,
                count(eflags), length(eflags), mlo, mhi, eflags, af, fibs, r.pixels))
    end

    # ---- persist the masks and the evidence maps -------------------------
    mfile = joinpath(OUT, "masks_$(tele)_$(chip).h5")
    h5open(mfile, "w") do h
        h["union_mask"] = UInt8.(U)
        h["cnt2d"] = cnt2d
        h["meanz_2dflat"] = meanz2d
        h["maxS_2dflat"] = maxS
        haskey(static, "dark_rate") && (h["det_dark"] = UInt8.(static["dark_rate"].det))
        haskey(static, "flat_im") && (h["det_flat"] = UInt8.(static["flat_im"].det))
        haskey(static, "dark_rate") && (h["smax_dark"] = Float32.(static["dark_rate"].smax))
        haskey(static, "flat_im") && (h["smax_flat"] = Float32.(static["flat_im"].smax))
        h["ar_badpix"] = UInt8.(arbad)
        h["epoch_mjd"] = [e[1] for e in epochs]
        h["epoch_kind"] = [e[2] for e in epochs]
        h["nepoch_frames"] = length(epochs)
    end
    P("    wrote $(basename(mfile))")
    out, (epochs = epochs, repfile = repfile)
end

# --------------------------------------------------------------------------

allm = mjds()
allregions = Region[]
meta = Dict{String, Any}()
for tele in ("apo", "lco")
    have = [m for m in allm
            if findexp(m, tele, "G", "domeflat") !== nothing &&
               findexp(m, tele, "G", "quartzflat") !== nothing]
    ml = spread(have, NMJD)
    P("=== $tele : $(length(have)) MJDs with dome+quartz on G; using $(length(ml)): $ml")
    meta["mjds_$tele"] = ml
    for chip in ("R", "G", "B")
        t0 = time()
        rs, _ = sweep(tele, chip, ml)
        append!(allregions, rs)
        P("    done in $(round(time()-t0,digits=1)) s\n")
    end
end
serialize(joinpath(OUT, "regions.jls"), (regions = allregions, meta = meta,
    thresh = THRESH, marginal = MARGINAL, min_area = MIN_AREA,
    scales = DEFAULT_SCALES, nmjd = NMJD))
P("TOTAL regions: $(length(allregions))")
close(log_io)
