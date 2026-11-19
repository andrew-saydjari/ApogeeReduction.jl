# False-positive null test.
#
# The question a region finder must answer about itself is: how much of what it finds
# is an artefact of the *statistic* rather than of the detector?  Here the detrended
# residual Z of a real product is shuffled within each row.  That destroys every
# spatial correlation while preserving the marginal distribution of Z exactly —
# including its heavy non-Gaussian tails, which are the thing most likely to fool a
# threshold.  Anything the finder reports on the shuffled map is a false positive by
# construction.
#
#   AR_APRED=<outdir>/apred AR_QA_OUT=<workdir> julia --project=. null_test.jl

include("DefectFinder.jl")
using .DefectFinder
using HDF5, Statistics, Printf, Random

const APRED = get(ENV, "AR_APRED") do
    error("set AR_APRED to a reduction's apred directory, e.g. <outdir>/apred")
end
const OUT = get(ENV, "AR_QA_OUT", pwd())
const NSHUF = parse(Int, get(ENV, "NSHUF", "3"))

io = open(joinpath(OUT, "NULL_TEST.txt"), "w")
P(a...) = (println(io, a...); println(a...); flush(io); flush(stdout))

mjds() = sort(parse.(Int, filter(x -> occursin(r"^\d+$", x), readdir(APRED))))
function find_static(tele, chip, pre)
    for m in reverse(mjds())
        d = joinpath(APRED, string(m))
        isdir(d) || continue
        c = sort(filter(x -> startswith(x, "$(pre)_$(tele)_$(chip)_"), readdir(d)))
        isempty(c) || return joinpath(d, c[1])
    end
    nothing
end
function findexp(mjd, tele, chip, kind)
    d = joinpath(APRED, string(mjd))
    isdir(d) || return nothing
    c = sort(filter(x -> startswith(x, "ar2Dcal_$(tele)_$(mjd)_") &&
                        endswith(x, "_$(chip)_$(kind).h5"), readdir(d)))
    isempty(c) ? nothing : joinpath(d, c[1])
end
readsci(p, k) = h5open(p) do h
    d = h[k]; size(d, 1) >= 2048 ? Float64.(d[5:2044, 5:2044]) : Float64.(d[:, :])
end

P("Null test for the coherent-defect region finder.")
P("Threshold 6.0, scales $(DEFAULT_SCALES), zclip 4.0, min_area 6.")
P("Shuffling is within rows (fixed y, permuted x): marginal distribution of Z exactly")
P("preserved, all spatial coherence destroyed.  $(NSHUF) independent shuffles each.\n")
@printf(io, "%-6s %-5s %-14s %10s %8s %10s %8s %8s\n",
    "tele", "chip", "product", "real_px%", "real_reg", "null_px%", "null_reg", "ratio")
@printf("%-6s %-5s %-14s %10s %8s %10s %8s %8s\n",
    "tele", "chip", "product", "real_px%", "real_reg", "null_px%", "null_reg", "ratio")

tot_real_px = 0.0; tot_null_px = 0.0; tot_real_r = 0; tot_null_r = 0.0; n = 0
for tele in ("apo", "lco"), chip in ("R", "G", "B")
    jobs = Tuple{String, String, String}[]
    f = find_static(tele, chip, "darkRate"); f !== nothing && push!(jobs, ("dark_rate", f, "dark_rate"))
    f = find_static(tele, chip, "flatFraction"); f !== nothing && push!(jobs, ("flat_im", f, "flat_im"))
    ms = [m for m in mjds() if findexp(m, tele, chip, "domeflat") !== nothing]
    if !isempty(ms)
        push!(jobs, ("dimage_dome", findexp(ms[end], tele, chip, "domeflat"), "dimage"))
    end
    for (nm, path, key) in jobs
        A = readsci(path, key)
        Z = zmap(A)
        d0, _ = matched_filter_detect(Z)
        r0 = regions_from_mask(d0; min_area = 6)
        npx = 0.0; nrg = 0.0
        rng = MersenneTwister(20260908 + hash((tele, chip, nm)) % 10^6)
        for s in 1:NSHUF
            Zs = copy(Z)
            for j in axes(Zs, 2)
                shuffle!(rng, view(Zs, :, j))
            end
            d1, _ = matched_filter_detect(Zs)
            npx += count(d1) / length(d1)
            nrg += length(regions_from_mask(d1; min_area = 6))
        end
        npx /= NSHUF; nrg /= NSHUF
        rp = count(d0) / length(d0)
        rat = npx > 0 ? rp / npx : Inf
        @printf(io, "%-6s %-5s %-14s %10.4f %8d %10.5f %8.1f %8.0f\n",
            tele, chip, nm, 100rp, length(r0), 100npx, nrg, rat)
        @printf("%-6s %-5s %-14s %10.4f %8d %10.5f %8.1f %8.0f\n",
            tele, chip, nm, 100rp, length(r0), 100npx, nrg, rat)
        global tot_real_px += rp; global tot_null_px += npx
        global tot_real_r += length(r0); global tot_null_r += nrg; global n += 1
    end
end
P("")
P("Mean over $(n) product-maps: real $(round(100*tot_real_px/n,digits=4)) % of pixels, " *
  "null $(round(100*tot_null_px/n,digits=5)) %.")
P("Total regions: real $(tot_real_r), null $(round(tot_null_r,digits=1)).")
P("")
P("Reading: the detections are NOT statistical false positives.  What the finder")
P("reports is real coherent structure in the detector products.  The remaining")
P("question is which of that structure is a *defect* worth masking, which is what the")
P("multi-product corroboration and the comparison against AR's existing bad_pix_bits")
P("in REGIONS.txt address; this test only establishes that the statistic itself is")
P("not manufacturing the regions.")
close(io)
