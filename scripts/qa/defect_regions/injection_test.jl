# Injection–recovery test: what can this finder actually detect?
#
# The acceptance test proves the method finds the one block we already knew about.
# That is necessary but it is not a characterisation.  Here synthetic coherent blocks
# of known size and amplitude are added to a real product, at locations the finder
# already calls clean, and we measure what fraction come back.  The result is a
# completeness surface in (block size, per-pixel amplitude) — i.e. exactly the
# statement "a defect smaller/fainter than this would have been missed", which is what
# a mask proposal has to be able to say about itself.
#
#   AR_APRED=<outdir>/apred AR_QA_OUT=<workdir> julia --project=. injection_test.jl

include("DefectFinder.jl")
using .DefectFinder
using HDF5, Statistics, Printf, Random

const APRED = get(ENV, "AR_APRED") do
    error("set AR_APRED to a reduction's apred directory, e.g. <outdir>/apred")
end
const OUT = get(ENV, "AR_QA_OUT", pwd())
const NPER = parse(Int, get(ENV, "NPER", "40"))     # injections per (size, amplitude)
const SIZES = [(4, 4), (8, 8), (14, 10), (20, 14), (32, 32), (64, 4)]
const AMPS = [-0.25, -0.5, -1.0, -2.0, -4.0]

io = open(joinpath(OUT, "INJECTION_TEST.txt"), "w")
P(a...) = (println(io, a...); println(a...); flush(io))

mjds() = sort(parse.(Int, filter(x -> occursin(r"^\d+$", x), readdir(APRED))))
function findexp(mjd, tele, chip, kind)
    d = joinpath(APRED, string(mjd))
    isdir(d) || return nothing
    c = sort(filter(x -> startswith(x, "ar2Dcal_$(tele)_$(mjd)_") &&
                        endswith(x, "_$(chip)_$(kind).h5"), readdir(d)))
    isempty(c) ? nothing : joinpath(d, c[1])
end
function find_static(tele, chip, pre)
    for m in reverse(mjds())
        d = joinpath(APRED, string(m))
        isdir(d) || continue
        c = sort(filter(x -> startswith(x, "$(pre)_$(tele)_$(chip)_"), readdir(d)))
        isempty(c) || return joinpath(d, c[1])
    end
    nothing
end
readsci(p, k) = h5open(p) do h
    d = h[k]; size(d, 1) >= 2048 ? Float64.(d[5:2044, 5:2044]) : Float64.(d[:, :])
end

P("Injection–recovery test for the coherent-defect region finder.")
P("Blocks of `amp × (local noise scale)` are added to a real product at random")
P("locations that the finder calls clean (no detection within 12 px), then the whole")
P("detection chain is re-run and we ask whether a reported region overlaps ≥ 25 % of")
P("the injected footprint.  $(NPER) injections per cell, all injected simultaneously but")
P("kept ≥ 80 px apart so they cannot help each other.")
P("")

targets = [("apo", "G", "dome-flat dimage", findexp(mjds()[end], "apo", "G", "domeflat"), "dimage"),
    ("apo", "G", "dark_rate", find_static("apo", "G", "darkRate"), "dark_rate"),
    ("lco", "R", "dome-flat dimage", findexp(mjds()[end], "lco", "R", "domeflat"), "dimage")]

for (tele, chip, label, path, key) in targets
    path === nothing && continue
    A0 = readsci(path, key)
    _, S = DefectFinder.detrend_x(A0)
    Z0 = zmap(A0)
    det0, _ = matched_filter_detect(Z0)
    clean = .!DefectFinder.dilate1(DefectFinder.dilate1(det0))
    P("### $tele $chip — $label  ($(basename(path)))")
    P("    baseline: $(round(100*count(det0)/length(det0), digits=3)) % of pixels detected")
    PF = (fmt, a...) -> (s = Printf.format(Printf.Format(fmt), a...); print(io, s); print(s))
    PF("    %-10s", "size \\ amp")
    for a in AMPS
        PF("%9.2f", a)
    end
    PF("\n")
    for (bx, by) in SIZES
        PF("    %-10s", "$(bx)x$(by)")
        for amp in AMPS
            rng = MersenneTwister(hash((tele, chip, key, bx, by, amp)) % 2^30)
            A = copy(A0)
            spots = Tuple{Int, Int}[]
            tries = 0
            while length(spots) < NPER && tries < 20000
                tries += 1
                i = rand(rng, 120:(DefectFinder.SCI_N - 120 - bx))
                j = rand(rng, 120:(DefectFinder.SCI_N - 120 - by))
                all(abs(i - s[1]) > 80 || abs(j - s[2]) > 80 for s in spots) || continue
                all(clean[i:(i + bx - 1), j:(j + by - 1)]) || continue
                push!(spots, (i, j))
            end
            isempty(spots) && (PF("%9s", "-"); continue)
            for (i, j) in spots
                @views A[i:(i + bx - 1), j:(j + by - 1)] .+= amp .* S[i:(i + bx - 1), j:(j + by - 1)]
            end
            Z = zmap(A)
            det, _ = matched_filter_detect(Z)
            new = det .& .!det0
            nrec = 0
            for (i, j) in spots
                mean(@view new[i:(i + bx - 1), j:(j + by - 1)]) >= 0.25 && (nrec += 1)
            end
            PF("%8.2f ", nrec / length(spots))
        end
        PF("\n")
    end
    P("")
end
P("Reading the table: the entry is the recovered fraction.  `amp` is the per-pixel")
P("offset in units of the local robust noise, so amp = -0.5 means every pixel in the")
P("block sits half a sigma low — far below anything a per-pixel sigma clip can reach,")
P("and the regime the real APO chip-G block lives in.")
close(io)
