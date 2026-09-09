## Census helper for arm_census.sh: tally arMADGICS per-spectrum flag columns
## straight out of the batch OUTPUT PRODUCTS, not out of the job log.
##
## WHY NOT THE LOG. The log route has four failure modes that all silently
## understate or overstate the truth, and three of them have already bitten:
##   - the sky-guard verdict is EXPOSURE-level but was printed once per TARGET
##     FIBER, so raw line counts overstated it ~130x (job 7001233: 479,570
##     lines, 3,676 distinct exposures). arMADGICS now suppresses that print
##     entirely, so grepping for it would report a silent ZERO;
##   - the AR-side census appends to the same file it greps, so on a resumed run
##     every pattern matched its own previous output (fixed separately, but the
##     class of bug is inherent to logs);
##   - a rotated, truncated or redirected log undercounts and looks exactly like
##     a clean run;
##   - any reword of a println silently zeroes a grep-based census.
## `ingestBit` and `skyBit` are already per-spectrum columns in the products.
## Read those. The log could only ever show FATAL ingestBit values anyway
## ("Skipping spectrum"); the informational bits were invisible.
##
## Usage:  julia --project=<AR base_dir> arm_census.jl <arm-raw-dir>
## Env:    ARM_CENSUS_PROCS  worker processes (default min(8, CPU threads))
##
## Prints "arM census: ..." lines on stdout. Exit 0 on success; exit 1 only when
## it could not count at all, having first said NOT COUNTED. It never reports a
## silent zero -- that failure mode is the reason this census exists.
##
## n.b. worker PROCESSES, not threads. HDF5.jl serialises its API behind a
## global lock, so a threaded version of this loop pins one core and did not
## finish 16,523 files in 22 minutes. Processes actually parallelise the reads.

using Distributed, Printf

const COLUMNS = ["ingestBit", "skyBit"]

# ---- arMADGICS ingestBit bit table (src/ingest.jl); fatal bits skip the solve.
# Bits 7/8 are the per-fiber throughput flag and are informational by default, so
# `ingestBit != 0` is NO LONGER equivalent to `RV_flag == 64`;
# `ingest_fatal(ingestBit)` is. Keep in sync with the arMADGICS README.
const INGEST_FATAL_BITS = 2^0 | 2^1 | 2^2 | 2^6
const INGEST_BIT_NAMES = Dict(
    0 => "runtime error (fatal)",
    1 => "flux all NaN/zero (fatal)",
    2 => "too few good pixels (fatal)",
    3 => "non-finite flux masked",
    4 => "non-finite/non-positive ivar masked",
    5 => "tiny-ivar pixels masked",
    6 => "starscale0 non-finite or <=0 (fatal)",
    7 => "AR says fiber throughput unusable",
    8 => "AR throughput flag ABSENT (unknown)")
const SKY_BIT_NAMES = Dict(
    0 => "sky fiber excluded upstream",
    1 => "sky fiber failed the z-cut",
    2 => "too few sky fibers (component skipped)",
    3 => "no sky fibers at all",
    4 => "kept fiber has non-positive median",
    5 => "non-finite sky decomposition dropped",
    6 => "sky fiber excluded by AR throughput flag")

say(s) = (println("arM census: ", s); flush(stdout))
fail(msg) = (say(msg); exit(1))

length(ARGS) == 1 ||
    fail("no arM raw directory given - arM diagnostics NOT counted (not zero, UNCOUNTED)")
raw = ARGS[1]
isdir(raw) ||
    fail("$raw is not a readable directory - arM diagnostics NOT counted (not zero, UNCOUNTED)")

# ---- enumerate the batch products -------------------------------------------
files = String[]
for f in 1:600
    d = joinpath(raw, lpad(f, 3, "0"))
    isdir(d) || continue
    for fn in readdir(d)
        endswith(fn, ".h5") && push!(files, joinpath(d, fn))
    end
end
isempty(files) &&
    fail("no batch products under $raw - arM diagnostics NOT counted (not zero, UNCOUNTED)")
say("$(length(files)) batch products to read")

# ---- batch_info gives (tele, mjd, expnum) per spectrum, in file order --------
# Format: linear_index, tele, mjd, expnum, adjfiberindx. `linear_index` is
# PER-FIBER and runs 1..n in file order, so a file named ..._batch_<start>.h5
# holding L rows covers that fiber's batch_info rows start .. start+L-1.
binfo = joinpath(raw, "batch_info.txt")
have_exposures = isfile(binfo)
fib_exp = [Tuple{UInt8, Int32, Int32}[] for _ in 1:600]
if have_exposures
    nbad = 0
    open(binfo) do io
        for ln in eachline(io)
            (isempty(ln) || startswith(ln, "#")) && continue
            p = split(ln, ",")
            length(p) < 5 && continue
            f = tryparse(Int, strip(p[5]))
            mj = tryparse(Int32, strip(p[3]))
            ex = tryparse(Int32, strip(p[4]))
            if f === nothing || mj === nothing || ex === nothing || !(1 <= f <= 600)
                nbad += 1
                continue
            end
            push!(fib_exp[f], (strip(p[2]) == "apo" ? UInt8(1) : UInt8(2), mj, ex))
        end
    end
    nbad > 0 && say("WARNING $nbad unparseable batch_info rows ignored")
else
    say("batch_info.txt missing - per-spectrum counts only, exposure counts NOT COUNTED")
end

# ---- read the two flag columns, in parallel ---------------------------------
nprocs_want = something(tryparse(Int, get(ENV, "ARM_CENSUS_PROCS", "")),
    min(8, Sys.CPU_THREADS))
nprocs_want = clamp(nprocs_want, 1, 32)
t0 = time()
if nprocs_want > 1
    addprocs(nprocs_want; exeflags = "--project=$(Base.active_project())")
end
@everywhere using HDF5

@everywhere function read_flags(fn)
    try
        m = match(r"fiber_(\d+)_batch_(\d+)\.h5$", basename(fn))
        m === nothing && return nothing
        fib = parse(Int, m.captures[1])
        st = parse(Int, m.captures[2])
        h5open(fn) do h
            (haskey(h, "ingestBit") && haskey(h, "skyBit")) || return nothing
            (fib, st, Int.(read(h["ingestBit"])), Int.(read(h["skyBit"])))
        end
    catch
        nothing
    end
end

results = pmap(read_flags, files; batch_size = 32)
elapsed = time() - t0
nprocs_want > 1 && rmprocs(workers())

nread = count(!isnothing, results)
nfail = length(files) - nread
nfail > 0 && say("WARNING $nfail of $(length(files)) batch products unreadable or " *
                 "missing the flag columns - those spectra NOT COUNTED")
nread == 0 &&
    fail("no batch product could be read - arM diagnostics NOT counted (not zero, UNCOUNTED)")

# ---- tally ------------------------------------------------------------------
# In a function, not at top level: a bare `for` at global scope puts every
# accumulator in soft scope and Julia treats the assignments as new locals.
function tally(results, fib_exp, have_exposures)
    ing_spec = Dict{Int, Int}()
    sky_spec = Dict{Int, Int}()
    ing_exp = Dict{Int, Set{Tuple{UInt8, Int32, Int32}}}()
    sky_exp = Dict{Int, Set{Tuple{UInt8, Int32, Int32}}}()
    allexp = Set{Tuple{UInt8, Int32, Int32}}()
    nspec = 0
    nmis = 0
    for r in results
        isnothing(r) && continue
        fib, st, ing, sky = r
        nspec += length(ing)
        for k in eachindex(ing)
            ing_spec[ing[k]] = get(ing_spec, ing[k], 0) + 1
            sky_spec[sky[k]] = get(sky_spec, sky[k], 0) + 1
            have_exposures || continue
            row = st + k - 1
            if 1 <= fib <= 600 && row <= length(fib_exp[fib])
                e = fib_exp[fib][row]
                push!(allexp, e)
                push!(get!(ing_exp, ing[k], Set{Tuple{UInt8, Int32, Int32}}()), e)
                push!(get!(sky_exp, sky[k], Set{Tuple{UInt8, Int32, Int32}}()), e)
            else
                nmis += 1
            end
        end
    end
    return ing_spec, sky_spec, ing_exp, sky_exp, allexp, nspec, nmis
end

ing_spec, sky_spec, ing_exp, sky_exp, allexp, nspec, nmis = tally(results, fib_exp, have_exposures)
nmis > 0 && say("WARNING $nmis spectra had no batch_info row - counted per-spectrum, " *
                "NOT COUNTED per-exposure")

say(@sprintf("read %d batch products, %d spectra%s in %.1f s with %d workers",
    nread, nspec,
    have_exposures ? ", $(length(allexp)) distinct exposures" : "",
    elapsed, nprocs_want))

function report(name, spec, expsets, names, fatalmask)
    say("--- $name ---")
    tot = sum(values(spec); init = 0)
    nz = sum(v for (k, v) in spec if k != 0; init = 0)
    say(@sprintf("%d of %d spectra have %s != 0 (%.5f)", nz, tot, name,
        tot == 0 ? 0.0 : nz / tot))
    if fatalmask != 0
        nf = sum(v for (k, v) in spec if (k & fatalmask) != 0; init = 0)
        say(@sprintf("%d of %d spectra are FATAL (never fitted, RV_flag=64) (%.5f)",
            nf, tot, tot == 0 ? 0.0 : nf / tot))
    end
    for k in sort(collect(keys(spec)))
        ne = haskey(expsets, k) ? @sprintf("%9d", length(expsets[k])) : "UNCOUNTED"
        bits = [b for b in 0:15 if (k >> b) & 1 == 1]
        lbl = isempty(bits) ? "no problems" :
              join([get(names, b, "bit $b") for b in bits], "; ")
        say(@sprintf("  %9d spectra  %s exposures  %s=%-4d  %s",
            spec[k], ne, name, k, lbl))
    end
end

report("ingestBit", ing_spec, ing_exp, INGEST_BIT_NAMES, INGEST_FATAL_BITS)
report("skyBit", sky_spec, sky_exp, SKY_BIT_NAMES, 0)
say("skyBit is EXPOSURE-level -- every spectrum of an exposure carries the same")
say("value, so the exposure column is the one with physical meaning; the spectrum")
say("column is inflated by the fiber multiplicity, which is exactly the distortion")
say("the old per-target-fiber log duplication produced.")
exit(0)
