#!/usr/bin/env julia
# h5diff_tree.jl — golden-diff comparison of two HDF5 output trees.
#
# Walks two directory trees (golden vs new), pairs files by relative path, and
# compares every HDF5 dataset and attribute. Per-dataset verdicts:
#   identical           : bit-for-bit equal (NaN == NaN, elementwise isequal)
#   within-rtol         : |a-b| <= atol + rtol*max(|a|,|b|) everywhere it differs
#                         (max abs/rel diff reported)
#   DIFFERS             : anything else, with localization (which indices along
#                         each axis differ; axes of length 300 are labeled
#                         "fiber", per the AR product schemas)
#   missing-in-{golden,new} : dataset or file present on only one side
#   ignored             : matched an ignore pattern (provenance metadata);
#                         still compared and its true status reported, but it
#                         never fails the diff
#
# Known AR product schemas (2026_05_01 outputs; see REFACTOR_PLAN.md v1 §0):
#   ar2D:       dimage/ivarimage/chisqimage/CRimage/last_unsaturated {2048,2560} + metadata group
#   ar2Dcal:    dimage/ivarimage/pix_bitmask {2048,2560}
#   ar1Dcal:    flux_1d/ivar_1d/mask_1d/... {300,2048}, relthrpt/bitmsk_relthrpt/fiberTypeList {300}
#   ar1Dunical: flux_1d/ivar_1d/mask_1d {300,8700}, extract_trace_coords {3,300,8700}, relthrpt {300,3}
#               (NOTE: mask_1d is a Bool GOOD mask here, an Int bad-bit mask in ar1Dcal)
#   traceMain:  trace_params/regularized_trace_params {3,300,2048}, trace_param_covs {3,3,300,2048}
#   waveCalFPI: chipWaveSoln {3,300,2048}, ditherParams {3,300}, linParams {5,300}, nlParams {4,300}
#   relFlux/domeFlux/quartzFlux: absthrpt/relthrpt/bitmsk_relthrpt {300}
#
# Usage:
#   julia --project=<this dir> h5diff_tree.jl <goldenDir> <newDir> \
#       [--out report.md] [--rtol 0] [--atol 0] [--ignore pat1,pat2] \
#       [--no-default-ignores] [--max-list 10]
#
# Exit code 0 when every non-ignored dataset is identical or within-rtol,
# 1 when anything DIFFERS or is missing, 2 on usage/IO errors.

using ArgParse
using HDF5
using Printf

# Provenance/metadata keys that legitimately differ between runs of the same
# code. Discovered by inspecting real 2026_05_01 outputs:
#   metadata/{git_branch, git_clean, git_commit}  — code provenance
#   metadata/trace_orig_param_fname               — absolute path into the run's outdir
# Patterns are matched as path suffixes against "group/dataset" paths and as
# "path@attrname" against attributes.
const DEFAULT_IGNORES = [
    "metadata/git_branch",
    "metadata/git_clean",
    "metadata/git_commit",
    "metadata/trace_orig_param_fname",
    "git_branch", "git_clean", "git_commit", # top-level variants (some products)
    # top-level dataset in ar1Dcal/ar2Dresidualscal embedding the absolute path
    # of the run's own outdir (always differs between two runs in different
    # outdirs, even for identical code+data)
    "trace_used_param_fname",
]

function parse_commandline()
    s = ArgParseSettings(prog = "h5diff_tree.jl")
    @add_arg_table s begin
        "goldenDir"
        help = "golden (baseline) output tree"
        required = true
        "newDir"
        help = "new (comparison) output tree"
        required = true
        "--out"
        help = "write the markdown report here (default: stdout only)"
        arg_type = String
        default = ""
        "--rtol"
        help = "relative tolerance for within-rtol verdict (0 = exact only)"
        arg_type = Float64
        default = 0.0
        "--atol"
        help = "absolute tolerance for within-rtol verdict"
        arg_type = Float64
        default = 0.0
        "--ignore"
        help = "comma-separated extra ignore path suffixes"
        arg_type = String
        default = ""
        "--no-default-ignores"
        help = "drop the built-in provenance ignore list"
        action = :store_true
        "--max-list"
        help = "max indices/paths listed per localization line"
        arg_type = Int
        default = 10
        "--glob-ext"
        help = "comma-separated file extensions to compare"
        arg_type = String
        default = ".h5,.hdf5,.jdat"
    end
    return parse_args(s)
end

# ---------------------------------------------------------------------------
# tree walking
# ---------------------------------------------------------------------------

"""Collect relpath => abspath for all files under `root` with an extension in
`exts`. Follows symlinks (AR outdirs symlink darkRate/flatFraction cals into
apred/<mjd>/)."""
function collect_files(root::String, exts::Vector{String})
    files = Dict{String, String}()
    root = rstrip(abspath(root), '/')
    for (dir, _, fnames) in walkdir(root; follow_symlinks = true)
        for f in fnames
            any(endswith(f, e) for e in exts) || continue
            ap = joinpath(dir, f)
            rp = relpath(ap, root)
            files[rp] = ap
        end
    end
    return files
end

# ---------------------------------------------------------------------------
# HDF5 object enumeration
# ---------------------------------------------------------------------------

"""Return (datasets, attrs) where datasets :: Vector{String} of all dataset
paths and attrs :: Dict{String,Any} keyed by "objpath@attrname"."""
function enumerate_h5(fid::HDF5.File)
    dsets = String[]
    attrs_out = Dict{String, Any}()
    function collect_attrs(obj, path)
        for aname in keys(HDF5.attrs(obj))
            attrs_out["$(path)@$(aname)"] = try
                read(HDF5.attrs(obj)[aname])
            catch e
                "<<unreadable attr: $(typeof(e))>>"
            end
        end
    end
    function walk(g, prefix)
        collect_attrs(g, prefix == "" ? "/" : prefix)
        for name in keys(g)
            path = prefix == "" ? name : "$(prefix)/$(name)"
            obj = g[name]
            if obj isa HDF5.Group
                walk(obj, path)
            elseif obj isa HDF5.Dataset
                push!(dsets, path)
                collect_attrs(obj, path)
            end
            close(obj)
        end
    end
    walk(fid, "")
    return dsets, attrs_out
end

ignored(path::String, patterns::Vector{String}) =
    any(p -> endswith(path, p), patterns)

# ---------------------------------------------------------------------------
# dataset comparison
# ---------------------------------------------------------------------------

struct DsetDiff
    verdict::Symbol       # :identical, :within_rtol, :differs, :missing_in_golden, :missing_in_new
    detail::String        # human-readable localization / stats
end

axis_label(n::Int) = n == 300 ? "fiber" : (n == 8700 ? "unipix" : (n == 2048 ? "xpix" : "axis"))

"""Summarize the set of differing CartesianIndices along each dimension."""
function localize(idxs, dims::Dims, maxlist::Int)
    parts = String[]
    nd = length(dims)
    for d in 1:nd
        vals = sort(unique(getindex.(Tuple.(idxs), d)))
        lab = axis_label(dims[d])
        shown = length(vals) <= maxlist ? join(vals, ",") :
                join(vals[1:maxlist], ",") * ",…"
        push!(parts, "dim$(d)[$(lab), n=$(dims[d])]: $(length(vals)) distinct ($(shown))")
    end
    return join(parts, "; ")
end

isnumericarray(x) = x isa AbstractArray{<:Number} || x isa Number

function compare_dataset(a, b, rtol::Float64, atol::Float64, maxlist::Int)
    # exact match first (isequal: NaN==NaN true, works for strings/any type)
    if isequal(a, b)
        return DsetDiff(:identical, "")
    end
    # shape / type mismatches
    if a isa AbstractArray && b isa AbstractArray && size(a) != size(b)
        return DsetDiff(:differs, "shape mismatch: $(size(a)) vs $(size(b))")
    end
    if !(isnumericarray(a) && isnumericarray(b))
        # non-numeric (strings etc.) that failed isequal
        if a isa AbstractArray && b isa AbstractArray
            idxs = findall(i -> !isequal(a[i], b[i]), CartesianIndices(a))
            n = length(idxs)
            ex = isempty(idxs) ? "" :
                 " e.g. $(Tuple(idxs[1])): $(repr(a[idxs[1]])) vs $(repr(b[idxs[1]]))"
            return DsetDiff(:differs,
                "non-numeric: $(n)/$(length(a)) elements differ;$(ex)")
        else
            return DsetDiff(:differs, "scalar: $(repr(a)) vs $(repr(b))")
        end
    end
    # numeric comparison
    A = a isa Number ? fill(float(a)) : float.(a)
    B = b isa Number ? fill(float(b)) : float.(b)
    dims = size(A)
    ci = CartesianIndices(A)
    diff_idx = findall(i -> !isequal(A[i], B[i]), ci)
    n = length(diff_idx)
    # NaN / Inf pattern changes are never within-rtol
    nan_mismatch = count(i -> isnan(A[i]) != isnan(B[i]), diff_idx)
    inf_mismatch = count(i -> isinf(A[i]) != isinf(B[i]), diff_idx)
    finite_idx = [i for i in diff_idx if isfinite(A[i]) && isfinite(B[i])]
    maxabs, maxrel, argmax_i = 0.0, 0.0, nothing
    within = (nan_mismatch == 0 && inf_mismatch == 0)
    for i in finite_idx
        ad = abs(A[i] - B[i])
        denom = max(abs(A[i]), abs(B[i]))
        rd = denom > 0 ? ad / denom : 0.0
        if ad > maxabs
            maxabs = ad
            argmax_i = i
        end
        maxrel = max(maxrel, rd)
        if ad > atol + rtol * denom
            within = false
        end
    end
    stats = @sprintf("%d/%d elements differ; max abs diff %.6g, max rel diff %.6g",
        n, length(A), maxabs, maxrel)
    if argmax_i !== nothing
        stats *= " at $(Tuple(argmax_i)) ($(A[argmax_i]) vs $(B[argmax_i]))"
    end
    if nan_mismatch > 0
        stats *= "; $(nan_mismatch) NaN-pattern mismatches"
    end
    if inf_mismatch > 0
        stats *= "; $(inf_mismatch) Inf-pattern mismatches"
    end
    if within && (rtol > 0 || atol > 0)
        return DsetDiff(:within_rtol, stats)
    end
    loc = (isempty(diff_idx) || length(dims) == 0) ? "" :
          "  \n    localization: " * localize(diff_idx, dims, maxlist)
    return DsetDiff(:differs, stats * loc)
end

# ---------------------------------------------------------------------------
# per-file comparison
# ---------------------------------------------------------------------------

struct FileReport
    relpath::String
    status::Symbol   # :identical, :within_rtol, :differs, :missing_in_golden, :missing_in_new, :error
    lines::Vector{String}          # markdown detail lines (non-identical items only)
    counts::Dict{Symbol, Int}      # verdict counts over datasets+attrs
end

function compare_file(relpath, goldpath, newpath, patterns, rtol, atol, maxlist)
    counts = Dict{Symbol, Int}()
    lines = String[]
    bump(v) = counts[v] = get(counts, v, 0) + 1
    local gd, ga, nd, na
    try
        h5open(goldpath, "r") do fg
            gd, ga = enumerate_h5(fg)
        end
        h5open(newpath, "r") do fn
            nd, na = enumerate_h5(fn)
        end
    catch e
        return FileReport(relpath, :error,
            ["- ERROR opening/enumerating: $(sprint(showerror, e))"], counts)
    end
    gset, nset = Set(gd), Set(nd)
    worst = :identical
    rank = Dict(:identical => 0, :ignored => 0, :within_rtol => 1,
        :differs => 2, :missing_in_new => 2, :missing_in_golden => 2, :error => 3)
    function record(verdict, path, detail; isattr = false)
        ign = ignored(path, patterns)
        tag = isattr ? "attr" : "dataset"
        if ign
            bump(:ignored)
            if verdict != :identical
                push!(lines, "- IGNORED ($(verdict)) $(tag) `$(path)`: $(detail)")
            end
        else
            bump(verdict)
            if verdict != :identical
                vstr = verdict == :differs ? "DIFFERS" : string(verdict)
                push!(lines, "- $(vstr) $(tag) `$(path)`" *
                             (detail == "" ? "" : ": $(detail)"))
                worst = rank[verdict] > rank[worst] ? verdict : worst
            end
        end
    end
    # datasets present on only one side
    for path in sort(collect(setdiff(gset, nset)))
        record(:missing_in_new, path, "present in golden only")
    end
    for path in sort(collect(setdiff(nset, gset)))
        record(:missing_in_golden, path, "present in new only")
    end
    # shared datasets
    shared = sort(collect(intersect(gset, nset)))
    h5open(goldpath, "r") do fg
        h5open(newpath, "r") do fn
            for path in shared
                d = try
                    a = read(fg[path])
                    b = read(fn[path])
                    compare_dataset(a, b, rtol, atol, maxlist)
                catch e
                    DsetDiff(:differs, "read/compare error: $(sprint(showerror, e))")
                end
                record(d.verdict, path, d.detail)
            end
        end
    end
    # attributes
    for key in sort(collect(union(keys(ga), keys(na))))
        if !haskey(na, key)
            record(:missing_in_new, key, "attr in golden only"; isattr = true)
        elseif !haskey(ga, key)
            record(:missing_in_golden, key, "attr in new only"; isattr = true)
        else
            d = compare_dataset(ga[key], na[key], rtol, atol, maxlist)
            record(d.verdict, key, d.detail; isattr = true)
        end
    end
    return FileReport(relpath, worst, lines, counts)
end

# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

function main()
    parg = parse_commandline()
    goldroot, newroot = parg["goldenDir"], parg["newDir"]
    for d in (goldroot, newroot)
        isdir(d) || (println(stderr, "Not a directory: $d"); exit(2))
    end
    patterns = parg["no-default-ignores"] ? String[] : copy(DEFAULT_IGNORES)
    if parg["ignore"] != ""
        append!(patterns, split(parg["ignore"], ","))
    end
    patterns = String.(patterns)
    exts = String.(split(parg["glob-ext"], ","))
    rtol, atol, maxlist = parg["rtol"], parg["atol"], parg["max-list"]

    gfiles = collect_files(goldroot, exts)
    nfiles = collect_files(newroot, exts)
    allrel = sort(collect(union(keys(gfiles), keys(nfiles))))

    reports = FileReport[]
    for rp in allrel
        if !haskey(nfiles, rp)
            push!(reports, FileReport(rp, :missing_in_new, String[], Dict()))
        elseif !haskey(gfiles, rp)
            push!(reports, FileReport(rp, :missing_in_golden, String[], Dict()))
        else
            push!(reports, compare_file(rp, gfiles[rp], nfiles[rp],
                patterns, rtol, atol, maxlist))
        end
    end

    # ---- markdown report ----
    io = IOBuffer()
    println(io, "# h5diff_tree report")
    println(io)
    println(io, "- golden: `$(abspath(goldroot))`")
    println(io, "- new:    `$(abspath(newroot))`")
    println(io, "- rtol=$(rtol), atol=$(atol); ignore patterns: ",
        join(["`$p`" for p in patterns], ", "))
    println(io, "- generated: $(string(Libc.strftime("%Y-%m-%d %H:%M:%S", time()))) on $(gethostname())")
    println(io)
    nfile = length(reports)
    by(s) = count(r -> r.status == s, reports)
    tot = Dict{Symbol, Int}()
    for r in reports, (k, v) in r.counts
        tot[k] = get(tot, k, 0) + v
    end
    println(io, "## Summary")
    println(io)
    println(io, "| | count |")
    println(io, "|---|---|")
    println(io, "| files compared | $(nfile) |")
    println(io, "| files identical | $(by(:identical)) |")
    println(io, "| files within-rtol | $(by(:within_rtol)) |")
    println(io, "| files with DIFFERS | $(by(:differs)) |")
    println(io, "| files missing in new | $(by(:missing_in_new)) |")
    println(io, "| files missing in golden | $(by(:missing_in_golden)) |")
    println(io, "| files with errors | $(by(:error)) |")
    println(io, "| datasets+attrs identical | $(get(tot, :identical, 0)) |")
    println(io, "| datasets+attrs within-rtol | $(get(tot, :within_rtol, 0)) |")
    println(io, "| datasets+attrs DIFFERS | $(get(tot, :differs, 0)) |")
    println(io, "| datasets+attrs missing (either side) | $(get(tot, :missing_in_new, 0) + get(tot, :missing_in_golden, 0)) |")
    println(io, "| ignored (provenance) | $(get(tot, :ignored, 0)) |")
    println(io)
    bad = [r for r in reports if r.status != :identical]
    if isempty(bad)
        println(io, "**VERDICT: all files identical (modulo ignore list).**")
    else
        pass = all(r -> r.status in (:identical, :within_rtol), reports)
        println(io, pass ?
                    "**VERDICT: all differences within tolerance.**" :
                    "**VERDICT: DIFFERENCES FOUND.**")
        println(io)
        println(io, "## Details (non-identical files)")
        for r in bad
            println(io)
            println(io, "### `$(r.relpath)` — $(r.status == :differs ? "DIFFERS" : r.status)")
            if r.status == :missing_in_new && isempty(r.lines)
                println(io, "- file present in golden only")
            elseif r.status == :missing_in_golden && isempty(r.lines)
                println(io, "- file present in new only")
            end
            for l in r.lines
                println(io, l)
            end
        end
    end
    report = String(take!(io))
    print(report)
    if parg["out"] != ""
        mkpath(dirname(abspath(parg["out"])))
        write(parg["out"], report)
        println(stderr, "\nreport written to $(parg["out"])")
    end
    ok = all(r -> r.status in (:identical, :within_rtol), reports)
    exit(ok ? 0 : 1)
end

main()
