# Night-level comparison of the two wavelength-chain orders on existing extracted 1D products
# (no re-extraction, no resampling): runs the legacy chain (the pipeline's calls, in its order) and
# the arc-first chain (src/wavecal_arcfirst.jl) on COPIES of a night's inputs, then compares the
# adopted per-exposure solutions with each other and with an archive refit (#420) output.
#
#   julia --project=. scripts/arcfirst/compare_wavecal_order.jl --mode legacy    --inputs IN --outroot OUT --nights N.csv
#   julia --project=. scripts/arcfirst/compare_wavecal_order.jl --mode arc_first --inputs IN --outroot OUT --nights N.csv
#   julia --project=. scripts/arcfirst/compare_wavecal_order.jl --mode compare   --inputs IN --outroot OUT --nights N.csv --refit_root R
#   julia --project=. scripts/arcfirst/compare_wavecal_order.jl --mode night_all --source PASS1_APRED --scratch /tmp/x \
#         --outroot OUT --nights N.csv --refit_root R --worker i --nworkers n
#     (validation job: per night, copy the inputs to node-local scratch, run legacy + arc_first + compare
#      there, keep only the comparison rows and the small night products under OUT, delete the scratch)
#
# IN holds apred/<mjd>/ (ar1Dcal, skyLinePeaks, fpiPeaks, arclampPeaks copies) and the almanac.
# Each mode works in OUT/<mode>/ whose apred/<mjd>/ holds SYMLINKS to the input copies (outputs are
# new files there; the input copies are read-only, so an "r+" open fails loudly).
using ArgParse, HDF5, JLD2, CSV, DataFrames, Statistics, Dates
using ApogeeReduction
const AR = ApogeeReduction

s = ArgParseSettings()
@add_arg_table! s begin
    "--mode"; arg_type = String; required = true
    "--inputs"; arg_type = String; required = true
    "--outroot"; arg_type = String; required = true
    "--nights"; arg_type = String; required = true
    "--refit_root"; arg_type = String; default = ""
    "--source"; arg_type = String; default = ""
    "--almanac"; arg_type = String; default = ""
    "--scratch"; arg_type = String; default = ""
    "--worker"; arg_type = Int; default = 1
    "--nworkers"; arg_type = Int; default = 1
    "--claim"; action = :store_true
    "--almanac_name"; arg_type = String; default = "allobs_57618_61230.h5"
    "--chip_transform_cal"; arg_type = String
    default = "/mnt/ceph/users/sdssv/work/asaydjari/cal_ref/chip_transform/chiptransform_pass1dr21_v2_allblock.h5"
    "--shape_template"; arg_type = String
    default = "/mnt/ceph/users/sdssv/work/asaydjari/2026_11_17/wavecal_refit/cal/template_shape_v1.h5"
end
const A_ = parse_args(s)
const NIGHTS = CSV.read(A_["nights"], DataFrame)
const IN = Ref(A_["inputs"])
const OUTROOT = Ref(A_["outroot"])

"Symlink unless it exists (parallel workers share the almanac link)."
function safe_symlink(target, l)
    (islink(l) || isfile(l)) && return
    try
        symlink(target, l)
    catch
        (islink(l) || isfile(l)) || rethrow()
    end
end

"Mode working directory with symlinks to the input copies of a night."
function setup_mode(mode, tele, mjd)
    od = joinpath(OUTROOT[], mode)
    md = joinpath(od, "apred", string(mjd))
    mkpath(md)
    mkpath(joinpath(od, "almanac"))
    mkpath(joinpath(od, "wavecal"))
    alm = joinpath(od, "almanac", A_["almanac_name"])
    safe_symlink(joinpath(IN[], A_["almanac_name"]), alm)
    src = joinpath(IN[], "apred", string(mjd))
    for f in readdir(src)
        occursin("_$(tele)_$(mjd)_", f) || continue
        l = joinpath(md, f)
        safe_symlink(joinpath(src, f), l)
    end
    return od, md, alm
end

files(md, pre, tele, mjd, suf) = sort(filter(f -> startswith(basename(f), "$(pre)_$(tele)_$(mjd)_") &&
                                                   endswith(f, suf), readdir(md; join = true)))

function run_legacy(tele, mjd)
    od, md, alm = setup_mode("legacy", tele, mjd)
    cal = AR.load_chip_transform_cal(A_["chip_transform_cal"])
    sky = files(md, "skyLinePeaks", tele, mjd, "_$(FIRST_CHIP)_object.h5")
    wc = [AR.get_and_save_sky_wavecal(f; checkpoint_mode = "commit_same", chip_transform_cal = cal) for f in sky]
    good = .!isnothing.(wc)
    wcs = convert(Vector{String}, wc[good])
    sky = sky[good]
    mjdl = fill(mjd, length(wcs))
    na = AR.skyline_medwavecal_skyline_dither(tele, mjd, mjdl, wcs, sky, "commit_same"; outdir = od * "/",
        almanac_file = alm)
    fpi = files(md, "fpiPeaks", tele, mjd, "_$(FIRST_CHIP)_arclamp.h5")
    if !isempty(fpi)
        AR.fpi_medwavecal_skyline_dither(mjd, fill(mjd, length(fpi)), mjdl, fpi, sky, na; verbose = false,
            checkpoint_mode = "commit_same", almanac_file = alm)
    end
    return "ok"
end

function run_arcfirst(tele, mjd)
    od, md, alm = setup_mode("arc_first", tele, mjd)
    cfg = AR.RefitConfig(chip_transform_file = A_["chip_transform_cal"], shape_template_file = A_["shape_template"],
        anchor_samepos = "none")
    exps = AR.arcfirst_exps_from_names(files(md, "ar1Dcal", tele, mjd, ".h5"), tele, mjd)
    fpi = files(md, "fpiPeaks", tele, mjd, "_$(FIRST_CHIP)_arclamp.h5")
    r = AR.arcfirst_night(tele, mjd; outdir = od, almanac_file = alm, exps, fpi_fnames = fpi, cfg,
        checkpoint_mode = "commit_exists")
    return "route=$(r.route) arc=$(r.arc) comb=$(r.comb) fpifib=$(r.fpifib) assemble=$(r.assemble) n=$(r.n_exp_written)/$(length(exps))"
end

"Legacy adopted per-exposure grid (the reinterp_spectra choice), its type, or nothing."
function legacy_grid(md, od, tele, mjd, e, kind)
    na = joinpath(od, "wavecal", "wavecalNightAve_$(tele)_$(mjd).h5")
    bwt = isfile(na) ? h5open(f -> haskey(f, "best_wave_type") ? read(f["best_wave_type"]) : "fpi", na, "r") : "fpi"
    for wt in (bwt == "sky" ? ["sky"] : ["fpi", "sky"])
        fn = joinpath(md, "waveCalNight$(wt)Dither_$(tele)_$(mjd)_$(AR._arc_expstr(e))_$(kind).h5")
        if isfile(fn)
            G = h5open(f -> read(f["chipWaveSoln"]), fn, "r")
            all(isfinite, G) && return G, "waveCalNight$(wt)Dither"
        end
    end
    return nothing, ""
end

"Per-chip comparison of grid a against reference b on the flat pixel grid (m/s): zero point, shape."
function cmp_grids(a, b)
    _, zp, sh, _ = AR.refit_vs_pass1(a[AR.REFIT_PXG, :, :], b[AR.REFIT_PXG, :, :])
    return zp, sh
end

function compare_rows()
    return DataFrame(tele = String[], mjd = Int[], era = String[], expnum = Int[], kind = String[],
        pair = String[], chip = String[], zp_ms = Float64[], shape_ms = Float64[], type_a = String[], type_b = String[],
        bits_a = Int64[])
end

function compare_night!(rows, tele, mjd)
        era = AR.refit_era(tele, mjd)
        odA = joinpath(OUTROOT[], "arc_first"); mdA = joinpath(odA, "apred", string(mjd))
        odL = joinpath(OUTROOT[], "legacy"); mdL = joinpath(odL, "apred", string(mjd))
        rf = joinpath(A_["refit_root"], "exp", "waveCalRefit_$(tele)_$(mjd).h5")
        rexp = isfile(rf) ? h5open(f -> Int.(read(f["exp_expnum"])), rf, "r") : Int[]
        isdir(mdA) || return rows
        exps = AR.arcfirst_exps_from_names(files(mdA, "ar1Dcal", tele, mjd, ".h5"), tele, mjd)
        for (e, kind) in exps
            kind == "object" || continue
            fa = AR.arcfirst_exp_name(joinpath(odA, "apred"), tele, mjd, e, kind)
            isfile(fa) || continue
            Ga, ta, ba = h5open(f -> (read(f["chipWaveSoln"]), String(read(f["wavecal_type"])),
                Int64(read(f["bitmsk_wavecal"]))), fa, "r")
            Gl, tl = legacy_grid(mdL, odL, tele, mjd, e, kind)
            refs = Any[("arcfirst_vs_legacy", Gl, tl)]
            if e in rexp
                Gr = AR.refit_exposure_wave(rf, e)
                tr = h5open(f -> String(read(f["wavecal_type"])[findfirst(==(e), rexp)]), rf, "r")
                push!(refs, ("arcfirst_vs_refit", Gr, tr))
            end
            for (pair, Gb, tb) in refs
                Gb === nothing && continue
                zp, sh = cmp_grids(Ga, Gb)
                for c in 1:3
                    push!(rows, (tele, mjd, era, e, kind, pair, CHIP_LIST[c], zp[c], sh[c], ta, tb, ba))
                end
            end
        end
    return rows
end

function summarize(rows)
    CSV.write(joinpath(OUTROOT[], "compare_exposures.csv"), rows)
    # night summary: median zero point, per-exposure scatter (std of the exposure zero points), median shape
    g = combine(groupby(rows, [:tele, :mjd, :era, :pair, :chip]),
        :zp_ms => (x -> median(x)) => :zp_med_ms, :zp_ms => (x -> length(x) > 1 ? std(x) : NaN) => :zp_exp_std_ms,
        :shape_ms => (x -> median(x)) => :shape_med_ms, nrow => :n_exp,
        :type_a => (x -> join(unique(x), "|")) => :types_arcfirst, :type_b => (x -> join(unique(x), "|")) => :types_ref)
    sort!(g, [:pair, :tele, :era, :mjd, :chip])
    CSV.write(joinpath(OUTROOT[], "compare_nights.csv"), g)
    show(stdout, g; allrows = true, allcols = true)
    println()
end

function run_compare()
    rows = compare_rows()
    for r in eachrow(NIGHTS)
        compare_night!(rows, String(r.tele), Int(r.mjd))
    end
    summarize(rows)
end

"Validation job: one night end to end in node-local scratch; keeps rows + small products."
function run_night_all(tele, mjd)
    keep = joinpath(A_["outroot"], "nights", "$(tele)_$(mjd)")
    isfile(joinpath(keep, "rows.csv")) && return "done"
    sc = joinpath(A_["scratch"], "$(tele)_$(mjd)")
    rm(sc; force = true, recursive = true)
    inp = joinpath(sc, "inputs"); md = joinpath(inp, "apred", string(mjd)); mkpath(md)
    src = joinpath(A_["source"], string(mjd))
    rx = Regex("^(ar1Dcal_$(tele)_$(mjd)_\\d{4}_[RGB]_(object|arclamp)|skyLinePeaks_$(tele)_$(mjd)_|fpiPeaks_$(tele)_$(mjd)_|arclampPeaks_$(tele)_$(mjd)_)")
    n = 0
    for f in readdir(src)
        occursin(rx, f) || continue
        cp(joinpath(src, f), joinpath(md, f)); n += 1    # plain read of the source; copies are made read-only
        chmod(joinpath(md, f), 0o444)
    end
    symlink(A_["almanac"], joinpath(inp, A_["almanac_name"]))
    IN[] = inp; OUTROOT[] = joinpath(sc, "runs")
    t0 = time(); sl = try run_legacy(tele, mjd) catch e; "error: " * sprint(showerror, e)[1:min(end, 200)] end
    t1 = time(); sa = try run_arcfirst(tele, mjd) catch e; "error: " * sprint(showerror, e)[1:min(end, 200)] end
    t2 = time()
    rows = compare_rows()
    try
        compare_night!(rows, tele, mjd)
    catch e
        @warn "compare failed" tele mjd exception = e
    end
    mkpath(keep)
    od = joinpath(OUTROOT[], "arc_first", "wavecal_arcfirst")
    isdir(od) && cp(od, joinpath(keep, "wavecal_arcfirst"); force = true)
    na = joinpath(OUTROOT[], "legacy", "wavecal", "wavecalNightAve_$(tele)_$(mjd).h5")
    isfile(na) && cp(na, joinpath(keep, basename(na)); force = true)
    open(joinpath(keep, "status.txt"), "w") do io
        println(io, "n_input_files=$(n)\tlegacy_s=$(round(t1 - t0, digits = 1))\tarcfirst_s=$(round(t2 - t1, digits = 1))")
        println(io, "legacy: ", sl); println(io, "arc_first: ", sa)
    end
    CSV.write(joinpath(keep, "rows.csv"), rows)
    rm(sc; force = true, recursive = true)
    return "legacy: $(sl) | arc_first: $(sa) | rows=$(nrow(rows))"
end

mode = A_["mode"]
if mode == "compare"
    run_compare()
elseif mode == "summarize"
    # concatenate the validation job's per-night rows and summarize them
    rows = compare_rows()
    for d in readdir(joinpath(OUTROOT[], "nights"); join = true)
        f = joinpath(d, "rows.csv")
        isfile(f) && filesize(f) > 0 && append!(rows, CSV.read(f, DataFrame; types = Dict(:type_a => String, :type_b => String)); promote = true)
    end
    summarize(rows)
else
    for (i, r) in enumerate(eachrow(NIGHTS))
        (i - 1) % A_["nworkers"] == A_["worker"] - 1 || continue
        if A_["claim"]      # validation job: nights claimed through lock directories (any node count/size)
            lk = joinpath(OUTROOT[], "locks_" * get(ENV, "SLURM_JOB_ID", "local"), "$(r.tele)_$(r.mjd)")
            mkpath(dirname(lk))
            try
                mkdir(lk)
            catch
                continue
            end
        end
        tele, mjd = String(r.tele), Int(r.mjd)
        t0 = time()
        st = try
            mode == "legacy" ? run_legacy(tele, mjd) : mode == "arc_first" ? run_arcfirst(tele, mjd) :
            mode == "night_all" ? run_night_all(tele, mjd) :
            error("unknown mode $(mode)")
        catch e
            @warn "night failed" tele mjd exception = (e, catch_backtrace())
            "error: " * sprint(showerror, e)[1:min(end, 200)]
        end
        println("$(now()) $(mode) $(tele) $(mjd) $(round(time() - t0, digits = 1)) s: $(st)")
        flush(stdout)
    end
end
