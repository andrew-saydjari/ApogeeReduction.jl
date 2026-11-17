# Rebuild the comb-measured inter-chip transform calibration from a whole reduction.
# See README.md in this folder for what it does, why, and how to run it before a release.
#
#   julia --project=<repo> scripts/cal/chip_transform/build_chip_transform.jl \
#       --root <reduction root with apred/ and wavecal/> --tele apo,lco --workdir <dir> \
#       --stage all --workers 6 --version <name> [--cover_from 56893] [--cover_to 99999]
#
# Stages (each resumable; a later stage reads the earlier stages' outputs in --workdir):
#   collect  stage 1: per night, the comb fit's per-fiber chip parameters (running the
#            production comb fit, into --workdir only, for nights whose comb file is missing
#            if --run_missing true)
#   derive   stage 2: blocks of time + per-fiber chip transform with robust statistics, QA plots,
#            a table of everything rejected
#   write    stage 3: the versioned calibration file with full provenance
#   check    stage 4: before/after sky (and arc) solutions on a sample of nights, held-out lines
#   all      1-4

using ArgParse, Distributed, Dates

function parse_commandline()
    s = ArgParseSettings(description = "Rebuild the FPI-comb chip-transform calibration")
    @add_arg_table s begin
        "--stage"
        help = "collect | derive | write | check | all"
        default = "all"
        "--root"
        help = "reduction root holding apred/<mjd>/ and wavecal/ (read-only unless --run_missing writes into --workdir)"
        required = true
        "--tele"
        help = "comma-separated telescopes"
        default = "apo,lco"
        "--workdir"
        help = "directory for every intermediate and output of the build"
        required = true
        "--workers"
        help = "local single-threaded workers (ignored under Slurm, where SlurmManager takes all tasks)"
        arg_type = Int
        default = 1
        "--mjd_min"
        arg_type = Int
        default = 0
        "--mjd_max"
        arg_type = Int
        default = 999999
        "--run_missing"
        help = "stage 1: run the production comb fit (into --workdir) for nights with FPI exposures but no comb file"
        arg_type = Bool
        default = false
        # ---- stage 2
        "--block_mode"
        help = "detect (from the data) | supplied (--boundaries) | supplied+detect | none (one block)"
        default = "detect"
        "--boundaries"
        help = "file with lines '<tele> <mjd>' = first MJD of a new block"
        default = ""
        "--min_block_nights"
        help = "minimum comb nights on each side of a detected boundary"
        arg_type = Int
        default = 20
        "--block_nsigma"
        help = "significance of a step in the night-median chip scale to open a new block"
        arg_type = Float64
        default = 5.0
        "--block_min_step"
        help = "minimum |step| in chip scale to open a new block (1e-5 ~ 3 m/s of per-chip shape)"
        arg_type = Float64
        default = 1.0e-5
        "--offset_event_min_step"
        help = "QA only: minimum |step| in chip offset (ximport units; 1 px = 4.9e-4) reported as a hardware event"
        arg_type = Float64
        default = 2.5e-4
        "--block_min_used_nights"
        help = "warn if a block has fewer usable nights"
        arg_type = Int
        default = 10
        "--fiber_min_peaks"
        help = "a fiber is used on a night only if the comb fit used at least this many peaks on every chip"
        arg_type = Int
        default = 100
        "--night_min_fibers"
        help = "a night is used only if at least this many fibers survive"
        arg_type = Int
        default = 200
        "--clip_nsigma"
        help = "robust clipping threshold for fiber and night outliers"
        arg_type = Float64
        default = 5.0
        "--fiber_model"
        help = "auto (split-half cross-validation) | per_fiber | poly"
        default = "auto"
        "--fiber_poly_degree"
        arg_type = Int
        default = 2
        # ---- stage 3
        "--version"
        help = "calibration version string (also the file name)"
        default = ""
        "--cover_from"
        help = "extend the first block of every telescope back to this MJD (e.g. the release start, to cover the plate era); 0 = no extension"
        arg_type = Int
        default = 0
        "--cover_to"
        help = "extend the last block forward to this MJD; 0 = no extension"
        arg_type = Int
        default = 0
        "--cover_note"
        help = "free text recorded in the provenance explaining the coverage choice (e.g. the validation behind --cover_from)"
        default = ""
        # ---- stage 4
        "--cal"
        help = "calibration file to check (default: the stage-3 output)"
        default = ""
        "--check_nights"
        help = "file with lines '<tele> <mjd> <category>'; default: --check_per_category nights per category"
        default = ""
        "--check_per_category"
        help = "nights per telescope and category (comb / nocomb), evenly spaced in MJD"
        arg_type = Int
        default = 40
        "--check_tag"
        default = "check"
        "--arc_lines_dir"
        help = "optional directory of arclines_<tele>_<mjd>.h5 identified arc line tables"
        default = ""
        "--max_sky_exp"
        arg_type = Int
        default = 4
    end
    return parse_args(s)
end

const ARGS_P = parse_commandline()
const P = ARGS_P
const TELES = split(P["tele"], ",")
const WORKDIR = abspath(P["workdir"])
const ROOT = abspath(P["root"])
const REPO = dirname(dirname(dirname(@__DIR__)))
mkpath(WORKDIR)

function add_workers()
    if haskey(ENV, "SLURM_NTASKS")
        @eval using SlurmClusterManager
        addprocs(SlurmManager(), exeflags = ["--project=$(REPO)"])
    elseif P["workers"] > 1
        addprocs(P["workers"], exeflags = ["--project=$(REPO)"],
            env = ["OPENBLAS_NUM_THREADS" => "1", "OMP_NUM_THREADS" => "1", "JULIA_NUM_THREADS" => "1"])
    end
end
add_workers()
@everywhere begin
    using LinearAlgebra
    BLAS.set_num_threads(1)
    include(joinpath($(@__DIR__), "chip_transform_stages.jl"))
    include(joinpath($(@__DIR__), "chip_transform_check.jl"))
end
using CairoMakie

stages = P["stage"] == "all" ? ["collect", "derive", "write", "check"] : [P["stage"]]
t_start = time()
logt(msg) = (println(Dates.format(now(), "HH:MM:SS"), "  [", round(Int, time() - t_start), " s] ", msg); flush(stdout))

function git_info()
    c = try; readchomp(`git -C $(REPO) rev-parse HEAD`); catch; "unknown"; end
    d = try; isempty(readchomp(`git -C $(REPO) status --porcelain --untracked-files=no`)) ? "clean" : "dirty"; catch; "unknown"; end
    b = try; readchomp(`git -C $(REPO) rev-parse --abbrev-ref HEAD`); catch; "unknown"; end
    return c, d, b
end

# ------------------------------------------------------------------------------ stage 1
if "collect" in stages
    for tele in TELES
        mjds = list_nights(ROOT, tele; mjd_min = P["mjd_min"], mjd_max = P["mjd_max"])
        logt("stage 1 $(tele): $(length(mjds)) nights with a nightly wavelength file")
        st = pmap(m -> (m, try; stage1_night(ROOT, tele, m, WORKDIR; run_missing = P["run_missing"]); "ok"; catch e; "ERROR " * sprint(showerror, e)[1:min(end, 200)]; end), mjds)
        errs = [s for s in st if s[2] != "ok"]
        isempty(errs) || (println("stage 1 errors:"); foreach(println, errs))
        S = stage1_gather(WORKDIR, tele, [m for (m, s) in st if s == "ok"])
        CSV.write(joinpath(WORKDIR, "stage1_$(tele)_nights.csv"), DataFrame(tele = tele, mjd = S.mjds, status = S.status, source = S.src))
        c = Dict{String, Int}(); for s in S.status; c[s] = get(c, s, 0) + 1; end
        logt("stage 1 $(tele): " * join(["$(k)=$(v)" for (k, v) in sort(collect(c))], ", "))
    end
end

# ------------------------------------------------------------------------------ stage 2
function stage2_params()
    Dict(k => P[k] for k in ("block_mode", "boundaries", "min_block_nights", "block_nsigma", "block_min_step",
        "offset_event_min_step", "block_min_used_nights", "fiber_min_peaks", "night_min_fibers", "clip_nsigma",
        "fiber_model", "fiber_poly_degree"))
end

function qa_plots(tele, S, D)
    set_theme!(merge(theme_black(), theme_latexfonts()))
    pdir = joinpath(WORKDIR, "qa"); mkpath(pdir)
    fig = Figure(size = (1500, 900))
    labs = ("chip-R scale  b_R - 1", "chip-B scale  b_B - 1", "chip-R offset a_R (px)", "chip-B offset a_B (px)")
    for q in 1:4
        ax = Axis(fig[q, 1], ylabel = labs[q], xlabel = q == 4 ? "MJD" : "")
        y = q <= 2 ? (D.ser[:, q] .- 1) .* 1e4 : D.ser[:, q] .* N_XPIX
        q <= 2 && (ax.ylabel = labs[q] * " (1e-4)")
        g = D.night_ok
        scatter!(ax, S.mjds[g], y[g], markersize = 4, color = :deepskyblue)
        for m in D.bmjd; vlines!(ax, [m], color = :orange, linewidth = 2); end
        for m in D.starts_offset_mjd; vlines!(ax, [m], color = :gray, linestyle = :dash); end
        if q <= 2
            for b in D.blocks
                lines!(ax, [b.data_mjd_min, b.data_mjd_max], fill((median(b.scale[:, q]) - 1) * 1e4, 2), color = :yellow, linewidth = 3)
            end
        end
    end
    Label(fig[0, 1], "$(tele): comb chip transform per night (fiber median). orange = block boundaries (scale), gray dashed = offset events (QA only)")
    save(joinpath(pdir, "qa_$(tele)_night_series.png"), fig)
    fig = Figure(size = (1500, 700))
    for (j, nm) in enumerate(("R", "B"))
        ax = Axis(fig[1, j], xlabel = "fiber index", ylabel = "b_$(nm) - 1 (1e-4)", title = "$(tele) chip $(nm): per block")
        for (k, b) in enumerate(D.blocks)
            scatter!(ax, 1:N_FIBERS, (b.scale[:, j] .- 1) .* 1e4, markersize = 4, label = "block $(k): $(b.data_mjd_min)-$(b.data_mjd_max) ($(b.n_nights) n, $(D.model[k]))")
        end
        j == 1 && axislegend(ax, position = :lb, labelsize = 10)
    end
    save(joinpath(pdir, "qa_$(tele)_per_fiber.png"), fig)
end

if "derive" in stages
    for tele in TELES
        st = CSV.read(joinpath(WORKDIR, "stage1_$(tele)_nights.csv"), DataFrame)
        S = stage1_gather(WORKDIR, tele, st.mjd)
        rejected = NamedTuple[]
        D = stage2_derive(S, tele; P = stage2_params(), rejected = rejected)
        h5open(joinpath(WORKDIR, "stage2_$(tele).h5"), "w") do f
            f["mjds"] = S.mjds; f["night_series"] = D.ser; f["night_ok"] = UInt8.(D.night_ok)
            f["block_starts_mjd"] = isempty(D.bmjd) ? Int[] : D.bmjd
            f["offset_event_mjd"] = isempty(D.starts_offset_mjd) ? Int[] : D.starts_offset_mjd
            for (k, b) in enumerate(D.blocks)
                g = create_group(f, "block_$(k)")
                g["data_mjd"] = [b.data_mjd_min, b.data_mjd_max]; g["n_nights"] = b.n_nights
                g["scale"] = b.scale; g["scale_err"] = b.scale_err; g["offset"] = b.offset; g["fiber_filled"] = UInt8.(b.fiber_filled)
                g["fiber_model"] = D.model[k]
            end
        end
        CSV.write(joinpath(WORKDIR, "stage2_$(tele)_rejected.csv"), DataFrame(rejected))
        CSV.write(joinpath(WORKDIR, "stage2_$(tele)_fiber_model_cv.csv"), DataFrame(D.cv))
        CSV.write(joinpath(WORKDIR, "stage2_$(tele)_blocks.csv"), DataFrame(
            [(tele = tele, block = k, data_mjd_min = b.data_mjd_min, data_mjd_max = b.data_mjd_max, n_nights = b.n_nights,
              fiber_model = D.model[k], bR_median = median(b.scale[:, 1]), bB_median = median(b.scale[:, 2]),
              bR_err_median = median(filter(isfinite, b.scale_err[:, 1])), bB_err_median = median(filter(isfinite, b.scale_err[:, 2])),
              aR_median_px = median(b.offset[:, 1]) * N_XPIX, aB_median_px = median(b.offset[:, 2]) * N_XPIX,
              n_fibers_filled = count(b.fiber_filled)) for (k, b) in enumerate(D.blocks)]))
        CSV.write(joinpath(WORKDIR, "stage2_$(tele)_night_series.csv"), DataFrame(tele = tele, mjd = S.mjds, used = D.night_ok,
            status = S.status, bR = D.ser[:, 1], bB = D.ser[:, 2], aR = D.ser[:, 3], aB = D.ser[:, 4]))
        qa_plots(tele, S, D)
        logt("stage 2 $(tele): $(count(D.night_ok)) nights used, $(length(D.blocks)) block(s) starting $(join([b.data_mjd_min for b in D.blocks], ", ")); " *
             "fiber model $(join(D.model, ",")); offset events (QA) at $(join(D.starts_offset_mjd, ", ")); $(length(rejected)) rejections")
    end
end

# ------------------------------------------------------------------------------ stage 3
version = isempty(P["version"]) ? "chiptransform_" * Dates.format(now(), "yyyymmdd") : P["version"]
calfile = isempty(P["cal"]) ? joinpath(WORKDIR, "$(version).h5") : abspath(P["cal"])
if "write" in stages
    blocks = Dict{String, Vector{ChipTransformBlock}}()
    for tele in TELES
        bl = ChipTransformBlock[]
        h5open(joinpath(WORKDIR, "stage2_$(tele).h5"), "r") do f
            k = 1
            while haskey(f, "block_$(k)")
                g = f["block_$(k)"]; dm = read(g["data_mjd"])
                push!(bl, ChipTransformBlock(dm[1], dm[2], dm[1], dm[2], read(g["n_nights"]), read(g["scale"]),
                    read(g["scale_err"]), read(g["offset"]), read(g["fiber_filled"]) .!= 0))
                k += 1
            end
        end
        blocks[tele] = set_coverage(bl; cover_from = P["cover_from"] > 0 ? P["cover_from"] : nothing,
            cover_to = P["cover_to"] > 0 ? P["cover_to"] : nothing)
    end
    commit, clean, branch = git_info()
    prov = Dict{String, String}(
        "created_utc" => Dates.format(now(UTC), "yyyy-mm-ddTHH:MM:SS"),
        "code_commit" => commit, "code_tree" => clean, "code_branch" => branch,
        "builder" => "scripts/cal/chip_transform/build_chip_transform.jl",
        "input_root" => ROOT, "telescopes" => join(TELES, ","),
        "mjd_min_arg" => string(P["mjd_min"]), "mjd_max_arg" => string(P["mjd_max"]),
        "cover_from" => string(P["cover_from"]), "cover_to" => string(P["cover_to"]), "cover_note" => P["cover_note"],
        "fixed_terms_in_pipeline" => "scale (b_R, b_B); offsets stored for QA only",
        "nlParams_layout" => "a_R, b_R, a_B, b_B (x' = a_c + b_c * ximport, ximport = (pixel - 1024) / 2048)")
    for (k, v) in stage2_params(); prov["stage2_" * k] = string(v); end
    for tele in TELES
        bl = blocks[tele]
        prov["$(tele)_blocks"] = join(["[$(b.mjd_start), $(b.mjd_end)] data $(b.data_mjd_min)-$(b.data_mjd_max) n=$(b.n_nights)" for b in bl], "; ")
        st = CSV.read(joinpath(WORKDIR, "stage1_$(tele)_nights.csv"), DataFrame)
        prov["$(tele)_stage1_counts"] = join(["$(r.status)=$(r.nrow)" for r in eachrow(combine(groupby(st, :status), nrow))], ", ")
        cv = CSV.read(joinpath(WORKDIR, "stage2_$(tele)_fiber_model_cv.csv"), DataFrame)
        prov["$(tele)_fiber_model"] = join(string.(cv.fiber_model), ",")
        prov["$(tele)_n_rejections"] = string(nrow(CSV.read(joinpath(WORKDIR, "stage2_$(tele)_rejected.csv"), DataFrame)))
    end
    write_chip_transform_cal(calfile, blocks, prov; calibration_version = version)
    logt("stage 3: wrote $(calfile)")
end

# ------------------------------------------------------------------------------ stage 4
if "check" in stages
    cal = load_chip_transform_cal(calfile)
    nights = if !isempty(P["check_nights"])
        rows = NamedTuple[]
        for l in eachline(P["check_nights"])
            s = split(strip(l)); (isempty(s) || startswith(s[1], "#")) && continue
            push!(rows, (tele = String(s[1]), mjd = parse(Int, s[2]), category = length(s) >= 3 ? String(s[3]) : "sample"))
        end
        DataFrame(rows)
    else
        rows = NamedTuple[]
        for tele in TELES
            s1 = CSV.read(joinpath(WORKDIR, "stage1_$(tele)_nights.csv"), DataFrame)
            for (cat, sel) in (("comb", startswith.(s1.status, "used")), ("nocomb", s1.status .== "no_fpi_exposures"))
                m = sort(s1.mjd[sel]); isempty(m) && continue
                m = filter(x -> chip_transform_block(cal, tele, x) !== nothing, m)
                isempty(m) && continue
                pick = unique(m[round.(Int, range(1, length(m), length = min(length(m), P["check_per_category"])))])
                append!(rows, [(tele = String(tele), mjd = x, category = cat) for x in pick])
            end
        end
        DataFrame(rows)
    end
    tag = P["check_tag"]
    logt("stage 4: $(nrow(nights)) nights, tag $(tag)")
    jobs = [(r.tele, r.mjd) for r in eachrow(nights)]
    calpath = calfile; root = ROOT; wd = WORKDIR; arcdir = P["arc_lines_dir"]; mx = P["max_sky_exp"]
    st4 = pmap(jobs) do (t, m)
        c = load_chip_transform_cal(calpath)
        t0 = time()
        s = try; string(check_night(root, t, m, c, wd, tag; arcdir = arcdir, max_exp = mx)); catch e; "ERROR " * sprint(showerror, e)[1:min(end, 200)]; end
        (t, m, s, time() - t0)
    end
    CSV.write(joinpath(WORKDIR, "stage4_$(tag)_status.csv"), DataFrame(tele = first.(st4), mjd = [s[2] for s in st4], status = [s[3] for s in st4], seconds = [s[4] for s in st4]))
    T = check_collect(WORKDIR, tag, nights)
    CSV.write(joinpath(WORKDIR, "stage4_$(tag)_per_night.csv"), T)
    if nrow(T) > 0
        S1 = check_summary(T; ref = "scale1"); CSV.write(joinpath(WORKDIR, "stage4_$(tag)_summary_vs_scale1.csv"), S1)
        S2 = check_summary(T; ref = "free_scale"); CSV.write(joinpath(WORKDIR, "stage4_$(tag)_summary_vs_free_scale.csv"), S2)
        F = check_flat_rv(T); CSV.write(joinpath(WORKDIR, "stage4_$(tag)_flat_rv.csv"), F)
        show(stdout, S1[S1.model .== "cal_scale", [:tele, :category, :tracer, :chip, :n_nights, :loo_rms_mA, :dloo_mA, :dloo_lo, :dloo_hi, :frac_improved, :dzp_ms, :shape_ms, :shape_ref_ms]]; allrows = true, allcols = true)
        println()
    end
    logt("stage 4 done")
end
