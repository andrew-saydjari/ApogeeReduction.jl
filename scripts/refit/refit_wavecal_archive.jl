# Archive-scale arc + FPI wavelength REFIT driver (design stages 3-6) on existing pass-1 products.
#
# Phases (run in order; each is resumable and skips finished work):
#   prep   write <outroot>/config.txt (switches, code commit, calibration fingerprints) and the
#          night list <outroot>/nights.csv; refuses to run into an outroot made with other switches.
#   1      per night: arc night fit (stage A) -> arc-anchored comb (stage B; nights without an arc
#          anchor are marked `deferred`).
#   2      per night: comb for deferred nights from the nearest arc-anchored neighbour cavity
#          (arcwin), FPI fibers (stage C), assembly per exposure (stage D).
#   3      single process: run table, cross-night comb checks (FSR continuity, cavity jumps),
#          counts; writes run_table.csv, crossnight_flags.csv, run_summary.txt.
# Workers claim nights through lock directories (atomic mkdir), so any number of workers on any
# number of nodes can run a phase concurrently; a failed night is logged and flagged, never fatal.
# pass1_dr21 is only ever opened read-only; outputs go only under --outroot.
#
# usage (from the repository root):
#   julia --project=. scripts/refit/refit_wavecal_archive.jl --phase prep --outroot DIR [switches]
#   julia --project=. scripts/refit/refit_wavecal_archive.jl --phase 1 --worker K --outroot DIR [switches]
#   julia --project=. scripts/refit/refit_wavecal_archive.jl --phase 2 --worker K --outroot DIR [switches]
#   julia --project=. scripts/refit/refit_wavecal_archive.jl --phase 3 --outroot DIR [switches]
# Variant reruns: give a new --outroot and --reuse_from BASE; products of stages that the changed
# switches do not invalidate are symlinked from BASE instead of recomputed.
using ApogeeReduction, ArgParse, CSV, DataFrames, HDF5, Dates, Printf, Statistics
const AR = ApogeeReduction
BLAS_THREADS = try
    using LinearAlgebra; LinearAlgebra.BLAS.set_num_threads(1); 1
catch
    0
end

function parse_cli()
    s = ArgParseSettings(description = "arc + FPI wavelength refit driver")
    @add_arg_table! s begin
        "--phase"; arg_type = String; required = true
        "--outroot"; arg_type = String; default = "/mnt/ceph/users/sdssv/work/asaydjari/2026_11_17/wavecal_refit/run1"
        "--pass1"; arg_type = String; default = "/mnt/ceph/users/sdssv/work/asaydjari/pass1_dr21"
        "--almanac"; arg_type = String; default = ""
        "--nights"; arg_type = String; default = ""; help = "CSV with columns tele,mjd (default: <outroot>/nights.csv; prep builds it from the pass-1 wavecal folder)"
        "--worker"; arg_type = Int; default = 1
        "--retry_failed"; action = :store_true
        "--reuse_from"; arg_type = String; default = ""
        "--max_nights"; arg_type = Int; default = 0; help = "prep: keep only the first N nights of the list (tests)"
        # switches (RefitConfig)
        "--neon"; arg_type = String; default = "none"
        "--neon_offset_file"; arg_type = String; default = ""
        "--neon_offset_set"; arg_type = String; default = "all"
        "--neon_line_offset_file"; arg_type = String; default = "/mnt/ceph/users/sdssv/work/asaydjari/2026_11_17/wavecal_refit/cal/neon_line_offsets_v1.csv"
        "--chip_transform_file"; arg_type = String; default = "/mnt/ceph/users/sdssv/work/asaydjari/cal_ref/chip_transform/chiptransform_pass1dr21_v2_allblock.h5"
        "--exp_terms"; arg_type = String; default = "chipzp_stretch_tilt"
        "--sky_fiber_term"; arg_type = String; default = "full"
        "--prefps_time"; arg_type = String; default = "regularised"
        "--sky_exp_terms"; arg_type = String; default = "chipzp_stretch_tilt"
        "--reg_weight_apo"; arg_type = Float64; default = 0.825
        "--reg_weight_lco"; arg_type = Float64; default = 0.65
        "--shape_template_file"; arg_type = String; default = "/mnt/ceph/users/sdssv/work/asaydjari/2026_11_17/wavecal_refit/cal/template_shape_v1.h5"
        "--shape_template_block"; arg_type = String; default = "all"
        "--plate_shape"; arg_type = String; default = "arc"; help = "plate-era shape, both telescopes: arc | template"
        "--fps_shape_apo"; arg_type = String; default = "template"
        "--fps_shape_lco"; arg_type = String; default = "template"
        "--anchor_block"; arg_type = String; default = "nearest_science"; help = "nearest_science | morning | any"
        "--anchor_prefer"; arg_type = String; default = "strict"; help = "nearest_science preference: strict | none"
        "--anchor_samepos"; arg_type = String; default = "pass1"; help = "same-position input: pass1 (later-stage solutions) | none"
        "--anchor_adjacency"; arg_type = String; default = "strict_then_nomove"
        "--anchor_max_candidates"; arg_type = Int; default = 3
        "--anchor_pattern_tol_px"; arg_type = Float64; default = 0.02
        "--anchor_lamp_tol_px"; arg_type = Float64; default = 0.1
        "--anchor_pairing"; arg_type = String; default = "block"
        "--anchor_block_h"; arg_type = Float64; default = 1.5
        "--arcwin_nights"; arg_type = Int; default = 15
        "--nsig"; arg_type = Float64; default = 8.0
        "--sky_anchor_window_h"; arg_type = Float64; default = 1.5
        "--pos_match_px"; arg_type = Float64; default = 0.25
        "--xcheck_threshold_ms"; arg_type = Float64; default = 60.0
        "--comb_arc_qa_ms"; arg_type = Float64; default = 150.0
        "--coeff_degree"; arg_type = Int; default = 9
    end
    return parse_args(s)
end

const A_ = parse_cli()
const OUT = A_["outroot"]
const PASS1 = A_["pass1"]
const APRED = joinpath(PASS1, "apred")
const WAVECAL = joinpath(PASS1, "wavecal")
const ALMANAC = isempty(A_["almanac"]) ? joinpath(PASS1, "almanac", "allobs_57618_61230.h5") : A_["almanac"]
const CFG = AR.RefitConfig(neon = A_["neon"], neon_offset_file = A_["neon_offset_file"],
    neon_offset_set = A_["neon_offset_set"], neon_line_offset_file = A_["neon_line_offset_file"],
    plate_shape = A_["plate_shape"],
    fps_shape_apo = A_["fps_shape_apo"], fps_shape_lco = A_["fps_shape_lco"],
    anchor_block = A_["anchor_block"], anchor_prefer = A_["anchor_prefer"], anchor_samepos = A_["anchor_samepos"],
    anchor_adjacency = A_["anchor_adjacency"], anchor_max_candidates = A_["anchor_max_candidates"],
    anchor_pattern_tol_px = A_["anchor_pattern_tol_px"], anchor_lamp_tol_px = A_["anchor_lamp_tol_px"], anchor_pairing = A_["anchor_pairing"], anchor_block_h = A_["anchor_block_h"],
    chip_transform_file = A_["chip_transform_file"],
    exp_terms = A_["exp_terms"], sky_fiber_term = A_["sky_fiber_term"], prefps_time = A_["prefps_time"],
    sky_exp_terms = A_["sky_exp_terms"], reg_weight_apo = A_["reg_weight_apo"], reg_weight_lco = A_["reg_weight_lco"],
    shape_template_file = A_["shape_template_file"], shape_template_block = A_["shape_template_block"],
    arcwin_nights = A_["arcwin_nights"], nsig = A_["nsig"], sky_anchor_window_h = A_["sky_anchor_window_h"],
    pos_match_px = A_["pos_match_px"], xcheck_threshold_ms = A_["xcheck_threshold_ms"],
    comb_arc_qa_ms = A_["comb_arc_qa_ms"], coeff_degree = A_["coeff_degree"])
AR.refit_check_config(CFG)
(abspath(OUT) == abspath(PASS1) || startswith(abspath(OUT), abspath(PASS1) * "/")) &&
    error("--outroot must not be inside the pass-1 tree")

const CODE_COMMIT = AR.git_commit
const CODE_CLEAN = AR.git_clean

# ------------------------------------------------------------------------------------------ config
function config_lines()
    d = AR.refit_config_dict(CFG)
    L = ["$(k)=$(d[k])" for k in sort(collect(keys(d)))]
    push!(L, "fingerprint_chip_transform_file=$(AR.refit_fingerprint(CFG.chip_transform_file))")
    push!(L, "fingerprint_neon_offset_file=$(AR.refit_fingerprint(CFG.neon_offset_file))")
    push!(L, "fingerprint_neon_line_offset_file=$(AR.refit_fingerprint(CFG.neon_line_offset_file))")
    push!(L, "fingerprint_shape_template_file=$(AR.refit_fingerprint(CFG.shape_template_file))")
    ll = AR.load_arclamp_linelist()
    push!(L, "linelist_file=$(ll.path)")
    push!(L, "fingerprint_linelist=$(ll.fingerprint)")
    push!(L, "pass1=$(PASS1)")
    push!(L, "almanac=$(ALMANAC)")
    return L
end
function check_or_write_config(write::Bool)
    fn = joinpath(OUT, "config.txt")
    mine = config_lines()
    if isfile(fn)
        theirs = filter(l -> !startswith(l, "#"), readlines(fn))
        if theirs != mine
            println("config mismatch with $(fn):")
            for l in setdiff(mine, theirs); println("  this run : ", l); end
            for l in setdiff(theirs, mine); println("  outroot  : ", l); end
            error("refusing to mix configurations in one outroot; use a new --outroot (with --reuse_from)")
        end
        return
    end
    write || error("$(fn) missing: run --phase prep first")
    mkpath(OUT)
    open(fn * ".tmp", "w") do io
        println(io, "# arc + FPI wavelength refit configuration; written $(now()) by refit_wavecal_archive.jl")
        println(io, "# code commit $(CODE_COMMIT) clean=$(CODE_CLEAN) (switch values below; see REFIT_PLAN.md)")
        foreach(l -> println(io, l), mine)
    end
    mv(fn * ".tmp", fn; force = true)
end

# ------------------------------------------------------------------------------------------ nights
function build_nights()
    rx = r"^wavecalNightAve_(apo|lco)_(\d+)\.h5$"
    rows = [(tele = String(m.captures[1]), mjd = parse(Int, m.captures[2]))
            for m in (match(rx, f) for f in readdir(WAVECAL)) if m !== nothing]
    df = DataFrame(rows)
    # most expensive first (FPS nights with a comb, newest first), then the plate era: workers
    # claim nights dynamically, so this keeps the tail of each phase short
    df.fps = [AR.refit_era(r.tele, r.mjd) == "fps" for r in eachrow(df)]
    sort!(df, [order(:fps, rev = true), order(:mjd, rev = true), :tele])
    select!(df, :tele, :mjd)
    A_["max_nights"] > 0 && (df = df[1:min(end, A_["max_nights"]), :])
    return df
end
function load_nights()
    fn = isempty(A_["nights"]) ? joinpath(OUT, "nights.csv") : A_["nights"]
    df = CSV.read(fn, DataFrame; types = Dict(:tele => String, :mjd => Int))
    return df
end

# ------------------------------------------------------------------------------------------ locks
function claim(phase, tele, mjd)
    lk = joinpath(OUT, "locks", "phase$(phase)", "$(tele)_$(mjd)")
    mkpath(dirname(lk))
    try
        mkdir(lk)
        return true
    catch
        return false
    end
end

# ------------------------------------------------------------------------------------------ reuse
const STAGE_PRODUCT = Dict("arc" => "arc", "comb" => "comb", "fpifib" => "fpifib", "assemble" => "assemble")
function invalidated_stages(base)
    fn = joinpath(base, "config.txt")
    isfile(fn) || error("--reuse_from $(base): no config.txt")
    theirs = Dict(split(l, "=", limit = 2)[1] => split(l, "=", limit = 2)[2]
                  for l in readlines(fn) if !startswith(l, "#") && occursin("=", l))
    mine = AR.refit_config_dict(CFG)
    inv = Set{String}()
    for (k, v) in mine
        get(theirs, k, "") == v && continue
        foreach(s -> push!(inv, s), split(get(AR.REFIT_SWITCH_STAGES, k, "arc,comb,fpifib,assemble"), ","))
    end
    # downstream closure
    order = ["arc", "comb", "fpifib", "assemble"]
    if !isempty(inv)
        first_bad = minimum(findfirst(==(s), order) for s in inv)
        foreach(s -> push!(inv, s), order[first_bad:end])
    end
    return inv
end
const REUSE = A_["reuse_from"]
const INVALID = isempty(REUSE) ? Set{String}() : invalidated_stages(REUSE)
function link_reused(stage, tele, mjd)
    isempty(REUSE) && return
    stage in INVALID && return
    src = AR.refit_path(REUSE, stage, tele, mjd)
    dst = AR.refit_path(OUT, stage, tele, mjd)
    if isfile(src) && !ispath(dst)
        mkpath(dirname(dst)); symlink(src, dst)
    end
    ss = AR.refit_status_path(REUSE, stage, tele, mjd)
    sd = AR.refit_status_path(OUT, stage, tele, mjd)
    if isfile(ss) && !ispath(sd)
        mkpath(dirname(sd)); cp(ss, sd)
    end
end

# ------------------------------------------------------------------------------------------ running a stage
terminal(st) = !(st in ("", "deferred")) && !(A_["retry_failed"] && st == "fail")
function run_stage(stage, tele, mjd, f; force_retry_deferred = false)
    link_reused(stage, tele, mjd)
    st = AR.refit_read_status(OUT, stage, tele, mjd)
    if isfile(AR.refit_path(OUT, stage, tele, mjd)) || (terminal(st) && !(force_retry_deferred && st == "deferred"))
        return (st == "" ? "done" : st), 0.0
    end
    t0 = time()
    res = try
        f()
    catch err
        msg = sprint(showerror, err)
        bt = sprint(Base.show_backtrace, catch_backtrace())
        println(stderr, "FAIL $(stage) $(tele) $(mjd): $(msg)\n", first(bt, 2000))
        AR.refit_write_status(OUT, stage, tele, mjd, "fail", first(msg, 400))
    end
    return String(res), time() - t0
end

function log_timing(phase, tele, mjd, recs)
    d = joinpath(OUT, "timing")
    mkpath(d)
    fn = joinpath(d, "phase$(phase)_w$(A_["worker"])_$(gethostname()).csv")
    new = !isfile(fn)
    open(fn, "a") do io
        new && println(io, "tele,mjd,stage,status,seconds,maxrss_mb,time")
        for (stage, st, dt) in recs
            println(io, "$(tele),$(mjd),$(stage),$(st),$(round(dt, digits = 2)),$(round(Sys.maxrss() / 2^20, digits = 0)),$(now())")
        end
    end
end

function phase12(phase)
    check_or_write_config(false)
    nights = load_nights()
    LL = Dict(t => AR.refit_linelist(CFG, t)[1] for t in ("apo", "lco"))
    n_done = 0
    t_start = time()
    for r in eachrow(nights)
        tele, mjd = String(r.tele), Int(r.mjd)
        claim(phase, tele, mjd) || continue
        foreach(st -> link_reused(st, tele, mjd), AR.REFIT_STAGES)   # products of valid stages from --reuse_from
        recs = Tuple{String, String, Float64}[]
        if phase == 1
            st, dt = run_stage("arc", tele, mjd, () -> AR.refit_night_arc(tele, mjd; apred_dir = APRED,
                almanac_file = ALMANAC, outroot = OUT, cfg = CFG))
            push!(recs, ("arc", st, dt))
            st, dt = run_stage("comb", tele, mjd, () -> AR.refit_night_comb(tele, mjd; apred_dir = APRED,
                wavecal_dir = WAVECAL, outroot = OUT, cfg = CFG, allow_arcwin = false, almanac_file = ALMANAC))
            push!(recs, ("comb", st, dt))
        else
            st, dt = run_stage("comb", tele, mjd, () -> AR.refit_night_comb(tele, mjd; apred_dir = APRED,
                wavecal_dir = WAVECAL, outroot = OUT, cfg = CFG, allow_arcwin = true, almanac_file = ALMANAC); force_retry_deferred = true)
            push!(recs, ("comb2", st, dt))
            st, dt = run_stage("fpifib", tele, mjd, () -> AR.refit_night_fpifib(tele, mjd; apred_dir = APRED,
                outroot = OUT, cfg = CFG))
            push!(recs, ("fpifib", st, dt))
            st, dt = run_stage("assemble", tele, mjd, () -> AR.refit_night_assemble(tele, mjd; apred_dir = APRED,
                wavecal_dir = WAVECAL, outroot = OUT, cfg = CFG, ll = LL[tele]))
            push!(recs, ("assemble", st, dt))
        end
        log_timing(phase, tele, mjd, recs)
        n_done += 1
        println(@sprintf("[w%d %s] phase %d %s %d  %s  (%d nights, %.1f min)", A_["worker"], gethostname(), phase,
            tele, mjd, join(["$(s)=$(x)($(round(d, digits = 1))s)" for (s, x, d) in recs], " "), n_done,
            (time() - t_start) / 60))
        flush(stdout)
    end
    println("WORKER $(A_["worker"]) PHASE $(phase) DONE ($(n_done) nights)")
end

# ------------------------------------------------------------------------------------------ phase 3
function phase3()
    check_or_write_config(false)
    nights = load_nights()
    rows = DataFrame()
    combrows = DataFrame()
    for r in eachrow(nights)
        tele, mjd = String(r.tele), Int(r.mjd)
        st = Dict(s => AR.refit_read_status(OUT, s, tele, mjd) for s in AR.REFIT_STAGES)
        row = Dict{Symbol, Any}(:tele => tele, :mjd => mjd, :era => (AR.refit_era(tele, mjd) == "plate" ? "plate" :
            (tele == "apo" ? (mjd < 60930 ? "fps1" : "fps2") : (mjd < 60700 ? "fps1" : "fps2"))))
        for s in AR.REFIT_STAGES; row[Symbol("status_" * s)] = st[s]; end
        ef = AR.refit_path(OUT, "assemble", tele, mjd)
        row[:n_exp] = 0
        if isfile(ef)
            h5open(ef, "r") do f
                wt = read(f["wavecal_type"]); bm = read(f["bitmsk_wavecal"])
                row[:n_exp] = length(wt)
                for w in unique(wt); row[Symbol("n_" * w)] = count(==(w), wt); end
                row[:n_pass1_copy] = count(b -> (b & AR.WAVE_PASS1_COPY) != 0, bm)
                row[:n_level_inconsistent] = count(b -> (b & AR.WAVE_LEVEL_INCONSISTENT) != 0, bm)
                row[:n_fpifib_unavailable] = count(b -> (b & AR.WAVE_FPIFIB_UNAVAILABLE) != 0, bm)
                row[:n_fps_no_comb] = count(b -> (b & AR.WAVE_FPS_NO_COMB) != 0, bm)
                z = read(f["pass1_zp_ms"])
                for c in 1:3; row[Symbol("med_dzp_pass1_c$c")] = AR._arc_nanmedian(z[c, :]); end
                x = read(f["xcheck_resid_ms"])
                row[:med_abs_xcheck_ms] = AR._arc_nanmedian(abs.(Float64.(x)))
            end
        end
        cf = AR.refit_path(OUT, "comb", tele, mjd)
        if isfile(cf)
            h5open(cf, "r") do f
                cr = Dict{Symbol, Any}(:tele => tele, :mjd => mjd, :anchor => String(read(f["metadata/anchor"])),
                    :L => read(f["fpi_cavity_size"]), :eps => read(f["fpi_m0_offset"]), :L_pass1 => read(f["pass1_cavity_size"]),
                    :anchor_kind => haskey(f, "anchor_kind") ? String(read(f["anchor_kind"])) : "",
                    :anchor_bits => haskey(f, "anchor_bits") ? Int(read(f["anchor_bits"])) : 0,
                    :n_slip => read(f["n_mode_slip"]), :n_order_changed => read(f["n_fiberchip_order_changed"]),
                    :rms_mA => read(f["comb_resid_rms_mA"]))
                lin = read(f["linParams"]); nl = read(f["nlParams"])
                cr[:lam_G_f150_p1024] = AR.refit_prod_wave(lin, nl, 150, 2, [1024.0])[1]
                a = read(f["chip_arc_minus_comb_ms"])
                for c in 1:3; cr[Symbol("arc_minus_comb_c$c")] = a[c]; row[Symbol("arc_minus_comb_c$c")] = a[c]; end
                row[:n_mode_slip] = cr[:n_slip]; row[:n_order_changed] = cr[:n_order_changed]; row[:comb_anchor] = cr[:anchor]
                push!(combrows, cr; cols = :union)
            end
        end
        push!(rows, row; cols = :union)
    end
    # cross-night comb checks: (ii) FSR continuity (lambda at a fixed fiber/pixel jumping by ~ an
    # integer FSR between consecutive comb nights) and (iii) cavity L against the trailing 15 nights (> 8 sigma)
    flags = DataFrame(tele = String[], mjd = Int[], check = String[], value = Float64[])
    for tele in ("apo", "lco")
        C = sort(combrows[combrows.tele .== tele, :], :mjd)
        nrow(C) < 2 && continue
        for i in 2:nrow(C)
            fsr = C.lam_G_f150_p1024[i]^2 / (2 * C.L[i])
            d = (C.lam_G_f150_p1024[i] - C.lam_G_f150_p1024[i - 1]) / fsr
            k = round(d)
            (k != 0 && abs(d - k) < 0.2) && push!(flags, (tele, C.mjd[i], "fsr_continuity", d))
            lo = max(1, i - 15)
            if i - lo >= 5
                prev = C.L[lo:(i - 1)]
                m = median(prev); s = max(1.4826 * median(abs.(prev .- m)), 1.0)
                abs(C.L[i] - m) > 8s && push!(flags, (tele, C.mjd[i], "cavity_jump_8sigma", (C.L[i] - m) / s))
            end
        end
    end
    CSV.write(joinpath(OUT, "run_table.csv"), rows)
    CSV.write(joinpath(OUT, "comb_table.csv"), combrows)
    CSV.write(joinpath(OUT, "crossnight_flags.csv"), flags)
    open(joinpath(OUT, "run_summary.txt"), "w") do io
        println(io, "arc + FPI refit run summary, $(now()); outroot $(OUT); code $(CODE_COMMIT)")
        for s in AR.REFIT_STAGES
            c = Dict{String, Int}()
            for x in rows[!, Symbol("status_" * s)]; c[x] = get(c, x, 0) + 1; end
            println(io, "status $(s): ", join(["$(k == "" ? "(none)" : k)=$(v)" for (k, v) in sort(collect(c))], " "))
        end
        cols = [n for n in names(rows) if startswith(n, "n_") && n != "n_exp"]
        for n in sort(cols)
            v = coalesce.(rows[!, n], 0)
            println(io, @sprintf("%-40s %8d", n, sum(v)))
        end
        println(io, "exposures: ", sum(rows.n_exp))
        println(io, "cross-night comb flags: ", nrow(flags))
    end
    print(read(joinpath(OUT, "run_summary.txt"), String))
end

# ------------------------------------------------------------------------------------------ main
ph = A_["phase"]
if ph == "prep"
    check_or_write_config(true)
    nf = joinpath(OUT, "nights.csv")
    isfile(nf) || CSV.write(nf, build_nights())
    rm(joinpath(OUT, "locks"); force = true, recursive = true)
    foreach(d -> mkpath(joinpath(OUT, d)), ("arc", "comb", "fpifib", "exp", "status", "timing", "logs"))
    println("prep done: $(nrow(load_nights())) nights; config $(joinpath(OUT, "config.txt"))")
elseif ph == "1" || ph == "2"
    phase12(parse(Int, ph))
elseif ph == "3"
    phase3()
else
    error("unknown --phase $(ph)")
end
