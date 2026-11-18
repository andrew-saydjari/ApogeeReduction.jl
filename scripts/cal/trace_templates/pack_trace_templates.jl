# pack_trace_templates.jl — pack the per-unit template builds and the epoch table into ONE
# trace-template calibration file (format: src/trace_template_cal.jl). Rerunnable: it reads
# whatever units exist, and with --allow_incomplete true writes a file marked incomplete that
# lists the missing units (lookup then returns nothing for them = template-free behaviour).
#
#   julia --project=<repo> pack_trace_templates.jl --epochs trace_epochs.csv \
#       --units_dir <dir with template_<tele>_<chip>_<epoch>.h5> --out <file.h5> --version <name> \
#       [--open_last_epoch true] [--allow_incomplete false]
#
# The epochs CSV has columns tele,epoch,mjd_start,mjd_end (inclusive; contiguous per telescope).
# Also writes <out without .h5>_units.csv (one row per expected unit: present or missing, and the
# build summary) next to the file. Inputs are opened read-only; nothing is installed anywhere.
using ArgParse, CSV, DataFrames, Dates, LibGit2
using ApogeeReduction
using ApogeeReduction: TraceEpoch, TraceTemplate, read_template_unit_file, write_trace_template_cal,
                       load_trace_template_cal, TRACE_EPOCH_OPEN_END
using SHA: sha256

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--epochs"
        help = "epoch table CSV: tele,epoch,mjd_start,mjd_end"
        arg_type = String
        required = true
        "--units_dir"
        help = "directory with the builder outputs template_<tele>_<chip>_<epoch>.h5 and _log.csv"
        arg_type = String
        required = true
        "--out"
        help = "calibration file to write"
        arg_type = String
        required = true
        "--version"
        help = "calibration_version string recorded in the file (e.g. tracetemplates_pass1dr21_v1)"
        arg_type = String
        required = true
        "--open_last_epoch"
        help = "serve every MJD after the last epoch's mjd_end with the last epoch (open-ended; nights after its newest input flat are flagged provisional). false: the table ends at mjd_end and later nights get no template. (my choice: true; needs AKS)"
        arg_type = Bool
        default = true
        "--unlit_positions"
        help = "template positions of the fibers dark in the epoch: interpolate (re-place, between their lit neighbours, the dark fibers the builder never fitted or that left their slot; faint lit fibers it fitted are kept; interpolate_unlit_template in src/trace_template_cal.jl) or builder (keep the builder's values; a unit whose dark fiber left its slot then fails validation and is reported missing). (my choice: interpolate; needs AKS)"
        arg_type = String
        default = "interpolate"
        "--allow_incomplete"
        help = "write even if some (tele, epoch, chip) units are missing (listed in the file and the units CSV)"
        arg_type = Bool
        default = false
    end
    return parse_args(s)
end
parg = parse_commandline()
parg["unlit_positions"] in ("interpolate", "builder") || error("--unlit_positions must be interpolate or builder")

E = CSV.read(parg["epochs"], DataFrame; types = Dict(:tele => String, :epoch => String))
epochs = Dict{String, Vector{TraceEpoch}}()
for tele in unique(E.tele)
    sub = sort(E[E.tele .== tele, :], :mjd_start)
    n = nrow(sub)
    epochs[tele] = [TraceEpoch(tele, sub.epoch[i], sub.mjd_start[i],
                        (i == n && parg["open_last_epoch"]) ? TRACE_EPOCH_OPEN_END : sub.mjd_end[i],
                        i == n && parg["open_last_epoch"]) for i in 1:n]
    ApogeeReduction.check_trace_epochs(tele, epochs[tele])
end

templates = TraceTemplate[]
rows = DataFrame()
for tele in sort(collect(keys(epochs))), ep in epochs[tele], chip in CHIP_LIST
    base = joinpath(parg["units_dir"], "template_$(tele)_$(chip)_$(ep.name)")
    present = isfile(base * ".h5")
    r = (tele = tele, chip = chip, epoch = ep.name, mjd_start = ep.mjd_start, mjd_end = ep.mjd_end,
        present = present, n_alllit = -1, n_candidates = -1, data_mjd_min = -1, data_mjd_max = -1,
        drift_pattern_rms_px = NaN, last_change_px = NaN, iterations = -1, reg_d_median = NaN,
        modal_missing_slots = "", unlit_repaired = "", invalid_reason = "", fingerprint = "",
        source_file = present ? base * ".h5" : "")
    if present
        t = read_template_unit_file(base * ".h5"; epoch = ep, log_csv = base * "_log.csv")
        rep = Int[]
        parg["unlit_positions"] == "interpolate" && ((t, rep) = ApogeeReduction.interpolate_unlit_template(t))
        bad = try
            ApogeeReduction.check_trace_template(t, ep); ""
        catch e
            sprint(showerror, e)
        end
        if bad != ""
            @warn "unit $(tele)_$(chip)_$(ep.name) fails validation and is left out: $(bad)"
            push!(rows, merge(r, (present = false, invalid_reason = bad)))
            continue
        end
        isempty(rep) || println("$(tele)_$(chip)_$(ep.name): dark fibers re-placed between lit neighbours: $(join(rep, ", "))")
        push!(templates, t)
        pv = t.provenance
        r = merge(r, (n_alllit = parse(Int, get(pv, "n_alllit", "-1")), n_candidates = parse(Int, get(pv, "n_candidates", "-1")),
            data_mjd_min = t.data_mjd_min, data_mjd_max = t.data_mjd_max,
            drift_pattern_rms_px = parse(Float64, get(pv, "drift_pattern_rms_px", "NaN")),
            last_change_px = parse(Float64, get(pv, "last_change_px", "NaN")),
            iterations = parse(Int, get(pv, "iterations", "-1")), reg_d_median = t.reg_d_median,
            modal_missing_slots = get(pv, "modal_missing_slots", ""), unlit_repaired = join(rep, ";"),
            fingerprint = t.fingerprint))
    end
    push!(rows, r)
end
nmiss = count(.!rows.present)
println("units present $(nrow(rows) - nmiss) of $(nrow(rows))", nmiss > 0 ? "; missing: " * join(["$(r.tele)_$(r.chip)_$(r.epoch)" for r in eachrow(rows) if !r.present], ", ") : "")

commit = try
    repo = LibGit2.GitRepo(pkgdir(ApogeeReduction))
    string(LibGit2.head_oid(repo)) * (LibGit2.isdirty(repo) ? "+dirty" : "")
catch
    "unknown"
end
prov = Dict{String, String}(
    "created_utc" => string(now(UTC)),
    "packer" => "scripts/cal/trace_templates/pack_trace_templates.jl",
    "packer_commit" => commit,
    "epochs_csv" => abspath(parg["epochs"]),
    "epochs_csv_sha256" => bytes2hex(sha256(read(parg["epochs"]))),
    "units_dir" => abspath(parg["units_dir"]),
    "open_last_epoch" => string(parg["open_last_epoch"]),
    "unlit_positions" => parg["unlit_positions"],
    "builder_code_sha256" => join(sort(unique([get(t.provenance, "code_sha256", "?") for t in templates])), ";"),
    "n_units_present" => string(length(templates)),
    "needs_trace_template_file" => "no: the builder reads only 2D quartz flats (README, Bootstrap)")
mkpath(dirname(abspath(parg["out"])))
write_trace_template_cal(parg["out"], epochs, templates, prov; calibration_version = parg["version"],
    allow_incomplete = parg["allow_incomplete"])
cal = load_trace_template_cal(parg["out"])          # read back through the pipeline's loader
unitcsv = replace(parg["out"], r"\.h5$" => "") * "_units.csv"
CSV.write(unitcsv, rows)
println("wrote $(parg["out"]) ($(filesize(parg["out"])) bytes): $(cal)")
println("content fingerprint $(cal.fingerprint); id $(ApogeeReduction.trace_template_cal_id(cal)); units table $(unitcsv)")
