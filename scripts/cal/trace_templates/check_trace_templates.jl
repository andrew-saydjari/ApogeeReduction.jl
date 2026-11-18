# check_trace_templates.jl — QA of a trace-template calibration file and of the unit builds it
# was packed from. Read-only on every input.
#
#   julia --project=<repo> check_trace_templates.jl --cal <packed.h5> --units_dir <builder outputs> \
#       --outdir <qa dir> [--reference <previous packed.h5>] [--require_complete false] \
#       [--max_dc_flag 0.3]
#
# Per unit (telescope, chip, epoch):
#   - coverage: present in the file, served MJD range, data range, all-lit / candidate flats,
#     modal dark slots, dark fibers re-placed at packing;
#   - drift within the epoch: rigid drift (range and NMAD of the per-flat registration shift d at
#     x = 1024, px; absorbed by the nightly registration), pattern drift (last third - first third,
#     px rms; what a template cannot absorb), per-flat residual about the template (resid_rms, px;
#     flats above 3x the unit median are listed);
#   - leave-nights-out: held-out fit_rms, NMAD and max of (held-out template - all-flats template)
#     centers on each held-out flat (px);
#   - no-template mode: NMAD / p99 / max of (no-template - template) centers per held-out flat (px).
# Held-out flats whose max difference exceeds --max_dc_flag px are listed for inspection
# (qa_flagged_flats.csv); README "QA" explains how they were diagnosed for the first build.
# With --reference: per common unit, the template change (common shift per band removed, lit
# fibers), px rms / max.
# Writes qa_units.csv, qa_flagged_flats.csv (and qa_vs_reference.csv) to --outdir and prints a summary.
using ArgParse, CSV, DataFrames, Statistics, Printf
using ApogeeReduction
using ApogeeReduction: load_trace_template_cal

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--cal"
        arg_type = String
        required = true
        "--units_dir"
        help = "builder outputs (template_*_log.csv, _loo.csv, _notemplate.csv)"
        arg_type = String
        required = true
        "--outdir"
        arg_type = String
        required = true
        "--reference"
        help = "previous calibration file to compare with"
        arg_type = String
        default = ""
        "--require_complete"
        arg_type = Bool
        default = false
        "--max_dc_flag"
        help = "list held-out flats whose largest center difference exceeds this (px)"
        arg_type = Float64
        default = 0.3
    end
    return parse_args(s)
end
parg = parse_commandline()
mkpath(parg["outdir"])
cal = load_trace_template_cal(parg["cal"]; require_complete = parg["require_complete"])
println(cal)
nmad(v) = (w = filter(isfinite, v); length(w) < 2 ? NaN : 1.4826 * median(abs.(w .- median(w))))
fmed(v) = (w = filter(isfinite, v); isempty(w) ? NaN : median(w))
fmax(v) = (w = filter(isfinite, v); isempty(w) ? NaN : maximum(w))
readcsv(p) = isfile(p) ? CSV.read(p, DataFrame) : nothing

units = DataFrame(); flagged = DataFrame()
for tele in sort(collect(keys(cal.epochs))), ep in cal.epochs[tele], chip in CHIP_LIST
    t = get(cal.templates, (tele, chip, ep.name), nothing)
    base = joinpath(parg["units_dir"], "template_$(tele)_$(chip)_$(ep.name)")
    lg = readcsv(base * "_log.csv"); lo = readcsv(base * "_loo.csv"); nt = readcsv(base * "_notemplate.csv")
    fin = isnothing(lg) ? nothing : lg[lg.stage .== "final", :]
    r = (tele = tele, chip = chip, epoch = ep.name, mjd_start = ep.mjd_start, mjd_end = ep.mjd_end,
        in_file = !isnothing(t), n_alllit = isnothing(t) ? -1 : length(t.input_mjd),
        n_candidates = isnothing(lg) ? -1 : count(lg.stage .== "census"),
        data_mjd_min = isnothing(t) ? -1 : t.data_mjd_min, data_mjd_max = isnothing(t) ? -1 : t.data_mjd_max,
        modal_dark = isnothing(t) ? "" : join(findall(.!t.lit), ";"),
        repaired = isnothing(t) ? "" : get(t.provenance, "unlit_repaired_fibers", ""),
        rigid_d_range_px = isnothing(fin) ? NaN : fmax(fin.d_center) - (-fmax(-fin.d_center)),
        rigid_d_nmad_px = isnothing(fin) ? NaN : nmad(fin.d_center),
        drift_pattern_rms_px = isnothing(t) ? NaN : parse(Float64, get(t.provenance, "drift_pattern_rms_px", "NaN")),
        resid_rms_med_px = isnothing(fin) ? NaN : fmed(fin.resid_rms),
        n_resid_gt3med = isnothing(fin) ? -1 : count(fin.resid_rms .> 3 * fmed(fin.resid_rms)),
        fit_rms_med_px = isnothing(fin) ? NaN : fmed(fin.fit_rms),
        loo_n = isnothing(lo) ? 0 : nrow(lo),
        loo_fit_rms_med_px = isnothing(lo) ? NaN : fmed(lo.fit_rms_heldout),
        loo_nmad_dc_med_px = isnothing(lo) ? NaN : fmed(lo.nmad_dc),
        loo_max_dc_px = isnothing(lo) ? NaN : fmax(lo.max_dc),
        loo_k_nonzero = isnothing(lo) ? -1 : count(lo.k .!= 0),
        nt_n = isnothing(nt) ? 0 : nrow(nt),
        nt_nmad_dc_med_px = isnothing(nt) ? NaN : fmed(nt.nmad_dc),
        nt_p99_dc_med_px = isnothing(nt) ? NaN : fmed(nt.p99_dc),
        nt_max_dc_px = isnothing(nt) ? NaN : fmax(nt.max_dc),
        nt_failed = isnothing(nt) ? -1 : count(nt.k_notemplate .!= 0))
    push!(units, r)
    for (df, kind) in ((lo, "leave_nights_out"), (nt, "no_template"))
        isnothing(df) && continue
        for x in eachrow(df)
            (isfinite(x.max_dc) && x.max_dc > parg["max_dc_flag"]) || continue
            push!(flagged, (tele = tele, chip = chip, epoch = ep.name, test = kind, mjd = x.mjd,
                max_dc_px = round(x.max_dc, digits = 3), nmad_dc_px = round(x.nmad_dc, digits = 4)))
        end
    end
    if !isnothing(fin)
        m = fmed(fin.resid_rms)
        for x in eachrow(fin)
            x.resid_rms > 3m || continue
            push!(flagged, (tele = tele, chip = chip, epoch = ep.name, test = "resid_rms>3x_median", mjd = x.mjd,
                max_dc_px = round(x.resid_rms, digits = 3), nmad_dc_px = NaN))
        end
    end
end
CSV.write(joinpath(parg["outdir"], "qa_units.csv"), units)
CSV.write(joinpath(parg["outdir"], "qa_flagged_flats.csv"), flagged)

@printf("%-4s %-1s %-6s %5s %9s %9s %9s %9s %9s %9s %9s\n", "tele", "c", "epoch", "nflat", "rigidNMAD",
    "drift", "LOOnmad", "LOOmax", "NTnmad", "NTp99", "NTmax")
for r in eachrow(units)
    @printf("%-4s %-1s %-6s %5d %9.4f %9.4f %9.4f %9.3f %9.4f %9.4f %9.3f%s\n", r.tele, r.chip, r.epoch, r.n_alllit,
        r.rigid_d_nmad_px, r.drift_pattern_rms_px, r.loo_nmad_dc_med_px, r.loo_max_dc_px, r.nt_nmad_dc_med_px,
        r.nt_p99_dc_med_px, r.nt_max_dc_px, r.in_file ? "" : "  MISSING")
end
println("flagged held-out flats (> $(parg["max_dc_flag"]) px) and residual outliers: $(nrow(flagged)); see qa_flagged_flats.csv")

if parg["reference"] != ""
    ref = load_trace_template_cal(parg["reference"])
    cmp = DataFrame()
    for (key, t) in cal.templates
        r0 = get(ref.templates, key, nothing)
        isnothing(r0) && continue
        size(r0.center) == size(t.center) || continue
        l = t.lit .& r0.lit
        d = t.center[:, l] .- r0.center[:, l]
        d = d .- median(d, dims = 2)
        push!(cmp, (tele = key[1], chip = key[2], epoch = key[3], same_fingerprint = t.fingerprint == r0.fingerprint,
            rms_px = sqrt(mean(d .^ 2)), max_px = maximum(abs.(d))))
    end
    CSV.write(joinpath(parg["outdir"], "qa_vs_reference.csv"), cmp)
    println("vs reference $(ref.calibration_version): $(nrow(cmp)) common units; max rms $(round(maximum(cmp.rms_px; init = 0.0), digits = 4)) px")
end
