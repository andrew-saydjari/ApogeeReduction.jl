# Own-night fallback levels for a night with no eligible flat (TRACE_ROBUST_V2 §5.2, §6), all from
# the night's own frames plus the epoch template; never another night's trace:
#   T1  the epoch template registered on another calibration frame of the night
#       (dome flat, UNe / ThAr / FPI arcs, internal flat)        -> TRACE_REGISTERED_CAL
#   T2  registered on the median of up to 8 science frames         -> TRACE_REGISTERED_SCI
#   T3  the template at the epoch-median registration              -> TRACE_TEMPLATE_ONLY (bad)
# Registration-only (`register_frame(...; mode = :register_only)`): label shift fixed at 0,
# |d - d_epoch| <= 3 px, no per-fiber delta. Each registered frame gets a QA record (as the flats
# do, PR A) under `<outdir>/trace_fallback/<mjd>/`, so a re-run reuses it and an audit can read
# why a frame was or was not used. Runs inside the per-(night, chip) `make_traceMain` task, after
# the 2D stage, and only on nights whose flats all failed (no change to staging or runlists).

"Default T1 frame kinds, in order of preference (V2 §5.2 measured registration quality: dome best)."
const TRACE_T1_KINDS = ["domeflat", "arc_UNe", "arc_ThAr", "arc_FPI", "internalflat"]
"Maximum number of science frames stacked for T2 (V2 §5.2)."
const TRACE_T2_NMAX = 8

"""
    parse_trace_fallback(s) -> (levels::Vector{Int}, fixed_file::Bool)

`--trace_fallback` values: `none` (no trace), `fixed_file` (the old survey-wide file),
`template` (= `T1,T2,T3`), or a comma list of own-night levels, e.g. `T1,T2` (stop before the
template-only level: such a night is then not extracted).
"""
function parse_trace_fallback(s::AbstractString)
    s = strip(s)
    s == "none" && return (Int[], false)
    s == "fixed_file" && return (Int[], true)
    s == "template" && return ([TRACE_LEVEL_REGISTERED_CAL, TRACE_LEVEL_REGISTERED_SCI, TRACE_LEVEL_TEMPLATE], false)
    lv = Int[]
    for t in split(s, ",")
        t = strip(t)
        t == "T1" ? push!(lv, TRACE_LEVEL_REGISTERED_CAL) :
        t == "T2" ? push!(lv, TRACE_LEVEL_REGISTERED_SCI) :
        t == "T3" ? push!(lv, TRACE_LEVEL_TEMPLATE) :
        throw(ArgumentError("--trace_fallback must be none, fixed_file, template or a comma list of T1,T2,T3; got \"$(s)\""))
    end
    issorted(lv) && allunique(lv) || throw(ArgumentError("--trace_fallback levels must be increasing and distinct: \"$(s)\""))
    return (lv, false)
end
"True if this `--trace_fallback` value uses the epoch templates (needs `--trace_template_cal`)."
trace_fallback_uses_templates(s::AbstractString) = !isempty(first(parse_trace_fallback(s)))

_toint(x) = x isa AbstractString ? (isempty(strip(x)) ? 0 : parse(Int, x)) : Int(x)

"""
    night_calibration_frames(outdir, tele, mjd, chip; almanac_file, kinds = TRACE_T1_KINDS)

The night's 2D calibration frames usable for T1 (`ar2Dcal_<tele>_<mjd>_<exp>_<chip>_<type>.h5` in
`apred/<mjd>/`), as `(kind, expnum, path)` sorted by `kinds` then exposure. Arc lamps are split by
the almanac lamp columns (row = exposure index): `arc_ThAr` (ThAr on), `arc_UNe` (UNe on), else
`arc_FPI`; without an almanac (or a row for the exposure) an arc is skipped.
"""
function night_calibration_frames(outdir, tele, mjd, chip; almanac_file = "", kinds = TRACE_T1_KINDS)
    dir = joinpath(outdir, "apred", "$(mjd)")
    isdir(dir) || return NamedTuple{(:kind, :expnum, :path), Tuple{String, Int, String}}[]
    lamps = Dict{Int, Tuple{Int, Int}}()
    if almanac_file != "" && isfile(almanac_file)
        try
            df = read_almanac_exp_df(almanac_file, tele, mjd)
            th = hasproperty(df, :lamp_thar) ? _toint.(df.lamp_thar) : zeros(Int, nrow(df))
            un = hasproperty(df, :lamp_une) ? _toint.(df.lamp_une) : zeros(Int, nrow(df))
            for i in 1:nrow(df)
                lamps[i] = (th[i], un[i])
            end
        catch e
            @warn "night_calibration_frames: could not read lamps from $(almanac_file) for $(tele) $(mjd): $(sprint(showerror, e))"
        end
    end
    rx = Regex("^ar2Dcal_$(tele)_$(mjd)_(\\d+)_$(chip)_([a-z]+)\\.h5\$")
    out = NamedTuple{(:kind, :expnum, :path), Tuple{String, Int, String}}[]
    for fn in readdir(dir)
        m = match(rx, fn)
        isnothing(m) && continue
        e = parse(Int, m.captures[1])
        ty = m.captures[2]
        kind = if ty == "arclamp"
            haskey(lamps, e) || continue
            lamps[e][1] == 1 ? "arc_ThAr" : lamps[e][2] == 1 ? "arc_UNe" : "arc_FPI"
        else
            ty
        end
        kind in kinds || continue
        push!(out, (kind = kind, expnum = e, path = joinpath(dir, fn)))
    end
    sort!(out, by = r -> (findfirst(==(r.kind), kinds), r.expnum))
    return out
end

"""
    night_science_frames(outdir, tele, mjd, chip; nmax = 8)

Up to `nmax` of the night's 2D science frames (`*_object.h5`), evenly spaced over the night.
"""
function night_science_frames(outdir, tele, mjd, chip; nmax = TRACE_T2_NMAX)
    dir = joinpath(outdir, "apred", "$(mjd)")
    isdir(dir) || return String[]
    rx = Regex("^ar2Dcal_$(tele)_$(mjd)_(\\d+)_$(chip)_object\\.h5\$")
    fs = sort([joinpath(dir, fn) for fn in readdir(dir) if !isnothing(match(rx, fn))])
    length(fs) <= nmax && return fs
    return fs[unique(round.(Int, range(1, length(fs), length = nmax)))]
end

function _load_frame(path)
    h5open(path, "r") do f
        d = f["dimage"][1:N_XPIX, 1:N_XPIX]
        iv = f["ivarimage"][1:N_XPIX, 1:N_XPIX]
        m = f["pix_bitmask"][1:N_XPIX, 1:N_XPIX]
        return d, iv, m
    end
end

"Band profiles of one 2D frame for registration (arcs: 65-column bands, no clipping)."
function frame_bands(path, kind; nband = 32)
    d, iv, m = _load_frame(path)
    good = ((m .& bad_pix_bits) .== 0) .& (iv .> 0) .& isfinite.(d)
    line = startswith(kind, "arc")
    wide = line || kind == "object"
    return slit_bands(d, 1 ./ sqrt.(iv), good; nb = nband, half = wide ? 32 : 8,
        nsigma = line ? Inf : 5)
end

trace_fallback_dir(outdir, mjd) = joinpath(outdir, "trace_fallback", "$(mjd)")
regcal_qa_path(outdir, tele, mjd, expnum, chip) = joinpath(trace_fallback_dir(outdir, mjd),
    "regcalTraceQA_$(tele)_$(mjd)_$(lpad(expnum, 4, '0'))_$(chip).h5")
regsci_qa_path(outdir, tele, mjd, chip) = joinpath(trace_fallback_dir(outdir, mjd),
    "regsciTraceQA_$(tele)_$(mjd)_$(chip).h5")

"QA record of a registration-only fit (T1/T2)."
function registration_qa(rf, tmpl, kind, src; gate_kw = (;))
    ok, reasons = registration_gate(rf; gate_kw...)
    qa = template_qa_fields(rf; tmpl = tmpl)
    qa["kind"] = String(kind)
    qa["source"] = String(src)
    qa["outcome"] = String(rf.outcome)
    qa["reg_gate_pass"] = ok
    qa["reg_gate_reasons"] = join(reasons, ";")
    qa["fit_rms"] = Float64(rf.fit_rms)
    qa["nlit"] = Int(rf.nlit)
    return qa
end

function _cached_qa(path, tmpl_id, checkpoint_mode)
    check_file(path, mode = checkpoint_mode) || return nothing
    qa = read_trace_qa(path)
    (isnothing(qa) || string(get(qa, "tmpl_id", "")) != tmpl_id) && return nothing
    return qa
end

"""
    trace_fallback_chain(outdir, tele, mjd, chip; levels, template_cal, almanac_file,
        data_path, checkpoint_mode, t1_kinds, gate_kw)

Run the own-night levels in `levels` (subset of T1, T2, T3) in order for a night whose flats are
all ineligible. Returns `nothing` (no template for this night, or no level produced a trace) or
a NamedTuple `(level, reg::TemplateRegistration, tmpl, source, qa, info)`; `info` holds the
candidate table of every frame tried (kind, exposure, gate, fit_rms), for the traceMain metadata.
T1 ranks passing frames by kind (`t1_kinds` order), then `fit_rms`, then the later exposure; the
kinds after the first one with a passing frame are not registered (they cannot win). QA records
go under `qa_root` (default `outdir`) `/trace_fallback/<mjd>/`.
"""
function trace_fallback_chain(outdir, tele, mjd, chip; levels, template_cal, almanac_file = "",
        data_path = "./data/", checkpoint_mode = "commit_same", t1_kinds = TRACE_T1_KINDS,
        t2_nmax = TRACE_T2_NMAX, gate_kw = (;), qa_root = outdir)
    tmpl = trace_template(template_cal, tele, chip, parse(Int, string(mjd)))
    if isnothing(tmpl)
        @warn "trace_fallback_chain: no epoch template serves $(tele) $(chip) $(mjd); own-night levels skipped"
        return nothing
    end
    tid = template_id(tmpl)
    P = nothing
    getP() = (isnothing(P) && (P = slit_profile(tele, chip; data_path = data_path)); P)
    info = Dict{String, Any}("trace_fb_cand_source" => String[], "trace_fb_cand_kind" => String[],
        "trace_fb_cand_pass" => Bool[], "trace_fb_cand_fit_rms" => Float64[],
        "trace_fb_cand_reasons" => String[])
    function note!(qa)
        push!(info["trace_fb_cand_source"], basename(string(qa["source"])))
        push!(info["trace_fb_cand_kind"], string(qa["kind"]))
        push!(info["trace_fb_cand_pass"], Bool(qa["reg_gate_pass"]))
        push!(info["trace_fb_cand_fit_rms"], Float64(qa["fit_rms"]))
        push!(info["trace_fb_cand_reasons"], string(qa["reg_gate_reasons"]))
    end
    for level in levels
        if level == TRACE_LEVEL_REGISTERED_CAL
            frames = night_calibration_frames(outdir, tele, mjd, chip; almanac_file = almanac_file,
                kinds = t1_kinds)
            qas = Dict{String, Any}[]
            for (i, fr) in enumerate(frames)
                # kinds are tried in order of preference: once a kind has a passing frame, the
                # less preferred kinds cannot win the ranking and are not registered
                if i > 1 && fr.kind != frames[i - 1].kind && any(q -> Bool(q["reg_gate_pass"]), qas)
                    break
                end
                qp = regcal_qa_path(qa_root, tele, mjd, fr.expnum, chip)
                qa = _cached_qa(qp, tid, checkpoint_mode)
                if isnothing(qa)
                    qa = try
                        xc, prof, err = frame_bands(fr.path, fr.kind)
                        rf = register_frame(getP(), tmpl, xc, prof, err; mode = :register_only)
                        registration_qa(rf, tmpl, fr.kind, fr.path; gate_kw = gate_kw)
                    catch e
                        Dict{String, Any}("kind" => fr.kind, "source" => fr.path,
                            "outcome" => "fit_error", "reg_gate_pass" => false,
                            "reg_gate_reasons" => "fit_error", "fit_rms" => NaN, "nlit" => 0,
                            "tmpl_id" => tid, "fit_error" => first(sprint(showerror, e), 300))
                    end
                    qa["expnum"] = fr.expnum
                    write_trace_qa(qp, qa)
                end
                note!(qa)
                push!(qas, qa)
            end
            ok = filter(q -> Bool(q["reg_gate_pass"]), qas)
            if !isempty(ok)
                # kind order, then fit_rms, then the later exposure (V2 §6.3)
                sort!(ok, by = q -> (something(findfirst(==(string(q["kind"])), t1_kinds), 99),
                    Float64(q["fit_rms"]), -Int(q["expnum"])))
                q = ok[1]
                return (level = level, reg = TemplateRegistration(Float64.(q["tmpl_reg"])), tmpl = tmpl,
                    source = "T1:$(q["kind"]):$(q["source"])", qa = q, info = info)
            end
        elseif level == TRACE_LEVEL_REGISTERED_SCI
            fs = night_science_frames(outdir, tele, mjd, chip; nmax = t2_nmax)
            isempty(fs) && continue
            qp = regsci_qa_path(qa_root, tele, mjd, chip)
            qa = _cached_qa(qp, tid, checkpoint_mode)
            if !isnothing(qa) && string(get(qa, "frames", "")) != join(basename.(fs), ";")
                qa = nothing
            end
            if isnothing(qa)
                qa = try
                    bs = [frame_bands(p, "object") for p in fs]
                    prof, err = stack_bands([b[2] for b in bs], [b[3] for b in bs])
                    rf = register_frame(getP(), tmpl, bs[1][1], prof, err; mode = :register_only)
                    registration_qa(rf, tmpl, "object_stack", "$(length(fs)) science frames"; gate_kw = gate_kw)
                catch e
                    Dict{String, Any}("kind" => "object_stack", "source" => "$(length(fs)) science frames",
                        "outcome" => "fit_error", "reg_gate_pass" => false, "reg_gate_reasons" => "fit_error",
                        "fit_rms" => NaN, "nlit" => 0, "tmpl_id" => tid,
                        "fit_error" => first(sprint(showerror, e), 300))
                end
                qa["frames"] = join(basename.(fs), ";")
                qa["expnum"] = -1
                write_trace_qa(qp, qa)
            end
            note!(qa)
            if Bool(qa["reg_gate_pass"])
                return (level = level, reg = TemplateRegistration(Float64.(qa["tmpl_reg"])), tmpl = tmpl,
                    source = "T2:object_stack:$(qa["frames"])", qa = qa, info = info)
            end
        elseif level == TRACE_LEVEL_TEMPLATE
            return (level = level, reg = epoch_median_registration(tmpl), tmpl = tmpl,
                source = "T3:template:$(tid)", qa = nothing, info = info)
        end
    end
    return nothing
end
