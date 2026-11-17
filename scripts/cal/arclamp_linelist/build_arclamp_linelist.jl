# Build the frozen arc-lamp line list `data/arclamp_linelist_v<N>.csv` used by the stage-1 arc
# wavelength fit (src/arclamp_wavecal.jl).
#
# usage:
#   julia --project=. scripts/cal/arclamp_linelist/build_arclamp_linelist.jl \
#       --empirical <empirical_linelist.csv> --labfield <labfield.csv> --lineset <lineset.csv> \
#       --fitsets <fitsets.csv> --status <line_status_pooled.csv> \
#       --out data/arclamp_linelist_v1.csv --version v1
#
# Inputs (see README.md in this folder for where each comes from):
#   empirical : recurring arc features found by the local-noise finder over the archive (one row
#               per lamp x chip x feature; columns lamp, chip_idx, wave, class, med_pix). Class "X"
#               rows (rejected features) are dropped; all others are kept, as association guards
#               when they have no laboratory match.
#   labfield  : laboratory lines per lamp (lamp, wave [vacuum A], species, spgroup), from the clean
#               tables (NIST SRD 161 / Redman 2014 Th, Redman 2011 U/Ar, NIST ASD Ne).
#   lineset   : the laboratory lines considered for fitting, with ids (lamp, wave, species, class,
#               drp_rederived, id).
#   fitsets   : per line id, the data-derived quality sets (clean, shape, clean_blend, shape_blend)
#               from the line-shape / stability classification.
#   status    : per line id, the pooled classification status.
#
# The output is deterministic given the inputs. Its header records the version, the date, the
# git commit of this repository and the sha256 of every input. The pipeline records the output's
# FNV-1a fingerprint in every product made with it.
#
# Selection columns (Bool) written; the fit uses one of them (`load_arclamp_linelist(; selection)`):
#   use_default      shape-clean lab Th/U/Ar lines, minus the DRP re-derived (blended) lines, minus Ne
#   use_with_blended use_default plus the DRP re-derived (blended) features, at lab wavelength
#   use_with_neon    use_default plus every lab-matched Ne feature
#   prelim_classB    the September lab class-B Th/U/Ar list
#   strict_clean     shape- and stability-clean with no significant offset

using CSV, DataFrames, Dates, Statistics

function parse_args(args)
    d = Dict{String, String}()
    i = 1
    while i <= length(args)
        startswith(args[i], "--") || error("unexpected argument $(args[i])")
        d[args[i][3:end]] = args[i + 1]
        i += 2
    end
    for k in ("empirical", "labfield", "lineset", "fitsets", "status", "out", "version")
        haskey(d, k) || error("missing --$(k)")
    end
    return d
end

sha256_of(path) = split(read(`sha256sum $(path)`, String))[1]

"""Match each empirical feature to the nearest laboratory line of its lamp, after removing a
per-(lamp, chip) median offset estimated from matches within 0.1 A. Returns (lab wave, lab
species group, separation) per feature (NaN / "" when no lab line lies within `tol`)."""
function match_lab(cat::DataFrame, labfield::DataFrame; tol = 0.06)
    labw = fill(NaN, nrow(cat))
    labsp = fill("", nrow(cat))
    sep = fill(NaN, nrow(cat))
    raw = fill(NaN, nrow(cat))
    for i in 1:nrow(cat)
        f = labfield[labfield.lamp .== cat.lamp[i], :]
        j = argmin(abs.(f.wave .- cat.wave[i]))
        raw[i] = cat.wave[i] - f.wave[j]
    end
    zp = Dict{Tuple{String, Int}, Float64}()
    for l in ("ThAr", "UNe"), c in 1:3
        m = (cat.lamp .== l) .& (cat.chip_idx .== c) .& (abs.(raw) .< 0.1) .& (cat.class .!= "X")
        zp[(l, c)] = count(m) >= 3 ? median(raw[m]) : 0.0
    end
    for i in 1:nrow(cat)
        f = labfield[labfield.lamp .== cat.lamp[i], :]
        w0 = cat.wave[i] - zp[(cat.lamp[i], cat.chip_idx[i])]
        j = argmin(abs.(f.wave .- w0))
        abs(f.wave[j] - w0) < tol || continue
        labw[i] = f.wave[j]
        labsp[i] = f.spgroup[j]
        sep[i] = f.wave[j] - w0
    end
    return labw, labsp, sep
end

function build(d)
    emp = CSV.read(d["empirical"], DataFrame)
    labfield = CSV.read(d["labfield"], DataFrame)
    lineset = CSV.read(d["lineset"], DataFrame)
    fitsets = CSV.read(d["fitsets"], DataFrame)
    status = CSV.read(d["status"], DataFrame)

    df = DataFrame(lamp = String.(emp.lamp), chip_idx = Int.(emp.chip_idx),
        emp_wave = Float64.(emp.wave), emp_class = String.(emp.class),
        med_pix = Float64.(emp.med_pix))
    df = df[df.emp_class .!= "X", :]
    labw, labsp, sep = match_lab(rename(copy(df), :emp_wave => :wave, :emp_class => :class), labfield)
    df.lab_id = zeros(Int, nrow(df))
    for i in 1:nrow(df)
        isfinite(labw[i]) || continue
        k = findfirst((lineset.lamp .== df.lamp[i]) .& (abs.(lineset.wave .- labw[i]) .< 0.002))
        k === nothing || (df.lab_id[i] = lineset.id[k])
    end
    df.lab_wave = [df.lab_id[i] > 0 ? labw[i] : NaN for i in 1:nrow(df)]
    # The DRP re-derived (blended) lines sit 40-210 mA from their lab wavelength, beyond the
    # 0.06 A lab-match tolerance, so the match above leaves their features unlabelled. Label the
    # nearest unlabelled feature of the same lamp within `drp_tol` A, so that they can be switched
    # on (`use_with_blended`) by a one-line selection change. They are never in `use_default`.
    drp_tol = 0.25
    for r in eachrow(lineset[Bool.(lineset.drp_rederived), :])
        any(==(r.id), df.lab_id) && continue
        cand = findall((df.lamp .== r.lamp) .& (df.lab_id .== 0) .& (abs.(df.emp_wave .- r.wave) .< drp_tol))
        isempty(cand) && continue
        i = cand[argmin(abs.(df.emp_wave[cand] .- r.wave))]
        df.lab_id[i] = r.id
        df.lab_wave[i] = r.wave
    end
    byid(col, default) = [df.lab_id[i] > 0 ? begin
        r = findfirst(==(df.lab_id[i]), col === :status ? status.id : lineset.id)
        r === nothing ? default : (col === :status ? status.status[r] : lineset[r, col])
    end : default for i in 1:nrow(df)]
    df.lab_species = String.(string.(byid(:species, "")))
    df.lab_class = String.(string.(byid(:class, "")))
    df.lab_status = String.(string.(byid(:status, "")))
    fs = Dict(r.id => r for r in eachrow(fitsets))
    getfs(i, k) = df.lab_id[i] > 0 && haskey(fs, df.lab_id[i]) ? Bool(fs[df.lab_id[i]][k]) : false
    df.store = df.lab_id .> 0
    df.shape_clean = [getfs(i, :shape) for i in 1:nrow(df)]
    df.shape_blend = [getfs(i, :shape_blend) for i in 1:nrow(df)]
    df.strict_clean = [getfs(i, :clean) for i in 1:nrow(df)]
    df.drp_rederived = Bool[df.lab_id[i] > 0 &&
                            Bool(lineset.drp_rederived[findfirst(==(df.lab_id[i]), lineset.id)])
                            for i in 1:nrow(df)]
    df.neon = [df.lab_id[i] > 0 && occursin("Ne", df.lab_species[i]) for i in 1:nrow(df)]
    df.prelim_classB = [df.lab_id[i] > 0 && df.lab_class[i] == "B" && !df.neon[i] &&
                        df.lab_species[i] != "U?" for i in 1:nrow(df)]
    df.use_default = df.shape_clean .& .!df.drp_rederived .& .!df.neon
    df.use_with_blended = df.use_default .| (df.drp_rederived .& df.store)
    df.use_with_neon = df.use_default .| (df.neon .& df.store)
    # association wavelength: the lab wavelength for the candidate fit lines (prelim class B,
    # shape-clean, shape+modelled blends), the empirical wavelength for everything else
    # (guards, Ne, the DRP re-derived features)
    fitcand = (df.prelim_classB .| df.shape_clean .| df.shape_blend) .& .!df.drp_rederived
    df.assoc_wave = [fitcand[i] ? df.lab_wave[i] : df.emp_wave[i] for i in 1:nrow(df)]
    df.row_id = 1:nrow(df)
    dup = [id for id in unique(df.lab_id[df.lab_id .> 0]) if count(==(id), df.lab_id) > 1]
    isempty(dup) || @warn "lab ids matched by more than one feature: $(dup)"
    cols = [:row_id, :lamp, :chip_idx, :assoc_wave, :med_pix, :emp_wave, :emp_class, :lab_id,
        :lab_wave, :lab_species, :lab_class, :lab_status, :store, :shape_clean, :shape_blend,
        :strict_clean, :prelim_classB, :drp_rederived, :neon, :use_default, :use_with_blended,
        :use_with_neon]
    df = df[:, cols]

    commit = try
        strip(read(`git -C $(@__DIR__) rev-parse HEAD`, String))
    catch
        "unknown"
    end
    hdr = String[]
    push!(hdr, "# APOGEE arc-lamp line list $(d["version"]) (ApogeeReduction.jl, stage-1 arc wavelength calibration)")
    push!(hdr, "# built $(Dates.format(now(UTC), "yyyy-mm-ddTHH:MM:SS"))Z by scripts/cal/arclamp_linelist/build_arclamp_linelist.jl at git commit $(commit)")
    push!(hdr, "# wavelengths: vacuum Angstrom. lab_wave from the clean laboratory tables (never the apogee_drp copies).")
    for k in ("empirical", "labfield", "lineset", "fitsets", "status")
        push!(hdr, "# input $(k): $(basename(d[k])) sha256 $(sha256_of(d[k]))")
    end
    push!(hdr, "# rows $(nrow(df)); lab-matched $(count(df.store)); use_default $(count(df.use_default)); " *
               "use_with_blended $(count(df.use_with_blended)); use_with_neon $(count(df.use_with_neon)); " *
               "prelim_classB $(count(df.prelim_classB)); strict_clean $(count(df.strict_clean))")
    tmp = d["out"] * ".tmp"
    open(tmp, "w") do io
        for h in hdr
            println(io, h)
        end
        CSV.write(io, df; append = false, writeheader = true)
    end
    mv(tmp, d["out"]; force = true)
    println(join(hdr, "\n"))
    return df
end

if abspath(PROGRAM_FILE) == @__FILE__
    build(parse_args(ARGS))
end
