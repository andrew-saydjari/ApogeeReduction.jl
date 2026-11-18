# Per-epoch trace templates: calibration file, validation, lookup by MJD, identity.
#
# The robust-trace design (2026_11_17/trace_robust/TRACE_ROBUST_V2.md, Sec. 3 and 7) places every
# fiber of a night from a per-epoch TEMPLATE of the slit (fiber centers, widths and relative
# throughput along x), registered to the night by a few rigid numbers. The templates and the
# epoch table are learned from the archive's all-lit quartz flats by
# `scripts/cal/trace_templates/` and live in ONE versioned HDF5 file for both telescopes and all
# three chips. The pipeline is pointed at it by path (`--trace_template_cal`), exactly like the
# chip-transform file; without it (`none`) the pipeline runs as in PR A (template-free).
#
# This file only holds the file format, its validation, the lookup and the identity strings. The
# slit model, the template-seeded fits and the fallback levels that USE a template are elsewhere.
#
# Conventions (checked against pass-1 traceMain files):
#   fiber f = 1..300 is the second axis of production `trace_params` (slot 1 at the bottom of the
#   chip); centers are detector rows on the scale of `trace_params[:, :, 2]`; widths are the
#   Gauss-Hermite profile width of `trace_params[:, :, 3]`; x is the detector column 1..2048 and
#   polynomial coefficients are low -> high in u = (x - 1024) / 1024. The template sits in a fixed
#   gauge (each input flat's rigid registration removed); a night's centers are
#   T_f(x) + d(u) + s(u) * (T_f(x) - 1024) [+ per-fiber residual].
#
# File layout (format version 1):
#   attributes  format_version, calibration_version, content_fingerprint, complete (UInt8)
#   /epochs/<tele>/name (String[n]), mjd_start, mjd_end (Int64[n]), open_ended (UInt8[n])
#   /templates/<tele>_<chip>_<epoch>/
#       xc [NB], center [NB,300], width [NB,300], center_poly [5,300], width_poly [3,300],
#       height_rel [300], lit (UInt8 [300]), input_mjd, input_exp (Int64 [n]),
#       input_d_center, input_s_center (Float64 [n])
#       attributes tele, chip, epoch, mjd_start, mjd_end, fingerprint
#       /provenance  String attributes of the unit's build
#   /provenance  String attributes of the packing (inputs, code commit, missing units, ...)

using HDF5
using Statistics: median

const TRACE_TEMPLATE_FORMAT_VERSION = 1
"mjd_end stored for an open-ended last epoch."
const TRACE_EPOCH_OPEN_END = 99999

struct TraceEpoch
    tele::String
    name::String
    mjd_start::Int
    mjd_end::Int
    open_ended::Bool
end

struct TraceTemplate
    tele::String
    chip::String
    epoch::String
    mjd_start::Int
    mjd_end::Int
    xc::Vector{Float64}
    center::Matrix{Float64}
    width::Matrix{Float64}
    center_poly::Matrix{Float64}
    width_poly::Matrix{Float64}
    height_rel::Vector{Float64}
    lit::Vector{Bool}
    input_mjd::Vector{Int}
    input_exp::Vector{Int}
    input_d_center::Vector{Float64}
    input_s_center::Vector{Float64}
    reg_d_median::Float64
    reg_s_median::Float64
    data_mjd_min::Int
    data_mjd_max::Int
    fingerprint::String
    provenance::Dict{String, String}
end

struct TraceTemplateCal
    path::String
    calibration_version::String
    fingerprint::String
    complete::Bool
    epochs::Dict{String, Vector{TraceEpoch}}
    templates::Dict{Tuple{String, String, String}, TraceTemplate}
    provenance::Dict{String, String}
end

function Base.show(io::IO, cal::TraceTemplateCal)
    ne = join(["$(t): $(length(e)) epoch(s)" for (t, e) in sort(collect(cal.epochs))], ", ")
    print(io, "TraceTemplateCal($(cal.calibration_version); $(ne); $(length(cal.templates)) template(s)",
        cal.complete ? "" : ", INCOMPLETE", "; $(cal.path))")
end
Base.show(io::IO, t::TraceTemplate) =
    print(io, "TraceTemplate($(t.tele) $(t.chip) $(t.epoch) [$(t.mjd_start), $(t.mjd_end)]; ",
        "$(length(t.input_mjd)) flats; $(t.fingerprint))")

# ------------------------------------------------------------------ content fingerprint
# FNV-1a 64 over a canonical serialization of the CONTENT (arrays as little-endian Float64/Int64,
# NaN canonicalized, strings as UTF-8 with a separator). Unlike a hash of the file bytes it does
# not change when identical contents are re-written (HDF5 stores object times).
const _FNV_OFFSET = 0xcbf29ce484222325
const _FNV_PRIME = 0x00000100000001b3

_fnv(h::UInt64, b::UInt8) = (h ⊻ UInt64(b)) * _FNV_PRIME
function _fnv(h::UInt64, bytes::AbstractVector{UInt8})
    for b in bytes
        h = _fnv(h, b)
    end
    return h
end
_fnv(h::UInt64, s::AbstractString) = _fnv(_fnv(h, Vector{UInt8}(codeunits(String(s)))), 0x1f)
_fnv(h::UInt64, x::Integer) = _fnv(h, reinterpret(UInt8, [htol(Int64(x))]))
function _fnv(h::UInt64, a::AbstractArray{<:AbstractFloat})
    v = [isnan(x) ? NaN : Float64(x) for x in vec(a)]
    h = _fnv(h, Int64(length(v)))
    return _fnv(h, reinterpret(UInt8, htol.(v)))
end
function _fnv(h::UInt64, a::AbstractArray{<:Union{Integer, Bool}})
    h = _fnv(h, Int64(length(a)))
    return _fnv(h, reinterpret(UInt8, htol.(Int64.(vec(a)))))
end
_hex(h::UInt64) = string(h, base = 16, pad = 16)

"Content fingerprint of one template unit (its arrays and the epoch row it serves)."
function template_fingerprint(tele, chip, epoch, mjd_start, mjd_end, xc, center, width, center_poly,
        width_poly, height_rel, lit, input_mjd, input_exp, input_d_center, input_s_center)
    h = _FNV_OFFSET
    for s in (tele, chip, epoch)
        h = _fnv(h, s)
    end
    h = _fnv(h, mjd_start); h = _fnv(h, mjd_end)
    for a in (xc, center, width, center_poly, width_poly, height_rel)
        h = _fnv(h, a)
    end
    h = _fnv(h, lit); h = _fnv(h, input_mjd); h = _fnv(h, input_exp)
    h = _fnv(h, input_d_center); h = _fnv(h, input_s_center)
    return _hex(h)
end
template_fingerprint(t::TraceTemplate) =
    template_fingerprint(t.tele, t.chip, t.epoch, t.mjd_start, t.mjd_end, t.xc, t.center, t.width,
        t.center_poly, t.width_poly, t.height_rel, t.lit, t.input_mjd, t.input_exp, t.input_d_center,
        t.input_s_center)

"Content fingerprint of a whole calibration: the epoch tables and the unit fingerprints."
function cal_content_fingerprint(epochs::Dict{String, Vector{TraceEpoch}}, templates)
    h = _FNV_OFFSET
    for tele in sort(collect(keys(epochs)))
        h = _fnv(h, tele)
        for e in epochs[tele]
            h = _fnv(h, e.name); h = _fnv(h, e.mjd_start); h = _fnv(h, e.mjd_end)
            h = _fnv(h, Int(e.open_ended))
        end
    end
    for t in sort(collect(templates), by = t -> (t.tele, t.chip, t.epoch))
        h = _fnv(h, t.fingerprint)
    end
    return _hex(h)
end

"""
    TraceTemplate(tele, chip, epoch::TraceEpoch, xc, center, width, center_poly, width_poly,
        height_rel, lit, input_mjd, input_exp, input_d_center, input_s_center, provenance)

Build a template unit; derives the epoch-median registration, the data range and the fingerprint.
"""
function TraceTemplate(tele, chip, epoch::TraceEpoch, xc, center, width, center_poly, width_poly,
        height_rel, lit, input_mjd, input_exp, input_d_center, input_s_center,
        provenance::Dict{String, String})
    fin(v) = filter(isfinite, v)
    rd = isempty(fin(input_d_center)) ? NaN : median(fin(input_d_center))
    rs = isempty(fin(input_s_center)) ? NaN : median(fin(input_s_center))
    im = Int.(collect(input_mjd))
    fp = template_fingerprint(tele, chip, epoch.name, epoch.mjd_start, epoch.mjd_end, xc, center,
        width, center_poly, width_poly, height_rel, Bool.(lit), im, Int.(input_exp),
        Float64.(input_d_center), Float64.(input_s_center))
    return TraceTemplate(String(tele), String(chip), epoch.name, epoch.mjd_start, epoch.mjd_end,
        Float64.(collect(xc)), Float64.(center), Float64.(width), Float64.(center_poly),
        Float64.(width_poly), Float64.(collect(height_rel)), Bool.(collect(lit)), im,
        Int.(collect(input_exp)), Float64.(collect(input_d_center)), Float64.(collect(input_s_center)),
        rd, rs, isempty(im) ? -1 : minimum(im), isempty(im) ? -1 : maximum(im), fp, provenance)
end

# ------------------------------------------------------------------ validation
"""
    check_trace_epochs(tele, epochs)

Epochs of one telescope must be sorted, non-overlapping and CONTIGUOUS (every MJD between the
first start and the last end is served by exactly one epoch), with unique names; only the last
may be open-ended.
"""
function check_trace_epochs(tele, ep::Vector{TraceEpoch})
    isempty(ep) && throw(ArgumentError("telescope $(tele): empty epoch table"))
    length(unique(e.name for e in ep)) == length(ep) ||
        throw(ArgumentError("telescope $(tele): duplicate epoch names"))
    for (i, e) in enumerate(ep)
        e.tele == tele || throw(ArgumentError("epoch $(e.name) belongs to $(e.tele), listed under $(tele)"))
        e.mjd_start <= e.mjd_end ||
            throw(ArgumentError("telescope $(tele) epoch $(e.name): mjd_start > mjd_end"))
        if i > 1
            p = ep[i - 1]
            e.mjd_start == p.mjd_end + 1 ||
                throw(ArgumentError("telescope $(tele): epochs $(p.name) [$(p.mjd_start), $(p.mjd_end)] and " *
                                    "$(e.name) [$(e.mjd_start), $(e.mjd_end)] are unsorted, overlap or leave a gap"))
        end
        e.open_ended && i != length(ep) &&
            throw(ArgumentError("telescope $(tele): only the last epoch may be open-ended ($(e.name))"))
    end
    return nothing
end

"Shape and finiteness checks of one template unit against the epoch it claims to serve."
function check_trace_template(t::TraceTemplate, e::TraceEpoch)
    id = "$(t.tele) $(t.chip) $(t.epoch)"
    t.chip in CHIP_LIST || throw(ArgumentError("$(id): unknown chip"))
    (t.mjd_start, t.mjd_end) == (e.mjd_start, e.mjd_end) ||
        throw(ArgumentError("$(id): template range [$(t.mjd_start), $(t.mjd_end)] differs from its epoch [$(e.mjd_start), $(e.mjd_end)]"))
    nb = length(t.xc)
    (nb >= 4 && issorted(t.xc) && all(isfinite, t.xc)) || throw(ArgumentError("$(id): bad xc"))
    size(t.center) == (nb, N_FIBERS) || throw(ArgumentError("$(id): center must be $(nb) x $(N_FIBERS)"))
    size(t.width) == (nb, N_FIBERS) || throw(ArgumentError("$(id): width must be $(nb) x $(N_FIBERS)"))
    size(t.center_poly, 2) == N_FIBERS || throw(ArgumentError("$(id): center_poly must have $(N_FIBERS) columns"))
    size(t.width_poly, 2) == N_FIBERS || throw(ArgumentError("$(id): width_poly must have $(N_FIBERS) columns"))
    length(t.height_rel) == N_FIBERS && length(t.lit) == N_FIBERS ||
        throw(ArgumentError("$(id): height_rel and lit need $(N_FIBERS) entries"))
    all(isfinite, t.center) && all(isfinite, t.width) && all(>(0), t.width) ||
        throw(ArgumentError("$(id): non-finite centers or non-positive widths"))
    all(isfinite, t.center_poly) && all(isfinite, t.width_poly) ||
        throw(ArgumentError("$(id): non-finite polynomial coefficients"))
    # fibers keep their order along the slit in every band
    all(b -> issorted(view(t.center, b, :)), 1:nb) ||
        throw(ArgumentError("$(id): template centers are not increasing with fiber index"))
    length(t.input_mjd) == length(t.input_exp) == length(t.input_d_center) == length(t.input_s_center) ||
        throw(ArgumentError("$(id): input lists of unequal length"))
    all(m -> t.mjd_start <= m <= t.mjd_end, t.input_mjd) ||
        throw(ArgumentError("$(id): an input flat lies outside the epoch"))
    return nothing
end

"Missing (tele, chip, epoch) units of an epoch table, sorted."
function missing_template_units(epochs::Dict{String, Vector{TraceEpoch}}, have)
    miss = Tuple{String, String, String}[]
    for tele in sort(collect(keys(epochs))), e in epochs[tele], chip in CHIP_LIST
        (tele, chip, e.name) in have || push!(miss, (tele, chip, e.name))
    end
    return miss
end

# ------------------------------------------------------------------ write / read
"""
    write_trace_template_cal(path, epochs, templates, provenance; calibration_version,
        allow_incomplete = false) -> path

Write the calibration file (format above) after validating it. With `allow_incomplete = false`
every (telescope, epoch, chip) must have a template; with `true` the missing units are listed in
the provenance (`missing_units`) and the file is marked incomplete. Written to a temporary name
and moved into place.
"""
function write_trace_template_cal(path::AbstractString, epochs::Dict{String, Vector{TraceEpoch}},
        templates::AbstractVector{TraceTemplate}, provenance::Dict{String, String};
        calibration_version::AbstractString, allow_incomplete::Bool = false)
    for (tele, ep) in epochs
        check_trace_epochs(tele, ep)
    end
    have = Set((t.tele, t.chip, t.epoch) for t in templates)
    length(have) == length(templates) || throw(ArgumentError("duplicate template units"))
    for t in templates
        ep = get(epochs, t.tele, TraceEpoch[])
        i = findfirst(e -> e.name == t.epoch, ep)
        isnothing(i) && throw(ArgumentError("template $(t.tele) $(t.chip) $(t.epoch): epoch not in the table"))
        check_trace_template(t, ep[i])
        template_fingerprint(t) == t.fingerprint ||
            throw(ArgumentError("template $(t.tele) $(t.chip) $(t.epoch): stale fingerprint"))
    end
    miss = missing_template_units(epochs, have)
    !isempty(miss) && !allow_incomplete &&
        throw(ArgumentError("$(length(miss)) template unit(s) missing: " *
                            join(["$(a)_$(b)_$(c)" for (a, b, c) in miss], ", ")))
    fp = cal_content_fingerprint(epochs, templates)
    tmp = path * ".tmp$(getpid())"
    h5open(tmp, "w") do f
        a = HDF5.attributes(f)
        a["format_version"] = TRACE_TEMPLATE_FORMAT_VERSION
        a["calibration_version"] = String(calibration_version)
        a["content_fingerprint"] = fp
        a["complete"] = UInt8(isempty(miss))
        ge = create_group(f, "epochs")
        for (tele, ep) in epochs
            g = create_group(ge, tele)
            g["name"] = [e.name for e in ep]
            g["mjd_start"] = Int64[e.mjd_start for e in ep]
            g["mjd_end"] = Int64[e.mjd_end for e in ep]
            g["open_ended"] = UInt8[e.open_ended for e in ep]
        end
        gt = create_group(f, "templates")
        for t in templates
            g = create_group(gt, "$(t.tele)_$(t.chip)_$(t.epoch)")
            g["xc"] = t.xc; g["center"] = t.center; g["width"] = t.width
            g["center_poly"] = t.center_poly; g["width_poly"] = t.width_poly
            g["height_rel"] = t.height_rel; g["lit"] = UInt8.(t.lit)
            g["input_mjd"] = Int64.(t.input_mjd); g["input_exp"] = Int64.(t.input_exp)
            g["input_d_center"] = t.input_d_center; g["input_s_center"] = t.input_s_center
            ga = HDF5.attributes(g)
            ga["tele"] = t.tele; ga["chip"] = t.chip; ga["epoch"] = t.epoch
            ga["mjd_start"] = t.mjd_start; ga["mjd_end"] = t.mjd_end; ga["fingerprint"] = t.fingerprint
            gp = create_group(g, "provenance")
            for (k, v) in t.provenance
                HDF5.attributes(gp)[k] = String(v)
            end
        end
        p = create_group(f, "provenance")
        prov = copy(provenance)
        prov["missing_units"] = isempty(miss) ? "none" : join(["$(a)_$(b)_$(c)" for (a, b, c) in miss], ";")
        for (k, v) in prov
            HDF5.attributes(p)[k] = String(v)
        end
    end
    mv(tmp, path; force = true)
    return path
end

_attrs_dict(obj) = Dict{String, String}(k => string(read(HDF5.attributes(obj)[k])) for k in keys(HDF5.attributes(obj)))

"""
    load_trace_template_cal(path; require_complete = false) -> TraceTemplateCal

Read and validate a calibration file written by `write_trace_template_cal`. Throws if the file
does not exist, has an unknown format version, fails validation, or its stored content
fingerprint does not match the recomputed one (a configured but unreadable calibration input is
a hard error, never a silent fallback to the template-free path). With `require_complete = true`
an incomplete file is also an error.
"""
function load_trace_template_cal(path::AbstractString; require_complete::Bool = false)
    isfile(path) || throw(ArgumentError("trace-template calibration file not found: $(path)"))
    h5open(path, "r") do f
        a = HDF5.attributes(f)
        haskey(a, "format_version") ||
            throw(ArgumentError("$(path): not a trace-template calibration file (no format_version)"))
        fv = read(a["format_version"])
        fv == TRACE_TEMPLATE_FORMAT_VERSION ||
            throw(ArgumentError("$(path): trace-template format version $(fv), expected $(TRACE_TEMPLATE_FORMAT_VERSION)"))
        version = string(read(a["calibration_version"]))
        stored_fp = string(read(a["content_fingerprint"]))
        epochs = Dict{String, Vector{TraceEpoch}}()
        for tele in keys(f["epochs"])
            g = f["epochs"][tele]
            nm = read(g["name"]); ms = read(g["mjd_start"]); me = read(g["mjd_end"]); oe = read(g["open_ended"])
            epochs[tele] = [TraceEpoch(tele, nm[i], ms[i], me[i], oe[i] != 0) for i in eachindex(nm)]
            check_trace_epochs(tele, epochs[tele])
        end
        templates = Dict{Tuple{String, String, String}, TraceTemplate}()
        for key in keys(f["templates"])
            g = f["templates"][key]
            ga = HDF5.attributes(g)
            tele = read(ga["tele"]); chip = read(ga["chip"]); ename = read(ga["epoch"])
            ep = get(epochs, tele, TraceEpoch[])
            i = findfirst(e -> e.name == ename, ep)
            isnothing(i) && throw(ArgumentError("$(path): template $(key) has no epoch in the table"))
            t = TraceTemplate(tele, chip, ep[i], read(g["xc"]), read(g["center"]), read(g["width"]),
                read(g["center_poly"]), read(g["width_poly"]), read(g["height_rel"]),
                read(g["lit"]) .!= 0, read(g["input_mjd"]), read(g["input_exp"]),
                read(g["input_d_center"]), read(g["input_s_center"]),
                haskey(g, "provenance") ? _attrs_dict(g["provenance"]) : Dict{String, String}())
            check_trace_template(t, ep[i])
            t.fingerprint == string(read(ga["fingerprint"])) ||
                throw(ArgumentError("$(path): template $(key) content does not match its stored fingerprint"))
            templates[(tele, chip, ename)] = t
        end
        fp = cal_content_fingerprint(epochs, values(templates))
        fp == stored_fp ||
            throw(ArgumentError("$(path): content fingerprint $(fp) does not match the stored $(stored_fp)"))
        miss = missing_template_units(epochs, keys(templates))
        complete = isempty(miss)
        complete == (read(a["complete"]) != 0) ||
            throw(ArgumentError("$(path): 'complete' attribute disagrees with the templates present"))
        require_complete && !complete &&
            throw(ArgumentError("$(path): incomplete trace-template calibration, missing " *
                                join(["$(x)_$(y)_$(z)" for (x, y, z) in miss], ", ")))
        prov = haskey(f, "provenance") ? _attrs_dict(f["provenance"]) : Dict{String, String}()
        return TraceTemplateCal(String(path), version, fp, complete, epochs, templates, prov)
    end
end

# ------------------------------------------------------------------ lookup and identity
"""
    trace_epoch(cal, tele, mjd) -> TraceEpoch or nothing

The epoch whose inclusive [mjd_start, mjd_end] contains `mjd`.
"""
function trace_epoch(cal::TraceTemplateCal, tele::AbstractString, mjd::Integer)
    for e in get(cal.epochs, String(tele), TraceEpoch[])
        e.mjd_start <= mjd <= e.mjd_end && return e
    end
    return nothing
end
trace_epoch(::Nothing, tele, mjd) = nothing

"""
    trace_template(cal, tele, chip, mjd) -> TraceTemplate or nothing

The template of the epoch serving `mjd`; `nothing` without a calibration, when no epoch serves the
night, or when that (telescope, chip, epoch) unit is missing from an incomplete file.
"""
function trace_template(cal::TraceTemplateCal, tele::AbstractString, chip::AbstractString, mjd::Integer)
    e = trace_epoch(cal, tele, mjd)
    isnothing(e) && return nothing
    return get(cal.templates, (String(tele), String(chip), e.name), nothing)
end
trace_template(::Nothing, tele, chip, mjd) = nothing

"""
    trace_template_provisional(cal, tele, chip, mjd) -> Bool

True when the night lies in an open-ended (last) epoch AFTER the newest flat that built its
template: an instrument change since then would not be known to the epoch table yet
(`TRACE_EPOCH_PROVISIONAL`).
"""
function trace_template_provisional(cal, tele::AbstractString, chip::AbstractString, mjd::Integer)
    t = trace_template(cal, tele, chip, mjd)
    isnothing(t) && return false
    e = trace_epoch(cal, tele, mjd)
    return e.open_ended && mjd > t.data_mjd_max
end

"""
    trace_template_id(cal, tele, chip, mjd) -> String

Identity of the template input of one (telescope, chip, night), recorded as `trace_template_id`
in the traceMain and per-flat QA metadata and compared under checkpointing: "none" without a
calibration or a serving template, else "<tele>_<chip>_<epoch>:<fingerprint>". Does not contain
the calibration version, so re-packing identical templates does not force recomputation.
"""
function trace_template_id(cal, tele::AbstractString, chip::AbstractString, mjd::Integer)
    t = trace_template(cal, tele, chip, mjd)
    isnothing(t) && return "none"
    return "$(t.tele)_$(t.chip)_$(t.epoch):$(t.fingerprint)"
end

"File-level identity: \"none\" or \"<calibration_version>:<content fingerprint>\"."
trace_template_cal_id(::Nothing) = "none"
trace_template_cal_id(cal::TraceTemplateCal) = "$(cal.calibration_version):$(cal.fingerprint)"

"""
    trace_template_metadata(cal, tele, chip, mjd) -> Dict{String, Any}

The provenance block written into traceMain files and per-flat QA records.
"""
function trace_template_metadata(cal, tele::AbstractString, chip::AbstractString, mjd::Integer)
    t = trace_template(cal, tele, chip, mjd)
    return Dict{String, Any}(
        "trace_template_id" => trace_template_id(cal, tele, chip, mjd),
        "trace_template_cal_id" => trace_template_cal_id(cal),
        "trace_template_cal" => isnothing(cal) ? "none" : cal.path,
        "trace_template_epoch" => isnothing(t) ? "" : t.epoch,
        "trace_template_provisional" => trace_template_provisional(cal, tele, chip, mjd))
end

"""
    template_needs_recompute(stored_metadata, cal, tele, chip, mjd) -> Bool

True when a checkpointed product's recorded `trace_template_id` (absent = "none", i.e. made
without templates) differs from the one the current configuration gives.
"""
template_needs_recompute(md::AbstractDict, cal, tele, chip, mjd) =
    string(get(md, "trace_template_id", "none")) != trace_template_id(cal, tele, chip, mjd)

_upoly(c::AbstractMatrix, xs) = (u = (Float64.(xs) .- 1024) ./ 1024;
    [sum(c[p, f] * ui^(p - 1) for p in axes(c, 1)) for ui in u, f in axes(c, 2)])

"Template centers at arbitrary columns `xs`: [length(xs), 300] (px, template frame)."
template_centers(t::TraceTemplate, xs::AbstractVector) = _upoly(t.center_poly, xs)
"Template Gauss-Hermite widths at arbitrary columns `xs`: [length(xs), 300] (px)."
template_widths(t::TraceTemplate, xs::AbstractVector) = _upoly(t.width_poly, xs)

# ------------------------------------------------------------------ builder output -> unit
"""
    read_template_unit_file(path; epoch::TraceEpoch, log_csv = nothing) -> TraceTemplate

Read one output of the template builder (`template_<tele>_<chip>_<epoch>.h5`, datasets xc,
center, width, center_poly_u, width_poly_u, height_rel, lit_in_epoch, input_mjd, input_exp). The
per-input registration (d, s at x = 1024) comes from the builder's log CSV (rows with
stage == "final"), if given; otherwise NaN. The builder's own [mjd_start, mjd_end] must lie inside
`epoch` (an open-ended epoch extends the served range past the builder's).
"""
function read_template_unit_file(path::AbstractString; epoch::TraceEpoch, log_csv = nothing)
    h5open(path, "r") do f
        prov = _attrs_dict(f)
        tele = prov["tele"]; chip = prov["chip"]
        tele == epoch.tele || throw(ArgumentError("$(path): telescope $(tele) != epoch's $(epoch.tele)"))
        prov["epoch"] == epoch.name || throw(ArgumentError("$(path): epoch $(prov["epoch"]) != $(epoch.name)"))
        m0 = parse(Int, prov["mjd_start"]); m1 = parse(Int, prov["mjd_end"])
        (epoch.mjd_start == m0 && (m1 == epoch.mjd_end || epoch.open_ended)) ||
            throw(ArgumentError("$(path): built for [$(m0), $(m1)], epoch table says [$(epoch.mjd_start), $(epoch.mjd_end)]"))
        im = Int.(read(f["input_mjd"])); ie = Int.(read(f["input_exp"]))
        dc = fill(NaN, length(im)); sc = fill(NaN, length(im))
        if !isnothing(log_csv) && isfile(log_csv)
            reg = _read_final_registration(log_csv)
            for i in eachindex(im)
                v = get(reg, (im[i], ie[i]), nothing)
                isnothing(v) || ((dc[i], sc[i]) = v)
            end
        end
        prov["source_file"] = String(path)
        prov["builder_mjd_range"] = "$(m0)-$(m1)"
        return TraceTemplate(tele, chip, epoch, read(f["xc"]), read(f["center"]), read(f["width"]),
            read(f["center_poly_u"]), read(f["width_poly_u"]), read(f["height_rel"]),
            read(f["lit_in_epoch"]) .!= 0, im, ie, dc, sc, prov)
    end
end

# minimal CSV reader for the builder log (no CSV.jl dependency in the hot module path)
function _read_final_registration(log_csv)
    lines = readlines(log_csv)
    hdr = split(lines[1], ',')
    ix(n) = findfirst(==(n), hdr)
    ist, im, ie, id, is = ix("stage"), ix("mjd"), ix("exp"), ix("d_center"), ix("s_center")
    out = Dict{Tuple{Int, Int}, Tuple{Float64, Float64}}()
    for l in lines[2:end]
        c = split(l, ',')
        length(c) >= length(hdr) || continue
        c[ist] == "final" || continue
        out[(parse(Int, c[im]), parse(Int, c[ie]))] = (parse(Float64, c[id]), parse(Float64, c[is]))
    end
    return out
end
