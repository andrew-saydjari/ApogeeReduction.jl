# make_sky_linelist.jl -- generate data/APOGEE_sky_linelist.csv, the sky (airglow)
# line list consumed by the wavelength calibration (src/skyline_peaks.jl).
#
# The list is GENERATED from cached literature source data (never hand-edited):
#
#   PRIMARY SOURCE (positions, Einstein A, labels, verified Lambda-doublet splitting):
#     Brooke et al. 2016, JQSRT 168, 142 (arXiv:1503.08420) -- the MoLLIST OH line list,
#     via the ExoMol conversion "MoLLIST-OH" (states + trans files, compiled by Yixin Wang;
#     Tennyson & Yurchenko 2016 format), retrieved 2026-09-12 from
#     https://www.exomol.com/db/OH/16O-1H/MoLLIST-OH/16O-1H__MoLLIST-OH.states.bz2
#     https://www.exomol.com/db/OH/16O-1H/MoLLIST-OH/16O-1H__MoLLIST-OH.trans.bz2
#     Positions rest on the Abrams et al. 1994 (ApJS 93, 351) FTS measurements; H-band
#     line-position accuracy ~0.0035 cm^-1 ~ 0.01 A (Oliva et al. 2013, A&A 555, A78).
#
#   CROSS-CHECKS (never sources of the adopted wavelengths):
#     Oliva et al. 2015, A&A 581, A47 -- GIANO empirical sky spectrum, per-doublet photon
#       fluxes + component positions (VizieR J/A+A/581/A47 table1.dat, retrieved 2026-09-12
#       from https://cdsarc.cds.unistra.fr/ftp/J/A+A/581/A47/table1.dat). Also supplies the
#       per-band empirical intensity scale for the model amplitude column.
#     Rousselot et al. 2000, A&A 354, 1134 -- centroid cross-check only (vacuum-converted
#       copy rousselot.vac.txt from the SDSS apogee_drp repo, branch daily, retrieved
#       2026-09-12; that repo's README_airglow.txt documents the air->vacuum conversion).
#       NEVER used for Q-branch component positions (its Q-branch splittings are wrong).
#     SDSS apogee_drp data/skylines/airglow.txt (branch daily, retrieved 2026-09-12) --
#       cross-check of feature reality and of the DRP's own USEWAVE choices; never a source.
#     stack_sky_peaks.csv -- OUR measured stacked-sky peak heights (median running-median-
#       scaled flux, all fibers x 4 object exposures, healthy night 60105, both telescopes;
#       pass-1 DR21 products). Copied here for reproducibility.
#
# All cached raw files live in data/linelist_sources/ with sha256 checksums in
# CHECKSUMS.sha256; this script verifies the checksums before using them.
#
# Wavelengths are VACUUM, in Angstroms, throughout.
#
# Selection cuts (candidate pool; encoded here, stated in the output header):
#   - X2Pi-X2Pi (Meinel) MAIN-branch (F' == F'') Lambda-doublets with exactly two
#     e/f components, both inside the APOGEE chip ranges [15125, 16960] A;
#   - doublet separation <= 2.5 A (the peak fitter uses a 21-pixel window);
#   - cold-atmosphere model intensity >= 0.002 of the brightest pool line
#     (Einstein-A based, T_rot = 200 K; see below);
#   - isolation: no other OH Meinel line (any branch) with component model intensity
#     >= 15% of the candidate's weaker component within 0.8 A of either component.
#
# Model amplitude column (model_rel_int): per component I = A * g_up * exp(-(E_up -
# E0_v')/kT) at T_rot = 200 K (component ratio of a Lambda-doublet is 1:1 to <= 0.2%:
# Einstein-A e/f equality <= 0.03% from these very A values; LTE sublevel-population
# equality <= 0.2%, T-independent to <= 1e-3), summed over the doublet, with each
# vibrational band scaled to the Oliva et al. 2015 empirical per-doublet fluxes
# (median flux/model ratio over matched LOW-J pool lines of that band), then normalized to
# the brightest pool line. Measured stacked-sky peak heights from our own data are
# carried in separate, clearly-labeled columns (meas_peak_apo / meas_peak_lco).
#
# wavecal_class assignments are FROZEN outcomes of the 2026-09 pass-1 wavelength-
# solution census and offline association validation (2,088 good-FPI tele-nights;
# theft rates 25-63% -> <=0.1% typical with class-A guards + the quartic rough model):
#   B = associate + fit + USE in wavelength solution (18 lines);
#   A = associate + fit, EXCLUDED from the solution (4 association-guard lines);
#   X = reviewed and EXCLUDED from association entirely (1 line: OH(5-3)R1(2.5)
#       16442.15 A, whose feature falls off the usable detector at both telescopes);
#   "" = candidate pool row, inert (not consumed by the pipeline).
# The script asserts that every classed line passes the pool cuts and that the 18
# class-B centroids agree with the previous production list at the mA level (any
# discrepancy is REPORTED, never silently adopted).
#
# Usage (from the repo root):  julia --project=. scripts/cal/make_sky_linelist.jl
# Output: data/APOGEE_sky_linelist.csv (+ verification report on stdout)

using CSV, DataFrames, Printf, Statistics, SHA

const REPO = dirname(dirname(@__DIR__))
const SRCDIR = joinpath(REPO, "data", "linelist_sources")
const OUTFILE = joinpath(REPO, "data", "APOGEE_sky_linelist.csv")
const GENDATE = "2026-09-12"

# physics constants
const KB_CM1 = 0.69503476        # Boltzmann constant, cm^-1 / K
const T_ROT = 200.0              # nominal airglow rotational temperature, K
const KT = KB_CM1 * T_ROT

# pool cuts
const WAV_MIN = 15125.0          # A, blue edge of chip B range (roughwave_dict)
const WAV_MAX = 16960.0          # A, red edge of chip R range
const MAX_SEP = 2.5              # A, maximum doublet separation for the pool
const MIN_REL_INT = 0.002        # of the brightest pool line
const ISOL_RADIUS = 0.8          # A, isolation search radius around each component
const ISOL_FRAC = 0.15           # contaminant intensity fraction that breaks isolation

# chip wavelength ranges (from data/roughwave_dict.jld2; APO and LCO agree to <15 A)
const CHIP_RANGES = [("B", 15125.0, 15834.0), ("G", 15834.0, 16454.0), ("R", 16454.0, 16960.0)]

# ---------------------------------------------------------------------------------
# 0. verify source checksums
# ---------------------------------------------------------------------------------
checks = Dict{String, String}()
for ln in eachline(joinpath(SRCDIR, "CHECKSUMS.sha256"))
    (startswith(ln, "#") || isempty(strip(ln))) && continue
    h, f = split(strip(ln))
    checks[f] = h
end
sources = ["16O-1H__MoLLIST-OH.states", "16O-1H__MoLLIST-OH.trans",
    "oliva2015_table1.dat", "rousselot_vac.txt", "drp_airglow.txt"]
raw = Dict{String, Vector{UInt8}}()
for name in sources
    raw[name] = read(joinpath(SRCDIR, name))
    got = bytes2hex(sha256(raw[name]))
    got == checks[name] ||
        error("checksum mismatch for $name: got $got, expected $(checks[name])")
end
println("source checksums verified ($(length(sources)) files)")

# ---------------------------------------------------------------------------------
# 1. parse MoLLIST states + trans; build Meinel Lambda-doublets in the window
# ---------------------------------------------------------------------------------
struct OHState
    E::Float64      # cm^-1
    g::Int          # total degeneracy
    J::Float64
    parity::String  # e / f
    v::Int
    F::Int          # 1 = F1 (Omega = 3/2), 2 = F2 (Omega = 1/2)
    state::String
end

states = Dict{Int, OHState}()
for ln in split(String(raw["16O-1H__MoLLIST-OH.states"]), '\n')
    isempty(strip(ln)) && continue
    p = split(ln)
    states[parse(Int, p[1])] = OHState(parse(Float64, p[2]), parse(Int, p[3]),
        parse(Float64, p[4]), String(p[5]), parse(Int, p[6]), parse(Int, p[7]), String(p[8]))
end

# lowest rotational level of each X-state vibrational manifold (for E_rot of the model)
E0v = Dict{Int, Float64}()
for s in values(states)
    s.state == "X(2PI)" || continue
    E0v[s.v] = min(get(E0v, s.v, Inf), s.E)
end

struct OHLine
    u::Int
    l::Int
    A::Float64
    nu::Float64     # cm^-1
    lam::Float64    # vacuum A
end

window_lines = OHLine[]   # all X-X lines in a padded window (for isolation checks)
for ln in split(String(raw["16O-1H__MoLLIST-OH.trans"]), '\n')
    isempty(strip(ln)) && continue
    p = split(ln)
    nu = parse(Float64, p[4])
    nu > 100 || continue
    lam = 1e8 / nu
    (WAV_MIN - 10 <= lam <= WAV_MAX + 10) || continue
    u, l = parse(Int, p[1]), parse(Int, p[2])
    (states[u].state == "X(2PI)" && states[l].state == "X(2PI)") || continue
    push!(window_lines, OHLine(u, l, parse(Float64, p[3]), nu, lam))
end
println("MoLLIST X-X lines in window: $(length(window_lines))")

branch_letter(dJ) = dJ == -1 ? "P" : dJ == 0 ? "Q" : dJ == 1 ? "R" : ""
# cold-model intensity of one component (arbitrary units, common scale inside a band)
model_comp_int(li) = li.A * states[li.u].g * exp(-(states[li.u].E - E0v[states[li.u].v]) / KT)

# group MAIN-branch lines into Lambda-doublets keyed by (v', v'', branch, F, J'')
doublet_groups = Dict{Tuple{Int, Int, String, Int, Float64}, Vector{OHLine}}()
for li in window_lines
    su, sl = states[li.u], states[li.l]
    su.F == sl.F || continue                       # main branches only
    br = branch_letter(su.J - sl.J)
    isempty(br) && continue
    key = (su.v, sl.v, br, su.F, sl.J)
    push!(get!(doublet_groups, key, OHLine[]), li)
end

# ---------------------------------------------------------------------------------
# 2. cross-check tables
# ---------------------------------------------------------------------------------
# Oliva et al. 2015 table1: fixed-width; one row per doublet
# cols: lambda1 (1-9), Iden1 (11-24), lambda2 (26-34), Iden2 (36-49), Flux (51-55)
struct OlivaRow
    lam::Vector{Float64}
    parity::Vector{String}
    key::Tuple{Int, Int, String, Int, Float64}
    flux::Float64
end
oliva = OlivaRow[]
let iden_re = r"\[(\d+)-(\d+)\]([PQR])([12])([ef])\((\d+\.\d)\)"
    for ln in split(String(raw["oliva2015_table1.dat"]), '\n')
        length(strip(ln)) > 50 || continue
        l1 = tryparse(Float64, strip(ln[1:9]))
        l2 = tryparse(Float64, strip(ln[26:34]))
        fx = tryparse(Float64, strip(ln[51:55]))
        m1 = match(iden_re, ln[11:24])
        m2 = match(iden_re, ln[36:49])
        (isnothing(l1) || isnothing(l2) || isnothing(fx) || isnothing(m1) || isnothing(m2)) &&
            continue
        key = (parse(Int, m1[1]), parse(Int, m1[2]), String(m1[3]), parse(Int, m1[4]),
            parse(Float64, m1[6]))
        push!(oliva, OlivaRow([l1, l2], [String(m1[5]), String(m2[5])], key, fx))
    end
end
println("Oliva 2015 doublets parsed: $(length(oliva))")

rousselot = Float64[]
for ln in split(String(raw["rousselot_vac.txt"]), '\n')
    v = tryparse(Float64, strip(ln))
    isnothing(v) || push!(rousselot, v)
end

struct DRPRow
    wave::Float64
    emission::Float64
    doublet::Int
    dbl_wsep::Float64
    usewave::Int
end
drp = DRPRow[]
for ln in split(String(raw["drp_airglow.txt"]), '\n')
    (startswith(strip(ln), "#") || isempty(strip(ln))) && continue
    p = split(ln)
    length(p) >= 10 || continue
    push!(drp, DRPRow(parse(Float64, p[4]), parse(Float64, p[5]), parse(Int, p[6]),
        parse(Float64, p[7]), parse(Int, p[10])))
end

# our measured stacked-sky peak heights (see header); columns tele,mjd,chip,feature,lam,peak_height
meas = CSV.read(joinpath(SRCDIR, "stack_sky_peaks.csv"), DataFrame; comment = "#")

# ---------------------------------------------------------------------------------
# 3. build the candidate pool
# ---------------------------------------------------------------------------------
function chip_of(lam)
    for (c, lo, hi) in CHIP_RANGES
        lo <= lam <= hi && return c
    end
    return ""
end

pool = DataFrame()
for (key, comps) in doublet_groups
    length(comps) == 2 || continue
    vu, vl, br, F, Jl = key
    # sort components by wavelength (ascending)
    sort!(comps, by = li -> li.lam)
    lam1, lam2 = comps[1].lam, comps[2].lam
    (WAV_MIN <= lam1 && lam2 <= WAV_MAX) || continue
    sep = lam2 - lam1
    sep <= MAX_SEP || continue
    ints = model_comp_int.(comps)
    cen = (lam1 + lam2) / 2
    push!(pool,
        (vu = vu, vl = vl, branch = br, F = F, Jlow = Jl,
            label = "OH($(vu)-$(vl))$(br)$(F)($(Jl))",
            chip = chip_of(cen), wave_cen_ang = cen,
            subwave_1_ang = lam1, subwave_2_ang = lam2, sep_ang = sep,
            parity_1 = states[comps[1].u].parity, parity_2 = states[comps[2].u].parity,
            einA_1 = comps[1].A, einA_2 = comps[2].A,
            Eup_1_cm1 = states[comps[1].u].E, Eup_2_cm1 = states[comps[2].u].E,
            u1 = comps[1].u, l1 = comps[1].l, u2 = comps[2].u, l2 = comps[2].l,
            model_int_raw = sum(ints), model_int_weak = minimum(ints)))
end

# per-band empirical scale from Oliva fluxes (median flux / cold-model unit).
# Bands with no Oliva counterpart in the window are not airglow-emitting bands
# (airglow OH is populated only to v' <= 9 by H + O3; MoLLIST also carries e.g.
# (13-9), (11-8) transitions that never appear in the sky) and are dropped.
oliva_by_key = Dict(o.key => o for o in oliva)
pool.oliva_flux = [haskey(oliva_by_key, (r.vu, r.vl, r.branch, r.F, r.Jlow)) ?
                   oliva_by_key[(r.vu, r.vl, r.branch, r.F, r.Jlow)].flux : NaN
                   for r in eachrow(pool)]
# scale from LOW-ROTATIONAL-EXCITATION matches only (E_rot <= 500 cm^-1): the cold
# 200 K LTE model reproduces Oliva to ~+-10% there, while the non-LTE hot-OH
# component makes high-J lines orders of magnitude brighter than any cold model.
# (All classed lines are low-J; high-J hot-OH lines simply fall out of the pool
# on the model intensity threshold, which is the intended conservative behavior.)
pool.Erot_up = [min(r.Eup_1_cm1, r.Eup_2_cm1) - E0v[r.vu] for r in eachrow(pool)]
band_scale = Dict{Tuple{Int, Int}, Float64}()
for b in unique([(r.vu, r.vl) for r in eachrow(pool)])
    msk = [(r.vu, r.vl) == b && !isnan(r.oliva_flux) && r.Erot_up <= 500.0
           for r in eachrow(pool)]
    if count(msk) > 0
        band_scale[b] = median(pool.oliva_flux[msk] ./ pool.model_int_raw[msk])
    end
end
println("band scales (Oliva flux units per cold-model unit; low-J matches only):")
for (b, s) in sort(collect(band_scale))
    n = count([(r.vu, r.vl) == b && !isnan(r.oliva_flux) && r.Erot_up <= 500.0
               for r in eachrow(pool)])
    @printf("  (%d-%d): %.3g  (n_match = %d)\n", b[1], b[2], s, n)
end
pool = pool[[haskey(band_scale, (r.vu, r.vl)) for r in eachrow(pool)], :]
pool.model_rel_int = [band_scale[(r.vu, r.vl)] * r.model_int_raw for r in eachrow(pool)]
pool.model_rel_int ./= maximum(pool.model_rel_int)
pool = pool[pool.model_rel_int .>= MIN_REL_INT, :]

# isolation: no other OH airglow line with band-scaled component intensity
# >= ISOL_FRAC of the candidate's weaker (band-scaled) component within
# ISOL_RADIUS of either component. Contaminants are drawn from ALL branches
# (satellites included) of the airglow-emitting (Oliva-matched) bands.
band_of_state(u) = states[u].v
function scaled_int(li)
    b = (states[li.u].v, states[li.l].v)
    haskey(band_scale, b) || return 0.0    # non-airglow band: never a contaminant
    return band_scale[b] * model_comp_int(li)
end
isol_ok = trues(nrow(pool))
for (i, r) in enumerate(eachrow(pool))
    thresh = ISOL_FRAC * band_scale[(r.vu, r.vl)] * r.model_int_weak
    for li in window_lines
        (li.u == r.u1 && li.l == r.l1) && continue
        (li.u == r.u2 && li.l == r.l2) && continue
        if (abs(li.lam - r.subwave_1_ang) < ISOL_RADIUS ||
            abs(li.lam - r.subwave_2_ang) < ISOL_RADIUS) && scaled_int(li) >= thresh
            isol_ok[i] = false
            break
        end
    end
end
pool = pool[isol_ok, :]
select!(pool, Not([:u1, :l1, :u2, :l2, :model_int_weak]))
sort!(pool, :wave_cen_ang)
# fixed 1:1 Lambda-doublet component ratio (see header rationale)
pool.weight_1 = fill(0.5, nrow(pool))
pool.weight_2 = fill(0.5, nrow(pool))
println("candidate pool after cuts: $(nrow(pool)) Lambda-doublets")

# ---------------------------------------------------------------------------------
# 4. frozen class assignments (2026-09 census + offline association validation)
# ---------------------------------------------------------------------------------
# key = (v', v'', branch, F, J'') => (class, expected centroid A, note)
CLASSMAP = Dict(
    # chip B, class B
    (3, 1, "P", 2, 1.5) => ("B", 15187.13, "healthy on 2088 good-FPI tele-nights"),
    (3, 1, "P", 1, 2.5) => ("B", 15240.96, "healthy"),
    (3, 1, "P", 2, 2.5) => ("B", 15287.79, "healthy"),
    (3, 1, "P", 2, 3.5) => ("B", 15395.34, "healthy"),
    (3, 1, "P", 1, 4.5) => ("B", 15432.17, "healthy"),
    (3, 1, "P", 1, 5.5) => ("B", 15540.39,
        "healthy but theft victim of 15546.13 (90 apo flip nights); guarded by class-A OH(4-2)R1(3.5)"),
    (4, 2, "R", 1, 2.5) => ("B", 15597.65, "healthy"),
    # chip G, class B
    (4, 2, "P", 2, 1.5) => ("B", 15972.61, "healthy"),
    (4, 2, "P", 1, 2.5) => ("B", 16030.82, "healthy"),
    (4, 2, "P", 2, 2.5) => ("B", 16079.75, "healthy"),
    (4, 2, "P", 1, 3.5) => ("B", 16128.61, "healthy"),
    (4, 2, "P", 1, 4.5) => ("B", 16235.36, "healthy"),
    (4, 2, "P", 1, 5.5) => ("B", 16351.29, "healthy"),
    # chip R, class B
    (5, 3, "R", 1, 1.5) => ("B", 16502.35, "healthy"),
    (5, 3, "R", 2, 0.5) => ("B", 16553.81, "healthy"),
    (5, 3, "Q", 1, 1.5) => ("B", 16692.37,
        "healthy when correctly associated; theft victim of the 16702.0/16703.2 pair (67 apo flip nights); guarded by class-A OH(5-3)Q2(1.5)"),
    (5, 3, "Q", 1, 2.5) => ("B", 16708.85, "healthy 0.81 A resolved doublet"),
    (5, 3, "P", 1, 2.5) => ("B", 16903.67,
        "blend watch: stable -0.013 A pull both telescopes (external non-Lambda blend nearby); keep B, monitor"),
    # class A association guards (fit, never in the wavelength solution)
    (4, 2, "R", 1, 3.5) => ("A", 15546.13,
        "guard: measured 2.05x (apo) / 1.36x (lco) brighter than its 15540.39 victim in stacked sky; below the old intensity gate"),
    (5, 3, "R", 2, 2.5) => ("A", 16414.74,
        "guard: the real feature the retired 16442.15 label always landed on (delta +27.41 A both telescopes); chip-G blue edge"),
    (5, 3, "Q", 2, 0.5) => ("A", 16689.18,
        "guard: 3.2 A blueward of the brightest chip-R line; DRP USELSF=1 USEWAVE=0; safe only with the quartic rough model"),
    (5, 3, "Q", 2, 1.5) => ("A", 16702.64,
        "guard: the 16692.37 thief; resolved 5.1-px pair at ~1/28 of its victim's model flux, hence guard not class B"),
    # class X: excluded from association entirely
    (5, 3, "R", 1, 2.5) => ("X", 16442.15,
        "feature at rough pixel ~3(lco)/13(apo), inside the 64-px segment trim at BOTH telescopes; label can only mis-associate (100% wrong when measured)"))

pool.wavecal_class = fill("", nrow(pool))
pool.class_note = fill("", nrow(pool))
found = Set{Tuple{Int, Int, String, Int, Float64}}()
for r in eachrow(pool)
    key = (r.vu, r.vl, r.branch, r.F, r.Jlow)
    if haskey(CLASSMAP, key)
        cls, expect, note = CLASSMAP[key]
        # identification is by quantum numbers; the expected value (the OLD list's
        # Rousselot-based centroid) can differ from the equal-weight Brooke centroid
        # by up to ~60 mA for the wider doublets -- see the section-6 report
        abs(r.wave_cen_ang - expect) < 0.10 ||
            error("classed line $(r.label): centroid $(r.wave_cen_ang) != expected $expect")
        r.wavecal_class = cls
        r.class_note = note
        push!(found, key)
    end
end
missing_keys = setdiff(keys(CLASSMAP), found)
isempty(missing_keys) ||
    error("classed lines missing from the candidate pool (check cuts): $missing_keys")
@printf("classes assigned: B = %d, A = %d, X = %d, pool-only = %d\n",
    count(pool.wavecal_class .== "B"), count(pool.wavecal_class .== "A"),
    count(pool.wavecal_class .== "X"), count(pool.wavecal_class .== ""))

# ---------------------------------------------------------------------------------
# 5. cross-checks + measured amplitudes
# ---------------------------------------------------------------------------------
# (a) component positions vs Oliva (classed rows must agree to < 0.03 A)
println("\ncomponent-position cross-check vs Oliva 2015 (classed rows):")
for r in eachrow(pool)
    isempty(r.wavecal_class) && continue
    o = get(oliva_by_key, (r.vu, r.vl, r.branch, r.F, r.Jlow), nothing)
    isnothing(o) && (println("  $(r.label): no Oliva match"); continue)
    dv = [minimum(abs.(o.lam .- r.subwave_1_ang)), minimum(abs.(o.lam .- r.subwave_2_ang))]
    @printf("  %-18s dlam = %+.1f / %+.1f mA\n", r.label,
        1e3 * (sort(o.lam)[1] - r.subwave_1_ang), 1e3 * (sort(o.lam)[2] - r.subwave_2_ang))
    maximum(dv) < 0.03 || error("$(r.label): Oliva position mismatch $(maximum(dv)) A")
end

# (b) Rousselot centroid cross-check (report only; his Q splittings are known-wrong)
pool.rousselot_cen_ang = [isempty(rousselot) ? NaN :
                          (d = abs.(rousselot .- r.wave_cen_ang);
    minimum(d) < 0.3 ? rousselot[argmin(d)] : NaN) for r in eachrow(pool)]

# (c) DRP airglow cross-check: nearest entry within 0.10 A of centroid or a component
function drp_usewave_of(r)
    for d in drp
        if abs(d.wave - r.wave_cen_ang) < 0.10 ||
           abs(d.wave - r.subwave_1_ang) < 0.10 || abs(d.wave - r.subwave_2_ang) < 0.10
            return d.usewave
        end
    end
    return -1
end
pool.drp_usewave = [drp_usewave_of(r) for r in eachrow(pool)]

# (d) measured stacked-sky peak heights (max component height within 0.7 A, night 60105)
function meas_peak(r, tele)
    m = meas[(meas.tele .== tele) .&
             (min.(abs.(meas.lam .- r.subwave_1_ang), abs.(meas.lam .- r.subwave_2_ang),
                 abs.(meas.lam .- r.wave_cen_ang)) .< 0.7), :]
    vals = filter(!isnan, m.peak_height)
    return isempty(vals) ? NaN : maximum(vals)
end
pool.meas_peak_apo = [meas_peak(r, "apo") for r in eachrow(pool)]
pool.meas_peak_lco = [meas_peak(r, "lco") for r in eachrow(pool)]

# ---------------------------------------------------------------------------------
# 6. verification vs the previous production list (report, never silently adopt)
# ---------------------------------------------------------------------------------
old = CSV.read(joinpath(REPO, "data", "APOGEE_lines.csv"), DataFrame)
parsevec(s) = parse.(Float64, filter(!isempty, split(replace(s[2:(end - 1)], r"\s+" => ","), ",")))
function old_effective_centroid(row)
    # replicate the old get_subline_params consumption: top-2 sublines by weight,
    # renormalized, wavelength-weighted mean (in A)
    w = parsevec(row.subline_I)
    lam = parsevec(row.subline_wav)
    si = sortperm(w, rev = true)[1:2]
    w2 = w[si] ./ sum(w[si])
    return sum(lam[si] .* w2) * 10
end
println("\nclass-B centroids vs previous production list (data/APOGEE_lines.csv):")
println("  line               new cen (A)    vs old wav (mA)  vs old consumed outwave (mA)")
for r in eachrow(pool)
    r.wavecal_class == "B" || continue
    d = abs.(old.wav .* 10 .- r.wave_cen_ang)
    i = argmin(d)
    d[i] < 0.3 || error("$(r.label): no production match within 0.3 A")
    @printf("  %-18s %.4f     %+8.1f         %+8.1f\n", r.label, r.wave_cen_ang,
        1e3 * (r.wave_cen_ang - old.wav[i] * 10),
        1e3 * (r.wave_cen_ang - old_effective_centroid(old[i, :])))
end

# ---------------------------------------------------------------------------------
# 7. write the output file
# ---------------------------------------------------------------------------------
r4(x) = round(x, digits = 4)
r5(x) = round(x, sigdigits = 5)
out = select(pool,
    :label, :chip, :wavecal_class,
    :wave_cen_ang => ByRow(r4) => :wave_cen_ang,
    :subwave_1_ang => ByRow(r4) => :subwave_1_ang,
    :subwave_2_ang => ByRow(r4) => :subwave_2_ang,
    :sep_ang => ByRow(r4) => :sep_ang,
    :parity_1, :parity_2, :weight_1, :weight_2,
    :einA_1 => ByRow(r5) => :einA_1, :einA_2 => ByRow(r5) => :einA_2,
    :Eup_1_cm1 => ByRow(r4) => :Eup_1_cm1, :Eup_2_cm1 => ByRow(r4) => :Eup_2_cm1,
    :model_rel_int => ByRow(r5) => :model_rel_int,
    :oliva_flux,
    :meas_peak_apo => ByRow(r5) => :meas_peak_apo,
    :meas_peak_lco => ByRow(r5) => :meas_peak_lco,
    :rousselot_cen_ang, :drp_usewave, :class_note)

hdr = """
# APOGEE_sky_linelist.csv -- sky (airglow) line list for the ApogeeReduction.jl
# wavelength calibration. GENERATED by scripts/cal/make_sky_linelist.jl on $(GENDATE);
# never edit by hand -- edit the generator and rerun.
#
# All wavelengths VACUUM, Angstroms. Sublines ascending (subwave_1_ang < subwave_2_ang);
# component parity labels (parity_1/2) are the UPPER-state e/f parity, matching the
# Oliva et al. 2015 doublet notation. wave_cen_ang = weighted centroid (weights below).
#
# SOURCES (cached in data/linelist_sources/, sha256-verified at generation):
#   positions/Einstein A/labels/splittings: Brooke et al. 2016, JQSRT 168, 142
#     (MoLLIST OH via ExoMol "MoLLIST-OH" states+trans, retrieved 2026-09-12;
#     H-band position accuracy ~0.01 A, Abrams et al. 1994 FTS lineage)
#   empirical per-band intensity scale + position cross-check: Oliva et al. 2015,
#     A&A 581, A47 (VizieR J/A+A/581/A47 table1, retrieved 2026-09-12)
#   centroid cross-check only: Rousselot et al. 2000, A&A 354, 1134 (vacuum-converted
#     copy from SDSS apogee_drp; Q-branch SPLITTINGS from this source are known-wrong
#     and were never used)
#   cross-check only: SDSS apogee_drp data/skylines/airglow.txt (drp_usewave column;
#     -1 = no DRP counterpart within 0.10 A)
#   meas_peak_apo/lco: OUR stacked-sky peak heights (median running-median-scaled flux,
#     healthy night MJD 60105, pass-1 DR21 products; NaN = not measured there) --
#     MEASURED amplitudes, vs model_rel_int which is the MODEL amplitude estimate.
#
# MODEL AMPLITUDES (model_rel_int): Einstein-A-based cold-atmosphere estimate at
# T_rot = 200 K: per component A * g_up * exp(-(E_up - E0_v')/kT), summed over the
# doublet, each vibrational band scaled to the Oliva 2015 empirical fluxes (median
# over low-J matches, E_rot <= 500 cm^-1; the cold model matches Oliva to ~+-10%
# there, while non-LTE hot-OH makes high-J lines far brighter than any cold model),
# normalized to the brightest pool line = 1. Component weights are FIXED 1:1 (weight_1 = weight_2
# = 0.5): the Lambda-doublet e/f component ratio is 1.000 +/- 0.002, temperature-
# independent to <= 1e-3 (Einstein-A e/f equality <= 0.03% in this very list; LTE
# sublevel-population equality <= 0.2%), and our 24-night stacked-sky measurements
# confirm ratios consistent with 1 at the 1.5-11% measurement-noise level.
#
# SELECTION CUTS (candidate pool): OH X2Pi-X2Pi Meinel MAIN-branch Lambda-doublets,
# both components in [$(WAV_MIN), $(WAV_MAX)] A, separation <= $(MAX_SEP) A,
# model_rel_int >= $(MIN_REL_INT), and isolation (no other OH Meinel line with
# component intensity >= $(Int(100 * ISOL_FRAC))% of the candidate's weaker component
# within $(ISOL_RADIUS) A of either component).
#
# wavecal_class (FROZEN outcomes of the 2026-09 pass-1 census + offline association
# validation; changing a class is a science decision, not a regeneration):
#   B  = associate + fit + USE in the wavelength solution (18 lines)
#   A  = associate + fit, EXCLUDED from the solution (4 association-guard lines).
#        Guards make the nearest-catalog-line association claim the thief features
#        that otherwise steal the labels of nearby class-B lines. Only safe together
#        with the quartic rough wavelength model (data/roughwave_dict.jld2).
#   X  = reviewed, EXCLUDED from association entirely (1 line)
#   "" = candidate pool row, inert (not consumed; available for future promotion)
#
# Consumed by src/skyline_peaks.jl (get_sky_peaks): class A+B rows are associated and
# fit; class recorded in the skyLinePeaks product; class-A rows are NaN-masked at
# wavelength-solution ingest (src/wavecal.jl ingest_skyLines_file).
"""

open(OUTFILE, "w") do io
    write(io, hdr)
    CSV.write(io, out; append = true, writeheader = true)
end
println("\nwrote $(OUTFILE): $(nrow(out)) rows")
