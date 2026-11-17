using HDF5, Random, Statistics
using ApogeeReduction: trace_peak_census, trace_identify_slit, fill_unlit_fibers_robust!,
                       trace_params_qa, trace_id_match, trace_gate, select_night_trace,
                       TraceCandidate, trace_fiber_bits, trace_bits_from_metadata,
                       write_trace_qa, read_trace_qa, trace_qa_path, make_traceMain,
                       trace_reference_positions, trace_level_bits, trace_bit_names,
                       night_science_mid_exposure, read_metadata, flat_trace_qa,
                       trace_light_qa, trace_bits_for_1d, trace_exposure_summary,
                       trace_extract, get_default_trace_hyperparams, bad_pix_bits,
                       TRACE_FIBER_UNLIT, TRACE_FIXED_FALLBACK, TRACE_UNGATED, TRACE_QA_MISSING,
                       TRACE_OTHER_FLAT, TRACE_TEMPLATE_ONLY, TRACE_BAD_BITS,
                       TRACE_LEVEL_OWN_FLAT, TRACE_LEVEL_FIXED_FALLBACK, TRACE_LEVEL_NONE,
                       TRACE_LEVEL_TEMPLATE, TRACE_GATE_NLIT_MIN

const DATA_DIR = joinpath(@__DIR__, "..", "data")

# a center-band profile: Gaussian fibers at `pos` (NaN = absent) with heights `h`
function synth_profile(pos, h; width = 1.0, noise = 3.0, rng = MersenneTwister(1), n = 2048)
    p = zeros(n)
    for (y0, hh) in zip(pos, h)
        isfinite(y0) || continue
        for y in max(1, floor(Int, y0 - 8)):min(n, ceil(Int, y0 + 8))
            p[y] += hh * exp(-0.5 * ((y - y0) / width)^2)
        end
    end
    p .+= noise .* randn(rng, n)
    return p, fill(noise, n)
end

# a per-column trace (n_x, 300, 3): smooth curved traces following the reference spacing
function synth_trace(R; nx = 2048, shift = 0.0, height = 1000.0)
    tp = zeros(nx, length(R), 3)
    for x in 1:nx
        u = (x - 1024) / 1024
        tp[x, :, 1] .= height
        tp[x, :, 2] .= R .+ shift .+ 0.8 * u .- 1.5 * u^2 .+ 0.0004 .* (R .- 1024) .* u
        tp[x, :, 3] .= 1.0 + 0.05 * u
    end
    return tp
end

# a 2048 x 2048 flat image[x, y] with Gaussian fibers along `tp` (heights from tp[:, :, 1])
function synth_image(tp; noise = 3.0, rng = MersenneTwister(11), width = 1.0)
    nx = size(tp, 1)
    img = noise .* randn(rng, nx, 2048)
    for x in 1:nx, f in axes(tp, 2)
        c = tp[x, f, 2]
        h = tp[x, f, 1]
        h == 0 && continue
        for y in max(1, floor(Int, c - 6)):min(2048, ceil(Int, c + 6))
            img[x, y] += h * exp(-0.5 * ((y - c) / width)^2)
        end
    end
    return img, fill(1 / noise^2, nx, 2048), zeros(Int, nx, 2048)
end

@testset "trace QA and robust identification" begin
    R = trace_reference_positions("apo", "G"; data_path = DATA_DIR)
    @test length(R) == 300 && issorted(R)
    sp = median(diff(R))

    @testset "census and whole-slit identification" begin
        for shift in (0.0, 1.7, -2.9, 3.3, 12.4)
            p, e = synth_profile(R .+ shift, fill(1000.0, 300))
            cen = trace_peak_census(p, e)
            @test length(cen.y) == 300
            @test cen.nslots == 300
            id = trace_identify_slit(cen.y, R; census_nslots = cen.nslots, census_slots = cen.slots)
            @test id.ok
            @test id.label == collect(1:300)       # labels never alias, even half a spacing off
            @test isapprox(id.a + (id.b - 1) * 1024, shift; atol = 0.3)
            @test id.margin > 2
        end
    end

    @testset "dark bundle: identified and placed from both sides" begin
        h = fill(1000.0, 300)
        h[151:180] .= 0.0                           # the APO 57969-58015 dark MTP bundle
        truth = R .+ 1.2
        p, e = synth_profile(truth, h)
        cen = trace_peak_census(p, e)
        @test length(cen.y) == 270
        id = trace_identify_slit(cen.y, R; census_nslots = cen.nslots, census_slots = cen.slots)
        @test id.ok
        @test count(id.matched) == 270
        @test !any(id.matched[151:180])
        @test maximum(abs.(id.positions[151:180] .- truth[151:180])) < 0.2
    end

    @testset "refusals are named, never silent relabels" begin
        # light beyond fiber 300: a slit relabelled by one fiber would explain it
        p, e = synth_profile(vcat(R, R[end] + sp), fill(1000.0, 301))
        cen = trace_peak_census(p, e)
        id = trace_identify_slit(cen.y, R; census_nslots = cen.nslots, census_slots = cen.slots)
        @test !id.ok
        @test id.reason == "id_light_outside_slit"
        # too few peaks
        id = trace_identify_slit(collect(100.0:7.0:150.0), R)
        @test !id.ok && id.reason == "id_too_few_peaks"
        # every outcome is either the right labelling or a refusal (edge-dark fibers, shifts
        # near half a spacing): never a wrong labelling accepted
        for shift in (-3.4, -3.0, 3.0, 3.4), ndark in (5, 20)
            h = fill(1000.0, 300)
            h[1:ndark] .= 0.0
            p, e = synth_profile(R .+ shift, h)
            cen = trace_peak_census(p, e)
            id = trace_identify_slit(cen.y, R; census_nslots = cen.nslots, census_slots = cen.slots)
            if id.ok
                @test all(id.label .== collect((ndark + 1):300))
            else
                @test id.reason in ("id_ambiguous", "id_light_outside_slit", "id_unmatched_light")
            end
        end
    end

    @testset "unlit fibers: robust fill vs two-neighbour extrapolation" begin
        tp = synth_trace(R; shift = 0.7)
        truth = copy(tp)
        # per-fiber pattern noise like real fits, then a dark bundle
        rng = MersenneTwister(3)
        tp[:, :, 2] .+= 0.004 .* randn(rng, size(tp, 1), 300)
        unlit = collect(151:180)
        tp[:, unlit, 2] .= NaN
        fill_unlit_fibers_robust!(tp, unlit, R)
        err_robust = maximum(abs.(tp[:, unlit, 2] .- truth[:, unlit, 2]))
        @test err_robust < 0.05
        @test all(tp[:, unlit, 1] .== 0)
        # production's rule: linear extrapolation from the two nearest lit fibers in index
        pr = copy(truth)
        pr[:, :, 2] .+= 0.004 .* randn(MersenneTwister(3), size(pr, 1), 300)
        for f in unlit
            n1, n2 = 150, 149
            pr[:, f, 2] .= pr[:, n1, 2] .+ (pr[:, n1, 2] .- pr[:, n2, 2]) .* (f - n1) ./ (n1 - n2)
        end
        err_prod = maximum(abs.(pr[:, unlit, 2] .- truth[:, unlit, 2]))
        @test err_prod > 5 * err_robust
    end

    @testset "trace metrics and gate" begin
        tp = synth_trace(R)
        tp[:, :, 2] .+= 0.003 .* randn(MersenneTwister(5), size(tp, 1), 300)
        q = trace_params_qa(tp)
        @test q.nlit == 300
        @test q.fit_rms < 0.01
        @test q.spacing_rms < 0.01
        @test q.spacing_nbad == 0
        # a block of fibers latching half a spacing off over part of the chip
        bad = copy(tp)
        bad[1300:end, 151:200, 2] .+= 0.5 * sp
        qb = trace_params_qa(bad)
        @test qb.spacing_nbad > 0
        # identification agreement of the central positions
        p, e = synth_profile(tp[1025, :, 2], fill(1000.0, 300))
        cen = trace_peak_census(p, e)
        id = trace_identify_slit(cen.y, R; census_nslots = cen.nslots, census_slots = cen.slots)
        @test trace_id_match(q.center_mid, q.lit, id) > 0.99
        shifted = copy(q.center_mid)
        shifted[1:100] .+= sp                       # a relabelled block
        @test trace_id_match(shifted, q.lit, id) < 0.7

        good = Dict{String, Any}("outcome" => "ok", "lamp_ok" => true, "has_trace" => true,
            "nlit" => 300, "fit_rms" => 0.005, "light_mismatch_frac" => 0.0,
            "id_match_frac" => 1.0, "id_ok" => true, "id_reason" => "ok", "sat_frac" => 0.0)
        @test trace_gate(good) == (2, String[])
        sat = merge(good, Dict("sat_frac" => 0.2))
        @test trace_gate(sat) == (1, ["saturated"])
        l, r = trace_gate(merge(good, Dict("fit_rms" => 0.3, "nlit" => 50)))
        @test l == 0 && "fit_rms" in r && "nlit" in r
        l, r = trace_gate(Dict{String, Any}("outcome" => "lamp_off", "lamp_ok" => false,
            "has_trace" => false, "id_ok" => false, "id_reason" => "id_too_few_peaks", "sat_frac" => 0.0))
        @test l == 0 && r[1] == "outcome:lamp_off" && "lamp_off" in r && "no_trace" in r
    end

    @testset "QA record of a flat image, light check" begin
        tp = synth_trace(R; shift = 1.3)
        tp[:, 151:180, 1] .= 0.0
        img, iv, pm = synth_image(tp)
        q = flat_trace_qa(img, iv, pm, "apo", "G"; trace_params = tp, refpos = R)
        @test q["outcome"] == "ok" && q["lamp_ok"] && q["id_ok"]
        @test q["nlit"] == 270 && q["id_nlit"] == 270
        @test q["light_mismatch_frac"] == 0.0
        @test q["light_offset_med"] < 0.1
        @test q["id_match_frac"] == 1.0
        @test isapprox(q["trace_shift_px"], 1.3 - 1.5 * 0.0; atol = 0.2)
        @test q["gate_level"] == 2 && q["gate_reasons"] == ""
        # the same light with a trace that twists half the slit onto the wrong place at x > 1500
        # (LCO 60702 G/R in 0.4.0): every fiber is smooth, the light check fails it
        bad = copy(tp)
        for x in 1500:2048
            bad[x, 151:300, 2] .+= 3.5 * (x - 1500) / 548
        end
        qb = flat_trace_qa(img, iv, pm, "apo", "G"; trace_params = bad, refpos = R)
        @test qb["light_mismatch_frac"] > 0.3
        @test qb["gate_level"] == 0 && occursin("light", qb["gate_reasons"])
        # lamp off: a named outcome with its metrics
        dark = 2.0 .* randn(MersenneTwister(2), 2048, 2048)
        qd = flat_trace_qa(dark, iv, pm, "apo", "G"; trace_params = nothing, outcome = "lamp_off",
            refpos = R)
        @test !qd["lamp_ok"] && qd["outcome"] == "lamp_off" && qd["gate_level"] == 0
        @test startswith(qd["gate_reasons"], "outcome:lamp_off")
    end

    @testset "robust trace_extract on a synthetic flat with a dark bundle" begin
        truth = synth_trace(R; shift = 2.0)
        truth[:, 151:180, 1] .= 0.0
        img, iv, pm = synth_image(truth; noise = 2.0)
        hp = get_default_trace_hyperparams("apo", "G", profile_path = DATA_DIR, plot_path = mktempdir())
        qa = Dict{String, Any}()
        tp, _ = trace_extract(img, iv, "apo", "99999", "0001", "G", hp..., "synthetic";
            good_pixels = (pm .& bad_pix_bits) .== 0, median_trace_pos_path = DATA_DIR,
            identify = "robust", qa = qa)
        @test size(tp) == (2048, 300, 3)
        @test qa["outcome"] == "ok" && qa["id_reason"] == "ok" && qa["n_unlit"] == 30
        xs = 100:1948
        lit = setdiff(1:300, 151:180)
        dl = abs.(tp[xs, lit, 2] .- truth[xs, lit, 2])
        @test median(dl) < 0.05
        du = abs.(tp[xs, 151:180, 2] .- truth[xs, 151:180, 2])
        @test maximum(du) < 0.3
        @test all(tp[:, 151:180, 1] .== 0)
        @test_throws ErrorException trace_extract(img, iv, "apo", "99999", "0001", "G", hp...,
            "synthetic"; median_trace_pos_path = DATA_DIR, identify = "template")
    end

    @testset "choice of the night's trace" begin
        qa(lvl, rms) = Dict{String, Any}("gate_level" => lvl, "fit_rms" => rms, "gate_reasons" => "")
        c = [TraceCandidate("quartz", 9, "/x/quartzTrace_apo_57969_0009_B.h5", qa(1, 0.03)),
            TraceCandidate("quartz", 4, "/x/quartzTrace_apo_57969_0004_B.h5", qa(2, 0.010)),
            TraceCandidate("quartz", 7, "/x/quartzTrace_apo_57969_0007_B.h5", qa(2, 0.008)),
            TraceCandidate("dome", 30, "/x/domeTrace_apo_57969_0030_B.h5", qa(2, 0.004))]
        ch, info = select_night_trace(c; mode = "last")
        @test ch.expnum == 9                       # the previous rule: last quartz by name
        ch, info = select_night_trace(c; mode = "gated")
        @test ch.expnum == 7                       # gated quartz, lowest fit_rms
        ch, _ = select_night_trace(c; mode = "gated", sci_mid = 5)
        @test ch.expnum == 4                       # closest to the science
        ch, _ = select_night_trace(reverse(c); mode = "gated", sci_mid = 5)
        @test ch.expnum == 4                       # independent of input order
        @test info["trace_cand_files"][1] == "quartzTrace_apo_57969_0007_B.h5"
        # only soft-passing quartz: still before nothing, but after any fully passing flat
        c2 = [c[1], TraceCandidate("dome", 30, "/x/domeTrace_apo_57969_0030_B.h5", qa(2, 0.004))]
        ch, _ = select_night_trace(c2; mode = "gated")
        @test ch.flat_type == "dome"
        ch, _ = select_night_trace([c[1]]; mode = "gated")
        @test ch.expnum == 9
        ch, _ = select_night_trace([TraceCandidate("quartz", 3, "/x/q3.h5", qa(0, 0.5)),
            TraceCandidate("quartz", 5, "/x/q5.h5", nothing)]; mode = "gated")
        @test isnothing(ch)                        # refused and QA-less flats are not eligible
        @test_throws ErrorException select_night_trace(c; mode = "best")
        @test night_science_mid_exposure((image_type = ["dark", "object", "object", "object", "arclamp"],)) == 3
        @test isnothing(night_science_mid_exposure((image_type = ["dark", "quartzflat"],)))
    end

    @testset "trace bits" begin
        tp = synth_trace(R)
        tp[:, [24, 151, 152], 1] .= 0.0
        b = trace_fiber_bits(tp, TRACE_LEVEL_OWN_FLAT)
        @test findall(x -> (x & TRACE_FIBER_UNLIT) != 0, b) == [24, 151, 152]
        @test count(!=(0), b) == 3
        b = trace_fiber_bits(tp, TRACE_LEVEL_FIXED_FALLBACK)
        @test all(b .== TRACE_FIXED_FALLBACK)      # another night's unlit fibers are not this night's
        b = trace_fiber_bits(tp, TRACE_LEVEL_OWN_FLAT; other_flat = true, ungated = true)
        @test all((b .& (TRACE_OTHER_FLAT | TRACE_UNGATED)) .== (TRACE_OTHER_FLAT | TRACE_UNGATED))
        @test trace_level_bits(TRACE_LEVEL_TEMPLATE) == TRACE_TEMPLATE_ONLY
        @test (TRACE_BAD_BITS & TRACE_FIBER_UNLIT) != 0 && (TRACE_BAD_BITS & TRACE_UNGATED) == 0
        @test trace_bit_names(TRACE_FIXED_FALLBACK | TRACE_UNGATED) ==
              ["TRACE_FIXED_FALLBACK", "TRACE_UNGATED"]
        b = trace_bits_from_metadata(tp, Dict{String, Any}("trace_type" => "quartz_fallback"))
        @test all(b .== (TRACE_FIXED_FALLBACK | TRACE_QA_MISSING))
        b = trace_bits_from_metadata(tp, Dict{String, Any}("trace_type" => "quartz"))
        @test b[24] == (TRACE_FIBER_UNLIT | TRACE_QA_MISSING) && b[1] == TRACE_QA_MISSING
    end

    @testset "trace bits in the 1D products" begin
        mktempdir() do d
            tp = synth_trace(R)
            tp[:, 30, 1] .= 0.0
            legacy = joinpath(d, "legacy.h5")
            h5open(f -> (f["regularized_trace_params"] = tp), legacy, "w")
            bits, md = trace_bits_for_1d(legacy, tp, Dict{String, Any}("trace_type" => "quartz_fallback"))
            @test all(bits .== (TRACE_FIXED_FALLBACK | TRACE_QA_MISSING))
            @test md["trace_level"] == TRACE_LEVEL_FIXED_FALLBACK
            @test md["trace_flags"] == (TRACE_FIXED_FALLBACK | TRACE_QA_MISSING)
            bits, md = trace_bits_for_1d(legacy, tp, Dict{String, Any}("trace_type" => "quartz"))
            @test bits[30] == (TRACE_FIBER_UNLIT | TRACE_QA_MISSING) && md["trace_flags"] == TRACE_QA_MISSING
            stored = joinpath(d, "new.h5")
            h5open(stored, "w") do f
                f["regularized_trace_params"] = tp
                f["bitmsk_trace"] = fill(TRACE_UNGATED, 300)
            end
            bits, md = trace_bits_for_1d(stored, tp, Dict{String, Any}("trace_type" => "quartz"))
            @test all(bits .== TRACE_UNGATED) && isempty(md)
        end
        st = zeros(Int, 3, 300)
        st[2, 7] = TRACE_FIBER_UNLIT
        st[3, :] .= TRACE_OTHER_FLAT
        s = trace_exposure_summary(st, [0, 0, 0])
        @test s["trace_flags"] == TRACE_OTHER_FLAT && s["trace_level_max"] == 0
        @test trace_exposure_summary(zeros(Int, 3, 300), [0, 4, -2])["trace_level_max"] == 4
    end

    @testset "QA records and traceMain" begin
        mktempdir() do outdir
            tele, mjd, chip = "apo", 57969, "B"
            qpath = trace_qa_path(outdir, "quartz", tele, mjd, "0004", chip)
            write_trace_qa(qpath, Dict{String, Any}("outcome" => "ok", "gate_level" => 2,
                "fit_rms" => 0.01, "gate_reasons" => "", "identify" => "production",
                "nlit" => 270))
            @test isfile(qpath)
            @test !any(endswith.(readdir(dirname(qpath)), r"\.tmp\d+"))
            @test read_trace_qa(qpath)["gate_level"] == 2
            @test isnothing(read_trace_qa(qpath * "nope"))
            write_trace_qa(trace_qa_path(outdir, "quartz", tele, mjd, "0009", chip),
                Dict{String, Any}("outcome" => "ok", "gate_level" => 0, "fit_rms" => 0.04,
                    "gate_reasons" => "fit_rms", "identify" => "production"))
            # per-flat traces (as make_traces_from_flats.jl writes them)
            tp = synth_trace(R)
            tp[:, 151:180, 1] .= 0.0
            fdir = joinpath(outdir, "quartz_flats", string(mjd))
            mkpath(fdir)
            for e in ("0004", "0009")
                h5open(joinpath(fdir, "quartzTrace_$(tele)_$(mjd)_$(e)_$(chip).h5"), "w") do f
                    f["trace_params"] = tp
                    f["regularized_trace_params"] = tp
                    g = create_group(f, "metadata")
                    g["git_commit"] = "test"
                end
            end
            fb = joinpath(outdir, "fallback.h5")
            h5open(fb, "w") do f
                f["trace_params"] = synth_trace(R; shift = 0.6)
                f["regularized_trace_params"] = synth_trace(R; shift = 0.6)
                create_group(f, "metadata")["git_commit"] = "fallback"
            end
            tm = joinpath(outdir, "apred", string(mjd), "traceMain_$(tele)_$(mjd)_$(chip).h5")
            # default rule: last quartz by name (0009), which failed the gate -> flagged
            lv = make_traceMain(outdir, tele, mjd, chip; fallback_file = fb, checkpoint_mode = "clobber")
            @test lv == TRACE_LEVEL_OWN_FLAT
            md = read_metadata(tm)
            @test endswith(md["trace_orig_param_fname"], "_0009_B.h5")
            @test md["trace_type"] == "quartz" && md["trace_level"] == TRACE_LEVEL_OWN_FLAT
            @test md["trace_select_rule"] == "last" && md["trace_gate_level"] == 0
            @test (md["trace_flags"] & TRACE_UNGATED) != 0
            bits = h5open(f -> read(f["bitmsk_trace"]), tm)
            @test count(x -> (x & TRACE_FIBER_UNLIT) != 0, bits) == 30
            @test h5open(f -> read(f["trace_params"]), tm) == tp
            # gated rule: 0004
            make_traceMain(outdir, tele, mjd, chip; select = "gated", fallback_file = fb,
                checkpoint_mode = "clobber")
            md = read_metadata(tm)
            @test endswith(md["trace_orig_param_fname"], "_0004_B.h5")
            @test md["trace_gate_level"] == 2 && md["trace_flags"] == 0
            @test md["trace_cand_gate_level"] == [2, 0]
            # no eligible flat: fixed file (flagged) or nothing
            other = joinpath(outdir, "apred", "57970")
            lv = make_traceMain(outdir, tele, 57970, chip; select = "gated", fallback_file = fb,
                checkpoint_mode = "clobber")
            @test lv == TRACE_LEVEL_FIXED_FALLBACK
            tm2 = joinpath(other, "traceMain_$(tele)_57970_$(chip).h5")
            md2 = read_metadata(tm2)
            @test md2["trace_type"] == "quartz_fallback" && md2["trace_found_match"] == false
            @test all(h5open(f -> read(f["bitmsk_trace"]), tm2) .== TRACE_FIXED_FALLBACK)
            lv = make_traceMain(outdir, tele, 57970, chip; select = "gated", fallback = "none",
                fallback_file = fb, checkpoint_mode = "clobber")
            @test lv == TRACE_LEVEL_NONE
            @test !isfile(tm2)
            @test_throws ErrorException make_traceMain(outdir, tele, 57970, chip;
                fallback_file = joinpath(outdir, "missing.h5"), checkpoint_mode = "clobber")
            @test_throws ErrorException make_traceMain(outdir, tele, mjd, chip; select = "best")
        end
    end
end
