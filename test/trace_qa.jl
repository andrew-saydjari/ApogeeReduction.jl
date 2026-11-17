using HDF5, Random, Statistics
using ApogeeReduction: trace_peak_census, trace_identify_slit, fill_unlit_fibers_robust!,
                       trace_params_qa, trace_id_match, trace_gate, select_night_trace,
                       TraceCandidate, trace_fiber_bits, trace_bits_from_metadata,
                       write_trace_qa, read_trace_qa, trace_qa_path, make_traceMain,
                       trace_reference_positions, trace_level_bits, trace_bit_names,
                       night_science_mid_exposure, read_metadata,
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
            "nlit" => 300, "fit_rms" => 0.005, "spacing_rms" => 0.004, "spacing_nbad" => 0,
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
