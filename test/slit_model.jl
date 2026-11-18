using HDF5, Random, Statistics
using ApogeeReduction: TraceEpoch, TraceTemplate, TemplateRegistration, slit_profile, slit_bands,
                       register_frame, registration_gate, registered_centers, registered_trace_params,
                       trace_extract_template, template_fiber_delta, epoch_median_registration,
                       parse_trace_fallback, trace_fallback_uses_templates, trace_fiber_bits,
                       make_traceMain, read_metadata, write_trace_template_cal, load_trace_template_cal,
                       night_calibration_frames, night_science_frames, stack_bands, template_qa_fields,
                       trace_reference_positions, write_trace_qa, trace_qa_path, flat_trace_qa,
                       TRACE_REGISTERED_CAL, TRACE_REGISTERED_SCI, TRACE_TEMPLATE_ONLY,
                       TRACE_FIBER_UNLIT, TRACE_FIBER_DELTA, TRACE_EPOCH_PROVISIONAL,
                       TRACE_LEVEL_OWN_FLAT, TRACE_LEVEL_REGISTERED_CAL, TRACE_LEVEL_REGISTERED_SCI,
                       TRACE_LEVEL_TEMPLATE, TRACE_LEVEL_NONE, TRACE_PER_FIBER_BITS, TraceBits

const SM_DATA = joinpath(@__DIR__, "..", "data")

# a template following the reference positions R (fiber 1..300 near the chip center), with
# curvature in x, given relative throughputs and epoch-dark fibers
function sm_template(R; tele = "apo", chip = "G", epoch = TraceEpoch(tele, "apoE9", 60000, 99999, true),
        href = ones(length(R)), dark = Int[], regd = 0.0)
    nf = length(R)
    cp = zeros(5, nf)
    cp[1, :] .= R
    cp[2, :] .= 0.8 .+ 0.0004 .* (R .- 1024)
    cp[3, :] .= -1.5
    wp = zeros(3, nf)
    wp[1, :] .= 1.0
    wp[2, :] .= 0.05
    xc = collect(32.0:64.0:2016.0)
    u = (xc .- 1024) ./ 1024
    center = [sum(cp[p, f] * ui^(p - 1) for p in 1:5) for ui in u, f in 1:nf]
    width = [sum(wp[p, f] * ui^(p - 1) for p in 1:3) for ui in u, f in 1:nf]
    lit = trues(nf)
    lit[dark] .= false
    TraceTemplate(tele, chip, epoch, xc, center, width, cp, wp, href, lit, [60000, 60010, 60020],
        [1, 2, 3], fill(regd, 3), zeros(3), Dict("builder" => "test"))
end

# trace parameters (nx, 300, 3) of the template registered by `reg`, columns 1:nx
function sm_truth(t, reg; nx = 2048, heights = fill(1000.0, size(t.center_poly, 2)), fiber_offset = nothing)
    tp = registered_trace_params(t, reg; nx = nx, heights = heights)
    isnothing(fiber_offset) || (tp[:, :, 2] .+= reshape(fiber_offset, 1, :))
    return tp
end

# image[x, y] with fibers along tp drawn with the production profiles (so the truth is exactly
# representable); `cols` = lit columns (arc-like frames)
const SM_P = Ref{Any}(nothing)
function sm_image(tp; noise = 3.0, rng = MersenneTwister(5), cols = axes(tp, 1))
    P = SM_P[]
    nx, nf, _ = size(tp)
    img = noise .* randn(rng, nx, 2048)
    for x in cols
        A, y0 = ApogeeReduction.slit_unit_profiles(P, tp[x, :, 2], tp[x, :, 3], 1:nf; half = 7)
        for f in 1:nf
            tp[x, f, 1] == 0 && continue
            for i in axes(A, 2)
                y = y0[f] + i - 1
                1 <= y <= 2048 && (img[x, y] += tp[x, f, 1] * A[f, i])
            end
        end
    end
    return img, fill(1 / noise^2, nx, 2048), zeros(Int, nx, 2048)
end

function sm_write_2d(path, img, iv, m)
    mkpath(dirname(path))
    h5open(path, "w") do f
        f["dimage"] = img
        f["ivarimage"] = iv
        f["pix_bitmask"] = m
    end
end

@testset "slit model (epoch templates)" begin
    R = trace_reference_positions("apo", "G"; data_path = SM_DATA)
    P = slit_profile("apo", "G"; data_path = SM_DATA)
    SM_P[] = P
    rng = MersenneTwister(3)
    href = 0.6 .+ 0.8 .* rand(rng, 300)        # a distinctive throughput fingerprint

    @testset "parse --trace_fallback" begin
        @test parse_trace_fallback("none") == (Int[], false)
        @test parse_trace_fallback("fixed_file") == (Int[], true)
        @test parse_trace_fallback("template") == ([1, 2, 3], false)
        @test parse_trace_fallback("T1,T2") == ([1, 2], false)
        @test parse_trace_fallback("T3") == ([3], false)
        @test_throws ArgumentError parse_trace_fallback("T2,T1")
        @test_throws ArgumentError parse_trace_fallback("nearest_night")
        @test trace_fallback_uses_templates("template") && !trace_fallback_uses_templates("none")
    end

    @testset "registration recovers shift, rotation, bow and scale" begin
        t = sm_template(R; href = href)
        truth = TemplateRegistration([2.3, 0.15, -0.1], [0.0006, 0.0])
        tp = sm_truth(t, truth; heights = 1000.0 .* href)
        img, iv, m = sm_image(tp)
        xc, prof, err = slit_bands(img, 1 ./ sqrt.(iv), trues(size(img)))
        rf = register_frame(P, t, xc, prof, err; mode = :full)
        @test rf.outcome == "ok" && rf.kbest == 0
        C = registered_centers(t, rf.reg, 1:64:2048)
        Ct = registered_centers(t, truth, 1:64:2048)
        @test maximum(abs.(C .- Ct)) < 0.05
        @test rf.fit_rms < 0.02 && rf.nlit == 300
        ok, reasons = registration_gate(rf)
        @test ok && isempty(reasons)
        # registration-only searches only +-3 px around the epoch-median shift
        rr = register_frame(P, t, xc, prof, err; mode = :register_only)
        @test rr.outcome == "ok" && maximum(abs.(registered_centers(t, rr.reg, 1:64:2048) .- Ct)) < 0.05
    end

    @testset "labels: a registration on the wrong alias is refused, never relabelled" begin
        t = sm_template(R; href = href)
        # the slit moved up by one fiber spacing; a search window too narrow to reach the true
        # shift (+-3 px) lands on the alias, where the throughput fingerprint says k = -1
        sp = median(diff(R))
        tp = sm_truth(t, TemplateRegistration([sp, 0, 0], [0.0, 0]); heights = 1000.0 .* href)
        img, iv, m = sm_image(tp)
        xc, prof, err = slit_bands(img, 1 ./ sqrt.(iv), trues(size(img)))
        rf = register_frame(P, t, xc, prof, err; mode = :full, dmax = 3.0)
        @test rf.kbest != 0
        @test rf.outcome == "id_template_k"
        tpx, _, info = trace_extract_template(img, iv, t, P; dmax = 3.0)
        @test isempty(tpx) && info.outcome == "id_template_k"
        # with the default window the shift is found and the labels are right
        rf = register_frame(P, t, xc, prof, err; mode = :full)
        @test rf.outcome == "ok" && rf.kbest == 0 && abs(rf.reg.d[1] - sp) < 0.05
    end

    @testset "T0: template-seeded per-column fit, unlit bundle, per-fiber delta" begin
        nx = 512                                    # a 512-column chip keeps the test fast
        t = sm_template(R; href = href)
        reg = TemplateRegistration([-1.7, 0.05, 0.0], [0.0003, 0.0])
        off = zeros(300)
        off[40] = 0.5                               # one fiber departs from the template
        h = 1000.0 .* href
        h[151:180] .= 0                             # a dark bundle on this flat
        tp = sm_truth(t, reg; nx = nx, heights = h, fiber_offset = off)
        img, iv, m = sm_image(tp)
        good = trues(size(img))
        good[1:40, 1:400] .= false               # masked detector area (corner)
        tpx, cv, info = trace_extract_template(img, iv, t, P; good_pixels = good)
        @test info.outcome == "ok"
        mk = info.masked
        @test size(mk) == (nx, 300) && any(mk[1:40, :]) && !any(mk[41:end, :])
        Rs = registered_centers(t, info.reg, 1:nx)
        f1 = findfirst(f -> mk[10, f], 1:300)
        @test abs(tpx[10, f1, 2] - Rs[10, f1]) < 0.05 && all(isnan, cv[10, f1, :, :])
        @test size(tpx) == (nx, 300, 3) && size(cv) == (nx, 300, 3, 3)
        lit = trues(300)
        lit[151:180] .= false
        @test info.lit == lit
        D = tpx[41:end, lit, 2] .- tp[41:end, lit, 2]
        @test median(abs.(D)) < 0.02
        @test quantile(vec(abs.(D)), 0.99) < 0.1
        # unlit fibers sit exactly at the registered template, height 0
        Rt = registered_centers(t, info.reg, 1:nx)
        @test maximum(abs.(tpx[:, 151:180, 2] .- Rt[:, 151:180])) < 1e-9
        @test all(tpx[:, 151:180, 1] .== 0)
        @test info.delta[40] > 0.3
        @test count(v -> isfinite(v) && v > 0.3, info.delta) == 1
        @test all(isnan, info.delta[151:180])
        # the same fit seeded at the bare template agrees on a flat this clean
        tpy, _, _ = trace_extract_template(img, iv, t, P; seed = :template)
        @test median(abs.(tpy[:, lit, 2] .- tpx[:, lit, 2])) < 0.01
        xq = ApogeeReduction.template_free_crosscheck(img, iv, tpx, P, R; lit = info.lit)
        @test xq["tmpl_xcheck_id_reason"] == "ok" && xq["tmpl_xcheck_label_agree"] == 1.0
        @test xq["tmpl_xcheck_med_abs_px"] < 0.01 && !xq["tmpl_xcheck_flag"]
        # a trace one fiber off in labels is flagged
        bad = copy(tpx); bad[:, 1:299, 2] .= tpx[:, 2:300, 2]
        @test ApogeeReduction.template_free_crosscheck(img, iv, bad, P, R; lit = info.lit)["tmpl_xcheck_flag"]
        q = template_qa_fields(info; tmpl = t)
        @test q["tmpl_n_fiber_delta"] == 1 && q["tmpl_outcome"] == "ok" && length(q["tmpl_reg"]) == 5
        bits = trace_fiber_bits(tpx, TRACE_LEVEL_OWN_FLAT; fiber_delta = info.delta, provisional = true)
        @test bits[40] & TRACE_FIBER_DELTA != 0 && count(b -> b & TRACE_FIBER_DELTA != 0, bits) == 1
        @test all(b -> b & TRACE_EPOCH_PROVISIONAL != 0, bits)
        @test all(b -> b & TRACE_FIBER_UNLIT != 0, bits[151:180])
        @test TRACE_PER_FIBER_BITS == TRACE_FIBER_UNLIT | TRACE_FIBER_DELTA
    end

    @testset "registration-only on an arc-like frame and on stacked science frames" begin
        t = sm_template(R; href = href)
        truth = TemplateRegistration([1.2, 0.0, 0.0], [0.0, 0.0])
        tp = sm_truth(t, truth; heights = 3000.0 .* href)
        # light only in a few columns per band (lines)
        img, iv, m = sm_image(tp; cols = vcat([c:(c + 2) for c in 10:40:2040]...))
        xc, prof, err = slit_bands(img, 1 ./ sqrt.(iv), trues(size(img)); half = 32, nsigma = Inf)
        rf = register_frame(P, t, xc, prof, err; mode = :register_only)
        @test rf.outcome == "ok"
        @test maximum(abs.(registered_centers(t, rf.reg, 1:64:2048) .- registered_centers(t, truth, 1:64:2048))) < 0.1
        # stack: median of band profiles of three noisy frames
        bs = map(1:3) do k
            im, ivk, _ = sm_image(sm_truth(t, truth; heights = 100.0 .* href); rng = MersenneTwister(k))
            slit_bands(im, 1 ./ sqrt.(ivk), trues(2048, 2048); half = 32)
        end
        p, e = stack_bands([b[2] for b in bs], [b[3] for b in bs])
        @test size(p) == (32, 2048) && all(isfinite, e[:, 500])
        rs = register_frame(P, t, bs[1][1], p, e; mode = :register_only)
        @test rs.outcome == "ok" && abs(rs.reg.d[1] - 1.2) < 0.1
        # pure noise does not pass the registration gate
        nimg = 3.0 .* randn(MersenneTwister(9), 2048, 2048)
        xn, pn, en = slit_bands(nimg, fill(3.0, 2048, 2048), trues(2048, 2048))
        rn = register_frame(P, t, xn, pn, en; mode = :register_only)
        @test !first(registration_gate(rn))
    end

    @testset "make_traceMain: own-night levels T1 > T2 > T3, never another night" begin
        mktempdir() do dir
            E = TraceEpoch
            eps = Dict("apo" => [E("apo", "apoE9", 60000, ApogeeReduction.TRACE_EPOCH_OPEN_END, true)])
            tm = [sm_template(R; tele = "apo", chip = c, epoch = eps["apo"][1], href = href,
                dark = [276], regd = 0.4) for c in ("R", "G", "B")]
            calp = joinpath(dir, "tt.h5")
            write_trace_template_cal(calp, eps, tm, Dict("code_commit" => "test");
                calibration_version = "test_v1")
            cal = load_trace_template_cal(calp)
            outdir = joinpath(dir, "out")
            t = tm[2]
            truth = TemplateRegistration([0.9, 0.0, 0.0], [0.0, 0.0])
            mjd = 60100
            # no flats at all: T3 (the template at the epoch-median registration)
            lv = make_traceMain(outdir, "apo", mjd, "G"; fallback = "template", template_cal = cal,
                data_path = SM_DATA)
            @test lv == TRACE_LEVEL_TEMPLATE
            tmf = joinpath(outdir, "apred", "$(mjd)", "traceMain_apo_$(mjd)_G.h5")
            md = read_metadata(tmf)
            @test md["trace_level"] == TRACE_LEVEL_TEMPLATE && md["trace_type"] == "template"
            @test md["trace_flags"] & TRACE_TEMPLATE_ONLY != 0
            @test md["trace_flags"] & TRACE_EPOCH_PROVISIONAL != 0     # after the newest input flat
            @test md["trace_template_id"] == ApogeeReduction.template_id(t)
            @test md["trace_reg"][1] ≈ 0.4
            bits = h5open(f -> read(f["bitmsk_trace"]), tmf, "r")
            @test eltype(bits) == TraceBits
            @test bits[276] & TRACE_FIBER_UNLIT != 0 && count(b -> b & TRACE_FIBER_UNLIT != 0, bits) == 1
            tpm = h5open(f -> read(f["regularized_trace_params"]), tmf, "r")
            @test maximum(abs.(tpm[1000:1050, :, 2] .- registered_centers(t, epoch_median_registration(t), 1000:1050))) < 1e-6
            # T1/T2 levels only: nothing serves the night -> no traceMain, not extracted
            @test make_traceMain(outdir, "apo", mjd, "G"; fallback = "T1,T2", template_cal = cal,
                data_path = SM_DATA) == TRACE_LEVEL_NONE
            @test !isfile(tmf)
            # science frames of the night: T2
            for (k, e) in enumerate((11, 12, 13))
                img, iv, m = sm_image(sm_truth(t, truth; heights = 200.0 .* href); rng = MersenneTwister(20 + k))
                sm_write_2d(joinpath(outdir, "apred", "$(mjd)", "ar2Dcal_apo_$(mjd)_$(lpad(e, 4, '0'))_G_object.h5"), img, iv, m)
            end
            @test length(night_science_frames(outdir, "apo", mjd, "G")) == 3
            lv = make_traceMain(outdir, "apo", mjd, "G"; fallback = "template", template_cal = cal,
                data_path = SM_DATA)
            @test lv == TRACE_LEVEL_REGISTERED_SCI
            md = read_metadata(tmf)
            @test md["trace_flags"] & TRACE_REGISTERED_SCI != 0 && md["trace_flags"] & TRACE_TEMPLATE_ONLY == 0
            @test abs(md["trace_reg"][1] - 0.9) < 0.1
            @test isfile(ApogeeReduction.regsci_qa_path(outdir, "apo", mjd, "G"))
            # a dome flat of the night (its trace not usable as T0): T1, preferred over T2
            img, iv, m = sm_image(sm_truth(t, truth; heights = 1000.0 .* href); rng = MersenneTwister(40))
            sm_write_2d(joinpath(outdir, "apred", "$(mjd)", "ar2Dcal_apo_$(mjd)_0005_G_domeflat.h5"), img, iv, m)
            # an arc without an almanac is skipped (its lamp is unknown)
            sm_write_2d(joinpath(outdir, "apred", "$(mjd)", "ar2Dcal_apo_$(mjd)_0006_G_arclamp.h5"), img, iv, m)
            fr = night_calibration_frames(outdir, "apo", mjd, "G")
            @test [f.kind for f in fr] == ["domeflat"]
            lv = make_traceMain(outdir, "apo", mjd, "G"; fallback = "template", template_cal = cal,
                data_path = SM_DATA)
            @test lv == TRACE_LEVEL_REGISTERED_CAL
            md = read_metadata(tmf)
            @test md["trace_flags"] & TRACE_REGISTERED_CAL != 0
            @test md["trace_type"] == "registered_cal" && md["trace_reg_source"] == "domeflat"
            @test abs(md["trace_reg"][1] - 0.9) < 0.05
            @test md["trace_fb_cand_kind"] == ["domeflat"] && md["trace_fb_cand_pass"] == [true]
            @test startswith(md["trace_orig_param_fname"], "T1:domeflat:")
            # checkpoint: the same decision is not redone (file untouched)
            make_traceMain(outdir, "apo", mjd, "G"; fallback = "template", template_cal = cal,
                data_path = SM_DATA, checkpoint_mode = "commit_exists")
            mt = mtime(tmf)
            sleep(1.1)
            make_traceMain(outdir, "apo", mjd, "G"; fallback = "template", template_cal = cal,
                data_path = SM_DATA, checkpoint_mode = "commit_exists")
            @test mtime(tmf) == mt
            # the template levels need the template file
            @test_throws ErrorException make_traceMain(outdir, "apo", mjd, "G"; fallback = "template",
                template_cal = nothing)
            # a night no epoch serves: own-night levels skipped, no trace
            @test make_traceMain(outdir, "apo", 59000, "G"; fallback = "template", template_cal = cal,
                data_path = SM_DATA) == TRACE_LEVEL_NONE
        end
    end
end
