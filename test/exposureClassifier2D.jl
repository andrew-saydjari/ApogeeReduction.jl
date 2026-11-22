# Tests of the 2D content forest's features, the approved commanded-object rule,
# and the translated v6 calibration-masking policy.
#
# Real-data tests (prototype feature reproduction, full-archive rule
# reproduction, forest scoring) run only where the data are reachable and are
# skipped with an @info otherwise:
#   AR_TEST_PASS1_APRED          pass-1 apred directory with ar2D files
#   AR_TEST_CLASSIFIER2D_PROTO   prototype folder (tables/oof_full.csv,
#                                tuning_object/tables/acted/acted_recommended.csv)
#   AR_TEST_CLASSIFIER2D_MODEL   a classifier_2d artifact (JLD2)
using ApogeeReduction: object_frame_rule, exposure_predicted_bad, exposure_check_category,
                       exposure_class_metadata, exposure_class_unknown_metadata,
                       read_exposure_type_check, exposure_class_verdict,
                       classifier_content_class, expected_content_class,
                       classifier2d_feature_names, classifier2d_chip_products,
                       read_classifier2d_chip_products, exposure_class_features_2d,
                       load_classifier2d_references, load_exposure_classifier_2d,
                       load_exposure_classifier, classify_exposure_2d, classifier2d_era,
                       safe_jldsave, EXPFLAG_PREDICTED_BAD, EXPFLAG_NOTRUN,
                       OBJECT_RULE_P_NONSKY_MIN, OBJECT_RULE_OH_REL_MAX,
                       OBJECT_RULE_LINE_Z_MAX, OBJECT_RULE_TR_DYN_GUARD,
                       CLASSIFIER2D_CLASSES
using HDF5, Printf
const _CSV = ApogeeReduction.CSV
const _DF = ApogeeReduction.DataFrame

const C2D_DATA = joinpath(@__DIR__, "data", "classifier_2d")
const C2D_REFDIR = joinpath(dirname(@__DIR__), "data", "exposure_classifier_2d")
const C2D_TEMPLATES = joinpath(C2D_REFDIR, "line_templates.h5")
const C2D_LIGHTREF = joinpath(C2D_REFDIR, "light_reference.csv")
const C2D_APRED = get(ENV, "AR_TEST_PASS1_APRED",
    "/mnt/ceph/users/sdssv/work/asaydjari/pass1_dr21/apred")
const C2D_PROTO = get(ENV, "AR_TEST_CLASSIFIER2D_PROTO",
    "/mnt/ceph/users/sdssv/work/asaydjari/2026_11_19/classifier_2d")
const C2D_MODEL = get(ENV, "AR_TEST_CLASSIFIER2D_MODEL",
    "/mnt/ceph/users/sdssv/work/asaydjari/cal_ref/exposure_classifier/exposure_classifier_2d_v1.jld2")

# equal to numerical precision: identical, both NaN, or within 1e-12 relative.
# (Bitwise identity does not hold across CPU types: the BLAS dot / SIMD sums in
# the profile means differ by an ulp between the cluster nodes that made the
# prototype tables and a workstation; measured max relative difference 1.4e-14.)
c2d_close(a, b) = (isnan(a) && isnan(b)) || a == b || isapprox(a, b; rtol = 1e-12, atol = 1e-13)

@testset "exposureClassifier2D" begin
    @testset "approved thresholds" begin
        # AKS 2026-11-21/22; changing one of these changes which science is dropped
        @test OBJECT_RULE_P_NONSKY_MIN == 0.8
        @test OBJECT_RULE_OH_REL_MAX == -0.6
        @test OBJECT_RULE_LINE_Z_MAX == 2.2
        @test OBJECT_RULE_TR_DYN_GUARD == 0.9
    end

    @testset "object_frame_rule on synthetic inputs" begin
        r = object_frame_rule(0.95, 0.0, 3.0, 0.5)
        @test r.act && r.reason == "forest" && r.forest && !r.no_sky_light && !r.guard
        r = object_frame_rule(0.1, -1.0, 1.5, 0.5)
        @test r.act && r.reason == "no_sky_light"
        r = object_frame_rule(0.9, -1.0, 1.5, 0.5)
        @test r.act && r.reason == "forest+no_sky_light"
        r = object_frame_rule(0.99, -1.0, 1.5, 1.2)
        @test !r.act && r.reason == "guard_exempt" && r.guard
        r = object_frame_rule(0.1, 0.0, 3.0, 0.5)
        @test !r.act && r.reason == "none" && r.inputs_finite
        # no-sky-light needs BOTH the level and the detection test
        @test !object_frame_rule(0.1, -1.0, 2.5, 0.5).act
        @test !object_frame_rule(0.1, -0.5, 1.5, 0.5).act
        # boundaries: >= for P and the guard, strict < for oh_rel and line_z
        @test object_frame_rule(0.8, 0.0, 3.0, 0.5).act
        @test !object_frame_rule(prevfloat(0.8), 0.0, 3.0, 0.5).act
        @test !object_frame_rule(0.1, -0.6, 1.0, 0.5).act
        @test object_frame_rule(0.1, prevfloat(-0.6), 1.0, 0.5).act
        @test !object_frame_rule(0.1, -1.0, 2.2, 0.5).act
        @test !object_frame_rule(0.9, 0.0, 3.0, 0.9).act
        @test object_frame_rule(0.9, 0.0, 3.0, prevfloat(0.9)).act
        # NaN / missing / nothing: that test does not fire
        r = object_frame_rule(NaN, -1.0, 1.5, 0.5)
        @test r.act && r.reason == "no_sky_light" && !r.inputs_finite
        @test !object_frame_rule(NaN, NaN, NaN, NaN).act
        @test !object_frame_rule(missing, missing, missing, missing).act
        @test !object_frame_rule(nothing, 0.0, 3.0, 0.0).act
        @test !object_frame_rule(0.1, NaN, 1.5, 0.5).act
        @test !object_frame_rule(0.1, -1.0, NaN, 0.5).act
        # NaN tr_dyn (<= 50 trace peaks, e.g. an unlit frame) does NOT exempt
        r = object_frame_rule(0.99, -2.0, 0.5, NaN)
        @test r.act && !r.guard && !r.inputs_finite
        # thresholds are keyword-overridable (for scans), defaults are the approved ones
        @test !object_frame_rule(0.85, 0.0, 3.0, 0.5; p_min = 0.9).act
    end

    @testset "object rule reproduces the approved archive list (compact fixture)" begin
        T = _CSV.read(joinpath(C2D_DATA, "object_rule_fixture.csv"), _DF;
            types = Dict(:tele => String))
        act = [object_frame_rule(1 - r.p_sky, r.oh_rel, r.line_z, r.tr_dyn).act for r in eachrow(T)]
        @test sum(T.acted_reference) == 437
        @test act == T.acted_reference
        @test any(isnan, T.tr_dyn)   # the fixture exercises the NaN path
    end

    @testset "object rule reproduces acted_recommended.csv (full archive)" begin
        oof = joinpath(C2D_PROTO, "tables", "oof_full.csv")
        ref = joinpath(C2D_PROTO, "tuning_object", "tables", "acted", "acted_recommended.csv")
        if isfile(oof) && isfile(ref)
            O = _CSV.read(oof, _DF; select = [:tele, :mjd, :expnum, :image_type, :p_full_sky,
                    :oh_rel, :n_G_line_z, :n_G_tr_dyn], types = Dict(:tele => String))
            O = O[O.image_type .== "object", :]
            a = Set((r.tele, r.mjd, r.expnum) for r in eachrow(O)
            if object_frame_rule(1 - r.p_full_sky, r.oh_rel, r.n_G_line_z, r.n_G_tr_dyn).act)
            R = _CSV.read(ref, _DF; select = [:tele, :mjd, :expnum], types = Dict(:tele => String))
            @test size(O, 1) == 153869
            @test length(a) == 437
            @test a == Set(zip(R.tele, R.mjd, R.expnum))
        else
            @info "skipping full-archive rule reproduction: prototype tables not found under $C2D_PROTO"
            @test_skip false
        end
    end

    @testset "translated calibration policy (v6 intent, 2D classes)" begin
        @test expected_content_class("arclamp_q0t1u0") == "thar"
        @test expected_content_class("internalflat_q0t0u0") == "ledflat"
        @test expected_content_class("arclamp_q1t1u1") === nothing
        @test classifier_content_class("twilightflat_faint") == "sky"
        @test classifier_content_class("object_q0t0u0") == "sky"
        @test classifier_content_class("une") == "une"
        res(pred, prob) = (pred = pred, prob = prob, decision = prob < 0.5 ? "unknown" : "ok")
        @test exposure_check_category("arclamp_q0t1u0", res("thar", 0.99)) == "ok"
        @test exposure_check_category("arclamp_q0t1u0", res("une", 0.99)) == "mislabel_candidate"
        @test exposure_check_category("internalflat_q0t0u0", res("dark", 0.9)) == "lamp_off_candidate"
        @test exposure_check_category("twilightflat_q0t0u0", res("sky", 0.9)) == "faint_twilight"
        @test exposure_check_category("twilightflat_q0t0u0", res("twilight", 0.9)) == "ok"
        @test exposure_check_category("dark_q0t0u0", res("dark_persist", 0.9)) == "mislabel_candidate"
        @test exposure_check_category("quartzflat_q1t0u0", res("quartzflat", 0.4)) == "unknown"
        @test exposure_check_category("arclamp_q1t1u1", res("thar", 0.9)) == "mislabel_candidate"
        # masking
        @test exposure_predicted_bad("dark_q0t0u0", "dark_persist", "ok")   # below flag_tau still masked
        @test !exposure_predicted_bad("dark_q0t0u0", "dark", "persistence_prior")
        @test exposure_predicted_bad("dark_q0t0u0", "dark", "unknown")
        @test exposure_predicted_bad("arclamp_q0t0u0", "thar", "ok")
        @test !exposure_predicted_bad("internalflat_q0t0u0", "ledflat", "ok")
        @test exposure_predicted_bad("quartzflat_q1t0u0", "domeflat", "mislabel_candidate")
        @test !exposure_predicted_bad("quartzflat_q1t0u0", "quartzflat", "ok")
        # legacy v6 verdict strings read the same way as before
        @test !exposure_predicted_bad("dark_q0t0u0", "dark_q0t0u0", "ok")
        @test exposure_predicted_bad("dark_q0t0u0", "dark_persist", "mislabel_candidate")
        @test !exposure_predicted_bad("internalflat_q0t0u0", "internalflat_q0t0u0", "ok")
    end

    @testset "commanded objects: EXPFLAG_PREDICTED_BAD follows the rule" begin
        bad = exposure_predicted_bad("object_q0t0u0", "fpi", "mislabel_candidate";
            p_sky = 0.05, oh_rel = 0.0, line_z = 3.0, tr_dyn = 0.5)
        @test bad
        @test !exposure_predicted_bad("object_q0t0u0", "fpi", "mislabel_candidate";
            p_sky = 0.05, oh_rel = 0.0, line_z = 3.0, tr_dyn = 1.1)          # guard
        @test exposure_predicted_bad("object_q0t0u0", "sky", "ok";
            p_sky = 0.9, oh_rel = -1.2, line_z = 1.0, tr_dyn = 0.4)          # no sky light
        @test !exposure_predicted_bad("object_q0t0u0", "sky", "ok";
            p_sky = 0.9, oh_rel = 0.0, line_z = 3.0, tr_dyn = 0.4)
        # anomalous lamp flags are still commanded objects
        @test exposure_predicted_bad("object_q0t0u1", "sky", "ok";
            p_sky = 0.9, oh_rel = -1.2, line_z = 1.0, tr_dyn = 0.4)
        # no evidence (e.g. a v6-era table) cannot fire the rule
        @test !exposure_predicted_bad("object_q0t0u0", "dark_q0t0u0", "mislabel_candidate")
        @test !exposure_predicted_bad("object_q0t0u0", "dark", "checkfail"; p_sky = 0.0)

        md = exposure_class_metadata("object_q0t0u0", "dark", 0.97, "mislabel_candidate";
            p_sky = 0.01, oh_rel = -1.3, line_z = 1.2, tr_dyn = 0.6)
        @test (md["exposure_flags"] & EXPFLAG_PREDICTED_BAD) != 0x00
        @test exposure_class_verdict(md) == :bad
        @test md["exp_class_object_rule"] == "forest+no_sky_light"
        @test md["exp_class_p_sky"] == 0.01 && md["exp_class_oh_rel"] == -1.3
        @test md["exp_class_line_z"] == 1.2 && md["exp_class_tr_dyn"] == 0.6
        @test md["exp_class_rule_forest"] && md["exp_class_rule_no_sky_light"]
        @test !md["exp_class_rule_guard"]
        g = exposure_class_metadata("object_q0t0u0", "sky", 0.99, "ok";
            p_sky = 0.99, oh_rel = 0.1, line_z = 3.5, tr_dyn = 1.3)
        @test g["exposure_flags"] == 0x00 && g["exp_class_object_rule"] == "none"
        c = exposure_class_metadata("quartzflat_q1t0u0", "quartzflat", 0.99, "ok"; p_sky = 0.0)
        @test c["exp_class_object_rule"] == "n/a" && !c["exp_class_rule_forest"]
        u = exposure_class_unknown_metadata()
        @test (u["exposure_flags"] & EXPFLAG_NOTRUN) != 0x00
        @test u["exp_class_object_rule"] == "n/a" && isnan(u["exp_class_p_sky"])
        # every block carries the same provenance keys
        @test issubset(["exp_class_p_sky", "exp_class_oh_rel", "exp_class_line_z",
                "exp_class_tr_dyn", "exp_class_object_rule"], keys(u))
    end

    @testset "per-MJD table round trip with and without rule columns" begin
        mktempdir() do dir
            p = joinpath(dir, "t_new.h5")
            safe_jldsave(p, Dict{String, Any}(); tele = ["apo", "apo"], mjd = [60000, 60000],
                expnum = [5, 6], labeled = ["object_q0t0u0", "object_q0t0u0"],
                pred = ["sky", "sky"], prob = [0.9, 0.95], flag = ["ok", "ok"],
                p_sky = [0.9, 0.95], oh_rel = [-1.0, 0.0], line_z = [1.5, 3.0],
                tr_dyn = [0.5, 0.5], object_rule = ["no_sky_light", "none"],
                exposure_flags = UInt8[0x01, 0x00])
            t = read_exposure_type_check(p)
            @test exposure_class_verdict(t[5]) == :bad
            @test exposure_class_verdict(t[6]) == :fine
            @test t[5]["exp_class_object_rule"] == "no_sky_light"
            # a v6-era table: no rule columns -> objects are never bad, cals unchanged
            p2 = joinpath(dir, "t_old.h5")
            safe_jldsave(p2, Dict{String, Any}(); tele = ["apo", "apo"], mjd = [60000, 60000],
                expnum = [1, 2], labeled = ["object_q0t0u0", "dark_q0t0u0"],
                pred = ["dark_q0t0u0", "dark_persist"], prob = [0.99, 0.9],
                flag = ["mislabel_candidate", "mislabel_candidate"],
                exposure_flags = UInt8[0x00, 0x01])
            t2 = read_exposure_type_check(p2)
            @test exposure_class_verdict(t2[1]) == :fine
            @test isnan(t2[1]["exp_class_p_sky"])
            @test exposure_class_verdict(t2[2]) == :bad
        end
    end

    @testset "missing references / artifact: check disabled, logged" begin
        @test (@test_logs (:warn, r"DISABLED") match_mode = :any load_classifier2d_references(
            "/definitely/not/templates.h5", C2D_LIGHTREF)) === nothing
        @test (@test_logs (:warn, r"DISABLED") match_mode = :any load_classifier2d_references(
            C2D_TEMPLATES, "/definitely/not/light.csv")) === nothing
        @test (@test_logs (:warn, r"DISABLED") match_mode = :any load_exposure_classifier_2d(
            "/definitely/not/model.jld2", C2D_TEMPLATES, C2D_LIGHTREF)) === nothing
        @test (@test_logs (:warn, r"DISABLED") match_mode = :any load_exposure_classifier(
            "/definitely/not/model.jld2", C2D_TEMPLATES, C2D_LIGHTREF)) === nothing
        # the shipped references load
        refs = load_classifier2d_references(C2D_TEMPLATES, C2D_LIGHTREF)
        @test refs !== nothing
        @test length(refs.templates) == 24 && all(length(v) == 2048 for v in values(refs.templates))
        @test sort(collect(keys(refs.light))) == [("apo", 1), ("apo", 2), ("lco", 1), ("lco", 2)]
        # a v6 artifact is rejected loudly
        mktempdir() do dir
            p = joinpath(dir, "v6like.jld2")
            ApogeeReduction.JLD2.jldsave(p; model = 1, classes = ["object_q0t0u0"],
                unknown_tau = 0.5, flag_tau = 0.7)
            @test_throws ErrorException load_exposure_classifier(p, C2D_TEMPLATES, C2D_LIGHTREF)
        end
    end

    @testset "features on synthetic images" begin
        names_ = classifier2d_feature_names()
        @test length(names_) == 175 && allunique(names_)
        @test classifier2d_era("apo", 59422) == 1 && classifier2d_era("apo", 59423) == 2
        @test classifier2d_era("lco", 59808) == 1 && classifier2d_era("lco", 59809) == 2
        rng = MersenneTwister(7)
        img = randn(rng, 2560, 2048) .* 5
        for j in 3:7:2048            # trace-like columns
            img[1:2048, j] .+= 200
        end
        img[100:110, 500] .= NaN
        lu = fill(30, 2560, 2048)
        chi = fill(30.0, 2560, 2048)
        cr = zeros(Int, 2560, 2048)
        before = copy(img)
        prods = Dict("R" => classifier2d_chip_products(img), "B" => classifier2d_chip_products(img),
            "G" => classifier2d_chip_products(img; last_unsaturated = lu, chisqimage = chi,
                CRimage = cr, ndiff_used = 30, quarters = true))
        @test isequal(img, before)   # inputs are not mutated
        @test eltype(prods["G"].ps) == Float32
        F = exposure_class_features_2d(prods, "apo", 60000, 31)
        @test Set(keys(F)) == Set(names_)
        @test all(isnan, (F["oh_rel"], F["t_oh_G"]))   # no references -> NaN
        @test F["n_read"] == 31 && F["n_ndiff"] == 30 && F["n_sat_frac"] == 0.0
        @test F["d_G_nanfrac"] > 0
        refs = load_classifier2d_references(C2D_TEMPLATES, C2D_LIGHTREF)
        F2 = exposure_class_features_2d(prods, "apo", 60000, 31; refs = refs)
        @test isfinite(F2["oh_rel"]) && isfinite(F2["t_thar_R"])
        # a blank (all-zero) frame: finite deployed block, no traces -> NaN tr_dyn
        z = zeros(2560, 2048)
        pz = Dict("R" => classifier2d_chip_products(z), "B" => classifier2d_chip_products(z),
            "G" => classifier2d_chip_products(z; quarters = true))
        Fz = exposure_class_features_2d(pz, "apo", 59226, 46; refs = refs)
        @test isnan(Fz["n_G_tr_dyn"]) && Fz["n_G_line_z"] == -3.0
        @test all(isfinite(Fz[k]) for k in names_ if startswith(k, "d_"))
        @test isnan(Fz["n_sat_frac"])   # no ramp arrays given
    end

    @testset "features reproduce the prototype on real exposures" begin
        X = _CSV.read(joinpath(C2D_DATA, "prototype_features.csv"), _DF;
            types = Dict(:tele => String, :image_type => String))
        have = isdir(C2D_APRED) && all(eachrow(X)) do r
            all(isfile(joinpath(C2D_APRED, string(r.mjd), @sprintf("ar2D_%s_%d_%04d_%s_%s.h5",
                r.tele, r.mjd, r.expnum, c, r.image_type))) for c in CHIP_LIST)
        end
        if !have
            @info "skipping prototype feature reproduction: ar2D files not found under $C2D_APRED"
            @test_skip false
        else
            refs = load_classifier2d_references(C2D_TEMPLATES, C2D_LIGHTREF)
            names_ = classifier2d_feature_names()
            clf = isfile(C2D_MODEL) ? load_exposure_classifier_2d(C2D_MODEL, C2D_TEMPLATES, C2D_LIGHTREF) :
                  nothing
            clf === nothing && @info "no classifier_2d artifact at $C2D_MODEL; scoring check skipped"
            @test size(X, 1) == 33
            nbad = 0
            for r in eachrow(X)
                fn = Dict(c => joinpath(C2D_APRED, string(r.mjd), @sprintf(
                    "ar2D_%s_%d_%04d_%s_%s.h5", r.tele, r.mjd, r.expnum, c, r.image_type))
                for c in CHIP_LIST)
                prods = Dict(c => read_classifier2d_chip_products(fn[c], c) for c in CHIP_LIST)
                @test prods["G"].nread_total == r.n_read
                F = exposure_class_features_2d(prods, r.tele, r.mjd, prods["G"].nread_total;
                    refs = refs)
                ok = all(c2d_close(F[k], Float64(r[k])) for k in names_)
                ok || (nbad += 1; @warn "prototype feature mismatch" r.tele r.mjd r.expnum)
                @test ok
                if clf !== nothing
                    Fp = Dict(k => Float64(r[k]) for k in names_)
                    @test classify_exposure_2d(clf, F).proba == classify_exposure_2d(clf, Fp).proba
                end
            end
            @test nbad == 0
        end
    end
end
