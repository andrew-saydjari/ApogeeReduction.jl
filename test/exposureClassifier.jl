using ApogeeReduction: exposure_class_metadata, exposure_class_unknown_metadata,
                       exposure_predicted_bad, exposure_class_label,
                       exposure_type_check_path, read_exposure_type_check,
                       exposure_class_metadata_for, safe_jldsave,
                       EXP_CLASS_BAD_UNKNOWN, EXP_CLASS_BAD_FALSE, EXP_CLASS_BAD_TRUE,
                       EXP_CLASS_STATUS_NOTRUN, EXP_CLASS_UNKNOWN_STR,
                       DEFAULT_EXP_CLASS_MODEL

@testset "exposureClassifier" begin
    @testset "default model artifact is pinned" begin
        # The default must name one specific trained forest. "Whatever the
        # newest file in that directory is" would leave a reduction unable to
        # say which model judged it.
        @test occursin("_v6.jld2", DEFAULT_EXP_CLASS_MODEL)
        @test isabspath(DEFAULT_EXP_CLASS_MODEL)
        # deliberately NOT asserted to exist: the test suite must pass on a
        # machine without the 175 MB artifact
    end

    @testset "exposure_class_label" begin
        @test exposure_class_label("arclamp", false, true, false) == "arclamp_q0t1u0"
        @test exposure_class_label("quartzflat", 1, 0, 0) == "quartzflat_q1t0u0"
        @test exposure_class_label("dark", "F", "F", "F") == "dark_q0t0u0"
        # an unrecognized lamp value must be visible as "?", not silently 0
        @test exposure_class_label("dark", "maybe", 0, 0) == "dark_q?t0u0"
    end

    @testset "tri-state verdict: unknown is not 'fine'" begin
        u = exposure_class_unknown_metadata()
        @test u["exp_class_predicted_bad"] == EXP_CLASS_BAD_UNKNOWN
        @test u["exp_class_predicted_bad"] == Int8(-1)
        @test u["exp_class_status"] == EXP_CLASS_STATUS_NOTRUN
        @test u["exp_class_pred"] == EXP_CLASS_UNKNOWN_STR
        @test u["exp_class_labeled"] == EXP_CLASS_UNKNOWN_STR
        @test isnan(u["exp_class_prob"])
        # the whole point: unknown must be distinguishable from good
        @test u["exp_class_predicted_bad"] != EXP_CLASS_BAD_FALSE
    end

    @testset "verdict policy" begin
        # healthy quartzflat
        g = exposure_class_metadata("quartzflat_q1t0u0", "quartzflat_q1t0u0", 1.0, "ok")
        @test g["exp_class_predicted_bad"] == EXP_CLASS_BAD_FALSE
        @test g["exp_class_status"] == "ok"

        # lamp-off quartzflat
        b = exposure_class_metadata(
            "quartzflat_q1t0u0", "dark_q0t0u0", 1.0, "lamp_off_candidate")
        @test b["exp_class_predicted_bad"] == EXP_CLASS_BAD_TRUE
        @test b["exp_class_pred"] == "dark_q0t0u0"

        # standard science frames are never masked, whatever the forest thinks
        o = exposure_class_metadata(
            "object_q0t0u0", "dark_q0t0u0", 1.0, "mislabel_candidate")
        @test o["exp_class_predicted_bad"] == EXP_CLASS_BAD_FALSE

        # ...but the exemption is on the EXACT label, not image_type == object.
        # An object frame with an anomalous lamp flag (lco 57802 209-212 in
        # DR21) is labeled object_q0t0u1 and DOES come back masked. Harmless
        # while the only consumer is the fiber-flat runlist builder, which
        # never sees an object frame; asserted here so the nuance is recorded
        # rather than rediscovered.
        o2 = exposure_class_metadata(
            "object_q0t0u1", "object_q0t0u0", 1.0, "mislabel_candidate")
        @test o2["exp_class_predicted_bad"] == EXP_CLASS_BAD_TRUE

        # a dark following a bright exposure is informational, not bad
        p = exposure_class_metadata(
            "dark_q0t0u0", "dark_q0t0u0", 1.0, "persistence_prior")
        @test p["exp_class_predicted_bad"] == EXP_CLASS_BAD_FALSE

        # a CRASHED check is a failure to form an opinion, never an adverse
        # opinion: it must map to unknown, not to bad
        c = exposure_class_metadata("checkfail", "checkfail", NaN, "checkfail")
        @test c["exp_class_predicted_bad"] == EXP_CLASS_BAD_UNKNOWN
        @test c["exp_class_status"] == EXP_CLASS_STATUS_NOTRUN
        @test !exposure_predicted_bad("quartzflat_q1t0u0", "checkfail", "checkfail")

        # exposures with no reduced 2D data are not masked here
        @test !exposure_predicted_bad("domeflat_q0t0u0", "", "nofiles")
        @test !exposure_predicted_bad("domeflat_q0t0u0", "", "unclassified")
    end

    @testset "check table round trip and missing-file degradation" begin
        mktempdir() do dir
            tele, mjd = "lco", "60078"
            mkpath(joinpath(dir, "apred", mjd))
            path = exposure_type_check_path(dir, tele, mjd)
            safe_jldsave(path, Dict{String, Any}();
                tele = [tele, tele], mjd = [parse(Int, mjd), parse(Int, mjd)],
                expnum = [9, 44],
                labeled = ["domeflat_q0t0u0", "domeflat_q0t0u0"],
                pred = ["domeflat_q0t0u0", "arclamp_q0t0u0"],
                prob = [1.0, 0.7933],
                flag = ["ok", "mislabel_candidate"],
                predicted_bad = Int8[0, 1])

            tbl = read_exposure_type_check(path)
            @test length(tbl) == 2
            @test tbl[9]["exp_class_predicted_bad"] == EXP_CLASS_BAD_FALSE
            @test tbl[44]["exp_class_predicted_bad"] == EXP_CLASS_BAD_TRUE

            # the accessor used by process_1D
            @test exposure_class_metadata_for(dir, tele, mjd, 44)["exp_class_predicted_bad"] ==
                  EXP_CLASS_BAD_TRUE
            # an exposure absent from the table -> explicit unknown
            @test exposure_class_metadata_for(dir, tele, mjd, 999)["exp_class_status"] ==
                  EXP_CLASS_STATUS_NOTRUN
            # a night with no table at all -> explicit unknown, no error
            @test exposure_class_metadata_for(dir, "apo", "12345", 1)["exp_class_predicted_bad"] ==
                  EXP_CLASS_BAD_UNKNOWN
        end
        # a missing file must degrade quietly to "no verdicts", never throw:
        # this is advisory metadata and must not be able to break a reduction
        @test isempty(read_exposure_type_check("/definitely/not/a/path.h5"))
    end
end
