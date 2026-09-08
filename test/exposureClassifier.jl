using ApogeeReduction: exposure_class_metadata, exposure_class_unknown_metadata,
                       exposure_predicted_bad, exposure_class_label,
                       exposure_type_check_path, read_exposure_type_check,
                       exposure_class_metadata_for, safe_jldsave,
                       exposure_class_verdict, exposure_flag_bits,
                       exposure_ok_for_science,
                       EXPFLAG_PREDICTED_BAD, EXPFLAG_ENGINEERING, EXPFLAG_NOTRUN,
                       EXPFLAG_NO_SCIENCE,
                       EXP_CLASS_STATUS_NOTRUN, EXP_CLASS_UNKNOWN_STR

@testset "exposureClassifier" begin
    @testset "shared exposure_flags bit schema" begin
        # These values are shared with the engineering-carton check (PR #397).
        # Renumbering silently corrupts every consumer of the byte, so pin them.
        @test EXPFLAG_PREDICTED_BAD == 0x01
        @test EXPFLAG_ENGINEERING == 0x02
        @test EXPFLAG_NOTRUN == 0x04
        # notrun is NOT a no-science bit: an unjudged frame is not known-bad,
        # and excluding everything we failed to look at would hide data loss
        @test EXPFLAG_NO_SCIENCE == 0x03
        @test (EXPFLAG_NO_SCIENCE & EXPFLAG_NOTRUN) == 0x00
        @test exposure_ok_for_science(EXPFLAG_NOTRUN)
        @test exposure_flag_bits(false, false) == 0x00
        @test exposure_flag_bits(true, false) == EXPFLAG_PREDICTED_BAD
        @test exposure_flag_bits(false, true) == EXPFLAG_ENGINEERING
        @test exposure_flag_bits(true, true) == EXPFLAG_NO_SCIENCE
        @test exposure_flag_bits(false, false; notrun = true) == EXPFLAG_NOTRUN
        @test exposure_flag_bits(false, true; notrun = true) ==
              (EXPFLAG_ENGINEERING | EXPFLAG_NOTRUN)
        # mutual exclusion is enforced, not merely documented: "no verdict"
        # and "adverse verdict" cannot both be true
        @test_throws ArgumentError exposure_flag_bits(true, false; notrun = true)
        @test exposure_ok_for_science(0x00)
        @test !exposure_ok_for_science(EXPFLAG_PREDICTED_BAD)
        @test !exposure_ok_for_science(EXPFLAG_ENGINEERING)
    end

    @testset "exposure_class_label" begin
        @test exposure_class_label("arclamp", false, true, false) == "arclamp_q0t1u0"
        @test exposure_class_label("quartzflat", 1, 0, 0) == "quartzflat_q1t0u0"
        @test exposure_class_label("dark", "F", "F", "F") == "dark_q0t0u0"
        # an unrecognized lamp value must be visible as "?", not silently 0
        @test exposure_class_label("dark", "maybe", 0, 0) == "dark_q?t0u0"
    end

    @testset "notrun is a BIT: zero means judged-and-fine" begin
        u = exposure_class_unknown_metadata()
        # the load-bearing assertion: "never judged" is a distinct VALUE of the
        # mask, readable from the byte alone
        @test (UInt8(u["exposure_flags"]) & EXPFLAG_NOTRUN) != 0x00
        @test u["exposure_flags"] != 0x00
        @test exposure_class_verdict(u) == :unknown
        # invariant: no verdict cannot also be an adverse verdict
        @test (UInt8(u["exposure_flags"]) & EXPFLAG_PREDICTED_BAD) == 0x00
        @test u["exp_class_status"] == EXP_CLASS_STATUS_NOTRUN
        @test u["exp_class_pred"] == EXP_CLASS_UNKNOWN_STR
        @test isnan(u["exp_class_prob"])

        # a judged-fine exposure is exactly zero, and that is now unambiguous
        g = exposure_class_metadata("quartzflat_q1t0u0", "quartzflat_q1t0u0", 1.0, "ok")
        @test g["exposure_flags"] == 0x00
        @test exposure_class_verdict(g) == :fine
        @test g["exposure_flags"] != u["exposure_flags"]

        # compatibility shim: a product predating exposure_flags has no such
        # field and must still read as unknown, not as fine
        @test exposure_class_verdict(Dict{String, Any}("cartid" => 19)) == :unknown
    end

    @testset "verdict policy" begin
        # healthy quartzflat
        g = exposure_class_metadata("quartzflat_q1t0u0", "quartzflat_q1t0u0", 1.0, "ok")
        @test exposure_class_verdict(g) == :fine
        @test g["exp_class_status"] == "ok"

        # lamp-off quartzflat
        b = exposure_class_metadata(
            "quartzflat_q1t0u0", "dark_q0t0u0", 1.0, "lamp_off_candidate")
        @test exposure_class_verdict(b) == :bad
        @test (b["exposure_flags"] & EXPFLAG_PREDICTED_BAD) != 0x00
        @test b["exp_class_pred"] == "dark_q0t0u0"

        # the engineering bit ORs in at this same call site (PR #397 hand-off)
        e = exposure_class_metadata("object_q0t0u0", "object_q0t0u0", 1.0, "ok";
            engineering = true)
        @test e["exposure_flags"] == EXPFLAG_ENGINEERING
        @test exposure_class_verdict(e) == :fine   # engineering is not "bad"
        @test !exposure_ok_for_science(e["exposure_flags"])

        # standard science frames are never masked, whatever the forest thinks
        o = exposure_class_metadata(
            "object_q0t0u0", "dark_q0t0u0", 1.0, "mislabel_candidate")
        @test exposure_class_verdict(o) == :fine

        # ...but the exemption is on the EXACT label, not image_type == object.
        # An object frame with an anomalous lamp flag (lco 57802 209-212 in
        # DR21) is labeled object_q0t0u1 and DOES come back masked. Harmless
        # while the only consumer is the fiber-flat runlist builder, which
        # never sees an object frame; asserted here so the nuance is recorded
        # rather than rediscovered.
        o2 = exposure_class_metadata(
            "object_q0t0u1", "object_q0t0u0", 1.0, "mislabel_candidate")
        @test exposure_class_verdict(o2) == :bad

        # a dark following a bright exposure is informational, not bad
        p = exposure_class_metadata(
            "dark_q0t0u0", "dark_q0t0u0", 1.0, "persistence_prior")
        @test exposure_class_verdict(p) == :fine

        # a CRASHED check is a failure to form an opinion, never an adverse
        # opinion: it must map to unknown, not to bad
        c = exposure_class_metadata("checkfail", "checkfail", NaN, "checkfail")
        @test exposure_class_verdict(c) == :unknown
        @test (UInt8(c["exposure_flags"]) & EXPFLAG_NOTRUN) != 0x00
        @test (UInt8(c["exposure_flags"]) & EXPFLAG_PREDICTED_BAD) == 0x00
        # checkfail shares the notrun BIT but keeps its own status STRING, so
        # the diagnostic distinction survives without spending a shared bit
        @test c["exp_class_status"] == "checkfail"
        @test exposure_class_unknown_metadata()["exp_class_status"] ==
              EXP_CLASS_STATUS_NOTRUN
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
                exposure_flags = UInt8[0x00, 0x01])

            tbl = read_exposure_type_check(path)
            @test length(tbl) == 2
            @test exposure_class_verdict(tbl[9]) == :fine
            @test exposure_class_verdict(tbl[44]) == :bad

            # the accessor used by process_1D
            @test exposure_class_verdict(
                exposure_class_metadata_for(dir, tele, mjd, 44)) == :bad
            # an exposure absent from the table -> explicit unknown
            @test exposure_class_metadata_for(dir, tele, mjd, 999)["exp_class_status"] ==
                  EXP_CLASS_STATUS_NOTRUN
            # a night with no table at all -> explicit unknown, no error
            @test exposure_class_verdict(
                exposure_class_metadata_for(dir, "apo", "12345", 1)) == :unknown
        end
        # a missing file must degrade quietly to "no verdicts", never throw:
        # this is advisory metadata and must not be able to break a reduction
        @test isempty(read_exposure_type_check("/definitely/not/a/path.h5"))
    end
end
