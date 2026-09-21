using ApogeeReduction: exposure_class_metadata, exposure_class_unknown_metadata,
                       exposure_predicted_bad, exposure_class_label,
                       exposure_type_check_path, read_exposure_type_check,
                       exposure_class_metadata_for, safe_jldsave,
                       exposure_class_verdict, exposure_flag_bits,
                       exposure_ok_for_science,
                       EXPFLAG_PREDICTED_BAD, EXPFLAG_ENGINEERING, EXPFLAG_NOTRUN,
                       EXPFLAG_NO_SCIENCE,
                       EXP_CLASS_STATUS_NOTRUN, EXP_CLASS_UNKNOWN_STR,
                       is_engineering_carton, exposure_is_engineering,
                       almanac_config_cartons, engineering_verdict,
                       exposure_engineering_from_almanac,
                       ENGINEERING_CARTON_PREFIXES, ENGINEERING_CARTON_PURITY,
                       ENGINEERING_CARTON_WARN_FRAC,
                       ENGINEERING_CHECK_IMAGE_TYPES,
                       ENGINEERING_FLAG_SCIENCELESS_POSITION_STARS,
                       ENGINEERING_FLAG_DESIGNLESS,
                       ENGINEERING_DESIGNLESS_DESIGN_ID,
                       is_designless_design, apply_designless_clause
using HDF5

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

@testset "engineering carton flag" begin
    @testset "carton name matching" begin
        # all four manual_fps_position_stars* variants present in the
        # 57618-61230 corpus must match the single prefix
        for c in ("manual_fps_position_stars", "manual_fps_position_stars_10",
            "manual_fps_position_stars_apogee_10",
            "manual_fps_position_stars_lco_apogee_10")
            @test is_engineering_carton(c)
            @test is_engineering_carton(uppercase(c))   # case-insensitive
            @test is_engineering_carton("  " * c * " ") # whitespace tolerant
        end
        # real science cartons must not match
        for c in ("mwm_snc_100pc", "bhm_csc_apogee", "manual_nsbh_apogee",
            "openfibertargets_nov2020_20", "", "science")
            @test !is_engineering_carton(c)
        end
        # guard against the prefix list being silently emptied
        @test !isempty(ENGINEERING_CARTON_PREFIXES)
        @test "manual_fps_position_stars" in ENGINEERING_CARTON_PREFIXES
    end

    @testset "purity rule" begin
        eng = "manual_fps_position_stars"
        sci = "mwm_snc_100pc"
        # the real-world case: 100% engineering
        r = exposure_is_engineering(fill(eng, 254))
        @test r.engineering
        @test r.frac == 1.0
        @test r.carton == eng
        @test r.nsci == 254
        @test r.basis == "science"
        # 100% science
        r = exposure_is_engineering(fill(sci, 254))
        @test !r.engineering
        @test r.frac == 0.0
        @test r.carton == ""
        # PURITY: anything short of 100% is NOT engineering (AKS 2026-09-08).
        # 299/300 fails; a majority fails; a single stray fiber fails.
        @test !exposure_is_engineering(vcat(fill(eng, 299), [sci])).engineering
        @test !exposure_is_engineering(vcat(fill(eng, 51), fill(sci, 49))).engineering
        @test !exposure_is_engineering(vcat(fill(eng, 50), fill(sci, 50))).engineering
        @test !exposure_is_engineering(vcat([eng], fill(sci, 299))).engineering
        # this is what keeps the 27 low-share manual_* cartons, and the single
        # 100%-pure manual_mwm_crosscalib_apogee config, out of the flag
        @test !exposure_is_engineering(vcat(fill(eng, 71), fill(sci, 29))).engineering
        # constants are the documented ones
        @test ENGINEERING_CARTON_PURITY == 1.0
        @test ENGINEERING_CARTON_WARN_FRAC == 0.5
        # MOSTLY-but-not-purely must WARN loudly (never seen in DR21, so it is
        # a real signal if it fires) and must NOT flag
        @test_logs (:warn, r"MOSTLY but not PURELY") match_mode=:any begin
            r = exposure_is_engineering(vcat(fill(eng, 90), fill(sci, 10)))
            @test !r.engineering
        end
        # below the warn threshold: no warning, no flag
        @test_logs min_level=Base.CoreLogging.Warn begin
            @test !exposure_is_engineering(vcat(fill(eng, 10), fill(sci, 90))).engineering
        end
        # no science fibers (plate era / no config) -> never engineering, no error
        r = exposure_is_engineering(String[])
        @test !r.engineering
        @test isnan(r.frac)
        @test r.nsci == 0
        # the dominant matching carton is reported
        r = exposure_is_engineering(vcat(fill(eng * "_10", 200), fill(eng, 50)))
        @test r.engineering
        @test r.carton == eng * "_10"
    end

    @testset "almanac reads" begin
        mktempdir() do dir
            path = joinpath(dir, "alm.h5")
            h5open(path, "w") do f
                # FPS-era engineering config
                g = create_group(f, "raw/apo/59625/fibers/3472")
                g["category"] = vcat(fill("science", 250), fill("sky_apogee", 40))
                g["firstcarton"] = vcat(fill("manual_fps_position_stars", 250),
                    fill("", 40))
                # FPS-era science config
                g2 = create_group(f, "raw/apo/59625/fibers/3500")
                g2["category"] = vcat(fill("science", 250), fill("sky_apogee", 40))
                g2["firstcarton"] = vcat(fill("mwm_snc_100pc", 250), fill("", 40))
                # plate-era fiber table: no firstcarton column at all
                g3 = create_group(f, "raw/apo/57674/fibers/8662")
                g3["category"] = fill("science", 200)
                g3["fiber_id"] = collect(1:200)
            end
            h5open(path, "r") do f
                @test length(almanac_config_cartons(f, "apo", "59625", 3472).science) == 250
                @test exposure_engineering_from_almanac(
                    f, "apo", "59625", 3472, "object").basis == "science"
                @test exposure_engineering_from_almanac(
                    f, "apo", "59625", 3472, "object").engineering
                @test !exposure_engineering_from_almanac(
                    f, "apo", "59625", 3500, "object").engineering
                # plate era: no firstcarton column -> empty, not an error
                @test almanac_config_cartons(f, "apo", "57674", 8662).science == String[]
                @test almanac_config_cartons(f, "apo", "57674", 8662).all == String[]
                @test !exposure_engineering_from_almanac(
                    f, "apo", "57674", 8662, "object").engineering
                # plate-era config_id sentinel
                @test almanac_config_cartons(f, "apo", "57674", -1).science == String[]
                @test !exposure_engineering_from_almanac(
                    f, "apo", "57674", -1, "object").engineering
                # missing configuration group
                @test almanac_config_cartons(f, "apo", "59625", 999999).science == String[]
                # calibration exposures riding an engineering config are NOT
                # engineering: a dark taken under config 3472 is still a dark
                for it in ("dark", "quartzflat", "domeflat", "arclamp", "internalflat")
                    @test !exposure_engineering_from_almanac(
                        f, "apo", "59625", 3472, it).engineering
                end
                @test ENGINEERING_CHECK_IMAGE_TYPES == ["object"]
            end
        end
    end

    @testset "clause 2: science-less position-stars configs" begin
        # AKS 2026-09-08: "if any fibers are manual_fps_position_stars* and no
        # science fibers, then reject."
        #
        # MEASURED over the full DR21 almanac: exactly 5 configurations
        # (apo 59558/105, 59558/106, 59560/121, 59560/122, 59561/133) carry
        # manual_fps_position_stars on 207-254 of 300 fibers with ZERO
        # category=="science" fibers — categories are "", bonus, open_fiber,
        # sky_boss. They back 35 object exposures, with ZERO overlap with
        # clause 1, so enabling this cannot perturb the 2,448.
        @test ENGINEERING_FLAG_SCIENCELESS_POSITION_STARS

        eng = "manual_fps_position_stars"
        sci = "mwm_snc_100pc"

        # The rule is "ANY position-stars fiber", not purity: a single one is
        # enough when there are no science fibers.
        r = engineering_verdict(String[], vcat([eng], fill(sci, 299)))
        @test r.engineering
        @test r.basis == "scienceless_position_stars"
        @test r.nsci == 0            # zero science fibers is what triggered it
        # ...but the SAME fiber content WITH science fibers goes to clause 1 and
        # is rejected by purity. This is the pair that pins the two clauses apart.
        @test !engineering_verdict(vcat([eng], fill(sci, 299)),
            vcat([eng], fill(sci, 299))).engineering

        # Not a general all-fibers fallback: a science-less config carrying some
        # OTHER carton must NOT be flagged, however uniform it is.
        r2 = engineering_verdict(String[], fill(sci, 300))
        @test !r2.engineering
        @test r2.basis == "none"
        # ...and no mostly-but-not-purely warning is raised on the clause-2 path
        @test_logs min_level=Base.CoreLogging.Warn begin
            @test !engineering_verdict(String[],
                vcat(fill(sci, 200), fill("mwm_bin_rv", 100))).engineering
        end

        # nothing at all -> not engineering, no error
        @test !engineering_verdict(String[], String[]).engineering
        @test engineering_verdict(String[], String[]).basis == "none"

        mktempdir() do dir
            path = joinpath(dir, "alm2.h5")
            h5open(path, "w") do f
                # the real apo 59558/105 shape: 245/300 position-stars, 0 science
                g = create_group(f, "raw/apo/59558/fibers/105")
                g["category"] = vcat(fill("", 245), fill("sky_boss", 40),
                    fill("bonus", 15))
                g["firstcarton"] = vcat(fill(eng, 245), fill("", 55))
                # same shape, real science program: must NOT be flagged
                g2 = create_group(f, "raw/apo/59558/fibers/106")
                g2["category"] = vcat(fill("", 245), fill("sky_boss", 55))
                g2["firstcarton"] = vcat(fill(sci, 245), fill("", 55))
            end
            h5open(path, "r") do f
                c = almanac_config_cartons(f, "apo", "59558", 105)
                @test c.science == String[]        # zero science fibers
                @test length(c.all) == 245         # blank cartons dropped
                v = exposure_engineering_from_almanac(f, "apo", "59558", 105, "object")
                @test v.engineering
                @test v.basis == "scienceless_position_stars"
                @test v.frac == 1.0                # 245/245 carton-bearing fibers
                @test !exposure_engineering_from_almanac(
                    f, "apo", "59558", 106, "object").engineering
                # clause 2 still respects the image-type restriction
                @test !exposure_engineering_from_almanac(
                    f, "apo", "59558", 105, "dark").engineering
            end
        end
    end

    @testset "clause 3: designless configurations (design_id == -999)" begin
        # AKS 2026-09-09: "implement the -999 clause".
        #
        # MEASURED over the full DR21 almanac: 17,311 rows carry design_id ==
        # -999, but 17,296 (99.91%) are CALIBRATION frames. Exactly 15 are
        # `object`: apo 59697 (2), 59765 (11), 60212 (2). The image-type
        # restriction is therefore load-bearing, not cosmetic — the tests below
        # assert it directly.
        @test ENGINEERING_FLAG_DESIGNLESS
        @test ENGINEERING_DESIGNLESS_DESIGN_ID == -999

        # -- the predicate is deliberately total --
        @test is_designless_design(-999)
        @test !is_designless_design(-1)         # almanac's OWN missing sentinel
        @test !is_designless_design(0)          # never occurs in the corpus
        @test !is_designless_design(388590)     # a real design
        @test !is_designless_design(nothing)    # column absent: unknown != absent
        @test !is_designless_design("-999")     # not an Integer
        @test !is_designless_design(missing)

        none = (engineering = false, frac = NaN, carton = "", nsci = 0, basis = "none")
        # -- clause 3 flags an otherwise-clean verdict --
        v = apply_designless_clause(none, -999)
        @test v.engineering
        @test v.basis == "designless"
        @test v.nsci == 0
        # -- and leaves everything else untouched --
        @test !apply_designless_clause(none, -1).engineering
        @test !apply_designless_clause(none, 388590).engineering
        @test !apply_designless_clause(none, nothing).engineering
        # -- a verdict already flagged keeps its MORE SPECIFIC basis --
        c1 = (engineering = true, frac = 1.0, carton = "manual_fps_position_stars",
            nsci = 250, basis = "science")
        @test apply_designless_clause(c1, -999).basis == "science"
        c2 = (engineering = true, frac = 1.0, carton = "manual_fps_position_stars",
            nsci = 0, basis = "scienceless_position_stars")
        @test apply_designless_clause(c2, -999).basis == "scienceless_position_stars"
        # -- idempotent: applying twice cannot change the answer (both call
        #    sites apply it after a cache lookup, so this must hold) --
        @test apply_designless_clause(apply_designless_clause(none, -999), -999).basis ==
              "designless"

        # -- end to end against a fixture shaped like the real thing --
        mktempdir() do dir
            path = joinpath(dir, "alm_designless.h5")
            h5open(path, "w") do f
                # apo 60212 config 10680: the twilight-ladder configuration.
                # 300 fibers, 0 science, 0 sky, 0 cartons, 298 blank + 2 bonus,
                # assigned = 0 — exactly as measured.
                g = create_group(f, "raw/apo/60212/fibers/10680")
                g["category"] = vcat(fill("", 298), fill("bonus", 2))
                g["firstcarton"] = fill("", 300)
                # config 10681 on the same night IS science: the control that
                # proves the emptiness above is real and not an ingest failure
                g2 = create_group(f, "raw/apo/60212/fibers/10681")
                g2["category"] = vcat(fill("science", 175), fill("sky_apogee", 91),
                    fill("standard_apogee", 15), fill("", 19))
                g2["firstcarton"] = vcat(fill("mwm_snc_100pc", 175), fill("", 125))
            end
            h5open(path, "r") do f
                # the carton clauses alone say nothing: no science, no cartons
                @test !exposure_engineering_from_almanac(
                    f, "apo", "60212", 10680, "object").engineering
                # clause 3 catches it
                v = exposure_engineering_from_almanac(
                    f, "apo", "60212", 10680, "object"; design_id = -999)
                @test v.engineering
                @test v.basis == "designless"

                # THE LOAD-BEARING ASSERTION: 17,296 calibration frames carry
                # design_id == -999 and must survive untouched. If this ever
                # fails, most FPS-era APO arclamps and darks are being thrown
                # away.
                for it in ("dark", "quartzflat", "domeflat", "arclamp",
                    "internalflat", "twilightflat")
                    @test !exposure_engineering_from_almanac(
                        f, "apo", "60212", 10680, it; design_id = -999).engineering
                end

                # a real design on the same night is untouched
                @test !exposure_engineering_from_almanac(
                    f, "apo", "60212", 10681, "object"; design_id = 382526).engineering
                # omitting design_id evaluates the carton clauses alone
                @test !exposure_engineering_from_almanac(
                    f, "apo", "60212", 10680, "object").engineering
            end
        end
    end

    @testset "exposure_flags bitmask" begin
        @test exposure_flag_bits(false, false) == 0x00
        @test exposure_flag_bits(true, false) == EXPFLAG_PREDICTED_BAD
        @test exposure_flag_bits(false, true) == EXPFLAG_ENGINEERING
        @test exposure_flag_bits(true, true) == (EXPFLAG_PREDICTED_BAD | EXPFLAG_ENGINEERING)
        # `extra` lets a future producer OR in a bit without any side
        # renumbering bits 0, 1 and 2.
        #
        # MERGE NOTE (#397 + #398): these two assertions originally used
        # `extra = 0x04`, chosen when 2^2 was UNASSIGNED. #398 has since claimed
        # 2^2 as EXPFLAG_NOTRUN (AKS: "please don't reuse 2^2"), and NOTRUN is
        # mutually exclusive with PREDICTED_BAD — so the old
        # `exposure_flag_bits(true, true; extra = 0x04) == 0x07` now asserts a
        # byte the packer is required to REJECT. The intent of the test is the
        # `extra` MECHANISM, not the literal 0x04, so it moves to the first
        # still-unassigned bit (2^3). The mechanism is tested exactly as before.
        @test exposure_flag_bits(false, false; extra = 0x08) == 0x08
        @test exposure_flag_bits(true, true; extra = 0x08) == 0x0b
        # ...and the invariant is enforced on the COMPOSED byte, so smuggling
        # the notrun bit in through `extra` cannot evade the mutual exclusion
        @test_throws ArgumentError exposure_flag_bits(true, false; extra = EXPFLAG_NOTRUN)
        @test EXPFLAG_PREDICTED_BAD == 0x01
        @test EXPFLAG_ENGINEERING == 0x02
        @test EXPFLAG_NOTRUN == 0x04
        # the guard predicate downstream science-sample assembly relies on.
        # If someone drops a bit out of EXPFLAG_NO_SCIENCE these fail.
        @test exposure_ok_for_science(0x00)
        @test !exposure_ok_for_science(EXPFLAG_PREDICTED_BAD)
        @test !exposure_ok_for_science(EXPFLAG_ENGINEERING)
        @test !exposure_ok_for_science(EXPFLAG_PREDICTED_BAD | EXPFLAG_ENGINEERING)
        @test EXPFLAG_NO_SCIENCE & EXPFLAG_PREDICTED_BAD != 0
        @test EXPFLAG_NO_SCIENCE & EXPFLAG_ENGINEERING != 0
    end
end
