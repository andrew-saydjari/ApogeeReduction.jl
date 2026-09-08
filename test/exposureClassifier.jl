using ApogeeReduction: is_engineering_carton, exposure_is_engineering,
                       almanac_config_cartons, engineering_verdict,
                       exposure_engineering_from_almanac,
                       exposure_flag_bits, exposure_ok_for_science,
                       ENGINEERING_CARTON_PREFIXES, ENGINEERING_CARTON_PURITY,
                       ENGINEERING_CARTON_WARN_FRAC,
                       ENGINEERING_CHECK_IMAGE_TYPES,
                       ENGINEERING_FLAG_SCIENCELESS_POSITION_STARS,
                       EXPFLAG_PREDICTED_BAD, EXPFLAG_ENGINEERING, EXPFLAG_NO_SCIENCE
using HDF5

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

    @testset "exposure_flags bitmask" begin
        @test exposure_flag_bits(false, false) == 0x00
        @test exposure_flag_bits(true, false) == EXPFLAG_PREDICTED_BAD
        @test exposure_flag_bits(false, true) == EXPFLAG_ENGINEERING
        @test exposure_flag_bits(true, true) == (EXPFLAG_PREDICTED_BAD | EXPFLAG_ENGINEERING)
        # `extra` lets the concurrent classifier-propagation work OR in a bit
        # (e.g. NOT-RUN) without either side renumbering bits 0 and 1
        @test exposure_flag_bits(false, false; extra = 0x04) == 0x04
        @test exposure_flag_bits(true, true; extra = 0x04) == 0x07
        @test EXPFLAG_PREDICTED_BAD == 0x01
        @test EXPFLAG_ENGINEERING == 0x02
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
