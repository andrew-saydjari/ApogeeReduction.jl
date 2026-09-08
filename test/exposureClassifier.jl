using ApogeeReduction: is_engineering_carton, exposure_is_engineering,
                       almanac_science_cartons, exposure_engineering_from_almanac,
                       exposure_flag_bits, exposure_ok_for_science,
                       ENGINEERING_CARTON_PREFIXES, ENGINEERING_CARTON_MIN_FRAC,
                       ENGINEERING_CHECK_IMAGE_TYPES, ENGINEERING_FALLBACK_ALL_FIBERS,
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

    @testset "majority rule" begin
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
        # strict majority: 51/100 flags, 50/100 does not
        @test exposure_is_engineering(vcat(fill(eng, 51), fill(sci, 49))).engineering
        @test !exposure_is_engineering(vcat(fill(eng, 50), fill(sci, 50))).engineering
        # a single stray engineering fiber must not condemn a science config
        @test !exposure_is_engineering(vcat([eng], fill(sci, 299))).engineering
        # threshold constant is the documented one
        @test ENGINEERING_CARTON_MIN_FRAC == 0.5
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
                @test length(almanac_science_cartons(f, "apo", "59625", 3472).cartons) == 250
                @test almanac_science_cartons(f, "apo", "59625", 3472).basis == "science"
                @test exposure_engineering_from_almanac(
                    f, "apo", "59625", 3472, "object").engineering
                @test !exposure_engineering_from_almanac(
                    f, "apo", "59625", 3500, "object").engineering
                # plate era: no firstcarton column -> empty, not an error
                @test almanac_science_cartons(f, "apo", "57674", 8662).cartons == String[]
                @test !exposure_engineering_from_almanac(
                    f, "apo", "57674", 8662, "object").engineering
                # plate-era config_id sentinel
                @test almanac_science_cartons(f, "apo", "57674", -1).cartons == String[]
                @test !exposure_engineering_from_almanac(
                    f, "apo", "57674", -1, "object").engineering
                # missing configuration group
                @test almanac_science_cartons(f, "apo", "59625", 999999).cartons == String[]
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

    @testset "no-science-fiber fallback" begin
        # MEASURED: 5 early-FPS configurations (apo 59558/105, 59558/106,
        # 59560/121, 59560/122, 59561/133) carry manual_fps_position_stars on
        # 207-254 of 300 fibers with ZERO category=="science" fibers — their
        # categories are "", bonus, open_fiber, sky_boss. Without the fallback
        # they back 35 unflagged engineering object exposures.
        # (ENGINEERING_FALLBACK_ALL_FIBERS is an addition beyond AKS's
        # science-fibers instruction; flip it to false to get the literal rule.)
        @test ENGINEERING_FALLBACK_ALL_FIBERS
        mktempdir() do dir
            path = joinpath(dir, "alm2.h5")
            h5open(path, "w") do f
                g = create_group(f, "raw/apo/59558/fibers/105")
                g["category"] = vcat(fill("", 245), fill("sky_boss", 40),
                    fill("bonus", 15))
                g["firstcarton"] = vcat(fill("manual_fps_position_stars", 245),
                    fill("", 55))
                # same shape but a real science program: must NOT be flagged
                g2 = create_group(f, "raw/apo/59558/fibers/106")
                g2["category"] = vcat(fill("", 245), fill("sky_boss", 55))
                g2["firstcarton"] = vcat(fill("mwm_snc_100pc", 245), fill("", 55))
            end
            h5open(path, "r") do f
                r = almanac_science_cartons(f, "apo", "59558", 105)
                @test r.basis == "all_fibers_fallback"
                @test length(r.cartons) == 245   # empty cartons dropped
                @test exposure_engineering_from_almanac(
                    f, "apo", "59558", 105, "object").engineering
                @test exposure_engineering_from_almanac(
                    f, "apo", "59558", 105, "object").basis == "all_fibers_fallback"
                @test !exposure_engineering_from_almanac(
                    f, "apo", "59558", 106, "object").engineering
                # with the fallback disabled the literal AKS rule applies and
                # neither is flagged (this is the one-line revert)
                @test almanac_science_cartons(f, "apo", "59558", 105;
                    fallback_all_fibers = false).cartons == String[]
            end
        end
    end

    @testset "exposure_flags bitmask" begin
        @test exposure_flag_bits(false, false) == 0x00
        @test exposure_flag_bits(true, false) == EXPFLAG_PREDICTED_BAD
        @test exposure_flag_bits(false, true) == EXPFLAG_ENGINEERING
        @test exposure_flag_bits(true, true) == (EXPFLAG_PREDICTED_BAD | EXPFLAG_ENGINEERING)
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
