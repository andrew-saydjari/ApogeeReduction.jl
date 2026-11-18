using HDF5

# arMADGICS.jl's selection predicates, reproduced here so the test pins the
# contract between the two packages (arMADGICS.jl src/ingest.jl):
#   fitted:  x[1:3] in accepted_fibtypes, accepted_fibtypes = ["sci", "tel"]
#   sky:     x[1:3] == "sky"   (so "skyB" counts as sky)
_madgics_fits(t) = t[1:3] in ["sci", "tel"]
_madgics_sky(t) = t[1:3] == "sky"

# Fiber table modelled on apo confSummaryS-27968 (MJD 61330), a manually built
# design with ZERO `science` fibers: 15 standard_apogee, 35 sky_apogee,
# 226 veto_location_boss, 8 veto_location_apogee, 2 veto_location_boss fibers
# that were assigned but not reached (written as almanac fork PR #10 writes them:
# category "offtarget", original_category kept, sdss_id withheld), 12 fibers
# with an empty category, and 2 "bonus" FPI stubs. (almanac writes APOGEE rows
# only; BOSS rows never reach the fiber table.)
function _veto_design_rows()
    cats = vcat(fill("standard_apogee", 15), fill("sky_apogee", 35),
        fill("veto_location_boss", 226), fill("veto_location_apogee", 8),
        fill("offtarget", 2), fill("", 12), fill("bonus", 2))
    @assert length(cats) == 300
    fiber_id = collect(1:300)
    fiber_type = [c == "bonus" ? "" : "APOGEE" for c in cats]
    orig = [c == "offtarget" ? "veto_location_boss" : "" for c in cats]
    carton = map(cats) do c
        c == "standard_apogee" ? "ops_std_apogee" :
        c == "sky_apogee" ? "ops_sky_apogee_best" :
        c in ("veto_location_boss", "offtarget") ? "ops_gaia_brightneighbors" :
        c == "veto_location_apogee" ? "ops_2mass_psc_brightneighbors" : ""
    end
    sdss_id = [c in ("standard_apogee", "veto_location_boss", "veto_location_apogee") ?
               100_000 + i : -1 for (i, c) in enumerate(cats)]
    return (; category = cats, fiber_id, fiber_type, original_category = orig,
        firstcarton = carton, sdss_id)
end

function _write_veto_almanac(path; tele = "apo", mjd = "61330", config_id = 27968)
    rows = _veto_design_rows()
    h5open(path, "w") do f
        base = "raw/$(tele)/$(mjd)"
        write(f, "$(base)/exposures/image_type", ["object", "object"])
        write(f, "$(base)/exposures/config_id", [config_id, config_id])
        write(f, "$(base)/exposures/exposure", [1, 2])
        for k in keys(rows)
            write(f, "$(base)/fibers/$(config_id)/$(k)", getfield(rows, k))
        end
    end
    return rows
end

@testset "almanac category -> AR.jl fiber type" begin
    m = ApogeeReduction.ALMANAC_CATEGORY_TO_FIBER_TYPE

    @testset "mapping" begin
        @test m["science"] == "sci"
        @test m["veto_location_boss"] == "sci"
        @test m["veto_location_apogee"] == "sci"
        @test m["standard_apogee"] == "tel"
        @test m["sky_apogee"] == "sky"
        @test m["sky_boss"] == "skyB"
        @test m["offtarget"] == "offtarget"
        # every mapped type survives arMADGICS.jl's `x[1:3]` without a BoundsError
        @test all(length(t) >= 3 for t in values(m))
        # the veto-location stars are fitted like science
        @test _madgics_fits(m["veto_location_boss"])
        @test _madgics_fits(m["veto_location_apogee"])
        # off-target fibers are neither fitted nor used as sky
        @test !_madgics_fits(m["offtarget"])
        @test !_madgics_sky(m["offtarget"])
        @test !_madgics_fits("fiberTypeFail") && !_madgics_sky("fiberTypeFail")
    end

    mktempdir() do dir
        path = joinpath(dir, "alm_veto.h5")
        rows = _write_veto_almanac(path)

        @testset "get_fibTargDict on a veto-location design" begin
            fibtarg, sdss = h5open(path, "r") do f
                ApogeeReduction.get_fibTargDict(f, "apo", "61330", 1)
            end
            types = [fibtarg[i] for i in 1:300]
            @test count(==("sci"), types) == 226 + 8
            @test count(==("tel"), types) == 15
            @test count(==("sky"), types) == 35
            @test count(==("offtarget"), types) == 2
            @test count(==("fpiguide"), types) == 2
            @test count(==("fiberTypeFail"), types) == 12   # empty category only

            # per-fiber agreement with the fixture, by fiber index
            apogee = findall(==("APOGEE"), rows.fiber_type)
            for j in apogee
                cat = rows.category[j]
                idx = ApogeeReduction.fiberID2fiberIndx(rows.fiber_id[j])
                if cat in ("veto_location_boss", "veto_location_apogee")
                    @test fibtarg[idx] == "sci"
                    @test sdss[idx] > 0
                elseif cat == "offtarget"
                    @test fibtarg[idx] == "offtarget"
                    @test sdss[idx] == -1
                end
            end

            # what arMADGICS.jl would do with this exposure
            @test count(_madgics_fits, types) == 226 + 8 + 15
            @test count(_madgics_sky, types) == 35
            off = findall(==("offtarget"), types)
            @test !any(_madgics_fits, types[off]) && !any(_madgics_sky, types[off])
        end

        @testset "veto-location design is not an engineering / science-less config" begin
            h5open(path, "r") do f
                c = almanac_config_cartons(f, "apo", "61330", 27968)
                # 226 + 8 veto-location fibers, plus the 2 off-target ones
                # counted by their original category (design intent)
                @test length(c.science) == 226 + 8 + 2
                @test Set(c.science) ==
                      Set(["ops_gaia_brightneighbors", "ops_2mass_psc_brightneighbors"])
                v = exposure_engineering_from_almanac(f, "apo", "61330", 27968, "object")
                @test !v.engineering
                @test v.basis == "science"     # clause 1, not science-less
                @test v.nsci == 236
            end
        end
    end

    @testset "science-less clause still keyed on science, veto counts as science" begin
        # Same position-stars carton on one fiber: with only veto-location
        # "science" fibers present, clause 1 (purity) applies and does NOT flag;
        # with no science-like fibers at all, clause 2 flags.
        eng = "manual_fps_position_stars"
        mktempdir() do dir
            path = joinpath(dir, "alm_eng.h5")
            h5open(path, "w") do f
                g = create_group(f, "raw/apo/61330/fibers/1")
                g["category"] = vcat(["sky_boss"], fill("veto_location_boss", 299))
                g["firstcarton"] = vcat([eng], fill("ops_gaia_brightneighbors", 299))
                g2 = create_group(f, "raw/apo/61330/fibers/2")
                g2["category"] = vcat(["sky_boss"], fill("", 299))
                g2["firstcarton"] = vcat([eng], fill("ops_gaia_brightneighbors", 299))
            end
            h5open(path, "r") do f
                v1 = exposure_engineering_from_almanac(f, "apo", "61330", 1, "object")
                @test !v1.engineering && v1.basis == "science"
                v2 = exposure_engineering_from_almanac(f, "apo", "61330", 2, "object")
                @test v2.engineering && v2.basis == "scienceless_position_stars"
            end
        end
    end
end
