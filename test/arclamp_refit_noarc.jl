# Nights without a usable arc fit: the within-night no-arc ladder (src/arclamp_refit_noarc.jl).
using Test, HDF5, Statistics
using ApogeeReduction
const ARN = ApogeeReduction

@testset "noarc ladder: config, bits, stage map" begin
    c = ARN.RefitConfig()
    @test c.noarc_route == "ladder" && c.noarc_plate_shape == "sky" && c.ignore_arcs == false
    @test ARN.refit_check_config(c)
    @test_throws ArgumentError ARN.refit_check_config(ARN.RefitConfig(noarc_route = "pass1_copy"))   # retired
    @test_throws ArgumentError ARN.refit_check_config(ARN.RefitConfig(noarc_plate_shape = "arc"))
    @test ARN.refit_check_config(ARN.RefitConfig(noarc_route = "template_sky", noarc_plate_shape = "template"))
    # no-arc switches invalidate no stage globally (the driver recomputes no-arc nights per night)
    for k in ("noarc_route", "noarc_plate_shape", "ignore_arcs")
        @test ARN.REFIT_SWITCH_STAGES[k] == "none"
        @test haskey(ARN.REFIT_SWITCH_LEGACY, k)
    end
    @test ARN.REFIT_SWITCH_LEGACY["ignore_arcs"] == "false"
    @test ARN.refit_config_dict(c)["ignore_arcs"] == "false"
    @test ARN.REFIT_STAGES == ("arc", "sky", "comb", "fpifib", "assemble")
    @test basename(ARN.refit_path("/r", "sky", "lco", 60656)) == "skyNightRefit_lco_60656.h5"
    @test basename(dirname(ARN.refit_path("/r", "sky", "lco", 60656))) == "sky"
    @test ARN.WAVE_BIT_NAMES[19] == "WAVE_RESERVED_18"
    @test ARN.WAVE_BIT_NAMES[40:42] == ["WAVE_NOARC_COMB_SKYSEED", "WAVE_NOARC_TEMPLATE_SKY", "WAVE_NOARC_SKY_ONLY"]
    @test :noarc_tmpl in ARN.REFIT_SKYSHIFT_ROUTES && :noarc_sky in ARN.REFIT_SKYSHIFT_ROUTES &&
          :noarc_nightave in ARN.REFIT_SKYSHIFT_ROUTES && !(:fpifib in ARN.REFIT_SKYSHIFT_ROUTES)
end

@testset "noarc ladder: legacy sky ingest masks class-A guards (#407)" begin
    mktempdir() do d
        fn = joinpath(d, "skyLinePeaks_apo_60000_0007_R_object.h5")
        M = fill(NaN, 3, 9, ARN.N_FIBERS)
        for l in 1:3, f in 1:ARN.N_FIBERS
            M[l, 1, f] = 300.0 * l
            M[l, 2, f] = 16000.0 + l
        end
        h5open(h -> (h["sky_line_mat_clean"] = M), fn, "w")
        x, w = ARN.ingest_skyLines_file(fn)
        @test all(isfinite, x) && all(isfinite, w)                 # pre-#407 file: unmasked
        h5open(h -> (h["line_class"] = [2, 1, 2]), fn, "r+")
        x, w = ARN.ingest_skyLines_file(fn)
        @test all(isnan, x[2, :]) && all(isnan, w[2, :]) && all(isfinite, x[[1, 3], :])
        h5open(fn, "r+") do h
            delete_object(h, "line_class"); h["line_class"] = [2, 2]
        end
        @test_throws DimensionMismatch ARN.ingest_skyLines_file(fn)
    end
end

# synthetic night: a smooth per-fiber solution lambda = P_f(a_c + x + delta_k), x = (p - 1024)/2048, chip offsets
# a_R/a_B as the legacy fit's LCO starting values, 8 OH lines per chip, exposures dithered by delta_k
const SYN_LIN = [16156.09, -579.085, -43.4533, 1.08138, 0.0239653]
const SYN_A = (-1.07295, 0.0, 1.076)
syn_lambda(f, c, p, dk) = evalpoly(SYN_A[c] + (p - 1024) / 2048 + dk, SYN_LIN .+ [0.002 * (f - 150); zeros(4)])
function syn_pixel(f, c, lam, dk)
    p = 1024.0
    for _ in 1:40
        h = 1e-3
        d = (syn_lambda(f, c, p + h, dk) - syn_lambda(f, c, p - h, dk)) / (2h)
        st = (syn_lambda(f, c, p, dk) - lam) / d
        p -= st
        abs(st) < 1e-10 && break
    end
    return p
end
function syn_write_peaks(dir, tele, mjd, e, dk; guard = false)
    for (c, chip) in enumerate(ARN.CHIP_LIST)
        lams = [syn_lambda(150, c, p, 0.0) for p in range(150.0, 1900.0, length = 8)]
        nl = guard ? 9 : 8
        M = fill(NaN, nl, 9, ARN.N_FIBERS)
        for f in 1:ARN.N_FIBERS
            for l in 1:8
                M[l, 1, f] = syn_pixel(f, c, lams[l], dk)
                M[l, 2, f] = lams[l]
                M[l, 9, f] = 100c + l
            end
            if guard   # a class-A guard: a real peak with a WRONG wavelength (would bias the fit if used)
                M[9, 1, f] = syn_pixel(f, c, lams[4] + 0.5, dk)
                M[9, 2, f] = lams[4] + 2.0
                M[9, 9, f] = 100c + 9
            end
        end
        h5open(joinpath(dir, "skyLinePeaks_$(tele)_$(mjd)_$(lpad(e, 4, '0'))_$(chip)_object.h5"), "w") do h
            h["sky_line_mat_clean"] = M
            guard && (h["line_class"] = [fill(2, 8); 1])
        end
    end
end

@testset "noarc ladder: within-night sky-only chain recovers a synthetic night" begin
    mktempdir() do d
        dks = [0.0, 0.3 / 2048, -0.2 / 2048]       # dithers of 0, +0.3, -0.2 px
        for (e, dk) in enumerate(dks)
            syn_write_peaks(d, "lco", 59000, e, dk; guard = true)
        end
        @test ARN.refit_sky_chain("lco", 59000, joinpath(d, "nope"), [1, 2, 3]) === nothing
        S = ARN.refit_sky_chain("lco", 59000, d, [1, 2, 3, 4]; workroot = d)   # exposure 4 has no peaks: skipped
        @test S !== nothing && S.expnum == [1, 2, 3] && all(S.ok)
        @test size(S.lin) == (ARN.N_FIBERS, 5) && size(S.dither) == (ARN.N_FIBERS, 2, 3)
        @test isempty(filter(x -> startswith(x, "skychain_"), readdir(d)))    # scratch directory removed
        # per-exposure solution vs truth (m/s), all fibers, every chip, 97 pixels
        err = Float64[]
        for (k, dk) in enumerate(dks), c in 1:3, f in 10:10:290
            px = collect(65.0:20.0:1985.0)
            lm = ARN.refit_sky_eval(S.lin, S.nl, view(S.dither, f, :, k), f, c, px)
            append!(err, 299792458.0 .* (lm .- [syn_lambda(f, c, p, dk) for p in px]) ./ lm)
        end
        @test median(abs.(err)) < 5.0
        @test maximum(abs.(err)) < 100.0
        # the night average (d = (0, 1)) sits at one measured frame: its per-exposure offset is the dither
        lm0 = ARN.refit_sky_eval(S.lin, S.nl, (0.0, 1.0), 150, 2, [1024.0])[1]
        @test abs(299792458.0 * (lm0 - syn_lambda(150, 2, 1024.0, 0.0)) / lm0) < 3000.0   # within ~0.6 px
    end
end

@testset "noarc ladder: sky stage status paths" begin
    mktempdir() do d
        cfg = ARN.RefitConfig()
        # no arc product at all -> a no-arc night; no exposures in the apred tree -> nosky
        @test ARN.refit_night_is_noarc(d, "lco", 59000, cfg)
        @test ARN.refit_night_sky("lco", 59000; apred_dir = joinpath(d, "apred"), outroot = d, cfg) == "nosky"
        @test ARN.refit_read_status(d, "sky", "lco", 59000) == "nosky"
        @test ARN.refit_night_is_noarc(d, "lco", 59000, ARN.RefitConfig(ignore_arcs = true))
    end
end

@testset "noarc ladder: template shift" begin
    # reference solution linear in pixel on every chip, template displacement zero: a pooled offset v (m/s) must
    # move the template model by v
    nfib = ARN.N_FIBERS
    ref_lin = repeat([16000.0 -300.0 0.0 0.0 0.0], nfib)
    ref_nl = repeat([-1.07 1.0 1.07 1.0], nfib)
    Tp = ARN.RefitTemplate("lco", "all", zeros(7, 3, nfib), ref_lin, ref_nl, "x")
    @test all(ARN.refit_template_ok(Tp))
    for v in (-2000.0, 35.0, 900.0)
        dl = ARN.refit_template_shift_px(Tp, [v, v, NaN])
        @test isfinite(dl[3]) && dl[3] == median(dl[1:2])
        for c in 1:2
            l0 = ARN.refit_template_noarc(Tp, zeros(3), 150, c, 1024.0)
            l1 = ARN.refit_template_noarc(Tp, dl, 150, c, 1024.0)
            @test abs(299792458.0 * (l1 - l0) / l0 - v) < 0.01 * abs(v) + 0.1
        end
    end
    @test all(isnan, ARN.refit_template_shift_px(Tp, [NaN, NaN, NaN]))
    T2 = copy(Tp.T); T2[3, 2, 7] = NaN
    @test ARN.refit_template_ok(ARN.RefitTemplate("lco", "all", T2, ref_lin, ref_nl, "x"))[7] == false
end
