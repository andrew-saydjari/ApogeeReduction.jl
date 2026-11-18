using HDF5

# The trace-template calibration must be rebuildable FROM SCRATCH: the builder reads only 2D quartz
# flats, never a trace, a traceMain or a template file. This renders a few synthetic all-lit quartz
# flats (Gauss-Hermite profiles from the package's profile tables, known centers, a dead fiber, a
# different rigid shift every night), runs build_trace_template.jl and pack_trace_templates.jl as
# subprocesses with no template file anywhere (AR_TRACE_TEMPLATE_CAL pointing at a path that does
# not exist), and checks that the packed template recovers the true slit up to the gauge.
@testset "trace templates bootstrap (no template file)" begin
    tdir = joinpath(pkgdir(ApogeeReduction), "scripts", "cal", "trace_templates")
    tele, chip = "apo", "G"
    pattern = 0.25 .* sin.((1:N_FIBERS) ./ 5)                 # per-fiber slit pattern (px)
    truth(f, ui, d) = 40.0 + 6.62 * (f - 1) + pattern[f] + 1.5 * ui^2 + 0.4 * ui + d
    # the flats are rendered in a subprocess in the package environment (TemplateFit's deps)
    render = """
    include(joinpath($(repr(tdir)), "TemplateFit.jl"))
    using .TemplateFit, HDF5, Random
    const TF = TemplateFit
    root, tele, chip = ARGS[1], "apo", "G"
    P = TF.load_profiles(tele, chip)
    pattern = 0.25 .* sin.((1:300) ./ 5)
    truth(f, ui, d) = 40.0 + 6.62 * (f - 1) + pattern[f] + 1.5 * ui^2 + 0.4 * ui + d
    thr = 0.6 .+ 0.4 .* cos.((1:300) ./ 13) .^ 2
    thr[90] = 0.0                                            # dead in the whole epoch
    rng = MersenneTwister(5)
    buf = zeros(11)
    for (mjd, d) in [(60001, 0.0), (60002, 0.35), (60003, -0.25), (60004, 0.6)]
        img = zeros(2048, 2048)
        for x in 1:2048, f in 1:300
            thr[f] == 0 && continue
            mu = truth(f, (x - 1024) / 1024, d); yc = round(Int, mu)
            y0, y1 = yc - 5, yc + 5
            (y0 < 1 || y1 > 2048) && continue
            TF.unitprof!(buf, P, f, mu, 0.9, y0, y1)
            @views img[x, y0:y1] .+= 3000 * thr[f] .* buf
        end
        img .+= 5 .* randn(rng, size(img))
        dd = joinpath(root, "apred", string(mjd)); mkpath(dd)
        h5open(joinpath(dd, "ar2Dcal_\$(tele)_\$(mjd)_0007_\$(chip)_quartzflat.h5"), "w") do f
            f["dimage"] = img
            f["ivarimage"] = fill(1 / 25, size(img))
            f["pix_bitmask"] = zeros(Int64, size(img))
        end
    end
    """
    mktempdir() do dir
        root = joinpath(dir, "reduction")
        env = ("AR_TRACE_TEMPLATE_CAL" => joinpath(dir, "nonexistent", "tracetemplates.h5"),
            "JULIA_NUM_THREADS" => "1", "OPENBLAS_NUM_THREADS" => "1",
            # Pkg.test points the load path at its sandbox; the scripts run in the package's own project
            "JULIA_LOAD_PATH" => "@:@stdlib", "JULIA_PROJECT" => pkgdir(ApogeeReduction))
        @test !isfile(env[1][2])
        jl = `$(Base.julia_cmd()) --project=$(pkgdir(ApogeeReduction))`
        rf = joinpath(dir, "render.jl"); write(rf, render)
        logr = joinpath(dir, "render.log")
        ok = success(pipeline(addenv(`$jl $rf $root`, env...); stdout = logr, stderr = logr))
        ok || println(read(logr, String))
        @test ok
        units = joinpath(dir, "units"); mkpath(units)
        logf = joinpath(dir, "build.log")
        ok = success(pipeline(addenv(`$jl $(joinpath(tdir, "build_trace_template.jl")) --root $root
            --tele $tele --chip $chip --epoch e1 --mjd_start 60000 --mjd_end 60010 --outdir $units
            --maxit 3 --validate false`, env...); stdout = logf, stderr = logf))
        ok || println(read(logf, String))
        @test ok
        @test isfile(joinpath(units, "template_apo_G_e1.h5"))
        epochs = joinpath(dir, "epochs.csv")
        write(epochs, "tele,epoch,mjd_start,mjd_end\napo,e1,60000,60010\n")
        calp = joinpath(dir, "tt_bootstrap.h5")
        logp = joinpath(dir, "pack.log")
        ok = success(pipeline(addenv(`$jl $(joinpath(tdir, "pack_trace_templates.jl")) --epochs $epochs
            --units_dir $units --out $calp --version bootstrap_test --allow_incomplete true`, env...);
            stdout = logp, stderr = logp))
        ok || println(read(logp, String))
        @test ok
        cal = ApogeeReduction.load_trace_template_cal(calp)
        @test !cal.complete                                    # R and B were not built
        @test startswith(cal.provenance["needs_trace_template_file"], "no")
        t = ApogeeReduction.trace_template(cal, "apo", "G", 60005)
        @test !isnothing(t) && length(t.input_mjd) == 4
        @test findall(.!t.lit) == [90]
        @test ApogeeReduction.trace_template_provisional(cal, "apo", "G", 60500)   # open-ended, after the data
        # the template is the true slit up to the gauge (a rigid shift and scale per band)
        T0 = [truth(f, (x - 1024) / 1024, 0.0) for x in t.xc, f in 1:N_FIBERS]
        res = Float64[]
        for b in eachindex(t.xc)
            l = t.lit
            X = hcat(ones(count(l)), T0[b, l] .- 1024)
            β = X \ (t.center[b, l] .- T0[b, l])
            append!(res, t.center[b, l] .- T0[b, l] .- X * β)
        end
        @test maximum(abs.(res)) < 0.05
        @test sqrt(sum(abs2, res) / length(res)) < 0.01
        # the dead fiber sits in its slot (interpolated between its neighbours)
        @test maximum(abs.(t.center[:, 90] .- (t.center[:, 89] .+ t.center[:, 91]) ./ 2)) < 0.3
    end
end
