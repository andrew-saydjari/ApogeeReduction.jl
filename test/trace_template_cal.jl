using HDF5

# Synthetic trace-template calibration: smooth, increasing fiber centers, one template per
# (telescope, chip, epoch).
function synthetic_trace_template(tele, chip, ep::ApogeeReduction.TraceEpoch; shift = 0.0, nin = 5)
    nb = 32
    xc = collect(range(32.5, 2016.5, length = nb))
    u = (xc .- 1024) ./ 1024
    base = 30.0 .+ 6.6 .* (0:(N_FIBERS - 1)) .+ shift
    cp = zeros(5, N_FIBERS); cp[1, :] .= base; cp[2, :] .= 0.8; cp[3, :] .= -0.3
    wp = zeros(3, N_FIBERS); wp[1, :] .= 0.9; wp[3, :] .= 0.05
    center = [sum(cp[p, f] * ui^(p - 1) for p in 1:5) for ui in u, f in 1:N_FIBERS]
    width = [sum(wp[p, f] * ui^(p - 1) for p in 1:3) for ui in u, f in 1:N_FIBERS]
    lit = trues(N_FIBERS); lit[[28, 90]] .= false
    m1 = ep.open_ended ? ep.mjd_start + 20 : ep.mjd_end
    inm = round.(Int, range(ep.mjd_start, m1, length = nin))
    ApogeeReduction.TraceTemplate(tele, chip, ep, xc, center, width, cp, wp, ones(N_FIBERS), lit,
        inm, fill(7, nin), 0.1 .* (1:nin), fill(1e-4, nin), Dict("builder" => "synthetic"))
end

function synthetic_trace_template_cal(path; version = "tt_v1", drop = nothing, open_last = true,
        shift_apo_G_e2 = 0.0)
    E = ApogeeReduction.TraceEpoch
    epochs = Dict(
        "apo" => [E("apo", "apoE1", 57618, 59999, false),
            E("apo", "apoE2", 60000, open_last ? ApogeeReduction.TRACE_EPOCH_OPEN_END : 61230, open_last)],
        "lco" => [E("lco", "lcoE1", 57800, open_last ? ApogeeReduction.TRACE_EPOCH_OPEN_END : 61230, open_last)])
    tm = ApogeeReduction.TraceTemplate[]
    for (tele, eps) in epochs, ep in eps, chip in CHIP_LIST
        (tele, chip, ep.name) == drop && continue
        sh = (tele, chip, ep.name) == ("apo", "G", "apoE2") ? shift_apo_G_e2 : 0.0
        push!(tm, synthetic_trace_template(tele, chip, ep; shift = sh))
    end
    ApogeeReduction.write_trace_template_cal(path, epochs, tm,
        Dict("code_commit" => "abc", "epoch_table" => "synthetic");
        calibration_version = version, allow_incomplete = !isnothing(drop))
end

@testset "trace_template_cal" begin
    AR = ApogeeReduction
    @testset "round trip, lookup, identity" begin
        mktempdir() do dir
            p = joinpath(dir, "tt.h5")
            synthetic_trace_template_cal(p)
            cal = AR.load_trace_template_cal(p; require_complete = true)
            @test cal.complete
            @test cal.calibration_version == "tt_v1"
            @test length(cal.templates) == 9
            @test cal.provenance["missing_units"] == "none"
            t = AR.trace_template(cal, "apo", "G", 60500)
            @test t.epoch == "apoE2" && t.chip == "G"
            @test size(t.center) == (32, N_FIBERS) && size(t.center_poly) == (5, N_FIBERS)
            @test t.reg_d_median ≈ 0.3
            @test AR.trace_template(cal, "apo", "G", 59999).epoch == "apoE1"
            @test isnothing(AR.trace_template(cal, "apo", "G", 57000))    # before the first epoch
            @test isnothing(AR.trace_template(cal, "xyz", "G", 60000))
            @test isnothing(AR.trace_template(nothing, "apo", "G", 60000))
            # polynomial evaluation reproduces the band table
            @test maximum(abs.(AR.template_centers(t, t.xc) .- t.center)) < 1e-9
            @test maximum(abs.(AR.template_widths(t, t.xc) .- t.width)) < 1e-9
            # provisional: open-ended last epoch, after the newest input flat
            @test !AR.trace_template_provisional(cal, "apo", "G", 60010)
            @test AR.trace_template_provisional(cal, "apo", "G", 61500)
            @test !AR.trace_template_provisional(cal, "apo", "G", 59000)
            @test !AR.trace_template_provisional(nothing, "apo", "G", 61500)
            id = AR.trace_template_id(cal, "apo", "G", 60500)
            @test startswith(id, "apo_G_apoE2:") && length(split(id, ":")[2]) == 16
            @test AR.trace_template_id(nothing, "apo", "G", 60500) == "none"
            @test AR.trace_template_id(cal, "apo", "G", 50000) == "none"
            @test AR.trace_template_cal_id(cal) == "tt_v1:" * cal.fingerprint
            md = AR.trace_template_metadata(cal, "apo", "G", 61500)
            @test md["trace_template_id"] == AR.trace_template_id(cal, "apo", "G", 61500)
            @test md["trace_template_epoch"] == "apoE2" && md["trace_template_provisional"] == true
            @test AR.trace_template_metadata(nothing, "apo", "G", 1)["trace_template_id"] == "none"
            # checkpoint decision
            @test AR.template_needs_recompute(Dict{String, Any}(), cal, "apo", "G", 60500)
            @test !AR.template_needs_recompute(Dict{String, Any}(), nothing, "apo", "G", 60500)
            @test !AR.template_needs_recompute(Dict("trace_template_id" => id), cal, "apo", "G", 60500)
        end
    end

    @testset "fingerprints: content, not bytes; per unit" begin
        mktempdir() do dir
            p1 = joinpath(dir, "a.h5"); p2 = joinpath(dir, "b.h5"); p3 = joinpath(dir, "c.h5")
            synthetic_trace_template_cal(p1)
            synthetic_trace_template_cal(p2; version = "tt_v2")            # same content, new version
            synthetic_trace_template_cal(p3; shift_apo_G_e2 = 0.01)       # one unit changed
            c1, c2, c3 = AR.load_trace_template_cal.((p1, p2, p3))
            @test c1.fingerprint == c2.fingerprint
            @test AR.trace_template_id(c1, "apo", "G", 60500) == AR.trace_template_id(c2, "apo", "G", 60500)
            @test c3.fingerprint != c1.fingerprint
            @test AR.trace_template_id(c3, "apo", "G", 60500) != AR.trace_template_id(c1, "apo", "G", 60500)
            # other units keep their identity: nights there need no recompute
            @test AR.trace_template_id(c3, "apo", "R", 60500) == AR.trace_template_id(c1, "apo", "R", 60500)
            @test AR.trace_template_id(c3, "apo", "G", 59000) == AR.trace_template_id(c1, "apo", "G", 59000)
            # tampering with a stored array is caught at load
            p4 = joinpath(dir, "d.h5"); cp(p1, p4)
            h5open(p4, "r+") do f
                c = read(f["templates/apo_G_apoE1/center"]); c[5, 5] += 0.5
                delete_object(f, "templates/apo_G_apoE1/center")
                f["templates/apo_G_apoE1/center"] = c
            end
            @test_throws ArgumentError AR.load_trace_template_cal(p4)
        end
    end

    @testset "coverage validation and incomplete files" begin
        mktempdir() do dir
            E = AR.TraceEpoch
            # missing unit: refused unless allow_incomplete; listed; lookup gives nothing there
            @test_throws ArgumentError begin
                epochs = Dict("apo" => [E("apo", "apoE1", 57618, 59999, false)])
                AR.write_trace_template_cal(joinpath(dir, "x.h5"), epochs,
                    [synthetic_trace_template("apo", "G", epochs["apo"][1])], Dict{String, String}();
                    calibration_version = "x")
            end
            p = joinpath(dir, "inc.h5")
            synthetic_trace_template_cal(p; drop = ("lco", "B", "lcoE1"))
            cal = AR.load_trace_template_cal(p)
            @test !cal.complete
            @test cal.provenance["missing_units"] == "lco_B_lcoE1"
            @test isnothing(AR.trace_template(cal, "lco", "B", 60000))
            @test AR.trace_template_id(cal, "lco", "B", 60000) == "none"
            @test !isnothing(AR.trace_template(cal, "lco", "R", 60000))
            @test_throws ArgumentError AR.load_trace_template_cal(p; require_complete = true)
            # gaps, overlaps, two open-ended epochs
            @test_throws ArgumentError AR.check_trace_epochs("apo",
                [E("apo", "a", 1, 10, false), E("apo", "b", 12, 20, false)])
            @test_throws ArgumentError AR.check_trace_epochs("apo",
                [E("apo", "a", 1, 10, false), E("apo", "b", 10, 20, false)])
            @test_throws ArgumentError AR.check_trace_epochs("apo",
                [E("apo", "a", 1, 10, true), E("apo", "b", 11, 20, true)])
            @test isnothing(AR.check_trace_epochs("apo", [E("apo", "a", 1, 10, false), E("apo", "b", 11, 20, true)]))
            # template range must equal its epoch's
            ep = E("apo", "a", 100, 200, false)
            t = synthetic_trace_template("apo", "G", ep)
            @test_throws ArgumentError AR.check_trace_template(t, E("apo", "a", 100, 300, false))
            # missing file / not a calibration file: hard errors
            @test_throws ArgumentError AR.load_trace_template_cal(joinpath(dir, "nope.h5"))
            q = joinpath(dir, "junk.h5"); h5open(q, "w") do f; f["x"] = 1; end
            @test_throws ArgumentError AR.load_trace_template_cal(q)
        end
    end

    @testset "builder unit file -> template" begin
        mktempdir() do dir
            ep = AR.TraceEpoch("apo", "apoE5", 59534, 59815, false)
            t0 = synthetic_trace_template("apo", "R", ep)
            p = joinpath(dir, "template_apo_R_apoE5.h5")
            h5open(p, "w") do f
                f["xc"] = t0.xc; f["center"] = t0.center; f["width"] = t0.width; f["height_rel"] = t0.height_rel
                f["center_poly_u"] = t0.center_poly; f["width_poly_u"] = t0.width_poly
                f["lit_in_epoch"] = Int.(t0.lit); f["input_mjd"] = t0.input_mjd; f["input_exp"] = t0.input_exp
                a = HDF5.attributes(f)
                a["tele"] = "apo"; a["chip"] = "R"; a["epoch"] = "apoE5"; a["mjd_start"] = 59534; a["mjd_end"] = 59815
                a["n_alllit"] = 5
            end
            logp = joinpath(dir, "template_apo_R_apoE5_log.csv")
            open(logp, "w") do io
                println(io, "stage,mjd,exp,nslots,nmissing,alllit,fit_rms,d_center,s_center,k,resid_rms,note")
                println(io, "census,$(t0.input_mjd[1]),7,300,2,true,NaN,NaN,NaN,0,NaN,28;90")
                for m in t0.input_mjd
                    println(io, "final,$(m),7,300,2,true,0.009,0.25,0.0001,0,0.02,")
                end
            end
            t = AR.read_template_unit_file(p; epoch = ep, log_csv = logp)
            @test t.center == t0.center && t.lit == t0.lit
            @test all(t.input_d_center .== 0.25) && t.reg_d_median == 0.25
            @test t.provenance["n_alllit"] == "5"
            # served range wider than the builder's only for an open-ended epoch
            @test_throws ArgumentError AR.read_template_unit_file(p;
                epoch = AR.TraceEpoch("apo", "apoE5", 59534, 99999, false))
            t2 = AR.read_template_unit_file(p; epoch = AR.TraceEpoch("apo", "apoE5", 59534, 99999, true))
            @test t2.mjd_end == 99999
        end
    end
end
