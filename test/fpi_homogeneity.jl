using HDF5

# Build a minimal corpus-almanac fixture (raw/ layout): per night an
# `exposures` table (column vectors, as almanac writes them) and one `fibers`
# table per configuration. `configs` maps config_id => FPI pair (a vector of
# bonus fiber_ids; empty vector = configuration with no bonus rows).
function _write_alm_night!(f, tele, mjd; configs)
    image_type = String[]
    config_id = Int[]
    exposure = Int[]
    # a calibration exposure that must be ignored
    push!(image_type, "domeflat")
    push!(config_id, 0)
    push!(exposure, 1)
    for (i, (cid, _)) in enumerate(configs)
        push!(image_type, "object")
        push!(config_id, cid)
        push!(exposure, i + 1)
    end
    base = "raw/$(tele)/$(mjd)"
    write(f, "$(base)/exposures/image_type", image_type)
    write(f, "$(base)/exposures/config_id", config_id)
    write(f, "$(base)/exposures/exposure", exposure)
    for (cid, pair) in configs
        fiber_id = Int[10, 20, 30]
        category = String["science", "sky_apogee", "standard_apogee"]
        for fid in pair
            push!(fiber_id, fid)
            push!(category, "bonus")
        end
        write(f, "$(base)/fibers/$(cid)/fiber_id", fiber_id)
        write(f, "$(base)/fibers/$(cid)/category", category)
    end
end

@testset "fpi_homogeneity" begin
    dir = mktempdir()

    @testset "homogeneous corpus" begin
        fname = joinpath(dir, "alm_homogeneous.h5")
        h5open(fname, "w") do f
            _write_alm_night!(f, "apo", 60000, configs = [1101 => [75, 225], 1102 => [75, 225]])
            _write_alm_night!(f, "apo", 60001, configs = [1103 => [75, 225]])
            # plate-era night: out of scope, must be skipped even with an object config
            _write_alm_night!(f, "apo", 59000, configs = [901 => [1, 2]])
        end
        results = ApogeeReduction.survey_fpi_homogeneity(fname)
        @test length(results) == 1   # no lco group in the fixture
        res = results[1]
        @test res.tele == "apo"
        @test res.dominant == [75, 225]
        @test res.n_configs == 3
        @test res.n_nights == 2
        @test isempty(res.deviations)
        @test isempty(res.missing_bonus)
        report = sprint(io -> ApogeeReduction.report_fpi_homogeneity(io, results))
        @test occursin("HOMOGENEOUS", report)
        @test !occursin("[NEW]", report)
    end

    @testset "deviating night is flagged, not rejected" begin
        fname = joinpath(dir, "alm_deviating.h5")
        h5open(fname, "w") do f
            _write_alm_night!(f, "lco", 60000,
                configs = [2101 => [82, 213], 2102 => [82, 213], 2103 => [82, 213]])
            _write_alm_night!(f, "lco", 59820, configs = [2001 => [142, 153]])
            _write_alm_night!(f, "lco", 60002, configs = [2104 => [82, 213]])
            # a config with no bonus rows is reported separately, never a deviation
            _write_alm_night!(f, "lco", 60003, configs = [2105 => Int[]])
        end
        results = ApogeeReduction.survey_fpi_homogeneity(fname)
        res = results[1]
        @test res.tele == "lco"
        @test res.dominant == [82, 213]
        @test res.n_configs == 5
        @test length(res.deviations) == 1
        @test res.deviations[1].mjd == 59820
        @test res.deviations[1].config_id == 2001
        @test res.deviations[1].pair == [142, 153]
        @test res.missing_bonus == [(60003, 2105)]

        # without a known-epochs file the deviation is NEW
        summary = ApogeeReduction.report_fpi_homogeneity(devnull, results)
        @test summary == (n_flagged_nights = 1, n_known = 0, n_new = 1)
        report = sprint(io -> ApogeeReduction.report_fpi_homogeneity(io, results))
        @test occursin("[NEW]", report)
        @test occursin("142/153", report)
        @test occursin("ACTION REQUIRED", report)

        # with the epoch annotated as known, it is flagged-but-known
        epochs_file = joinpath(dir, "known_epochs.txt")
        write(epochs_file,
            """
            # comment line
            lco 59810 59850 re-fibering epoch, kept in DR21
            """)
        known = ApogeeReduction.parse_known_epochs(epochs_file)
        summary = ApogeeReduction.report_fpi_homogeneity(devnull, results, known_epochs = known)
        @test summary == (n_flagged_nights = 1, n_known = 1, n_new = 0)
        report = sprint(io -> ApogeeReduction.report_fpi_homogeneity(
            io, results, known_epochs = known))
        @test occursin("[KNOWN]", report)
        @test !occursin("[NEW]", report)
        @test occursin("re-fibering epoch", report)
    end

    @testset "parse_known_epochs" begin
        @test isempty(ApogeeReduction.parse_known_epochs(""))
        epochs_file = joinpath(dir, "epochs.txt")
        write(epochs_file,
            """
            # full-line comment

            lco 59810 59850 re-fibering epoch
            apo 60100 60101   # inline comment, no label
            """)
        epochs = ApogeeReduction.parse_known_epochs(epochs_file)
        @test length(epochs) == 2
        @test epochs[1] == (tele = "lco", mjd_lo = 59810, mjd_hi = 59850,
            label = "re-fibering epoch")
        @test epochs[2].tele == "apo"
        @test ApogeeReduction.known_epoch_label(epochs, "lco", 59820) == "re-fibering epoch"
        # a labelless epoch still annotates
        @test ApogeeReduction.known_epoch_label(epochs, "apo", 60100) !== nothing
        @test isnothing(ApogeeReduction.known_epoch_label(epochs, "lco", 59809))
        @test isnothing(ApogeeReduction.known_epoch_label(epochs, "apo", 59820))

        bad_file = joinpath(dir, "epochs_bad.txt")
        write(bad_file, "lco 59810\n")
        @test_throws ErrorException ApogeeReduction.parse_known_epochs(bad_file)
        @test_throws ErrorException ApogeeReduction.parse_known_epochs(
            joinpath(dir, "does_not_exist.txt"))
    end

    # the shipped corpus known-epochs file must stay parseable
    @testset "shipped metadata/fpi_known_epochs.txt parses" begin
        shipped = joinpath(@__DIR__, "..", "metadata", "fpi_known_epochs.txt")
        epochs = ApogeeReduction.parse_known_epochs(shipped)
        @test any(e -> e.tele == "lco" && e.mjd_lo == 59810 && e.mjd_hi == 59850, epochs)
    end
end
