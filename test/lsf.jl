# Tests for src/lsf.jl (get_lsf_matrix / build_sparse_lsf_mat).
#
# These tests require the regularized FPI LSF parameter file, which lives on
# ceph and is not shipped with the repo. They are skipped cleanly when the
# file is absent (e.g. in CI); set AR_LSF_PARAMS_FILE to point at a local copy
# to run them elsewhere.
#
# CONTRACT pinned here: the matrix returned by get_lsf_matrix is
# ROW-NORMALIZED — every row sums to exactly 1 (the LSF matrix applied to a
# ones vector returns a ones vector; see the norm_factor construction in
# build_sparse_lsf_mat). Downstream consumers (e.g. arMADGICS
# sample_starCont.jl) rely on this and apply no further normalization.

using SparseArrays: nnz, nonzeros, findnz

lsf_params_file = get(ENV, "AR_LSF_PARAMS_FILE",
    "/mnt/ceph/users/sdssv/work/asaydjari/2026_04_27/fpiLSFparams_REGULARIZED_apo_60861.h5")

# Vendored reference numbers, generated 2026-08-31 by evaluating kmckinnon's
# lsf_tools.jl (/mnt/ceph/users/sdssv/work/asaydjari/2026_04_27/lsf_tools.jl)
# fib_lam_lsf with force_mode_0 = true on
# fpiLSFparams_REGULARIZED_apo_60861.h5, fiberindx 150, with
# dlam_bin_edges = outwave_pad_bin_edges .- center_wave restricted to
# |dlam| < 5 Angstroms (outwave_pad_bin_edges as defined in get_lsf_matrix
# with n_pad = 50, fstep = 0).
REFERENCE_LSF_PROFILES = [
    (center_wave = 15500.0,
        ref = (n_edges_kept = 47, first_kept_edge_idx = 2042,
            sum_w = 1.0, argmax = 24,
            values = [22 => 0.13675000480521143,
                23 => 0.27154668784036134,
                24 => 0.29114951352717505,
                25 => 0.16727285113555457,
                26 => 0.046984710116369351])),
    (center_wave = 16000.0,
        ref = (n_edges_kept = 45, first_kept_edge_idx = 4341,
            sum_w = 1.0, argmax = 23,
            values = [21 => 0.13648092244138726,
                22 => 0.25317519084447332,
                23 => 0.27626874695375081,
                24 => 0.16910324617302966,
                25 => 0.052499931007110805])),
    (center_wave = 16500.0,
        ref = (n_edges_kept = 44, first_kept_edge_idx = 6569,
            sum_w = 1.0, argmax = 22,
            values = [20 => 0.093969687819554171,
                21 => 0.22356325287953491,
                22 => 0.29750237046257327,
                23 => 0.21913838844876926,
                24 => 0.084573640778203146]))
]

@testset "lsf" begin
    if !isfile(lsf_params_file)
        @warn "LSF parameter file not found; skipping LSF matrix tests" lsf_params_file
    else
        adjfibindx = 150 # apo fiber index 150

        # grid definitions matching get_lsf_matrix (fstep = 0)
        n_pad = 50
        outwave_log_step = 6.0e-6
        len_outwave = 8700
        outwave_log_min = 4.17825
        outwave = 10 .^ range(start = outwave_log_min, step = outwave_log_step,
            length = len_outwave)
        outwave_pad_bin_edges = 10 .^ range(
            start = outwave_log_min - n_pad * outwave_log_step - 0.5 * outwave_log_step,
            step = outwave_log_step, length = len_outwave + 2 * n_pad + 1)
        modelwave = 15000:(1 / 100):17000

        @testset "fib_lam_lsf vs lsf_tools.jl reference evaluation" begin
            # Reference numbers generated 2026-08-31 by evaluating kmckinnon's
            # /mnt/ceph/users/sdssv/work/asaydjari/2026_04_27/lsf_tools.jl
            # (fib_lam_lsf, force_mode_0 = true) on this same parameter file,
            # apo fiberindx 150, with bin edges = outwave_pad_bin_edges .- center
            # restricted to |dlam| < 5 Angstroms. See REFERENCE_VALUES below.
            params = ApogeeReduction.get_fib_lsf_params(150, lsf_params_file)
            for (; center_wave, ref) in REFERENCE_LSF_PROFILES
                dlam = outwave_pad_bin_edges .- center_wave
                keep = findall(abs.(dlam) .< 5)
                w = ApogeeReduction.fib_lam_lsf(center_wave, dlam[keep], params...;
                    force_mode_0 = true)
                @test length(keep) == ref.n_edges_kept
                @test keep[1] == ref.first_kept_edge_idx
                @test isapprox(sum(w), ref.sum_w; rtol = 1e-8)
                @test argmax(w) == ref.argmax
                for (k, v) in ref.values
                    @test isapprox(w[k], v; rtol = 1e-8)
                end
            end
        end

        @testset "build_sparse_lsf_mat assembly vs direct evaluation" begin
            # Independent dense reimplementation of the assembly +
            # normalization chain on a small problem, using only fib_lam_lsf.
            params = ApogeeReduction.get_fib_lsf_params(150, lsf_params_file)
            small_edges = 10 .^ range(start = log10(15990.0), step = outwave_log_step,
                length = 301) # 300 pixels around 16000 A
            small_modelwave = 15995.0:0.01:16005.0
            min_lsf_weight = 1e-15

            expected = zeros(300, length(small_modelwave))
            for (j, center_wave) in enumerate(small_modelwave)
                dlam = small_edges .- center_wave
                keep = findall(abs.(dlam) .< 5)
                length(keep) < 2 && continue
                w = ApogeeReduction.fib_lam_lsf(center_wave, dlam[keep], params...;
                    force_mode_0 = true)
                w[w .<= min_lsf_weight] .= 0
                s = sum(w)
                iszero(s) && continue
                w ./= s
                expected[keep[begin:(end - 1)], j] .= w
            end
            rowsums = vec(sum(expected, dims = 2))
            for i in 1:300
                if rowsums[i] > 0
                    expected[i, :] ./= rowsums[i]
                end
            end

            M_small = ApogeeReduction.build_sparse_lsf_mat(
                small_modelwave, small_edges, 150, lsf_params_file;
                force_mode_0 = true, min_lsf_weight = min_lsf_weight, n_pad = 0)
            @test size(M_small) == size(expected)
            @test isapprox(Matrix(M_small), expected; rtol = 1e-8)
        end

        @testset "get_lsf_matrix full-matrix contract" begin
            M = ApogeeReduction.get_lsf_matrix(adjfibindx, lsf_params_file)

            # shape: 8700 output pixels x 200001 model wavelengths
            @test size(M) == (len_outwave, length(modelwave))

            # non-negativity (strictly positive stored values after the
            # min_lsf_weight threshold)
            @test all(>(0), nonzeros(M))

            # ROW-NORMALIZATION contract: every row sums to 1
            row_sums = vec(sum(M, dims = 2))
            @test all(isapprox.(row_sums, 1.0; atol = 1e-8))
            # equivalently: LSF applied to a constant vector is that constant
            @test maximum(abs.(M * ones(length(modelwave)) .- 1.0)) < 1e-8

            # sparsity pattern sanity: every column's support lies within the
            # +-5 A evaluation window of its modelwave center, and the matrix
            # is genuinely sparse
            @test 0 < nnz(M) < 0.01 * length(M)
            rows, cols, _ = findnz(M)
            @test maximum(abs.(outwave[rows] .- modelwave[cols])) < 5.0
            # every output pixel row is covered
            @test all(>(0), row_sums) # implied by approx 1, but explicit
            @test length(unique(rows)) == len_outwave
        end
    end
end
