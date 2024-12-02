using StatsBase: iqr

bad_dark_pix_bits = 2^2 + 2^4 + 2^5;
bad_flat_pix_bits = 2^6;
# most multiread CR detections are bad for other reasons
bad_cr_pix_bits = 2^7 + 2^8; # could probably drop 2^7 at least in the future (happily correct 1 read CRs)
bad_chi2_pix_bits = 2^9;
bad_pix_bits = bad_dark_pix_bits + bad_flat_pix_bits + bad_cr_pix_bits + bad_chi2_pix_bits;

function initalize_git(git_dir)
    git_commit = LibGit2.head(git_dir)
    git_repo = LibGit2.GitRepo(git_dir)
    git_head = LibGit2.head(git_repo)
    git_branch = LibGit2.shortname(git_head)
    println("Running on branch: $git_branch, commit: $git_commit")
    flush(stdout)
    return git_branch, git_commit
end

"""
True if x is NaN or zero.
"""
function isnanorzero(x)
    return isnan(x) || iszero(x)
end

"""
Mean of the non-NaN, non-zero elements of x.
"""
nanzeromean(x) =
    if all(isnanorzero, x)
        NaN
    else
        mean(filter(!isnanorzero, x))
    end
nanzeromean(x, y) = mapslices(nanzeromean, x, dims = y)

"""
Sum of the non-NaN elements of x.
"""
nansum(x) = sum(filter(!isnan, x))
nansum(x, y) = mapslices(nansum, x, dims = y)

"""
Sum of the non-NaN, non-zero elements of x.
"""
nanzerosum(x) =
    if all(isnanorzero, x)
        NaN
    else
        sum(filter(!isnanorzero, x))
    end
nanzerosum(x, y) = mapslices(nanzerosum, x, dims = y)

"""
Median of the non-NaN, non-zero elements of x.
"""
nanzeromedian(x) =
    if all(isnanorzero, x)
        NaN
    else
        median(filter(!isnanorzero, x))
    end
nanzeromedian(x, y) = mapslices(nanzeromedian, x, dims = y)

"""
Median of the non-NaN elements of x.
"""
nanmedian(x) =
    if all(isnan, x)
        NaN
    else
        median(filter(!isnan, x))
    end
nanmedian(x, y) = mapslices(nanmedian, x, dims = y)

"""
The scaled interquartile range of the non-NaN, non-zero elements of x.
This gives σ for a Gaussian.
"""
nanzeroiqr(x) =
    if all(isnanorzero, x)
        NaN
    else
        iqr(filter(!isnanorzero, x)) / 1.34896
    end
nanzeroiqr(x, y) = mapslices(nanzeroiqr, x, dims = y)

function grow_msk2d(msk; rad = 1)
    (sx, sy) = size(msk)
    msknew = zeros(Bool, (sx, sy))
    for i in 1:sx
        for j in 1:sy
            srngx = maximum([1, i - rad]):minimum([sx, i + rad])
            srngy = maximum([1, j - rad]):minimum([sy, j + rad])
            msknew[i, j] = any(msk[srngx, srngy])
        end
    end
    return msknew
end

"""
Generate a design matrix for 2D sines + cosines + X + Y basis
"""
function gen_design_mat(nx, ny, fx, fy, X, Y)
    n_points = nx * ny
    n_basis = 2 * (fx + 1) * (fy + 1) - 2 + 3 #+ ny  # Total number of basis functions (sin + cos + constant)
    design_matrix = zeros(Float64, n_points, n_basis)

    col = 1
    for kx in 0:fx
        for ky in 0:fy
            if !((kx == 0) & (ky == 0))
                # Sine basis functions
                design_matrix[:, col] = vec(sin.(kx .* X + ky .* Y))
                col += 1

                # Cosine basis functions
                design_matrix[:, col] = vec(cos.(kx .* X + ky .* Y))
                col += 1
            end
        end
    end
    design_matrix[:, col] .= 1
    col += 1
    design_matrix[:, col] .= vec(X ./ maximum(X))
    col += 1
    design_matrix[:, col] .= vec(Y ./ maximum(Y))
    col += 1
    return design_matrix
end
