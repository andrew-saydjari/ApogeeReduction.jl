using StatsBase: iqr
using Jackknife
using JLD2

# used to record git branch and commit in safe_jldsave
using LibGit2
function initalize_git(git_dir)
    git_commit = LibGit2.head(git_dir)
    git_repo = LibGit2.GitRepo(git_dir)
    git_head = LibGit2.head(git_repo)
    git_branch = LibGit2.shortname(git_head)
    println("Running on branch: $git_branch, commit: $git_commit")
    flush(stdout)
    return git_branch, git_commit
end
# this will be reexecuted each time utils.jl is included somewhere, this is not inherently a problem
# but it is a symptom of the fact that the include situation is a bit tangled
git_branch, git_commit = initalize_git("./")

# ENV["SLACK_CHANNEL"] = "C08B7FKMP16" #apogee-reduction-jl
if !haskey(ENV, "SLACK_CHANNEL")
    ENV["SLACK_CHANNEL"] = "C07KQ7BJY5P" #apogee-reduction-jl-dev
end

# bad_dark_pix_bits = 2^2 + 2^4 #+ 2^5; temporarily remove 2^5 from badlist for now
bad_dark_pix_bits = 2^4
bad_flat_pix_bits = 2^6;
# most multiread CR detections are bad for other reasons
bad_cr_pix_bits = 2^7 + 2^8; # could probably drop 2^7 at least in the future (happily correct 1 read CRs)
bad_chi2_pix_bits = 2^9;
bad_pix_bits = bad_dark_pix_bits + bad_flat_pix_bits + bad_cr_pix_bits + bad_chi2_pix_bits;

# flags for 1d flux extraction
bad_1d_failed_extract = 2^10;
bad_1d_no_good_pix = 2^11;
bad_1d_neff = 2^12;

function isnanorzero(x)
    return isnan(x) | iszero(x)
end

nanzeromean(x) =
    if all(isnanorzero, x)
        NaN
    else
        mean(filter(!isnanorzero, x))
    end
nanzeromean(x, y) = mapslices(nanzeromean, x, dims = y)

nansum(x) = sum(filter(!isnan, x))
nansum(x, y) = mapslices(nansum, x, dims = y)

nanzerosum(x) =
    if all(isnanorzero, x)
        NaN
    else
        sum(filter(!isnanorzero, x))
    end
nanzerosum(x, y) = mapslices(nanzerosum, x, dims = y)

nanzeromedian(x) =
    if all(isnanorzero, x)
        NaN
    else
        median(filter(!isnanorzero, x))
    end
nanzeromedian(x, y) = mapslices(nanzeromedian, x, dims = y)

#@Kevin can we remove?
nanmedian(x) =
    if all(isnan, x)
        NaN
    else
        median(filter(!isnan, x))
    end
nanmedian(x, y) = mapslices(nanmedian, x, dims = y)

"Returns 1 for unit normal"
nanzeroiqr(x) =
    if all(isnanorzero, x)
        NaN
    else
        iqr(filter(!isnanorzero, x)) / 1.34896
    end
nanzeroiqr(x, y) = mapslices(nanzeroiqr, x, dims = y)

function log10n(x)
    if x <= 0
        return NaN
    else
        return log10(x)
    end
end

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

function getChipIndx(chip)
    if chip == "a"
        return 1
    elseif chip == "b"
        return 2
    elseif chip == "c"
        return 3
    end
end

function convert_to_int(x)
    if isnan(x)
        return 0
    else
        return convert(Int, x)
    end
end

function get_last_ind(x)
    if isnothing(x)
        return []
    else
        return x[end, :]
    end
end

function jack_std(x)
    y = filter(!isnanorzero, x)
    ly = length(y)
    if ly > 2
        return Jackknife.estimate(std, y)
    else
        return Inf
    end
end

normal_pdf(Δ, σ) = exp(-0.5 * Δ^2 / σ^2) / √(2π) / σ

"""
This function is a wrapper around JLD2.jldsave with a couple of extra features:
- It checks if the types of the values to be saved will result in a hard-to-read HDF5 file and warns if so.
- It converts BitArrays to Array{Bool} if necessary. This means that the saved data will be 8x
  larger (Bools are 1 byte), even when read back into Julia.
- It records the git branch and commit in the saved file.
"""
function safe_jldsave(filename; kwargs...)
    to_save = Dict{Symbol, Any}()
    for (k, v) in kwargs
        # convert BitArray to Array{Bool} if necessary
        if v isa BitArray
            to_save[k] = convert(Array{Bool}, v)
        else
            to_save[k] = v

            # if the value will result in a hard-to-read HDF5 file, warn
            t = if isa(v, Array)
                eltype(v)
            else
                typeof(v)
            end
            if !(t in [Bool, Int, Int64, Int32, Int16, Int8, UInt, UInt64, UInt32,
                UInt16, UInt8, Float64, Float32, String])
                #throw(ArgumentError("When saving to JLD, only types Strings and standard numerical types are supported. Type $t, which is being used for key $k, will result in a hard-to-read HDF5 file."))
                @warn "When saving to JLD, only types Strings and standard numerical types are supported. Type $t, which is being used for key $k, will result in a hard-to-read HDF5 file."
            end
        end
    end

    # record git branch and commit
    to_save[:git_branch] = git_branch
    to_save[:git_commit] = git_commit

    JLD2.jldsave(filename; to_save...)
end
