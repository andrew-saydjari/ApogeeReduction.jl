# Stage 1 of the 2D content-forest rebuild: extract the REFERENCE-FREE features
# and the per-chip spectral profiles for every manifest exposure from reduced
# ar2D images, through the pipeline's own code
# (ApogeeReduction.read_classifier2d_chip_products / exposure_class_features_2d),
# so training and inference features are the same function by construction.
#
# Output HDF5 (row-aligned with the manifest):
#   tele, mjd, expnum, image_type, ok (1 = all three chips read, 0 = missing/failed),
#   nread_total (metadata/nread_total of chip G; provenance only, not a feature),
#   features (nfeat x n, the
#   reference-free features of classifier2d_feature_names(); t_* / *_rel are NaN
#   here and are filled by assemble_features_2d.jl), feature_names,
#   PS (2048 x 3 x n Float32 spectral profiles, chips R,G,B: the input of
#   build_line_templates_2d.jl and of the template correlations).
#
# Crash-resilient like PR #408's extract_training_features.jl: chunks of --chunk
# exposures are checkpointed to <outfile>.parts/part_<i>.h5 and a restart skips
# finished parts. ar2D files are opened read-only.
#
# Run: julia --project=<repo> -p <N> scripts/classifier/extract_features_2d.jl \
#        --manifest <manifest.csv> --apred_dir <.../apred> --outfile <features_2d.h5>
# The full archive (~315k exposures) is a Slurm-sized job (~1.3 s per exposure).
using Distributed
@everywhere using ApogeeReduction
@everywhere using ApogeeReduction: read_classifier2d_chip_products, exposure_class_features_2d,
                                   classifier2d_feature_names
@everywhere using HDF5, Printf
using ArgParse, DataFrames, CSV

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--manifest"
        required = true
        help = "CSV with tele (or obs), mjd, expnum, image_type"
        arg_type = String
        "--apred_dir"
        required = true
        help = "apred directory holding <mjd>/ar2D_<tele>_<mjd>_<expnum>_<chip>_<image_type>.h5"
        arg_type = String
        "--outfile"
        required = true
        arg_type = String
        "--chunk"
        default = 200
        arg_type = Int
    end
    return parse_args(s)
end
parg = parse_commandline()

@everywhere function extract_one(apred, tele, mjd, expnum, imtype)
    names_ = classifier2d_feature_names()
    fail = (ok = false, nread = -1, feats = fill(NaN, length(names_)),
        PS = zeros(Float32, N_XPIX, 3))
    paths = Dict(chip => joinpath(apred, string(mjd),
                     @sprintf("ar2D_%s_%d_%04d_%s_%s.h5", tele, mjd, expnum, chip, imtype))
    for chip in CHIP_LIST)
    all(isfile, values(paths)) || return fail   # missing ar2D: recorded as ok = 0
    try
        prods = Dict(chip => read_classifier2d_chip_products(paths[chip], chip)
        for chip in CHIP_LIST)
        nr = prods["G"].nread_total   # provenance only; not a v2 feature
        F = exposure_class_features_2d(prods, tele, mjd)
        PS = hcat([prods[c].ps for c in CHIP_LIST]...)
        (ok = true, nread = nr, feats = Float64[F[k] for k in names_], PS = PS)
    catch e
        @warn "feature extraction failed" tele mjd expnum exception = e
        fail
    end
end

man = CSV.read(parg["manifest"], DataFrame)
"obs" in names(man) && !("tele" in names(man)) && rename!(man, :obs => :tele)
n = nrow(man)
names_ = classifier2d_feature_names()
nf = length(names_)
partd = parg["outfile"] * ".parts"
mkpath(partd)
nparts = cld(n, parg["chunk"])
println("extracting 2D features for $n exposures in $nparts parts on $(nworkers()) workers")
t0 = time()
for p in 1:nparts
    pf = joinpath(partd, @sprintf("part_%05d.h5", p))
    isfile(pf) && continue
    lo, hi = (p - 1) * parg["chunk"] + 1, min(p * parg["chunk"], n)
    res = pmap(eachrow(man[lo:hi, :]); batch_size = 2) do r
        extract_one(parg["apred_dir"], String(r.tele), r.mjd, r.expnum, String(r.image_type))
    end
    tmp = pf * ".tmp"
    h5open(tmp, "w") do f
        f["lo"] = lo
        f["hi"] = hi
        f["ok"] = UInt8.(getindex.(res, :ok))
        f["nread_total"] = Int.(getindex.(res, :nread))
        f["features"] = hcat(getindex.(res, :feats)...)
        f["PS"] = cat(getindex.(res, :PS)...; dims = 3)
    end
    mv(tmp, pf; force = true)
    @printf("part %d/%d rows %d:%d ok %d  %.1f min\n", p, nparts, lo, hi,
        sum(getindex.(res, :ok)), (time() - t0) / 60)
    flush(stdout)
end

F = fill(NaN, nf, n)
PS = zeros(Float32, N_XPIX, 3, n)
ok = zeros(UInt8, n)
nread = fill(-1, n)
for p in 1:nparts
    h5open(joinpath(partd, @sprintf("part_%05d.h5", p)), "r") do f
        lo, hi = read(f["lo"]), read(f["hi"])
        F[:, lo:hi] = read(f["features"])
        PS[:, :, lo:hi] = read(f["PS"])
        ok[lo:hi] = read(f["ok"])
        nread[lo:hi] = read(f["nread_total"])
    end
end
h5open(parg["outfile"], "w") do f
    f["tele"] = String.(man.tele)
    f["mjd"] = collect(man.mjd)
    f["expnum"] = collect(man.expnum)
    f["image_type"] = String.(man.image_type)
    f["ok"] = ok
    f["nread_total"] = nread
    f["features"] = F
    f["feature_names"] = names_
    f["PS"] = PS
end
println("wrote ", parg["outfile"], ": ", sum(ok), " / ", n, " exposures with all three chips")
