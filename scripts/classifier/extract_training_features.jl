# Extract exposure-type classifier features for every manifest exposure from
# reduced ar2D images, via ApogeeReduction.exposure_class_features itself (the
# deployed inference path), so training and pipeline features are bit-identical
# by construction. NaN features + ok=0 are recorded for exposures whose chip
# files are missing (e.g. never reduced).
#
# Crash-resilient: work is checkpointed in chunks of --chunk exposures to
# <outfile>.parts/part_<i>.h5; a restart skips finished parts. Worker deaths
# (a libhdf5 segfault on a bad/interrupted read takes the whole process down
# and try/catch cannot save it — observed once per ~5k exposures on ceph) are
# retried once on another worker, then recorded as ok=0.
#
# Run: julia --project=<repo> -p <N> scripts/classifier/extract_training_features.jl \
#        --manifest <manifest.csv> --apred_dir <.../apred> --outfile <features.h5>
using Distributed
@everywhere using ApogeeReduction
@everywhere using ApogeeReduction: exposure_class_features
@everywhere using HDF5, Printf

using ArgParse, DataFrames, CSV, Printf

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--manifest"
        required = true
        help = "training manifest CSV (from build_training_manifest.jl)"
        arg_type = String
        "--apred_dir"
        required = true
        help = "apred directory holding <mjd>/ar2D_<obs>_<mjd>_<expnum>_<chip>_<image_type>.h5"
        arg_type = String
        "--outfile"
        required = true
        help = "output features HDF5"
        arg_type = String
        "--chunk"
        default = 500
        help = "exposures per checkpointed chunk"
        arg_type = Int
    end
    return parse_args(s)
end
parg = parse_commandline()

@everywhere begin
    const CHIPS = ("R", "G", "B")
    const NCHIPFEAT = ApogeeReduction.CLASSIFIER_NCHIPFEAT

    function exposure_features(apred, obs, mjd, expnum, imtype)
        feats = fill(NaN, NCHIPFEAT, 3)
        ok = true
        for (ci, chip) in enumerate(CHIPS)
            path = joinpath(apred, string(mjd),
                @sprintf("ar2D_%s_%d_%04d_%s_%s.h5", obs, mjd, expnum, chip, imtype))
            try
                feats[:, ci] = exposure_class_features(h5read(path, "dimage"))
            catch
                ok = false
            end
        end
        (feats = feats, ok = ok)
    end
end

df = CSV.read(parg["manifest"], DataFrame)
nexp = nrow(df)
apred = parg["apred_dir"]
chunk = parg["chunk"]
nparts = cld(nexp, chunk)
partd = parg["outfile"] * ".parts"
mkpath(partd)
println("extracting features for $nexp exposures in $nparts parts on ",
    nworkers(), " workers")

nf = ApogeeReduction.CLASSIFIER_NCHIPFEAT
sentinel(_) = (feats = fill(NaN, nf, 3), ok = false)

t0 = time()
for p in 1:nparts
    partfile = joinpath(partd, @sprintf("part_%04d.h5", p))
    isfile(partfile) && continue
    lo, hi = (p - 1) * chunk + 1, min(p * chunk, nexp)
    results = pmap(eachrow(df[lo:hi, :]); batch_size = 4, retry_delays = zeros(1),
        on_error = sentinel) do r
        exposure_features(apred, r.obs, r.mjd, r.expnum, r.image_type)
    end
    F = fill(NaN, nf, 3, length(results))
    for i in eachindex(results)
        F[:, :, i] = results[i].feats
    end
    tmp = partfile * ".tmp"
    h5open(tmp, "w") do f
        f["lo"] = lo
        f["hi"] = hi
        f["features"] = F
        f["ok"] = UInt8.(getindex.(results, :ok))
    end
    mv(tmp, partfile)
    @printf("part %d/%d done (rows %d:%d, fail %d) %.1f min elapsed, workers %d\n",
        p, nparts, lo, hi, sum(.!getindex.(results, :ok)), (time() - t0) / 60,
        nworkers())
    flush(stdout)
end

# ---- assemble ----
F = fill(NaN, nf, 3, nexp)
okv = zeros(UInt8, nexp)
for p in 1:nparts
    h5open(joinpath(partd, @sprintf("part_%04d.h5", p)), "r") do f
        F[:, :, read(f["lo"]):read(f["hi"])] = read(f["features"])
        okv[read(f["lo"]):read(f["hi"])] = read(f["ok"])
    end
end
println("total failures: ", sum(okv .== 0), " / ", nexp)

featnames = vcat(["q" * string(q) for q in ApogeeReduction.CLASSIFIER_QUANTS],
    ["mad", "nanfrac", "spec_med", "spec_maxrat", "spec_npk5", "spec_npk20",
        "spec_conc20", "spec_acpk", "fib_med", "fib_contrast", "fib_npk2", "fib_acpk"])

h5open(parg["outfile"], "w") do f
    f["features"] = F
    f["feature_names"] = featnames
    f["chips"] = collect(String.(CHIPS))
    f["ok"] = okv
    f["observatory"] = String.(df.obs)
    f["mjd"] = collect(df.mjd)
    f["exposure"] = collect(df.expnum)
    f["image_type"] = String.(df.image_type)
    f["labeled"] = String.(df.labeled)
    f["role"] = String.(df.role)
    f["truth_source"] = String.(df.truth_source)
end
println("wrote ", parg["outfile"])
