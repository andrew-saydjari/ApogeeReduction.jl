# Extract exposure-type classifier features for every manifest exposure from
# reduced ar2D images, via ApogeeReduction.exposure_class_features itself (the
# deployed inference path), so training and pipeline features are bit-identical
# by construction. NaN features + ok=0 are recorded for exposures whose chip
# files are missing (e.g. never reduced).
#
# Run: julia --project=<repo> -p <N> scripts/classifier/extract_training_features.jl \
#        --manifest <manifest.csv> --apred_dir <.../apred> --outfile <features.h5>
using Distributed
@everywhere using Pkg
@everywhere using ApogeeReduction
@everywhere using ApogeeReduction: exposure_class_features
@everywhere using HDF5, Printf

using ArgParse, DataFrames, CSV, ProgressMeter

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
println("extracting features for ", nrow(df), " exposures on ", nworkers(), " workers")
apred = parg["apred_dir"]

t0 = time()
results = @showprogress pmap(eachrow(df); batch_size = 8) do r
    exposure_features(apred, r.obs, r.mjd, r.expnum, r.image_type)
end
println("done in ", round((time() - t0) / 60, digits = 1), " min; failures: ",
    sum(.!getindex.(results, :ok)))

nexp = nrow(df)
F = zeros(Float64, ApogeeReduction.CLASSIFIER_NCHIPFEAT, 3, nexp)
for i in 1:nexp
    F[:, :, i] = results[i].feats
end
featnames = vcat(["q" * string(q) for q in ApogeeReduction.CLASSIFIER_QUANTS],
    ["mad", "nanfrac", "spec_med", "spec_maxrat", "spec_npk5", "spec_npk20",
        "spec_conc20", "spec_acpk", "fib_med", "fib_contrast", "fib_npk2", "fib_acpk"])

h5open(parg["outfile"], "w") do f
    f["features"] = F
    f["feature_names"] = featnames
    f["chips"] = collect(String.(CHIPS))
    f["ok"] = UInt8.(getindex.(results, :ok))
    f["observatory"] = String.(df.obs)
    f["mjd"] = collect(df.mjd)
    f["exposure"] = collect(df.expnum)
    f["image_type"] = String.(df.image_type)
    f["labeled"] = String.(df.labeled)
    f["role"] = String.(df.role)
    f["truth_source"] = String.(df.truth_source)
end
println("wrote ", parg["outfile"])
