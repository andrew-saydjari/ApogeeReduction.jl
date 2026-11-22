# Stage 3 of the 2D content-forest rebuild: add the reference-dependent features
# (12 template correlations t_*, and oh_rel / cont_rel / trq90_rel) to the
# reference-free features of extract_features_2d.jl, and write the training
# feature table: one CSV row per exposure with all three chips (ok = 1), columns
# tele, mjd, expnum, image_type + classifier2d_feature_names() (n_read from the
# ar2D metadata/nread_total).
#
# Uses ApogeeReduction.classifier2d_add_reference_features!, the same function
# the pipeline uses, on the same Float32 spectral profiles.
#
# Run: JULIA_NUM_THREADS=n julia --project=<repo> scripts/classifier/assemble_features_2d.jl \
#        --features <features_2d.h5> --templates <line_templates.h5> \
#        --light_reference <light_reference.csv> --outfile <features_2d.csv>
using ArgParse, CSV, DataFrames, HDF5
using ApogeeReduction
using ApogeeReduction: classifier2d_add_reference_features!, classifier2d_feature_names,
                       load_classifier2d_references

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--features"
        required = true
        arg_type = String
        "--templates"
        required = true
        arg_type = String
        "--light_reference"
        required = true
        arg_type = String
        "--outfile"
        required = true
        arg_type = String
    end
    return parse_args(s)
end
parg = parse_commandline()
refs = load_classifier2d_references(parg["templates"], parg["light_reference"])
refs === nothing && error("references missing")
names_ = classifier2d_feature_names()

tele, mjd, expnum, imtype, ok, nread, X, fn, PS = h5open(parg["features"], "r") do f
    read(f["tele"]), read(f["mjd"]), read(f["expnum"]), read(f["image_type"]),
    read(f["ok"]) .== 1, read(f["nread_total"]), read(f["features"]),
    read(f["feature_names"]), read(f["PS"])
end
fn == names_ || error("feature_names in $(parg["features"]) differ from classifier2d_feature_names()")
rows = findall(ok)
out = fill(NaN, length(names_), length(rows))
Threads.@threads for j in eachindex(rows)
    i = rows[j]
    F = Dict{String, Float64}(zip(names_, X[:, i]))
    F["n_read"] = Float64(nread[i])
    classifier2d_add_reference_features!(F,
        Dict(c => PS[:, ci, i] for (ci, c) in enumerate(CHIP_LIST)), tele[i], mjd[i], refs)
    out[:, j] = [F[k] for k in names_]
end
D = DataFrame(tele = tele[rows], mjd = mjd[rows], expnum = expnum[rows], image_type = imtype[rows])
for (k, nm) in enumerate(names_)
    D[!, nm] = out[k, :]
end
CSV.write(parg["outfile"], D)
println("wrote ", parg["outfile"], ": ", nrow(D), " exposures x ", length(names_), " features")
