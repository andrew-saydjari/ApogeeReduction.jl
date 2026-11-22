# Rebuild the 2D content forest's LIGHT REFERENCE from scratch
# (data/exposure_classifier_2d/light_reference.csv).
#
# The reference is, per (telescope, era), the median over a set of normal
# science exposures of three chip-G metrics of the ar2D image: the OH line flux
# `oh`, the continuum `cont` and the trace-height q90 `tr_q90`
# (ApogeeReduction.classifier2d_profile_features, the EXPOSURE_LIGHT 2D metric).
# The forest features oh_rel / cont_rel / trq90_rel and the object-frame rule's
# no-sky-light test are log10(x / reference).
#
# Era = plate (1) / FPS (2) with the FPS start of CLASSIFIER2D_ERA2_START
# (APO 59423, LCO 59809).
#
# Reference set: --reference_list, a CSV of tele, mjd, expnum. The shipped file was
# built from metadata/exposure_classifier_2d_light_reference_list.csv: the 2,995
# random release-eclass NORMAL commanded-object exposures of the 2026-11-18
# EXPOSURE_LIGHT study (eclass NORMAL = measured post-extraction throughput and
# RV scatter consistent with the catalog; see that study), of which those with a
# readable chip-G ar2D are used. Any other set of demonstrably normal science
# frames is a valid from-scratch choice; the medians then move by the sampling
# noise of the median.
#
# NOTE on precision: the reference uses the Float64 profiles straight from the
# image (as EXPOSURE_LIGHT did), NOT the Float32-rounded profiles the forest
# features use. NaN and Inf pixels are set to 0.
#
# Run: julia --project=<repo> [-p N] scripts/classifier/build_light_reference_2d.jl \
#        --reference_list <list.csv> --apred_dir <.../apred> --outfile <light_reference.csv>
using Distributed
@everywhere using ApogeeReduction
@everywhere using ApogeeReduction: classifier2d_profile_features, classifier2d_era
@everywhere using HDF5, Printf, Statistics
using ArgParse, CSV, DataFrames

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--reference_list"
        required = true
        help = "CSV with tele, mjd, expnum of the normal reference exposures (commanded objects)"
        arg_type = String
        "--apred_dir"
        required = true
        arg_type = String
        "--outfile"
        required = true
        arg_type = String
        "--per_exposure"
        default = ""
        help = "optional CSV of the per-exposure metrics"
        arg_type = String
    end
    return parse_args(s)
end
parg = parse_commandline()

@everywhere function light_metrics(apred, tele, mjd, expnum)
    p = joinpath(apred, string(mjd), @sprintf("ar2D_%s_%d_%04d_G_object.h5", tele, mjd, expnum))
    isfile(p) || return (ok = false, oh = NaN, cont = NaN, tr_q90 = NaN)
    try
        d = h5open(p, "r") do f
            read(f["dimage"])
        end
        s = d[1:N_XPIX, :]
        s[.!isfinite.(s)] .= 0
        v = classifier2d_profile_features(vec(mean(s, dims = 2)), vec(mean(s, dims = 1)))
        # PROFNAMES: 1 oh, 2 cont, 14 tr_q90
        (ok = isfinite(v[14]), oh = v[1], cont = v[2], tr_q90 = v[14])
    catch e
        @warn "light reference: could not read $p" exception = e
        (ok = false, oh = NaN, cont = NaN, tr_q90 = NaN)
    end
end

L = CSV.read(parg["reference_list"], DataFrame)
"obs" in names(L) && !("tele" in names(L)) && rename!(L, :obs => :tele)
res = pmap(r -> light_metrics(parg["apred_dir"], String(r.tele), r.mjd, r.expnum), eachrow(L);
    batch_size = 4)
L.ok = getindex.(res, :ok)
L.oh = getindex.(res, :oh)
L.cont = getindex.(res, :cont)
L.tr_q90 = getindex.(res, :tr_q90)
L.era = [classifier2d_era(r.tele, r.mjd) for r in eachrow(L)]
isempty(parg["per_exposure"]) || CSV.write(parg["per_exposure"], L)
G = L[L.ok, :]
R = combine(groupby(G, [:tele, :era]), nrow => :n_ref, :oh => median => :oh,
    :cont => median => :cont, :tr_q90 => median => :tr_q90)
sort!(R, [:tele, :era])
CSV.write(parg["outfile"], R)
println("reference exposures: ", nrow(L), " listed, ", nrow(G), " used")
show(R, allrows = true)
println("\nwrote ", parg["outfile"])
