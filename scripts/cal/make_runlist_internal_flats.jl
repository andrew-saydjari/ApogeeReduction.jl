using Pkg;
Pkg.instantiate();
using JLD2, ArgParse, DataFrames, HDF5
using ApogeeReduction: safe_jldsave, read_almanac_exp_df

## Parse command line arguments
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--tele"
        required = false
        help = "telescope name (apo or lco)"
        arg_type = String
        default = "both"
        "--almanac_file"
        required = true
        help = "path to the almanac file"
        arg_type = String
        default = ""
        "--output"
        required = true
        help = "path to output runlist file"
        arg_type = String
        default = ""
    end
    return parse_args(s)
end

parg = parse_commandline()

mjdexp_list = Int[]
expid_list = String[]
dfindx_list = Int[]
tele_list = String[]
f = h5open(parg["almanac_file"])
tele2do = if parg["tele"] == "both"
    keys(f)
else
    [parg["tele"]]
end
for tele in tele2do
    mjd_list = keys(f[tele])
    for tstmjd in mjd_list
        tstmjd_int = parse(Int, tstmjd)
        df = read_almanac_exp_df(f, tele, tstmjd)
        good_exp = (df.image_type .== "internalflat") .& (df.n_read .> 3) .& (df.chip_flags .== 7)
        dfindx_list_loc = findall(good_exp)
        for dfindx in dfindx_list_loc
            push!(mjdexp_list, tstmjd_int)
            push!(expid_list, df.exposure_string[dfindx])
            push!(dfindx_list, dfindx)
            push!(tele_list, tele)
        end
    end
end

safe_jldsave(parg["output"], Dict{String, Any}(); tele = tele_list, mjd = mjdexp_list, expid = expid_list, dfindx = dfindx_list)
