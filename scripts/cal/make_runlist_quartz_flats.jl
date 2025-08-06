using Pkg;
Pkg.instantiate();
using JLD2, ArgParse, DataFrames, HDF5
using ApogeeReduction: safe_jldsave, read_almanac_exp_df

## Parse command line arguments
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--tele"
        required = true
        help = "telescope name (apo or lco)"
        arg_type = String
        default = ""
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
expid_list = Int[]
f = h5open(parg["almanac_file"])
mjd_list = keys(f[parg["tele"]])
for tstmjd in mjd_list
    df = read_almanac_exp_df(f, parg["tele"], tstmjd)
    expindx_list = findall((df.imagetyp .== "QuartzFlat"))
    for expindx in expindx_list
        push!(mjdexp_list, parse(Int, tstmjd))
        push!(expid_list, expindx)
    end
end

safe_jldsave(parg["output"]; mjd = mjdexp_list, expid = expid_list, no_metadata = true)
