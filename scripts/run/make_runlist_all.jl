# makes a runlist of all the exposures that should be run for a given night.
# called by run_all.sh
using Pkg;
Pkg.instantiate();
using HDF5, JLD2, ArgParse, DataFrames
using ApogeeReduction: read_almanac_exp_df, safe_jldsave # for safe_jldsave

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
    good_exp = (df.nreadInt .> 3) .|
               ((df.imagetyp .== "DomeFlat") .& (df.observatory .== "apo")) .|
               ((df.imagetyp .== "QuartzFlat") .& (df.nreadInt .== 3))
    expindx_list = findall(good_exp)
    for expindx in expindx_list
        push!(mjdexp_list, parse(Int, tstmjd))
        push!(expid_list, expindx)
    end
end

safe_jldsave(parg["output"]; mjd = mjdexp_list, expid = expid_list)
