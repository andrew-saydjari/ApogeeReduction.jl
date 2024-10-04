using Pkg;
Pkg.instantiate();
using HDF5, JLD2, ArgParse, DataFrames

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

mjd = []
expid = []
f = h5open(parg["almanac_file"])
mjd_list = keys(f[parg["tele"]])
for tstmjd in mjd_list
    df = DataFrame(read(f[parg["tele"] * "/$(tstmjd)/exposures"]))
    df.nreadInt = parse.(Int, df.nread)
    expindx_list = findall((df.nreadInt .> 3))
    for expindx in expindx_list
        if expindx > 1
            push!(mjd, parse(Int, tstmjd))
            push!(expid, expindx)
        end
    end
end

jldsave(parg["output"]; mjd = mjd, expid = expid)
