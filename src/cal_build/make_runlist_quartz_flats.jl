using Pkg;
Pkg.instantiate();
using JLD2, ArgParse, DataFrames
include("../utils.jl") # for safe_jldsave

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

darks_mjd = Int[]
darks_expid = Int[]
f = h5open(parg["almanac_file"])
mjd_list = keys(f[parg["tele"]])
for tstmjd in mjd_list
    df = DataFrame(read(f[parg["tele"] * "/$(tstmjd)/exposures"]))
    df.nreadInt = parse.(Int, df.nread)
    expindx_list = findall((df.imagetyp .== "QuartzFlat"))
    #    if f[parg["tele"]] == 'lco'
    #        expindx_list = findall((df.imagetyp .== "QuartzFlat") .& (df.nreadInt .== 3))
    #    elseif f[parg["tele"]] == 'apo'
    #        expindx_list = findall((df.imagetyp .== "QuartzFlat") .& (df.nreadInt .== 10))
    #    end
    for expindx in expindx_list
        push!(darks_mjd, parse(Int, tstmjd))
        push!(darks_expid, expindx)
    end
end

safe_jldsave(parg["output"]; mjd = darks_mjd, expid = darks_expid)
