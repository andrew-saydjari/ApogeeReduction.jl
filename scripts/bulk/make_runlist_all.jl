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
expid_list = Int[]
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
        df = read_almanac_exp_df(f, tele, tstmjd)
        good_exp = (df.nreadInt .> 3) .|
                   ((df.imagetyp .== "DomeFlat") .& (df.observatory .== "apo")) .|
                   ((df.imagetyp .== "QuartzFlat") .& (df.nreadInt .== 3))
        dfindx_list_loc = findall(good_exp)
        for dfindx in dfindx_list_loc
            push!(mjdexp_list, parse(Int, tstmjd))
            push!(expid_list, long_expid_to_short(df.exposure[dfindx]))
            push!(dfindx_list, dfindx)
            push!(tele_list, tele)
        end
    end
end

safe_jldsave(parg["output"]; tele = tele_list, mjd = mjdexp_list, expid = expid_list, dfindx = dfindx_list)
