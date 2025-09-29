# makes a runlist of all the exposures that should be run for a given night.
# called by run_all.sh
using Pkg;
Pkg.instantiate();
using HDF5, JLD2, ArgParse, DataFrames
using ApogeeReduction: read_almanac_exp_df, safe_jldsave, long_expid_to_short

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

# Check if almanac file is empty and exit with specific code
if length(keys(f)) == 0
    println("WARNING: No exposures found for the specified parameters")
    println("Exiting with code 16 to indicate empty runlist")
    exit(16)
end

for tele in tele2do
    mjd_list = keys(f[tele])
    for tstmjd in mjd_list
        tstmjd_int = parse(Int, tstmjd)
        df = read_almanac_exp_df(f, tele, tstmjd)
        good_exp = (df.n_read .> 3) .|
                   ((df.imagetyp .== "QuartzFlat") .& (df.n_read .== 3))
        dfindx_list_loc = findall(good_exp)
        for dfindx in dfindx_list_loc
            push!(mjdexp_list, tstmjd_int)
            push!(expid_list, long_expid_to_short(tstmjd_int,df.exposure_int[dfindx]))
            push!(dfindx_list, dfindx)
            push!(tele_list, tele)
        end
    end
end

# Save the data even if empty
safe_jldsave(parg["output"], Dict{String, Any}(); tele = tele_list, mjd = mjdexp_list, expid = expid_list, dfindx = dfindx_list)

# Check if list is empty and exit with specific code
if length(mjdexp_list) == 0
    println("WARNING: No exposures found for the specified parameters")
    println("Exiting with code 16 to indicate empty runlist")
    exit(16)
end
