#=
Automated FPI fiber-index homogeneity check (task #54).

Surveys the corpus of FPI guide fiber indices per telescope and tests whether
it is homogeneous. Design (AKS, 2026-09-08):

  "survey the corpus of fpiguide fiber indices per telescope and test whether
   it is homogeneous. Any night whose FPI fibers differ from that telescope's
   dominant pair is flagged for review/exclusion. This generalises past the
   three known nights and needs no constant."

The FPI pair per configuration comes from the almanac `bonus` rows, the same
derivation `get_fpi_fiberIDs_from_almanac` uses in production. The check
FLAGS inhomogeneity for review; it never hard-rejects (exit status is 0 even
when nights are flagged). The known example is the LCO re-fibering epoch
MJD 59810-59850 (pair 142/153 instead of the dominant 82/213), ruled
trustworthy and kept in DR21 -- listed in metadata/fpi_known_epochs.txt so it
shows up [KNOWN] rather than [NEW]. Only [NEW] deviations need human review.

Scope: FPS era only (apo mjd > 59423, lco mjd > 59808). The plate era has no
confSummary and no `bonus` stubs, so the FPI concept does not apply the same
way there and those nights are skipped by construction.

Usage:
  julia --project=. scripts/verify/check_fpi_homogeneity.jl \
      --almanac_file /path/to/allobs.h5 \
      --known_epochs metadata/fpi_known_epochs.txt

This belongs in the automated rejection checks for a delivered dataset, run
against the corpus almanac of the run.
=#

using ArgParse
using ApogeeReduction: survey_fpi_homogeneity, report_fpi_homogeneity, parse_known_epochs

## Parse command line arguments
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--almanac_file"
        required = true
        help = "path to the corpus almanac HDF5 file (raw/ layout)"
        arg_type = String
        "--known_epochs"
        required = false
        help = "file of known/accepted deviation epochs (`tele mjd_lo mjd_hi [label]` per line); pass \"\" for none"
        arg_type = String
        default = joinpath(dirname(@__DIR__), "..", "metadata", "fpi_known_epochs.txt")
        "--tele"
        required = false
        help = "comma-separated telescopes to survey"
        arg_type = String
        default = "apo,lco"
    end
    return parse_args(s)
end

parg = parse_commandline()

known_epochs = if isempty(parg["known_epochs"])
    parse_known_epochs("")
else
    parse_known_epochs(parg["known_epochs"])
end

teles = String.(split(parg["tele"], ","))
results = survey_fpi_homogeneity(parg["almanac_file"]; teles = teles)
summary = report_fpi_homogeneity(stdout, results; known_epochs = known_epochs)

# Inhomogeneity is a flag for review, never a hard error.
exit(0)
