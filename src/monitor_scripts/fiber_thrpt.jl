using Pkg;
Pkg.instantiate();
using JLD2, ProgressMeter, ArgParse, SlackThreads, Glob, StatsBase

src_dir = "../"
include(src_dir * "/utils.jl")
include(src_dir * "/makie_plotutils.jl")

## Parse command line arguments
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--outdir"
        required = true
        help = "outdir of trace extraction runs"
        arg_type = String
        default = ""
    end
    return parse_args(s)
end

parg = parse_commandline()

## This should be deprecated in favor of the new relFluxing
## however, one might want to use it as a scaffold for comparisons. Revisit TODO

#this just dumps everything we have in outdir into a throughput file/plot
#runs on both quartz and dome flats and both telescopes
function summarize_fiber_thrpt(flat_type, tele)
    fname_list = sort(glob("$(flat_type)_flats/$(flat_type)Trace_$(tele)_*_a*", parg["outdir"]))
    sjd5 = map(x -> parse(Int, split(x, "_")[end - 2]), fname_list)
    expid = map(x -> parse(Int, last(split(x, "_")[end - 1], 4)), fname_list)
    traceid = 1:300
    adjfiberindx = traceid .+ 300 * (tele == "lco")
    fiberid = 301 .- traceid

    thrpt_mat = zeros(300, 3, length(fname_list))
    @showprogress for (findx, fname) in enumerate(fname_list),
        (cindx, chip) in enumerate(["a", "b", "c"])

        fname_loc = replace(fname, "_a_" => "_$(chip)_")
        trace_params = load(fname_loc, "trace_params")
        thrpt_mat[:, cindx, findx] .= dropdims(nanzeromedian(trace_params[:, :, 1], 1), dims = 1)
    end
    fname_out = joinpath(parg["outdir"], "monitor", "$(flat_type)_thrpt_summary_$(tele).jld2")
    if !ispath(dirname(fname_out))
        mkpath(dirname(fname_out))
    end

    safe_jldsave(fname_out;
        adjfiberindx = adjfiberindx,
        traceid = traceid,
        fiberid = fiberid,
        sjd5 = sjd5,
        expid = expid,
        thrpt_mat
    )

    dat = thrpt_mat[:, 1, :]
    dat ./= nanzeromedian(dat, 1)

    msklow = (dat .< (1 .- 5 * nanzeroiqr(dat, 1)))
    dat[msklow] .= NaN

    fig = Figure(size = (600, 600), fontsize = 22)
    ax = Axis(fig[1, 1], title = "$tele Chip a", xlabel = "Time", ylabel = "Fiber TRACEID")
    ce = heatmap!(ax, dat',
        colormap = :linear_bgy_10_95_c74_n256,
        nan_color = :red        # colorrange=(vmin,vmax)
    )

    Colorbar(fig[1, 2], ce, width = 20, height = Relative(1.0))
    colgap!(fig.layout, 1, 20)  # Space between image 1 and colorbar 1
    colsize!(fig.layout, 1, Aspect(1, size(dat, 2) / size(dat, 1)))

    resize_to_layout!(fig)
    framePath = joinpath(parg["outdir"], "monitor", "$(flat_type)_thrpt_summary_$(tele).png")
    save(framePath, fig, px_per_unit = 3)
    return framePath
end

thread = SlackThread();
thread("Fiber throughput summary for $(parg["outdir"])")

# Run for both flat types and telescopes
for tele in ["apo", "lco"]
    for flat_type in ["dome", "quartz"]
        framePath = summarize_fiber_thrpt(flat_type, tele)
        thread("Throughput to date using $(flat_type) flats on $(tele)", framePath)
    end
end
