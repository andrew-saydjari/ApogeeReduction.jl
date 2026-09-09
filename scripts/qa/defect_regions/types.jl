# Result types shared by run_sweep.jl / run_trace.jl (which write them) and
# analyze.jl / make_plots.jl (which read them back).  They live at Main scope so that
# `Serialization` can round-trip them between scripts.

struct ProductEvidence
    name::String
    medz::Float64      # median detrended residual z inside the region
    maxS::Float64      # signed most-extreme matched-filter significance
    detfrac::Float64   # fraction of region pixels this product itself detected
end

struct Region
    tele::String
    chip::String
    id::Int
    x0::Int; x1::Int; y0::Int; y1::Int
    area::Int
    evidence::Vector{ProductEvidence}
    nepoch_det::Int
    nepoch::Int
    mjd_lo::Int
    mjd_hi::Int
    epoch_flags::Vector{Bool}
    arflag_frac::Float64     # fraction already carrying AR bad_pix_bits
    fibers::Vector{Int}      # fibre indices whose trace crosses the region
    pixels::Vector{Tuple{Int, Int}}   # science coords
end

struct TraceRegion
    tele::String; chip::String; id::Int
    x0::Int; x1::Int; f0::Int; f1::Int
    area::Int
    medz::Float64
    maxS::Float64
    nepoch_det::Int; nepoch::Int
    mjd_lo::Int; mjd_hi::Int
    class::String
    fibers::Vector{Int}
end
