# this has the plotting setup used for plotting with Makie

using CairoMakie, ColorSchemes
sas_prefix = "https://data.sdss5.org/sas/sdsswork/users/"

set_theme!(theme_black());

function nanify(x, msk)
    out = zeros(eltype(x), length(msk))
    out[msk] .= x
    out[.!msk] .= NaN
    return out
end
