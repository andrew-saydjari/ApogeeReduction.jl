# this has the plotting setup used for plotting with PythonPlot
# it should be eliminated once we switch to Makie

ENV["JULIA_CONDAPKG_BACKEND"] = "MicroMamba";
using CondaPkg;
using PythonCall;
pyimport("sys").stdout = pytextio(stdout);
pyimport("sys").stderr = pytextio(stderr);
import PythonPlot;
const plt = PythonPlot.pyplot;
plt.matplotlib.style.use("dark_background");
cc = pyimport("colorcet");
mplcolors = pyimport("matplotlib.colors");
mpltk = pyimport("mpl_toolkits.axes_grid1");
mplcm = pyimport("matplotlib.cm");
mplani = pyimport("matplotlib.animation");

PythonPlot.matplotlib.rcParams["figure.max_open_warning"] = 40

sas_prefix = "https://data.sdss5.org/sas/sdsswork/users/"

function nanify(x, msk)
    out = zeros(eltype(x), length(msk))
    out[msk] .= x
    out[.!msk] .= NaN
    return out
end

function nanify(x)
    msk = isnanorzero.(x)
    return nanify(x, msk)
end
