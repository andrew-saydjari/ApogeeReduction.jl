ENV["JULIA_CONDAPKG_BACKEND"] = "MicroMamba";
using CondaPkg; using PythonCall
pyimport("sys").stdout = pytextio(stdout); pyimport("sys").stderr = pytextio(stderr);
import PythonPlot; const plt = PythonPlot.pyplot; 
plt.matplotlib.style.use("dark_background"); cc=pyimport("colorcet");
mplcolors=pyimport("matplotlib.colors"); mpltk=pyimport("mpl_toolkits.axes_grid1");
mplcm = pyimport("matplotlib.cm"); mplani=pyimport("matplotlib.animation");

sas_prefix = "https://data.sdss5.org/sas/sdsswork/"