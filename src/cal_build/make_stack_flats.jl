using JLD2, ProgressMeter, ArgParse, SlackThreads, Glob, StatsBase

src_dir = "../"
include(src_dir * "/fileNameHandling.jl")
include(src_dir * "/utils.jl")
include(src_dir * "/plotutils.jl")

## Parse command line arguments
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--chip"
        required = false
        help = "chip name (a, b, c)"
        arg_type = String
        "--tele"
        required = true
        help = "telescope name (apo or lco)"
        arg_type = String
        default = ""
        "--mjd-start"
        required = true
        help = "start mjd"
        arg_type = Int
        default = 0
        "--mjd-end"
        required = true
        help = "end mjd"
        arg_type = Int
        default = 0
        "--flat_dir"
        required = true
        help = "directory where 2D extractions of flats are stored"
        arg_type = String
        default = ""
        "--runlist"
        required = true
        help = "path name to hdf5 file with keys specifying list of exposures to run"
        arg_type = String
        default = ""
        "--dark_path"
        required = true
        help = "directory where 2D extractions of flats are stored"
        arg_type = String
        default = ""
        "--dark-mjd-start"
        required = true
        help = "start mjd"
        arg_type = Int
        default = 0
        "--dark-mjd-end"
        required = true
        help = "end mjd"
        arg_type = Int
        default = 0
    end
    return parse_args(s)
end

parg = parse_commandline()

chip = parg["chip"] # Python plotting issues prevented this from looping, so just do with three calls

# make summary plots and gifs and send to slack in outer loop
thread = SlackThread();
thread("Thread for flat stack for $(parg["tele"]) $(chip) from $(parg["mjd-start"]) to $(parg["mjd-end"])")

bad_dark_pix_bits = 2^2 + 2^4 + 2^5;
flat_frac_cut = 0.2
pcut_flat = 0.2
fx, fy = 10, 10    # Number of frequencies in x and y

dirNamePlots = parg["flat_dir"] * "plots/"
if !ispath(dirNamePlots)
    mkpath(dirNamePlots)
end

cmap_g = PythonPlot.get_cmap("cet_gouldian")
cmap_g.set_bad("r")

# Load in the exact set of exposures
mjd = load(parg["runlist"], "mjd")
expid = load(parg["runlist"], "expid");
flist = get_cal_file.(
    Ref(parg["flat_dir"]), Ref(parg["tele"]), mjd, expid, Ref(chip), Ref("INTERNALFLAT"))

# this is dumb and should probably be replaced by ivar weighting
nx, ny = 2040, 2040
flat_im = zeros(nx, ny)
model_im = zeros(nx, ny)

## this is dumb, Andy Casey will replace with a working 2d FINUFFT
x = range(0, stop = π, length = nx)
y = range(0, stop = π, length = ny)
X, Y = [x[i] for i in 1:nx, j in 1:ny], [y[j] for i in 1:nx, j in 1:ny];

design_matrix = gen_design_mat(nx,ny,fx,fy,X,Y)

# fig = PythonPlot.figure(figsize = (8, 8), dpi = 300)
# ax = fig.add_subplot(1, 1, 1)

# divider = mpltk.make_axes_locatable(ax)
# cax = divider.append_axes("right", size = "5%", pad = 0.05)
# cbar = plt.colorbar(
#     mplcm.ScalarMappable(
#         norm = mplcolors.Normalize(vmin = -0.2, vmax = 0.2), cmap = "cet_bkr"),
#     cax = cax,
#     orientation = "vertical")
# im_lst = []

#terrible hard code that needs to be replaced by cal look up system
f = jldopen("$(parg["dark_path"])"*"darkRate_$(parg["tele"])_$(chip)_$(parg["dark-mjd-start"])_$(parg["dark-mjd-end"]).jld2")
dark_pix_bitmask = f["dark_pix_bitmask"]
close(f)
bad_pix_dark = (dark_pix_bitmask[5:2044,5:2044] .& bad_dark_pix_bits .!= 0);

@showprogress for (indx, fname) in enumerate(flist)
    sname = split(fname, "_")
    tele, mjdloc, chiploc, expidloc = sname[(end - 4):(end - 1)]
    floc = jldopen(fname)
    temp_im = floc["dimage"]
    close(floc)

    bmat = temp_im[5:2044,5:2044]
    b = bmat[:];
    ref_med = nanzeromedian(b)

    bad_pix = copy(bad_pix_dark)
    bad_pix .|= (bmat./ref_med.<pcut_flat)
    good_pix = .!grow_msk2d(bad_pix; rad=3)[:]

    cvec = design_matrix[good_pix,:] \ b[good_pix]
    modImage = reshape(design_matrix*cvec,nx,ny);

    flat_im .+= bmat./modImage
    model_im .+= modImage

    # img = ax.imshow(temp_im',
    #     vmin = -0.2,
    #     vmax = 0.2,
    #     interpolation = "none",
    #     cmap = "cet_bkr",
    #     origin = "lower",
    #     aspect = "auto"
    # )

    # ttl = plt.text(
    #     0.5, 1.01, "Tele: $(tele), MJD: $(mjd), Chip: $(chiploc) Expid: $(expid)",
    #     ha = "center", va = "bottom", transform = ax.transAxes)

    # push!(im_lst, [img, ttl])
end
# PythonPlot.plotclose(fig)
# ani = mplani.ArtistAnimation(
#     fig, im_lst, interval = 300, repeat_delay = 300, blit = false)
# vidPath = dirNamePlots *
#           "darkStack_$(parg["tele"])_$(chip)_$(parg["mjd-start"])_$(parg["mjd-end"]).mp4"
# ani.save(vidPath)
# ani = nothing

flat_im ./= length(flist)
model_im ./= length(flist)

pix_bit_mask = zeros(Int, 2040, 2040)
pix_bit_mask .|= (flat_im .< flat_frac_cut) * 2^6 # too low response in flat

# save dark_pix_bitmask and dark_rate (electron per read)
jldsave(
    parg["flat_dir"] *
    "flats/flatFraction_$(parg["tele"])_$(chip)_$(parg["mjd-start"])_$(parg["mjd-end"]).jld2";
    flat_im = flat_im,
    flat_pix_bitmask = pix_bit_mask,
    model_im = model_im)

# Figures for QA
vmin, vmax = percentile(model_im[:],[2,98])

fig = PythonPlot.figure(figsize = (8, 8), dpi = 300)
ax = fig.add_subplot(1, 1, 1)
img = ax.imshow(flat_im',
    vmin = 0.95,
    vmax = 1.05,
    interpolation = "none",
    cmap = cmap_g,
    origin = "lower",
    aspect = "auto"
)

plt.text(0.5,
    1.01,
    "Tele: $(parg["tele"]), MJD Range: $(parg["mjd-start"])-$(parg["mjd-end"]), Chip: $(chip)",
    ha = "center",
    va = "bottom",
    transform = ax.transAxes)

divider = mpltk.make_axes_locatable(ax)
cax = divider.append_axes("right", size = "5%", pad = 0.05)
cbar = plt.colorbar(img, cax = cax, orientation = "vertical")
flatPath = dirNamePlots *
           "flatFraction_$(parg["tele"])_$(chip)_$(parg["mjd-start"])_$(parg["mjd-end"]).png"
fig.savefig(flatPath, bbox_inches = "tight", pad_inches = 0.1)
thread("Flat stack for $(parg["tele"]) $(chip) from $(parg["mjd-start"]) to $(parg["mjd-end"]) done.")
PythonPlot.plotclose(fig)
thread("Here is the final flat image", flatPath)

# dark_im[pix_bit_mask .& 2^3 .== 0] .= 0;
# dark_im[pix_bit_mask .& bad_pix_bits .!= 0] .= NaN;

# totNum = length(pix_bit_mask[1:2048, 1:2048])
# badVec = pix_bit_mask[1:2048, 1:2048] .& bad_pix_bits .!= 0
# corrVec = (pix_bit_mask[1:2048, 1:2048] .& 2^3 .!= 0) .& .!badVec
# notCorVec = (pix_bit_mask[1:2048, 1:2048] .& 2^3 .== 0)

# fracBad = count(badVec) / totNum
# fracCorr = count(corrVec) / totNum
# fracNotCorr = count(notCorVec) / totNum

# fig = PythonPlot.figure(figsize = (8, 8), dpi = 300)
# ax = fig.add_subplot(1, 1, 1)
# img = ax.imshow(dark_im',
#     vmin = -0.2,
#     vmax = 0.2,
#     interpolation = "none",
#     cmap = cmap_w,
#     origin = "lower",
#     aspect = "auto"
# )

# plt.text(0.5,
#     1.01,
#     "Tele: $(parg["tele"]), MJD Range: $(parg["mjd-start"])-$(parg["mjd-end"]), Chip: $(chip)\n Bad: $(round(100*fracBad,digits=2))% Corrected: $(round(100*fracCorr,digits=2))% NoCorrection: $(round(100*fracNotCorr,digits=2))%",
#     ha = "center",
#     va = "bottom",
#     transform = ax.transAxes
# )

# divider = mpltk.make_axes_locatable(ax)
# cax = divider.append_axes("right", size = "5%", pad = 0.05)
# cbar = plt.colorbar(img, cax = cax, orientation = "vertical")
# maskPath = dirNamePlots *
#            "darkRateMask_$(parg["tele"])_$(chip)_$(parg["mjd-start"])_$(parg["mjd-end"]).png"
# fig.savefig(maskPath, bbox_inches = "tight", pad_inches = 0.1)
# thread("Here is the final dark mask image", maskPath)
# PythonPlot.plotclose(fig)

# len_vid = stat(vidPath).size
# # if the length is bigger than 1 gigabyte, we need to upload the link to slack
# if len_vid > 1e9 # 1 GB
#     vidSasPath = replace(abspath(vidPath), r".*users" => sas_prefix)
#     thread("Here is the video of all of the frames included in the stack: $vidSasPath")
# else
#     thread("Here is the video of all of the frames included in the stack", vidPath)
# end

# fig = PythonPlot.figure(figsize = (8, 8), dpi = 300)
# ax = fig.add_subplot(1, 1, 1)

# divider = mpltk.make_axes_locatable(ax)
# cax = divider.append_axes("right", size = "5%", pad = 0.05)
# cbar = plt.colorbar(
#     mplcm.ScalarMappable(
#         norm = mplcolors.Normalize(vmin = -0.2, vmax = 0.2), cmap = "cet_bkr"),
#     cax = cax,
#     orientation = "vertical")
# im_lst = []
# @showprogress for (indx, fname) in enumerate(flist)
#     sname = split(fname, "_")
#     tele, mjd, chiploc, expid = sname[(end - 4):(end - 1)]
#     f = jldopen(fname)
#     temp_im = f["dimage"]
#     close(f)
#     ref_val_vec[indx] = nanzeromedian(temp_im[1:2048, 1:2048])
#     temp_im[1:2048, 1:2048] .-= ref_val_vec[indx] .-dark_im

#     img = ax.imshow(temp_im[1:2048, 1:2048]',
#         vmin = -0.2,
#         vmax = 0.2,
#         interpolation = "none",
#         cmap = "cet_bkr",
#         origin = "lower",
#         aspect = "auto"
#     )

#     ttl = plt.text(
#         0.5, 1.01, "Tele: $(tele), MJD: $(mjd), Chip: $(chiploc) Expid: $(expid)",
#         ha = "center", va = "bottom", transform = ax.transAxes)

#     push!(im_lst, [img, ttl])
# end
# PythonPlot.plotclose(fig)
# ani = mplani.ArtistAnimation(
#     fig, im_lst, interval = 300, repeat_delay = 300, blit = false)
# vidPath = dirNamePlots *
#           "darkStackRes_$(parg["tele"])_$(chip)_$(parg["mjd-start"])_$(parg["mjd-end"]).mp4"
# ani.save(vidPath)
# ani = nothing

# len_vid = stat(vidPath).size
# # if the length is bigger than 1 gigabyte, we need to upload the link to slack
# if len_vid > 1e9 # 1 GB
#     vidSasPath = replace(abspath(vidPath), r".*users" => sas_prefix)
#     thread("Here is the video of all of the residuals for frames included in the stack: $vidSasPath")
# else
#     thread("Here is the video of all of the residuals for frames included in the stack", vidPath)
# end