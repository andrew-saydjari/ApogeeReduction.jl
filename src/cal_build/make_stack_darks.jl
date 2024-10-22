using JLD2, ProgressMeter, ArgParse, SlackThreads, Glob, StatsBase

src_dir = "../"
include(src_dir * "/utils.jl")
include(src_dir * "/fileNameHandling.jl")
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
        "--dark_dir"
        required = true
        help = "directory where 2D extractions of darks are stored"
        arg_type = String
        default = ""
        "--runlist"
        required = true
        help = "path name to hdf5 file with keys specifying list of exposures to run"
        arg_type = String
        default = ""
    end
    return parse_args(s)
end

parg = parse_commandline()

chip = parg["chip"] # Python plotting issues prevented this from looping, so just do with three calls

# make summary plots and gifs and send to slack in outer loop
thread = SlackThread();
thread("DARK stack for $(parg["tele"]) $(chip) from $(parg["mjd-start"]) to $(parg["mjd-end"])")

bad_pix_bits = 2^2 + 2^4 + 2^5;
sig_measure = 0
sig_bad_lower = 5
sig_bad_upper = 7

dirNamePlots = parg["dark_dir"] * "plots/"
if !ispath(dirNamePlots)
    mkpath(dirNamePlots)
end

cmap_w = PythonPlot.get_cmap("cet_bkr")
cmap_w.set_bad("w")

# Load in the exact set of exposures
mjd = load(parg["runlist"], "mjd")
expid = load(parg["runlist"], "expid");
flist = get_cal_file.(
    Ref(parg["dark_dir"]), Ref(parg["tele"]), mjd, expid, Ref(chip), Ref("DARK"))

# this is dumb and should probably be replaced by ivar weighting
ref_val_vec = zeros(length(flist))
dark_im = zeros(2560, 2048)

@showprogress for (indx, fname) in enumerate(flist)
    sname = split(fname, "_")
    tele, mjd, chiploc, expid = sname[(end - 4):(end - 1)]
    f = jldopen(fname)
    temp_im = f["dimage"]
    close(f)
    ref_val_vec[indx] = nanzeromedian(temp_im[1:2048, 1:2048])
    temp_im[1:2048, 1:2048] .-= ref_val_vec[indx]
    dark_im .+= temp_im
end

dark_im ./= length(flist)
cen_dark = nanzeromedian(dark_im)

dat = dark_im[1:2048, 1:2048][:]
sig_est = nanzeroiqr(dat)

pix_bit_mask = zeros(Int, 2560, 2048)
pix_bit_mask[2049:end, :] .|= 2^0 # refArray
pix_bit_mask[1:4, :] .|= 2^1 # refPixels
pix_bit_mask[2045:2048, :] .|= 2^1 # refPixels
pix_bit_mask[:, 1:4] .|= 2^1 # refPixels
pix_bit_mask[:, (end - 3):end] .|= 2^1 # refPixels
pix_bit_mask[:, end] .|= 2^2 # bad refPixels (this is hard coded, only for chipA and APO for now)
pix_bit_mask .|= (abs.(dark_im) .> sig_measure .* sig_est) * 2^3
pix_bit_mask .|= (dark_im .< -sig_bad_lower .* sig_est) * 2^4
pix_bit_mask .|= (dark_im .> sig_bad_upper .* sig_est) * 2^5

dat = dark_im[1:2048, 1:2048][pix_bit_mask[1:2048, 1:2048] .& bad_pix_bits .== 0]
sig_after = nanzeroiqr(dat)

# save dark_pix_bitmask and dark_rate (electron per read)
jldsave(
    parg["dark_dir"] *
    "darks/darkRate_$(parg["tele"])_$(chip)_$(parg["mjd-start"])_$(parg["mjd-end"]).jld2";
    dark_rate = dark_im,
    dark_pix_bitmask = pix_bit_mask,
    ref_val_vec = ref_val_vec,
    cen_dar = cen_dark,
    sig_est = sig_est,
    sig_after = sig_after)

# Figures for QA
fig = PythonPlot.figure(figsize = (8, 8), dpi = 300)
ax = fig.add_subplot(1, 1, 1)
img = ax.imshow(dark_im',
    vmin = -0.2,
    vmax = 0.2,
    interpolation = "none",
    cmap = "cet_bkr",
    origin = "lower",
    aspect = "auto"
)

plt.text(0.5,
    1.01,
    "Tele: $(parg["tele"]), MJD Range: $(parg["mjd-start"])-$(parg["mjd-end"]), Chip: $(chip)\n Scatter: $(round(sig_est,digits=4)) e-/read",
    ha = "center",
    va = "bottom",
    transform = ax.transAxes)

divider = mpltk.make_axes_locatable(ax)
cax = divider.append_axes("right", size = "5%", pad = 0.05)
cbar = plt.colorbar(img, cax = cax, orientation = "vertical")
ratePath = dirNamePlots *
           "darkRate_$(parg["tele"])_$(chip)_$(parg["mjd-start"])_$(parg["mjd-end"]).png"
fig.savefig(ratePath, bbox_inches = "tight", pad_inches = 0.1)
thread("Dark stack for $(parg["tele"]) $(chip) from $(parg["mjd-start"]) to $(parg["mjd-end"]) done.")
PythonPlot.plotclose(fig)
thread("Here is the final dark rate image", ratePath)

dark_im_msk = copy(dark_im)
dark_im_msk[pix_bit_mask .& 2^3 .== 0] .= 0;
dark_im_msk[pix_bit_mask .& bad_pix_bits .!= 0] .= NaN;

totNum = length(pix_bit_mask[1:2048, 1:2048])
badVec = pix_bit_mask[1:2048, 1:2048] .& bad_pix_bits .!= 0
corrVec = (pix_bit_mask[1:2048, 1:2048] .& 2^3 .!= 0) .& .!badVec
notCorVec = (pix_bit_mask[1:2048, 1:2048] .& 2^3 .== 0)

fracBad = count(badVec) / totNum
fracCorr = count(corrVec) / totNum
fracNotCorr = count(notCorVec) / totNum

fig = PythonPlot.figure(figsize = (8, 8), dpi = 300)
ax = fig.add_subplot(1, 1, 1)
img = ax.imshow(dark_im_msk',
    vmin = -0.2,
    vmax = 0.2,
    interpolation = "none",
    cmap = cmap_w,
    origin = "lower",
    aspect = "auto"
)

plt.text(0.5,
    1.01,
    "Tele: $(parg["tele"]), MJD Range: $(parg["mjd-start"])-$(parg["mjd-end"]), Chip: $(chip)\n Bad: $(round(100*fracBad,digits=2))% Corrected: $(round(100*fracCorr,digits=2))% NoCorrection: $(round(100*fracNotCorr,digits=2))%",
    ha = "center",
    va = "bottom",
    transform = ax.transAxes
)

divider = mpltk.make_axes_locatable(ax)
cax = divider.append_axes("right", size = "5%", pad = 0.05)
cbar = plt.colorbar(img, cax = cax, orientation = "vertical")
maskPath = dirNamePlots *
           "darkRateMask_$(parg["tele"])_$(chip)_$(parg["mjd-start"])_$(parg["mjd-end"]).png"
fig.savefig(maskPath, bbox_inches = "tight", pad_inches = 0.1)
thread("Here is the final dark mask image", maskPath)
PythonPlot.plotclose(fig)

fig = PythonPlot.figure(figsize = (8, 8), dpi = 300)
ax = fig.add_subplot(1, 1, 1)

divider = mpltk.make_axes_locatable(ax)
cax = divider.append_axes("right", size = "5%", pad = 0.05)
cbar = plt.colorbar(
    mplcm.ScalarMappable(
        norm = mplcolors.Normalize(vmin = -0.2, vmax = 0.2), cmap = "cet_bkr"),
    cax = cax,
    orientation = "vertical")
im_lst = []
@showprogress for (indx, fname) in enumerate(flist)
    sname = split(fname, "_")
    tele, mjd, chiploc, expid = sname[(end - 4):(end - 1)]
    f = jldopen(fname)
    temp_im = f["dimage"]
    close(f)
    temp_im[1:2048, 1:2048] .-= ref_val_vec[indx]

    imgloc = ax.imshow(temp_im',
        vmin = -0.2,
        vmax = 0.2,
        interpolation = "none",
        cmap = "cet_bkr",
        origin = "lower",
        aspect = "auto"
    )

    ttl = plt.text(
        0.5, 1.01, "Tele: $(tele), MJD: $(mjd), Chip: $(chiploc) Expid: $(expid)",
        ha = "center", va = "bottom", transform = ax.transAxes)

    push!(im_lst, [imgloc, ttl])
end
PythonPlot.plotclose(fig)
ani = mplani.ArtistAnimation(
    fig, im_lst, interval = 300, repeat_delay = 300, blit = false)
vidPath = dirNamePlots *
          "darkStack_$(parg["tele"])_$(chip)_$(parg["mjd-start"])_$(parg["mjd-end"]).mp4"
ani.save(vidPath)
ani = nothing

len_vid = stat(vidPath).size
# if the length is bigger than 1 gigabyte, we need to upload the link to slack
if len_vid > 1e9 # 1 GB
    vidSasPath = replace(abspath(vidPath), r".*users" => sas_prefix)
    thread("Here is the video of all of the frames included in the stack: $vidSasPath")
else
    thread("Here is the video of all of the frames included in the stack", vidPath)
end

fig = PythonPlot.figure(figsize = (8, 8), dpi = 300)
ax = fig.add_subplot(1, 1, 1)

divider = mpltk.make_axes_locatable(ax)
cax = divider.append_axes("right", size = "5%", pad = 0.05)
cbar = plt.colorbar(
    mplcm.ScalarMappable(
        norm = mplcolors.Normalize(vmin = -0.2, vmax = 0.2), cmap = "cet_bkr"),
    cax = cax,
    orientation = "vertical")
im_lst = []
@showprogress for (indx, fname) in enumerate(flist)
    sname = split(fname, "_")
    tele, mjd, chiploc, expid = sname[(end - 4):(end - 1)]
    f = jldopen(fname)
    temp_im = f["dimage"]
    close(f)
    ref_val_vec[indx] = nanzeromedian(temp_im[1:2048, 1:2048])
    temp_im[1:2048, 1:2048] .-= ref_val_vec[indx]
    temp_im .-= dark_im
    imgloc = ax.imshow(temp_im',
        vmin = -0.2,
        vmax = 0.2,
        interpolation = "none",
        cmap = "cet_bkr",
        origin = "lower",
        aspect = "auto"
    )
    ttl = plt.text(
        0.5, 1.01, "Tele: $(tele), MJD: $(mjd), Chip: $(chiploc) Expid: $(expid)",
        ha = "center", va = "bottom", transform = ax.transAxes)

    push!(im_lst, [imgloc, ttl])
end
PythonPlot.plotclose(fig)
ani = mplani.ArtistAnimation(
    fig, im_lst, interval = 300, repeat_delay = 300, blit = false)
vidPath = dirNamePlots *
          "darkStackRes_$(parg["tele"])_$(chip)_$(parg["mjd-start"])_$(parg["mjd-end"]).mp4"
ani.save(vidPath)
ani = nothing

len_vid = stat(vidPath).size
# if the length is bigger than 1 gigabyte, we need to upload the link to slack
if len_vid > 1e9 # 1 GB
    vidSasPath = replace(abspath(vidPath), r".*users" => sas_prefix)
    thread("Here is the video of all of the residuals for frames included in the stack: $vidSasPath")
else
    thread("Here is the video of all of the residuals for frames included in the stack",
        vidPath)
end
