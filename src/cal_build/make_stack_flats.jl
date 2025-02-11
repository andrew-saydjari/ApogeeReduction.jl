using JLD2, ProgressMeter, ArgParse, SlackThreads, Glob, StatsBase

src_dir = "../"
include(src_dir * "/fileNameHandling.jl")
include(src_dir * "/utils.jl")
include(src_dir * "/makie_plotutils.jl")
include(src_dir * "/ap3D.jl")
include(src_dir * "/ap2Dcal.jl")

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
        "--caldir_darks"
        required = true
        help = "outdir where to look for the dark cals"
        arg_type = String
        default = ""
    end
    return parse_args(s)
end

parg = parse_commandline()

chip = parg["chip"]

# make summary plots and gifs and send to slack in outer loop
thread = SlackThread();
thread("FLAT stack for $(parg["tele"]) $(chip) from $(parg["mjd-start"]) to $(parg["mjd-end"])")

bad_dark_pix_bits = 2^2 + 2^4 + 2^5;
bad_flat_pix_bits = 2^6;
flat_frac_cut = 0.2
pcut_flat = 0.2
fx, fy = 10, 10    # Number of frequencies in x and y

#hard coded for now
gainReadCalDir = "/uufs/chpc.utah.edu/common/home/u6039752/scratch1/working/2025_02_03/"
gainMatDict = load_gain_maps(gainReadCalDir, parg["tele"], chip)

dirNamePlots = parg["flat_dir"] * "plots/"
if !ispath(dirNamePlots)
    mkpath(dirNamePlots)
end

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

design_matrix = gen_design_mat(nx, ny, fx, fy, X, Y)

#this is not going to be scalable, but current the backslash solver is too slow
flat_im_mat = zeros(2040, 2040, length(flist))
# this is using the mjd of the first exposure, which works if we are processing an ultra dark run
# but we might want to reconsider moving this inside the loop if we decide to use nightly flats moving forward
darkFlist = sort(glob("darkRate*.jld2", parg["caldir_darks"] * "darks/"))
df_dark = cal2df(darkFlist)
calPath, calFlag = get_cal_path(df_dark, parg["tele"], mjd[1], chip)

f = jldopen(calPath)
dark_pix_bitmask = f["dark_pix_bitmask"]
close(f)
bad_pix_dark = (dark_pix_bitmask[5:2044, 5:2044] .& bad_dark_pix_bits .!= 0);

desc = "Stacking flats for $(parg["tele"]) $(chip) from $(parg["mjd-start"]) to $(parg["mjd-end"])"
@showprogress desc=desc for (indx, fname) in enumerate(flist)
    sname = split(fname, "_")
    tele, mjdloc, chiploc, expidloc = sname[(end - 4):(end - 1)]
    floc = jldopen(fname)
    temp_im = floc["dimage"]
    close(floc)

    bmat = temp_im[5:2044, 5:2044]
    bmatg = bmat .* gainMatDict[chip][5:2044, 5:2044]
    b = bmatg[:]
    ref_med = nanzeromedian(b)

    # if less than 2 counts, then all pixels bad (internal flat lamps off at LCO for example)
    if ref_med < 2
        flat_im_mat[:, :, indx] .= 0
        @warn "Flat image has less than 2 e-/read on average for $(fname)"
    else
        bad_pix = copy(bad_pix_dark)
        bad_pix .|= (bmat ./ ref_med .< pcut_flat)
        good_pix = .!grow_msk2d(bad_pix; rad = 3)[:]
        # TODO mask out the low gain islands so they don't distort the illumination fit

        cvec = design_matrix[good_pix, :] \ b[good_pix]
        modImage = reshape(design_matrix * cvec, nx, ny)

        flat_im_mat[:, :, indx] = bmat ./ modImage
        flat_im_mat[:, :, indx] ./= nanzeromedian(flat_im_mat[:, :, indx]) # handles explicit gain modeling
        flat_im .+= flat_im_mat[:, :, indx]
        model_im .+= modImage
    end
end

# divide by the number of indices in the last dimension of flat_im_mat that are not all zero
nflats = count(x -> sum(x) != 0, eachslice(flat_im_mat, dims = 3))
flat_im ./= nflats
model_im ./= nflats

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
if nflats == 0
    thread("No flats with enough flux found for $(parg["tele"]) $(chip) from $(parg["mjd-start"]) to $(parg["mjd-end"])")
else
    vmin, vmax = percentile(model_im[:], [2, 98])

    let # flat image
        fig = Figure(size = (1200, 800), fontsize = 24)
        ax = Axis(fig[1, 1])
        hm = heatmap!(ax, flat_im,
            colormap = :linear_kbgyw_5_98_c62_n256,
            colorrange = (0.92, 1.08),
            interpolate = false
        )

        text!(ax,
            0.5, 1.05,
            text = "Tele: $(parg["tele"]), MJD Range: $(parg["mjd-start"])-$(parg["mjd-end"]), Chip: $(chip)",
            align = (:center, :bottom),
            space = :relative
        )

        Colorbar(fig[1, 2], hm, width = 20, height = Relative(1.0))
        colgap!(fig.layout, 1, 20)  # Set spacing between image 1 and colorbar 1
        data_aspect = diff(hm[1][])[1] / (diff(hm[2][])[1])
        colsize!(fig.layout, 1, Aspect(1, data_aspect))
        resize_to_layout!(fig)

        flatPath = dirNamePlots *
                   "flatFraction_$(parg["tele"])_$(chip)_$(parg["mjd-start"])_$(parg["mjd-end"]).png"
        save(flatPath, fig, px_per_unit = 3)
        thread("Flat stack for $(parg["tele"]) $(chip) from $(parg["mjd-start"]) to $(parg["mjd-end"]) done.")
        thread("Here is the final flat image", flatPath)
    end

    let # flat image with bad pixels masked
        flat_im_msk = copy(flat_im)
        flat_im_msk[pix_bit_mask .& bad_flat_pix_bits .!= 0] .= NaN

        totNum = length(pix_bit_mask)
        badVec = pix_bit_mask .& bad_flat_pix_bits .!= 0
        fracBad = count(badVec) / totNum

        fig = Figure(size = (1200, 800), fontsize = 24)
        ax = Axis(fig[1, 1])
        hm = heatmap!(ax, flat_im_msk,
            colormap = :linear_kbgyw_5_98_c62_n256,
            colorrange = (0.95, 1.05),
            interpolate = false,
            nan_color = :red
        )

        text!(ax,
            0.5, 1.05,
            text = "Tele: $(parg["tele"]), MJD Range: $(parg["mjd-start"])-$(parg["mjd-end"]), Chip: $(chip)\nBad: $(round(100*fracBad,digits=2))%",
            align = (:center, :bottom),
            space = :relative
        )

        Colorbar(fig[1, 2], hm, width = 20, height = Relative(1.0))
        colgap!(fig.layout, 1, 20)  # Set spacing between image 1 and colorbar 1
        data_aspect = diff(hm[1][])[1] / (diff(hm[2][])[1])
        colsize!(fig.layout, 1, Aspect(1, data_aspect))
        resize_to_layout!(fig)

        maskPath = dirNamePlots *
                   "flatMask_$(parg["tele"])_$(chip)_$(parg["mjd-start"])_$(parg["mjd-end"]).png"
        save(maskPath, fig, px_per_unit = 3)
        thread("Here is the final flat mask image", maskPath)
    end

    let # video of each frame
        if length(flist) < 5 # if there are less than 5 frames, don't make a video
            # Post each frame individually since there are few frames
            for (i, fname) in enumerate(flist)
                sname = split(fname, "_")
                teleloc, mjdloc, chiploc, expidloc = sname[(end - 4):(end - 1)]

                fig = Figure(size = (1200, 800), fontsize = 24)
                ax = Axis(fig[1, 1])
                hm = heatmap!(ax, flat_im_mat[:, :, i],
                    colormap = :linear_kbgyw_5_98_c62_n256,
                    colorrange = (0.95, 1.05),
                    interpolate = false
                )

                text!(ax,
                    0.5, 1.05,
                    text = "Tele: $(teleloc), MJD: $(mjdloc), Chip: $(chiploc) Expid: $(expidloc)",
                    align = (:center, :bottom),
                    space = :relative
                )

                Colorbar(fig[1, 2], hm, width = 20, height = Relative(1.0))
                colgap!(fig.layout, 1, 20)  # Set spacing between image 1 and colorbar 1
                data_aspect = diff(hm[1][])[1] / (diff(hm[2][])[1])
                colsize!(fig.layout, 1, Aspect(1, data_aspect))
                resize_to_layout!(fig)

                framePath = dirNamePlots *
                            "flatFrame$(i)_$(parg["tele"])_$(chip)_$(parg["mjd-start"])_$(parg["mjd-end"]).png"
                save(framePath, fig, px_per_unit = 3)
                thread("Frame $i of flat stack", framePath)
            end
        else
            # video not implemented
        end
    end

    let
        if length(flist) < 5
            for i in eachindex(flist)
                fig = Figure(size = (1200, 800), fontsize = 24)
                ax = Axis(fig[1, 1])
                hm = heatmap!(ax, (flat_im_mat[:, :, i] .- flat_im),
                    colormap = :diverging_bkr_55_10_c35_n256,
                    colorrange = (-0.02, 0.02),
                    interpolate = false
                )

                Colorbar(fig[1, 2], hm, width = 20, height = Relative(1.0))
                colgap!(fig.layout, 1, 20)  # Set spacing between image 1 and colorbar 1
                data_aspect = diff(hm[1][])[1] / (diff(hm[2][])[1])
                colsize!(fig.layout, 1, Aspect(1, data_aspect))
                resize_to_layout!(fig)

                sname = split(flist[i], "_")
                teleloc, mjdloc, chiploc, expidloc = sname[(end - 4):(end - 1)]

                text!(ax,
                    0.5, 1.05,
                    text = "Tele: $(teleloc), MJD: $(mjdloc), Chip: $(chiploc) Expid: $(expidloc)",
                    align = (:center, :bottom),
                    space = :relative
                )

                framePath = dirNamePlots *
                            "flatStackRes_$(parg["tele"])_$(chip)_$(parg["mjd-start"])_$(parg["mjd-end"])_frame$i.png"
                save(framePath, fig, px_per_unit = 3)
                thread("Frame $i residual", framePath)
            end
        else
            # video not implemented
        end
    end
end
