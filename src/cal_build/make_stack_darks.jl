using JLD2, ProgressMeter, ArgParse, SlackThreads

src_dir = "../"
include(src_dir*"src/utils.jl")

## Parse command line arguments
function parse_commandline()
    s=ArgParseSettings()
    @add_arg_table s begin
        "--output"
            required = true
            help = "path to output runlist file"
            arg_type = String
            default = ""
        "--chip"
            required = true
            help = "chip name (a, b, c, or all), default is all"
            arg_type = String
            default = "all"
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
    end
    return parse_args(s)
end

parg = parse_commandline()

chip_list = if parg["chip"] == "all" ["a","b","c"] else [parg["chip"]] end

# make summary plots and gifs and send to slack in outer loop
thread = SlackThread();
thread("Starting dark stack for $(parg["tele"]) from $(parg["mjd-start"]) to $(parg["mjd-end"])")

for chip in chip_list
    flist_all = glob("ap2D_$(parg["tele"])*_$(chip)_*DARK.jld2",parg["dark_dir"]*"ap2D/");
    mjd_flist = map(x->parse(Int,split(x,"_")[end-1]),flist)
    mskMJD = parg["mjd-start"] .<= mjd_flist .<= parg["mjd-end"]
    flist  = flist_all[mskMJD]

    # this is dumb and should probably be replaced by ivar weighting
    ref_val_vec = zeros(length(flist))
    dark_im = zeros(2560,2048)
    @showprogress for (indx,fname) in enumerate(flist)
        f = jldopen(fname)
        temp_im = f["dimage"]
        close(f)
        ref_val_vec[indx] = nanzeromedian(temp_im[1:2048,1:2048])
        temp_im[1:2048,1:2048] .-= ref_val_vec[indx]
        dark_im .+= temp_im
    end
    dark_im ./= length(flist);

    cen_dark = nanzeromedian(dark_im)

    dat = dark_im[1:2048,1:2048][:];
    sig_est = nanzeroiqr(dat)

    bad_pix_bits = 2^2 + 2^4 + 2^5;
    pix_bit_mask = zeros(Int,2560,2048);
    pix_bit_mask[2049:end,:] .|=2^0 # refArray
    pix_bit_mask[1:4,:] .|=2^1 # refPixels
    pix_bit_mask[2045:2048,:] .|=2^1 # refPixels
    pix_bit_mask[:,1:4] .|=2^1 # refPixels
    pix_bit_mask[:,end-3:end] .|=2^1 # refPixels
    pix_bit_mask[:,end] .|=2^2 # bad refPixels (this is hard coded, only for chipA and APO for now)
    pix_bit_mask .|= (abs.(dark_im).>sig_measure.*sig_est)*2^3
    pix_bit_mask .|= (dark_im.<-sig_bad_lower.*sig_est)*2^4
    pix_bit_mask .|= (dark_im.>sig_bad_upper.*sig_est)*2^5;

    dat = dark_im[1:2048,1:2048][pix_bit_mask[1:2048,1:2048] .& bad_pix_bits .== 0];
    sig_after = nanzeroiqr(dat)

    # save dark_pix_bitmask and dark_rate (electron per read)
    jldsave(parg["dark_dir"]*"darks/darkRate_$(tele)_$(chip)_$(parg["mjd-start"])_$(parg["mjd-end"])"; dark_rate=dark_im, dark_pix_bitmask=pix_bit_mask, ref_val_vec=ref_val_vec, cen_dar=cen_dark, sig_est=sig_est, sig_after=sig_after)

    thread("Dark stack for $(parg["tele"]) $(chip) from $(parg["mjd-start"]) to $(parg["mjd-end"]) done.")
end