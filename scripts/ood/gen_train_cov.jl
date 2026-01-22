using InteractiveUtils
versioninfo(); flush(stdout);
@time "Package activation" begin
    import Pkg
    Pkg.instantiate()
    Pkg.precompile() # no need for Pkg.activate("./") because of invocation w/ environment
    using Distributed, ArgParse, TimerOutputs
end

## Parse command line arguments
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--tele"
        required = false
        help = "telescope name (apo or lco)"
        arg_type = String
        default = ""
        "--outdir"
        required = true
        help = "output directory"
        arg_type = String
        default = ""
        "--ood_dir"
        required = true
        help = "directory where ood data is stored"
        arg_type = String
        default = ""
        "--runlist"
        required = true
        help = "path name to hdf5 file with keys specifying list of exposures to run"
        arg_type = String
        default = ""
        "--runname"
        required = true
        help = "name of the run (specifically almanac file)"
        arg_type = String
        default = "test"
    end
    return parse_args(s)
end

parg = parse_commandline()

proj_path = dirname(Base.active_project()) * "/"
include(joinpath(proj_path, "src/makie_plotutils.jl"))

if "SLURM_NTASKS" in keys(ENV)
    using SlurmClusterManager
    addprocs(SlurmManager(), exeflags = ["--project=$proj_path"])
else
    addprocs(64)
end

using JLD2
@everywhere begin
    using ParallelDataTransfer
end

@passobj 1 workers() parg # make it available to all workers
@passobj 1 workers() proj_path

tele_list = if parg["runlist"] != ""
    load(parg["runlist"], "tele")
else
    [parg["tele"]]
end
unique_teles = unique(tele_list)
mskTele = tele_list .== parg["tele"];

mjd_list = if parg["runlist"] != ""
    load(parg["runlist"], "mjd")
else
    [parg["mjd"]]
end
unique_mjds = unique(mjd_list[mskTele]);

dfindx_list = if parg["runlist"] != ""
    load(parg["runlist"], "dfindx")
else
    [parg["dfindx"]]
end;

# predefined choices
@everywhere begin
    nsub = 10
    tpix = 16
    plot_pdir = parg["ood_dir"]
end

@everywhere begin
    exposure_types = ["fpi", "thar", "une"] # dark, internalflat, domeflat, quartzflat
    function get_1d_name_FPI(expid,df,mjd,outdir)
        if (df.image_type[expid] == "arclamp") & (df.lamp_thar[expid] == 0) & (df.lamp_une[expid] == 0)
            return outdir * "/apred/$(mjd)/" * get_1d_name(expid, df, cal = true) *".h5"
        else
            return nothing
        end
    end

    function get_1d_name_thar(expid,df,mjd,outdir)
        if (df.image_type[expid] == "arclamp") & (df.lamp_thar[expid] == 1)
            return outdir * "/apred/$(mjd)/" * get_1d_name(expid, df, cal = true) *".h5"
        else
            return nothing
        end
    end

    function get_1d_name_une(expid,df,mjd,outdir)
        if (df.image_type[expid] == "arclamp") & (df.lamp_une[expid] == 1)
            return outdir * "/apred/$(mjd)/" * get_1d_name(expid, df, cal = true) *".h5"
        else
            return nothing
        end
    end

    function get_1d_name_dark(expid,df,mjd,outdir)
        if (df.image_type[expid] == "dark")
            return outdir * "/apred/$(mjd)/" * get_1d_name(expid, df, cal = true) * ".h5"
        else
            return nothing
        end
    end

    function get_1d_name_domeflat(expid,df,mjd,outdir)
        if (df.image_type[expid] == "domeflat")
            return outdir * "/apred/$(mjd)/" * get_1d_name(expid, df, cal = true) * ".h5"
        else
            return nothing
        end
    end
    exposure_1dname_functions = [get_1d_name_FPI, get_1d_name_thar, get_1d_name_une] # get_1d_name_dark, get_1d_name_domeflat]
end

# predefined functions
@everywhere begin
    using ApogeeReduction
    import ApogeeReduction: read_almanac_exp_df, get_1d_name, bad_pix_bits
    using LowRankOps, LinearAlgebra, StatsBase
    using HDF5, JLD2, DelimitedFiles, ProgressMeter
    
    yranges = []
    xranges = []
    for i in 1:4
        push!(xranges, 128 + 512*(i-1) : 128 + tpix - 1 + 512*(i-1))
        push!(yranges, 256 + 512*(i-1) : 256 + tpix - 1 + 512*(i-1))
    end

    function get_testsubimage(fname)
        f = h5open(fname, "r")
        dimage = zeros(Float64, (length(xranges[1]), sum(length.(yranges))))
        pix_bitmask = zeros(Int64,  (length(xranges[1]), sum(length.(yranges))))
        for i in 1:4
            dimage[:,(1:tpix) .+ tpix*(i-1)] .= f["dimage"][xranges[i],yranges[i]]
            pix_bitmask[:,(1:tpix) .+ tpix*(i-1)] .= f["pix_bitmask"][xranges[i],yranges[i]]
        end
        mask_good = (pix_bitmask .& bad_pix_bits) .== 0;
        close(f)
        return dimage[mask_good], mask_good
    end
    
    function wood_precomp_mult(matList)
        Ainv = matList[1]
        V = matList[2]
        AinvV = Ainv*V
        return [(AinvV)*inv(I+V'*(AinvV))]
    end

    function wood_fxn_mult(matList,precompList,x)
        Ainv = matList[1]
        V = matList[2]
        arg1 = precompList[1]
        return Ainv*(x - V*(arg1'*x))
    end

    function get_chi2(fname, mat2D_cen, Vmat)
        f = h5open(fname, "r")
        dimage = zeros(Float64, (length(xranges[1]), sum(length.(yranges))))
        ivarimage = zeros(Float64, (length(xranges[1]), sum(length.(yranges))))
        pix_bitmask = zeros(Int64,  (length(xranges[1]), sum(length.(yranges))))
        for i in 1:4
            dimage[:,(1:tpix) .+ tpix*(i-1)] .= f["dimage"][xranges[i],yranges[i]]
            ivarimage[:,(1:tpix) .+ tpix*(i-1)] .= f["ivarimage"][xranges[i],yranges[i]]
            pix_bitmask[:,(1:tpix) .+ tpix*(i-1)] .= f["pix_bitmask"][xranges[i],yranges[i]]
        end
        close(f)
        mask_good = (pix_bitmask .& bad_pix_bits) .== 0;
        dimage_cen = (dimage .- mat2D_cen)[mask_good]
        npix = count(mask_good)
        Ainv = Diagonal(ivarimage[mask_good])
        Vmat_loc = Vmat[mask_good[:],:]
        Ctotinv = LowRankMultMat([Ainv,Vmat_loc],wood_precomp_mult,wood_fxn_mult);
        chi2 = dimage_cen'*(Ctotinv*dimage_cen)/npix

        # flux spread metric
        spread = std(dimage_cen)
        chi_fluxscatter = (spread - p50_spread)
        if chi_fluxscatter < 0
            chi_fluxscatter /= (p50_spread - p16_spread)
        else
            chi_fluxscatter /= (p84_spread - p50_spread)
        end
        return chi2, chi_fluxscatter
    end
end

# need an outer loop over the exposure types
for (exposure_type, exposure_1dname_function) in zip(exposure_types, exposure_1dname_functions)
    println("Making bad list for $(exposure_type) exposures"); flush(stdout);
    for chip in CHIP_LIST
        println("... on chip $(chip)"); flush(stdout);
        list1Dexp = []
        list1Dexp_all = []
        for mjd in unique_mjds
            df = read_almanac_exp_df(
                joinpath(parg["outdir"], "almanac/$(parg["runname"]).h5"), parg["tele"], mjd)
            exposure_1dname_function_partial(expid) = exposure_1dname_function(expid, df, mjd, parg["outdir"]);
            if any(df.image_type .== "object")
                mskMJD = (mjd_list .== mjd) .& mskTele
                local1D = exposure_1dname_function_partial.(dfindx_list[mskMJD])
                push!(list1Dexp, filter(!isnothing, local1D))
                push!(list1Dexp_all, filter(!isnothing, local1D))
            else
                mskMJD = (mjd_list .== mjd) .& mskTele
                local1D = exposure_1dname_function_partial.(dfindx_list[mskMJD])
                push!(list1Dexp_all, filter(!isnothing, local1D))
            end
        end

        train1Dfiles = convert(Vector{String}, vcat(list1Dexp...));
        train2Dfiles = replace.(train1Dfiles, "ar1Dcal" => "ar2Dcal");
        all1Dfiles = convert(Vector{String}, vcat(list1Dexp_all...));
        all2Dfiles = replace.(all1Dfiles, "ar1Dcal" => "ar2Dcal");

        oodfname = "$(exposure_type)_roughprior_$(chip)_$(nsub).h5"

        mat2D = zeros(length(xranges[1]), sum(length.(yranges)),length(train2Dfiles));
        cnt2D = zeros(length(xranges[1]), sum(length.(yranges)),length(train2Dfiles));
        
        pout = @showprogress pmap(get_testsubimage, train2Dfiles)
        
        for (i,(dimage,mask_good)) in enumerate(pout)
            mat2D[mask_good,i] .= dimage
            cnt2D[mask_good,i] .+= 1
        end

        # rough mask based on flux scatter
        spread_vec = dropdims(std(mat2D, dims=(1,2)), dims=(1,2))
        p16_spread, p50_spread, p84_spread = percentile(spread_vec, [16, 50, 84])
        chi_fluxscatter = (spread_vec.-p50_spread)
        msk_neg = (chi_fluxscatter .< 0)
        msk_pos = (chi_fluxscatter .> 0)
        chi_fluxscatter[msk_neg] ./= (p50_spread.-p16_spread)
        chi_fluxscatter[msk_pos] ./= (p84_spread.-p50_spread)

        # this seems to work fine for fpi, arclamps, darks
        # struggles with domeflats (TODO)
        msk_indist = -2.5 .< chi_fluxscatter .< 30;

        # get centers of distribution
        mat2D_cnt = dropdims(sum(cnt2D[:,:,msk_indist],dims=3),dims=3)
        mat2D_cen = dropdims(sum(mat2D[:,:,msk_indist],dims=3),dims=3)./(mat2D_cnt .+ (mat2D_cnt.==0));
        mat2Dvec = reshape(mat2D,:,length(train2Dfiles))[:,msk_indist]
        mat2Dvec_cnt = reshape(cnt2D,:,length(train2Dfiles))[:,msk_indist];
        mat2Dvec_cen = reshape(mat2D_cen,:);

        # build covariance matrix
        Covprior = (mat2Dvec .- mat2Dvec_cen) * (mat2Dvec .- mat2Dvec_cen)'
        norm_weights = (mat2Dvec_cnt * mat2Dvec_cnt')
        Covprior ./= (norm_weights .+ (norm_weights.==0));

        # SVD for low-rank approximation
        SF = svd(Covprior);
        EVEC = zeros(size(mat2Dvec,1),size(SF.U,2))
        EVEC.=SF.U;
        Vmat = EVEC[:,1:nsub]*Diagonal(sqrt.(SF.S[1:nsub]));

        # save results
        h5write(oodfname, "Vmat", Vmat)
        h5write(oodfname, "mat2D_cen", mat2D_cen);
        h5write(oodfname, "p16_spread", p16_spread);
        h5write(oodfname, "p50_spread", p50_spread);
        h5write(oodfname, "p84_spread", p84_spread);

        @everywhere begin
            Vmat = load($oodfname, "Vmat")
            mat2D_cen = load($oodfname, "mat2D_cen")
            p16_spread = load($oodfname, "p16_spread")
            p50_spread = load($oodfname, "p50_spread")
            p84_spread = load($oodfname, "p84_spread")
        end

        @everywhere get_chi2_partial(fname) = get_chi2(fname, mat2D_cen, Vmat);
        pout = @showprogress pmap(get_chi2_partial, all2Dfiles);
        chi2out = zeros(length(all2Dfiles));
        chifluxscatterout = zeros(length(all2Dfiles));
        for (i,(chi2, chi_fluxscatter)) in enumerate(pout)
            chi2out[i] = chi2
            chifluxscatterout[i] = chi_fluxscatter
        end
        pchi2 = sortperm(chi2out,rev=true);
        pscatter = sortperm(chifluxscatterout,rev=true);
        msk_good_chifluxscatter = (-2.5 .< chifluxscatterout .< 30);
        msk_good_chi2 = (0.3 .< chi2out .< 10);
        msk_good = msk_good_chifluxscatter .& msk_good_chi2;
        msk_bad = .!msk_good;

        # need to save
        bad_fname = "$(exposure_type)_$(chip)_bad.csv"
        # need tele, mjd, expid, chip
        tele_vec = map(x->split(basename(x), "_")[2], all2Dfiles[msk_bad]);
        mjd_vec = map(x->parse(Int, split(basename(x), "_")[3]), all2Dfiles[msk_bad]);
        expid_vec = map(x->parse(Int, split(basename(x), "_")[4]), all2Dfiles[msk_bad]);
        chip_vec = map(x->split(basename(x), "_")[5], all2Dfiles[msk_bad]);
        writedlm(bad_fname, hcat(tele_vec, mjd_vec, expid_vec, chip_vec));

        # need to do plots
        badinds = findall(msk_bad);
        just_goodinds = []

        # Helper to compute safe ranges
        function safe_range(start, stop, n)
            # Clamp range to [1, n] and handle possible negatives/missing
            start = max(1, min(n, start))
            stop = max(1, min(n, stop))
            if start > stop
                Int[]
            else
                collect(start:stop)
            end
        end

        # Just_goodinds block, robust to bounds
        npscatter = length(pscatter)
        npchi2 = length(pchi2)

        fst_indx = findfirst(chifluxscatterout[pscatter] .< 30)
        if !isnothing(fst_indx)
            inds = safe_range(fst_indx, min(fst_indx+4, npscatter), npscatter)
            push!(just_goodinds, pscatter[inds])
        end

        fst_indx = findlast(-2.5 .< chifluxscatterout[pscatter])
        if !isnothing(fst_indx)
            start_idx = fst_indx+1
            inds = safe_range(start_idx, min(start_idx+4, npscatter), npscatter)
            if !isempty(inds)
                push!(just_goodinds, pscatter[inds])
            end
        end

        fst_indx = findfirst(chi2out[pchi2] .< 10)
        if !isnothing(fst_indx)
            inds = safe_range(fst_indx, min(fst_indx+4, npchi2), npchi2)
            push!(just_goodinds, pchi2[inds])
        end

        fst_indx = findlast(chi2out[pchi2] .< 10)
        if !isnothing(fst_indx)
            start_idx = fst_indx+1
            inds = safe_range(start_idx, min(start_idx+4, npchi2), npchi2)
            if !isempty(inds)
                push!(just_goodinds, pchi2[inds])
            end
        end

        just_goodinds = unique(vcat(just_goodinds...));

        center_indx = pchi2[length(pchi2)÷2]
        center_fname = all2Dfiles[center_indx];
        center_dimage = load(center_fname, "dimage")
        center_pix_bitmask = load(center_fname, "pix_bitmask")
        center_mask_good = (center_pix_bitmask .& bad_pix_bits) .== 0;
        center_dimage[.!(center_mask_good)] .= NaN
        vmin, vmax = percentile(center_dimage[center_mask_good], [1, 99])

        for (ind_list, class_name) in zip([just_goodinds, badinds], ["just_good", "bad"])
            plot_dir = joinpath(plot_pdir, "$(exposure_type)_$(chip)_$(class_name)");
            mkpath(plot_dir);
            for indx in ind_list
                fname = all2Dfiles[indx];
                dimage = load(fname, "dimage");
                pix_bitmask = load(fname, "pix_bitmask")
                mask_good = (pix_bitmask .& bad_pix_bits) .== 0;
                dimage[.!(mask_good)] .= NaN

                fig = Figure(size=(800,800),fontsize=30)
                ax = Axis(fig[1,1], title="$(class_name): $(basename(fname))\nChi2: $(round(chi2out[indx], digits=2)) $(msk_good_chi2[indx] ? "Pass" : "Fail"), ChiFluxScatter: $(round(chifluxscatterout[indx], digits=2)) $(msk_good_chifluxscatter[indx] ? "Pass" : "Fail")")
                hm = heatmap!(ax, dimage,
                    colormap=ColorSchemes.cmr_chroma,
                    colorrange=(vmin, vmax),
                    nan_color=(93/300, 101/300, 106/300))

                ax = Axis(fig[2,1],)
                hm = heatmap!(ax, center_dimage,
                    colormap=ColorSchemes.cmr_chroma,
                    colorrange=(vmin, vmax),
                    nan_color=(93/300, 101/300, 106/300))
                data_aspect = size(center_dimage,1) / size(center_dimage,2)
                colsize!(fig.layout, 1, Aspect(1, data_aspect))

                ax = Axis(fig[1,2],)
                hm = heatmap!(ax, dimage[1024:1536,1024:1536],
                    colormap=ColorSchemes.cmr_chroma,
                    colorrange=(vmin, vmax),
                    nan_color=(93/300, 101/300, 106/300))

                ax = Axis(fig[2,2],)
                hm = heatmap!(ax, center_dimage[1024:1536,1024:1536],
                    colormap=ColorSchemes.cmr_chroma,
                    colorrange=(vmin, vmax),
                    nan_color=(93/300, 101/300, 106/300))
                data_aspect = size(center_dimage[1024:1536,1024:1536],1) / size(center_dimage[1024:1536,1024:1536],2)
                colsize!(fig.layout, 2, Aspect(1, data_aspect))
                
                Colorbar(fig[1:2, 3], hm, width = 10, height = Relative(1.0))
                resize_to_layout!(fig)
                save(joinpath(plot_dir, "$(basename(fname))_$(class_name).png"), fig)
            end
        end
    end
end