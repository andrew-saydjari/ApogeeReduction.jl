
function linear_loss_fit(x, y; wporder=2, returnL2only=false)
    A = positional_poly_mat(x,porder=wporder)
    linparam = A\y
    if returnL2only
        return sum((y .- A*linparam).^2)
    else
        return y .- A*linparam, linparam
    end
end

function positional_poly_mat(x; porder=2)
    if porder == 2
        return [ones(length(x)) x x.^2]
    elseif porder == 3
        return [ones(length(x)) x x.^2 x.^3]
    elseif porder == 4
        return [ones(length(x)) x x.^2 x.^3 x.^4]
    end
end

# make sure to pass this a copy so you can modify in place
function nonlinear_loss_fit!(chipPolyParams,inparams,x,y,chipInt;wporder=2,cporder=1,returnL2only=false)
    params2ChipPolyParams!(chipPolyParams,inparams,cporder)
    xt = zeros(Float64,length(x))
    for i = 1:3
        msk = chipInt .== i
        xt[msk] .= transform_x_chips(x[msk], chipPolyParams[i,:])
    end
    return linear_loss_fit(xt, y, wporder=wporder, returnL2only=returnL2only)
end

function transform_x_chips(x,chipPolyParams)
    porder = length(chipPolyParams)-1
    if porder == 1
        return chipPolyParams[1] .+ chipPolyParams[2].*x
    elseif porder == 2
        return chipPolyParams[1] .+ chipPolyParams[2].*x .+ chipPolyParams[3].*x.^2
    elseif porder == 3
        return chipPolyParams[1] .+ chipPolyParams[2].*x .+ chipPolyParams[3].*x.^2 .+ chipPolyParams[4].*x.^3
    end
end

function params2ChipPolyParams!(chipPolyParams,inparams,cporder)
    chipPolyParams[1,:] .= inparams[1:(cporder+1),:]
    chipPolyParams[3,:] .= inparams[(cporder+2):(2*(cporder+1)),:]
    return nothing
end

function ChipPolyParams2Params(chipPolyParams)
    return vcat(chipPolyParams[1,:],chipPolyParams[3,:])
end

# Sky line wavecal
function get_and_save_sky_wavecal(fname; cporder = 1, wporder = 2)
    # initial guess for the (low-order)chip polynomial parameters
    chipPolyParams0 = [
        -1.070 1
        0 1
        1.076 1]
    outname = replace(replace(fname,"skyLine_peaks"=>"wavecal_skyline"),"_a_"=>"_")
    sky_line_uxlst, sky_line_fwlst, sky_line_chipInt = ingest_skyLines_exp(fname)
    linParams, nlParams, resid_vec = get_sky_wavecal(sky_line_uxlst, sky_line_fwlst, sky_line_chipInt, chipPolyParams0; cporder = cporder, wporder = wporder)

    chipWaveSoln = zeros(Float64,2048,300,3);
    x = 1:2048
    ximport = (x.-1024)./2048
    for chip in ["a","b","c"]
        chipIndx = getChipIndx(chip)
        for fibIndx in 1:300
            params2ChipPolyParams!(chipPolyParams0,nlParams[fibIndx,:],cporder)
            xt = transform_x_chips(ximport,chipPolyParams0[chipIndx,:])
            Ax = positional_poly_mat(xt,porder=2)
            yt = Ax*linParams[fibIndx,:]
            chipWaveSoln[:,fibIndx,chipIndx] .= yt
        end
    end

    jldsave(outname; linParams = linParams, nlParams = nlParams, resid_vec = resid_vec, chipWaveSoln = chipWaveSoln)
end

function get_sky_wavecal(sky_line_uxlst, sky_line_fwlst, sky_line_chipInt, chipPolyParams0; cporder = 1, wporder = 2)
    linParams = zeros(Float64,300,wporder+1)
    nlParams = zeros(Float64,300,2*(cporder+1))
    resid_vec = zeros(Float64,300,size(sky_line_uxlst,1))
    fill!(resid_vec,NaN)
    for i = 1:300
        xv = sky_line_uxlst[:,i]
        yv = sky_line_fwlst[:,i]
        chipIntv = sky_line_chipInt[:,i]
        msk = .!isnan.(xv)
        chipPolyParams = copy(chipPolyParams0)
        inparams = ChipPolyParams2Params(chipPolyParams)
        nonlinear_loss_fit_partial(inparams) = nonlinear_loss_fit!(chipPolyParams, inparams, xv[msk], yv[msk], chipIntv[msk]; wporder=wporder,cporder=cporder,returnL2only=true)
        res = optimize(nonlinear_loss_fit_partial, inparams, LBFGS(), Optim.Options(show_trace=false))
        nlParamsOpt = Optim.minimizer(res)
        linResid, linParamsOpt = nonlinear_loss_fit!(chipPolyParams,nlParamsOpt,xv[msk],yv[msk],chipIntv[msk];wporder=wporder,cporder=cporder,returnL2only=false)
        nlParams[i,:] = nlParamsOpt
        linParams[i,:] = linParamsOpt
        resid_vec[i,msk] = linResid
    end
    return linParams, nlParams, resid_vec
end

function ingest_skyLines_file(fileName)
    # Read in sky line peaks data
    f = h5open(fileName, "r+")
    sky_line_mat_clean = try
        read(f["sky_line_mat_clean"])
    catch
        println(fileName)
        read(f["sky_line_mat_clean"])
    end
    close(f)
    sky_line_xlst = (sky_line_mat_clean[:,1,:] .- 1024)./2048
    sky_line_wlst = sky_line_mat_clean[:,2,:]
    return sky_line_xlst, sky_line_wlst
end

# this takes in a filename and replaces the chip index (make "a" default approx)
function ingest_skyLines_exp(fname)
    sky_line_uxlst = Matrix{Float64}[]
    sky_line_fwlst = Matrix{Float64}[]
    sky_line_chipInt = Matrix{Int}[]
    for chip in ["a","b","c"]
        fnameloc = replace(fname,"_a_" => "_$(chip)_")
        if isfile(fnameloc)
            sky_line_xlst, sky_line_wlst = ingest_skyLines_file(fnameloc)
            chipIndx = getChipIndx(chip)
            push!(sky_line_uxlst,sky_line_xlst)
            push!(sky_line_fwlst,sky_line_wlst)
            push!(sky_line_chipInt,chipIndx*ones(Int,size(sky_line_wlst)))
        else
            push!(sky_line_uxlst,[])
            push!(sky_line_fwlst,[])
            push!(sky_line_chipInt,[])
        end
    end
    sky_line_uxlst = vcat(sky_line_uxlst...);
    sky_line_fwlst = vcat(sky_line_fwlst...);
    sky_line_chipInt = vcat(sky_line_chipInt...);
    
    msk_large_scatter = dropdims(nanzeroiqr(sky_line_uxlst,2) .> 0.002,dims=2)
    sky_line_uxlst[msk_large_scatter,:] .= NaN
    sky_line_fwlst[msk_large_scatter,:] .= NaN
    # dims are num_sky_lines x num_fibers
    return sky_line_uxlst, sky_line_fwlst, sky_line_chipInt
end