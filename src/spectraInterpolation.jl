## Spectra Interpolation Module

function find_yinx(x::AbstractVector,y::AbstractVector)
    lx = length(x)
    ly = length(y)
    out = zeros(Int,ly)
    @inbounds for i=1:length(y)
        ind = searchsortedfirst(x,y[i])
        if ind == 1
            out[i] = ind
        elseif ind > lx
            out[i] = ind-1
        else
            if abs(x[ind]-y[i]) < abs(x[ind-1]-y[i])
                out[i] = ind
            else
                out[i] = ind-1
            end
        end
    end
    return out
end

function returnWeights_inv(obsCoordall::AbstractVector,obsBitMsk::Vector{Int},pixindx::AbstractVector,targVal::Float64,cindx::Int;
        kernsize::Int=4,linFallBack::Bool=true)
    obslen = length(obsCoordall)
    if obslen == 0 #when is this happening and why?
        return zeros(Int,2*kernsize), NaN*ones(2*kernsize)
    end
    diffwav = diff(obsCoordall[maximum([1,(cindx-1)]):minimum([(cindx+1),obslen])])
    diffpixind = diff(pixindx[maximum([1,(cindx-1)]):minimum([(cindx+1),obslen])])
    pscale = minimum(abs.(diffwav./diffpixind))
    offset = (obsCoordall[cindx].-targVal)/pscale
    cbit = obsBitMsk[cindx]
    if (cbit .& 2^1)!=0
        cchip = 1
    elseif (cbit .& 2^2)!=0
        cchip = 2
    elseif (cbit .& 2^3)!=0
        cchip = 3
    else
        cchip = 0
        # println("NO CHIP?")
    end

    indvec = (-kernsize:kernsize) .+ cindx
    offvec = (-kernsize:kernsize) .+ offset
    msk = (1 .<= indvec .<= obslen) # within bounds range
    msk .&= (-kernsize .<= offvec .<= kernsize) # within kernel bounds
    indvecr = indvec[msk]

    mskb = ((obsBitMsk[indvecr] .& 2^cchip).!=0) # same chip mask
    mskb .&= ((obsBitMsk[indvecr] .& 2^4).==0) #bad pix mask
    mskb .&= (-kernsize .<= (pixindx[indvecr].- pixindx[cindx]) .<= kernsize)
    indvecrr = indvecr[mskb]
    
    if offset == 0
        return cindx, 1.0
    end
    if (count(mskb) >= 2*kernsize)
        if (maximum(diff((pixindx[indvecr].- pixindx[cindx])[mskb]))==1)
            wvec = Interpolations.lanczos.(offvec[msk][mskb],kernsize)
            return indvecrr, wvec
        end
    end
    if linFallBack & (count(mskb) > 1)
        pixOffset = pixindx[indvecrr].- pixindx[cindx]
        obsOffset = offvec[msk][mskb]
        lindx = findlast(obsOffset.<=0)
        rindx = findfirst(obsOffset.>=0)
        if !isnothing(lindx) & !isnothing(rindx)
            if (pixOffset[lindx]>=-1) & (pixOffset[rindx]<=1)
                offlst = offvec[msk][mskb]
                totoff = offlst[rindx].-offlst[lindx]
                return [indvecrr[lindx],indvecrr[rindx]],[1-abs(offlst[lindx])/totoff, 1-abs(offlst[rindx])/totoff]
            end
        end
    end
    return zeros(Int,2*kernsize), NaN*ones(2*kernsize)
end

function generateInterpMatrix_sparse_inv(waveobs::AbstractVector,obsBitMsk::Vector{Int},wavemod::AbstractVector,pixindx::AbstractVector;kernsize::Int=4,linFallBack::Bool=true)
    obslen = length(waveobs)
    modlen = length(wavemod)
    cindx = find_yinx(waveobs,wavemod)
    row, col, val = Int[], Int[], Float64[]
    for (modind, modval) in enumerate(wavemod)
        indxvec, wvec = returnWeights_inv(waveobs,obsBitMsk,pixindx,modval,cindx[modind],kernsize=kernsize,linFallBack=linFallBack)
        if !isnan(wvec[1]) .& (wvec[1].!=1.0)
            wvec ./= sum(wvec)
            nz = (wvec.!=0)
            push!(row, (modind.*ones(Int,length(indxvec[nz])))...)
            push!(col, indxvec[nz]...)
            push!(val, wvec[nz]...)
        end
    end
    return sparse(row,col,val,modlen,obslen)
end