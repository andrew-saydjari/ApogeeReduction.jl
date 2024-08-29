# Handling the 3D data cube

function apz2cube(fname)
    f = FITS(fname)
    hdr = read_header(f[2])
    avg_dcounts = read(f[2])
    cubedat = zeros(Float32,size(avg_dcounts)...,length(f)-2) #use float bc of NaN
    cubedat[:,:,1] .= read(f[3])
    for i=2:(length(f)-2)
        cubedat[:,:,i] .= read(f[i+2]) .+ avg_dcounts .+ cubedat[:,:,i-1]
    end
    close(f)
    return cubedat, hdr
end

# feeling against this sort of subtraction, but we do see why one might want to do it in the darks
function refcorr(dcubedat)
    # subtracts reference array with proper read ordering
    dcubedat_out = copy(dcubedat[1:2048,1:2048,..])
    dcubedat_out[1:512,:,..].-= dcubedat[2049:end,..]
    dcubedat_out[513:1024,:,..].-= dcubedat[end:-1:2049,..]
    dcubedat_out[1025:1536,:,..].-= dcubedat[2049:end,..]
    dcubedat_out[1537:2048,:,..].-= dcubedat[end:-1:2049,..]
    return dcubedat_out
end

function refarray_zpt!(dcubedat)
    mean_ref_vec = mean(dcubedat[2049:end,:,..],dims=(1,2))
    dcubedat .-= mean_ref_vec
    mean_ref_val = mean(dcubedat[2049:end,..])
    dcubedat .-= mean_ref_val
    return
end

function vert_ref_edge_corr_amp!(dcubedat_out)
    dcubedat_out[1:512,..] .-= mean([mean(dcubedat_out[1:512,1:4,..]),mean(dcubedat_out[1:512,end-3:end-1,..])])
    dcubedat_out[513:1024,..] .-= mean([mean(dcubedat_out[513:1024,1:4,..]),mean(dcubedat_out[513:1024,end-3:end-1,..])])
    dcubedat_out[1025:1536,..] .-= mean([mean(dcubedat_out[1025:1536,1:4,..]),mean(dcubedat_out[1025:1536,end-3:end-1,..])])
    dcubedat_out[1537:2048,..] .-= mean([mean(dcubedat_out[1537:2048,1:4,..]),mean(dcubedat_out[1537:2048,end-3:end-1,..])])
    dcubedat_out[2049:end,..] .-= mean([mean(dcubedat_out[2049:end,1:4,..]),mean(dcubedat_out[2049:end,end-3:end-1,..])])
    return
end

function dcs(dcubedat,gainMat,readVarMat;firstind=1)
    dimage = gainMat.*(dcubedat[:,:,end].-dcubedat[:,:,firstind])
    # bad to use measured flux as the photon noise
    ivarimage = 1 ./(2 .*readVarMat .+ abs.(dimage))
    return dimage, ivarimage, nothing # no chisq, just mirroring sutr_tb
end

function sutr_tb(dcubedat,gainMat,readVarMat;firstind=1,good_diffs=nothing)
    # Last editted by Kevin McKinnon on Aug 20, 2024 
    # based on Tim Brandt SUTR python code (https://github.com/t-brandt/fitramp)
    # datacube has shape (npix_x,npix_y,n_reads)
    # read_var_mat has shape (npix_x,npix_y)
    # good_diffs is boolean and has shape (npix_x,npix_y,n_reads-1)
	        
    # assumes all images are sequential (ie separated by one read time)

    dimages = gainMat.*(dcubedat[:,:,firstind+1:end]-dcubedat[:,:,firstind:end-1])

    if isnothing(good_diffs)
        good_diffs = ones(Bool, size(dimages))
    end

    n_reads = size(dimages, 3)+1
    ndiffs = n_reads-1
    read_times = range(start=1,stop=n_reads,step=1)
    mean_t = read_times
    tau = read_times
    N = ones(Float64, n_reads)

    delta_t = mean_t[2:end] .- mean_t[1:end-1]

    alpha_readnoise = (1 ./ N[1:end-1] .+ 1 ./ N[2:end]) ./ delta_t.^2
    beta_readnoise = -1 ./ N[2:end-1] ./ (delta_t[2:end] .* delta_t[1:end-1])

    alpha_phnoise = (tau[1:end-1] .+ tau[2:end] .- 2 .* mean_t[1:end-1]) ./ delta_t.^2
    beta_phnoise = (mean_t[2:end-1] .- tau[2:end-1]) ./ (delta_t[2:end] .* delta_t[1:end-1])

    final_countrates = zeros(Float64, size(dcubedat, 1), size(dcubedat, 2))
    final_vars = zeros(Float64, size(dcubedat, 1), size(dcubedat, 2))
    final_chisqs = zeros(Float64, size(dcubedat, 1), size(dcubedat, 2))

    #slice the datacube to analyze each row sequentially to keep runtime down
    for s_ind in 1:size(dimages, 1)
        diffs = transpose(dimages[s_ind, :, :]) # shape = (ndiffs, npix)
        read_var = readVarMat[s_ind, :] # shape = (npix)
	    read_var = reshape(read_var,(1,size(read_var,1)))
        diffs2use = transpose(good_diffs[s_ind, :, :]) # shape = (ndiffs, npix)

        n_repeat = 2 # DO NOT CHANGE THIS!
        for r_ind in 1:n_repeat
            #repeat the process after getting a good first guess
	        #to remove a bias (as discussed in the paper)
            if r_ind == 1
                #initialize first guess of countrates with mean
                countrateguess = sum(diffs .* diffs2use, dims=1) ./ sum(diffs2use, dims=1)
		        countrateguess .*= (countrateguess .> 0)
            else
                countrateguess = final_countrates[s_ind, :] .* (final_countrates[s_ind, :] .> 0)
		        countrateguess = reshape(countrateguess, (1,size(countrateguess,1)))
            end
	    
            # Elements of the covariance matrix
            alpha = countrateguess .* alpha_phnoise .+ read_var .* alpha_readnoise
            beta = countrateguess .* beta_phnoise .+ read_var .* beta_readnoise

            # rescale the covariance matrix to a determinant of order 1 to
            # avoid possible overflow/underflow.  The uncertainty and chi
            # squared value will need to be scaled back later.
            scale = exp.(mean(log.(alpha), dims=1))

            alpha ./= scale
            beta ./= scale

            ndiffs, npix = size(alpha)
            # Mask resultant differences that should be ignored.  This is half
            # of what we need to do to mask these resultant differences; the
            # rest comes later.

            d = diffs .* diffs2use
            beta .= beta .* diffs2use[2:end, :] .* diffs2use[1:end-1, :]

            # All definitions and formulas here are in the paper.
            theta = ones(Float64, ndiffs + 1, npix)
            theta[2, :] .= alpha[1, :]
            for i in 3:ndiffs + 1
                theta[i, :] .= alpha[i - 1, :] .* theta[i - 1, :] .- beta[i - 2, :].^2 .* theta[i - 2, :]
            end

            phi = ones(Float64, ndiffs + 1, npix)
            phi[ndiffs, :] .= alpha[ndiffs, :]
            for i in ndiffs - 1:-1:1
                phi[i, :] .= alpha[i, :] .* phi[i + 1, :] .- beta[i, :].^2 .* phi[i + 2, :]
            end

            sgn = ones(Float64, ndiffs, npix)
            sgn[1:2:end, :] .= -1

            Phi = zeros(Float64, ndiffs, npix)
            for i in ndiffs - 1:-1:1
                Phi[i, :] .= Phi[i + 1, :] .* beta[i, :] .+ sgn[i + 1, :] .* beta[i, :] .* phi[i + 2, :]
            end

            # This one is defined later in the paper and is used for jump
            # detection and pedestal fitting.

            PhiD = zeros(Float64, ndiffs, npix)
            for i in ndiffs - 1:-1:1
                PhiD[i, :] .= (PhiD[i + 1, :] .+ sgn[i + 1, :] .* d[i + 1, :] .* phi[i + 2, :]) .* beta[i, :]
            end

            Theta = zeros(Float64, ndiffs, npix)
            Theta[1, :] .= -theta[1, :]
            for i in 2:ndiffs-1
                Theta[i, :] .= Theta[i - 1, :] .* beta[i - 1, :] .+ sgn[i, :] .* theta[i, :]
            end

            ThetaD = zeros(Float64, ndiffs + 1, npix)
            ThetaD[2, :] .= -d[1, :] .* theta[1, :]
            for i in 2:ndiffs-1
                ThetaD[i + 1, :] .= beta[i - 1, :] .* ThetaD[i, :] .+ sgn[i, :] .* d[i, :] .* theta[i, :]
            end

            beta_extended = ones(Float64, ndiffs, npix)
            beta_extended[2:end, :] .= beta

            # C' and B' in the paper

            theta_ndiffs = theta[ndiffs + 1, :]
            theta_ndiffs = reshape(theta_ndiffs,(1,size(theta_ndiffs,1)))
            dC = sgn ./ theta_ndiffs .* (phi[2:end, :] .* Theta .+ theta[1:end-1, :] .* Phi)
            dC .*= diffs2use

            dB = sgn ./ theta_ndiffs .* (phi[2:end, :] .* ThetaD[2:end, :] .+ theta[1:end-1, :] .* PhiD)

            # {\cal A}, {\cal B}, {\cal C} in the paper
            A = 2 * sum(d .* sgn ./ theta_ndiffs .* beta_extended .* phi[2:end, :] .* ThetaD[1:end-1, :], dims=1)
            A .+= sum(d.^2 .* theta[1:end-1, :] .* phi[2:end, :] ./ theta_ndiffs, dims=1)

            B = sum(d .* dC, dims=1)
            C = sum(dC, dims=1)

            countrate = B ./ C
            chisq = (A .- B.^2 ./ C) ./ scale
            var = scale ./ C
            weights = dC ./ C

            #use first countrate measurement to improve alpha,beta definitions
            #and extract better, unbiased countrate
            final_countrates[s_ind,:] .= countrate[begin,:]
            final_vars[s_ind,:] .= var[begin,:]
	        final_chisqs[s_ind,:] = chisq[begin,:]
	end
    end
    
    #to be similar to CDS outputs, we want to multiply by the number of differences
    count_mean = ndiffs .* final_countrates
    count_var = ndiffs.^2 .* final_vars
    
    return count_mean,count_var.^(-1),final_chisqs
end


