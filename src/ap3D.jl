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

function refcorr(dcubedat)
    # subtracts reference array with proper read ordering
    dcubedat_out = copy(dcubedat[1:2048,1:2048,..])
    dcubedat_out[1:512,:,..].-= dcubedat[2049:end,..]
    dcubedat_out[513:1024,:,..].-= dcubedat[end:-1:2049,..]
    dcubedat_out[1025:1536,:,..].-= dcubedat[2049:end,..]
    dcubedat_out[1537:2048,:,..].-= dcubedat[end:-1:2049,..]
    return dcubedat_out
end

function vert_ref_edge_corr(dcubedat_out)
    # change to inplace after validate
    # choosing NOT to do per quadrant without any clear cause to do so
    # seeing frequencies indicating that we need SIRS.jl
    dcubedat_out_v = copy(dcubedat_out)
    top = dropdims(mean(dcubedat_out[:,1:4,:],dims=(1,2)),dims=(1,2))
    bot = dropdims(mean(dcubedat_out[:,end-3:end,:],dims=(1,2)),dims=(1,2))
    dcubedat_out_v .-= reshape((top+bot)./2,(1,1,:))
    return dcubedat_out_v
end

# Need to revisit with SIRS.jl, this currently appears to add more bias than it removes
function horz_ref_edge_corr(dcubedat_out_v)
    # change to inplace after validate
    # choosing NOT to do per quadrant (or flipping?) without any clear cause to do so
    dcubedat_out_vh = copy(dcubedat_out_v)
    horzbias = dropdims(median(vcat(dcubedat_out_v[1:4,:,:],dcubedat_out_v[end-3:end,:,:]),dims=1),dims=1);
    dcubedat_out_vh .-= horzbias
    return dcubedat_out_vh
end

function dcs(dcubedat,gainMat,readVarMat;firstind=2)
    dimage = gainMat.*(dcubedat[:,:,end].-dcubedat[:,:,firstind])
    # bad to use measured flux as the photon noise
    ivarimage = 1 ./(2 .*readVarMat .+ abs.(dimage))
    return dimage, ivarimage
end