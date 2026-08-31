using Polynomials: Polynomial, roots
using SparseArrays: spzeros, sparse, nnz
using JLD2: jldopen

# for a given fiber index, return the LSF
# we may want this to be date dependent later
# could do that by making lsf_mjd a range and checking inclusion
# fstep runs fsteprng = 0.5:(-0.1):(-0.4) for 1:10 for the subpixel shifts, 6 being no shift
function get_lsf_matrix(adjfibindx, lsfpath; n_pad = 50, fstep = 0.0)
    fiberindx = adjFiberIndx2FiberIndx(adjfibindx)

    #define the parameters for the usual logUniWaveAPOGEE
    outwave_log_step = 6.0e-6
    len_outwave = 8700
    outwave_log_min = 4.17825 + fstep*outwave_log_step
    outwave_log_max = outwave_log_min + outwave_log_step*(len_outwave-1)
    # need to pad the outwave for the LSF to function as intended
    # outwave_pad[n_pad+1:end-n_pad] == outwave
    outwave_pad = 10 .^ range(start = outwave_log_min-n_pad*outwave_log_step, 
                            step = outwave_log_step, 
                            length = len_outwave + 2*n_pad)

    #now get the pixel edges that we need to integrate over
    outwave_bin_edges = 10 .^ ((outwave_log_min-0.5*outwave_log_step):outwave_log_step:(outwave_log_max+0.5*outwave_log_step))
    outwave_pad_bin_edges = 10 .^ range(start = outwave_log_min-n_pad*outwave_log_step-0.5*outwave_log_step, step = outwave_log_step, length = len_outwave + 2*n_pad + 1)

    #define the modelwave centers that the LSF is evaluated at
    min_wave, max_wave, dlam = 15000, 17000, 1/100
    modelwave = min_wave:dlam:max_wave

    return build_sparse_lsf_mat(modelwave,outwave_pad_bin_edges,fiberindx,lsfpath;
                                force_mode_0=true,min_lsf_weight=1e-15,n_pad=n_pad);
end

function build_sparse_lsf_mat(modelwave,outwave_bin_edges,fiberindx,lsfpath;
    force_mode_0=true,min_lsf_weight=1e-15,n_pad=0)
    # outwave_bin_edges is usually defined by logUniWaveAPOGEE 
    # but gives the EDGES of the pixel bins to integrate over (i.e. ~8701 pixels in length)
    # tele,fibIndx give the telescope and fiber index 
    # modelwave set the central wavelengths where the LSF is evaluated 

    # Returns sparse_lsf, a ~8700×200001 sparse matrix of the LSF evaluated at all the modelwave centers
    # that bins onto the 8700 outwave pixels

    #read in necessary LSF parameters for the chosen telescope and fiber index
    freq_bounds,freq_zeropoint,
    freq_scale,lsf_width_coeffs,
    gh_height_coeffs,throughput_coeffs,
    background_coeffs = get_fib_lsf_params(fiberindx,lsfpath)
    GH_order = size(gh_height_coeffs,1)

    nrows = length(outwave_bin_edges) - 1
    ncols = length(modelwave)
    row = Int[]
    col = Int[]
    val = Float64[]
    sizehint!(row, 40 * ncols)
    sizehint!(col, 40 * ncols)
    sizehint!(val, 40 * ncols)

    for j in 1:ncols
        center_wave = modelwave[j]
        curr_dlam_bin_edges = outwave_bin_edges .- center_wave
        curr_keep = findall(abs.(curr_dlam_bin_edges) .< 5)
        if length(curr_keep) < 2
            continue
        end
        
        curr_weights = fib_lam_lsf(center_wave,curr_dlam_bin_edges[curr_keep],
                                    freq_bounds,freq_zeropoint,freq_scale,
                                    lsf_width_coeffs,gh_height_coeffs,
                                    throughput_coeffs,background_coeffs,
                                    force_mode_0=force_mode_0)
        good_weights = curr_weights .> min_lsf_weight
        curr_weights[.!good_weights] .= 0
        s = sum(curr_weights)
        iszero(s) && continue
        curr_weights ./= s
        rw = curr_keep[begin:end-1][good_weights]
        vw = curr_weights[good_weights]
        nn = length(rw)
        append!(row, rw)
        append!(col, fill(j, nn))
        append!(val, vw)
    end
    sparse_lsf = sparse(row, col, val, nrows, ncols)

    # dot with a ones vector to get the normalizing factor 
    # (i.e. LSF applied to ones vector should give back ones vector)
    norm_factor = sparse(1:(length(outwave_bin_edges)-1), 
                        1:(length(outwave_bin_edges)-1), 
                        1 ./ (sparse_lsf * ones(length(modelwave))))
    if n_pad > 0
        return norm_factor[n_pad+1:end-n_pad,n_pad+1:end-n_pad] * sparse_lsf[n_pad+1:end-n_pad,:]
    else
        return norm_factor * sparse_lsf
    end
end

function get_fib_lsf_params(fiberindx, lsfpath)

    f = jldopen(lsfpath)
    fpi_cavity_size = f["fpi_cavity_size"]
    fpi_m0_offset = f["fpi_m0_offset"]
    freq_zeropoint = f["freq_zeropoint"]
    freq_scale = f["freq_scale"]
    core_HWHM_zeropoint = f["core_HWHM_zeropoint"]
    core_HWHM_slope = f["core_HWHM_slope"]
    throughput_x_scale = f["throughput_x_scale"]
    throughput_y_scale = f["throughput_y_scale"]
    background_scale = f["background_scale"]
    background_zeropoint = f["background_zeropoint"]

    n_fnames = f["n_fnames_used"]
    n_lsf_width_coeffs = f["porder_lsf_widths"]+1
    n_throughput_coeffs = f["porder_throughput_coeffs"]+1
    GH_order = f["GH_order"]
    n_gh_height_coeffs = f["porder_gh_height_coeffs"]+1
    n_background_coeffs = f["porder_background_coeffs"]+1

    #size (2 (lower,upper bound), N_CHIPS)
    freq_bounds = f["freq_bounds"][:,:,fiberindx]
    #size (n_fnames,N_CHIPS)
    throughput_mults = f["throughput_mults"][:,:,fiberindx]
    #size (n_throughput_coeffs,N_CHIPS)
    lsf_width_coeffs = f["lsf_width_coeffs"][:,:,fiberindx]
    #size (n_throughput_coeffs,N_CHIPS)
    throughput_coeffs = f["throughput_coeffs"][:,:,fiberindx]
    #size (GH_order,n_gh_height_coeffs,N_CHIPS)
    gh_height_coeffs = f["gh_height_coeffs"][:,:,:,fiberindx]
    #size (n_fnames,n_background_coeffs,N_CHIPS)
    background_coeffs = f["background_coeffs"][:,:,:,fiberindx]
    close(f)

    return freq_bounds,freq_zeropoint,freq_scale,lsf_width_coeffs,gh_height_coeffs,throughput_coeffs,background_coeffs
end

function fib_lam_lsf(lam,dlam_bin_edges,
    freq_bounds,freq_zeropoint,freq_scale,
    lsf_width_coeffs,gh_height_coeffs,throughput_coeffs,background_coeffs;
    force_mode_0 = false, speed_of_light=299792.458) #3e5 km/s
    #lam,dlam_bin_edges in Angstroms
    freq = (1e-9) * (speed_of_light * 1e3) / (lam / 1e10)
    dfreqs_bin_edges = (1e-9) .* (speed_of_light * 1e3) ./ ((dlam_bin_edges .+ lam) ./ 1e10) .- freq

    return fib_freq_lsf(freq,dfreqs_bin_edges,
        freq_bounds,freq_zeropoint,freq_scale,
        lsf_width_coeffs,gh_height_coeffs,throughput_coeffs,background_coeffs,force_mode_0=force_mode_0)
end

function fib_freq_lsf(freq,dfreqs_bin_edges,
        freq_bounds,freq_zeropoint,freq_scale,
        lsf_width_coeffs,gh_height_coeffs,throughput_coeffs,background_coeffs;
        force_mode_0 = false)
    #freq in GHz, to calculate the LSF at
    #dfreqs_bin_edges in GHz, list of distances away from freq to measure the LSF weights at
    #(freq_bounds,freq_zeropoint,freq_scale,lsf_width_coeffs,gh_height_coeffs)
    #                           are outputs from get_fib_lsf_params function
    GH_order = size(gh_height_coeffs,1)

    weights = zeros(length(dfreqs_bin_edges)-1)
    chip_weights = freq_interp_weights(freq,freq_bounds)
    #preallocate
    curr_gh_weights = zeros(length(dfreqs_bin_edges)-1)
    z_vals = zeros(length(dfreqs_bin_edges))

    scaled_freq = (freq .- freq_zeropoint) ./ freq_scale
    for chipIndx in 1:N_CHIPS
    if chip_weights[chipIndx] == 0
    continue
    end

    @views curr_lsf_width = Polynomial(lsf_width_coeffs[:,chipIndx])(scaled_freq)
    if force_mode_0
    root_poly = _gauss_hermite_return_poly(0+2)
    for gh_ind in 1:GH_order
    @views curr_gh_height = Polynomial(gh_height_coeffs[gh_ind,:,chipIndx])(scaled_freq)
    root_poly += curr_gh_height * _gauss_hermite_return_poly(gh_ind+2)
    end

    root_vals = roots(root_poly)
    good_roots = (imag.(root_vals) .== 0)
    real_roots = real.(root_vals[good_roots])
    mode_freq = real_roots[argmin(abs.(real_roots))] * curr_lsf_width
    else
    mode_freq = 0.0
    end

    # @views z_vals .= (dfreqs_bin_edges .- mode_freq) ./ curr_lsf_width
    @views z_vals .= (dfreqs_bin_edges .+ mode_freq) ./ curr_lsf_width

    @views curr_gh_weights .= diff(int_gauss_hermite_term(z_vals, 0; mean = 0.0, width = 1.0))
    if GH_order > 0
    for gh_ind in 1:GH_order
    @views curr_gh_height = Polynomial(gh_height_coeffs[gh_ind,:,chipIndx])(scaled_freq)
    @views curr_gh_weights .+= curr_gh_height .* 
                        diff(int_gauss_hermite_term(z_vals, gh_ind; mean = 0.0, width = 1.0))
    end
    end

    # println(freq," ",chipIndx," ",curr_lsf_width," ",sum(curr_gh_weights)," ",)

    #USER NEEDS TO ENSURE THAT THE dfreqs LIST IS LARGE ENOUGH TO CONTAIN THE FULL INTEGRAL        
    @views weights .+= chip_weights[chipIndx] .* curr_gh_weights ./ sum(curr_gh_weights)
    end

    return weights
end

function _gauss_hermite_return_poly(n)
    if n == 1
        Polynomial([1])
    elseif n == 2
        Polynomial([0,1])
    elseif n == 3   
        Polynomial([-1,0,1])
    elseif n == 4
        Polynomial([0,-3,0,1])
    elseif n == 5
        Polynomial([3,0,-6,0,1])
    elseif n == 6
        Polynomial([0,15,0,-10,0,1])
    elseif n == 7
        Polynomial([-15,0,45,0,-15,0,1])
    elseif n == 8
        Polynomial([0,-105,0,105,0,-21,0,1])
    elseif n == 9
        Polynomial([105,0,-420,0,210,0,-28,0,1])
    elseif n == 10
        Polynomial([0,945,0,-1260,0,378,0,-36,0,1])
    elseif n == 11
        Polynomial([-945,0,4725,0,-3150,0,630,0,-45,0,1])
    else
        throw(ArgumentError("n must be between 0 and 11"))
    end
end

function freq_interp_weights(freq,freq_bounds)
    if freq <= freq_bounds[2,1]
        #then on first chip
        weights = [1,0,0]
    elseif (freq >= freq_bounds[1,2]) & (freq <= freq_bounds[2,2])
        #then on middle chip
        weights = [0,1,0]
    elseif freq >= freq_bounds[1,3]
        #then on last chip
        weights = [0,0,1]
    elseif (freq > freq_bounds[2,1]) & (freq < freq_bounds[1,2])
        #then between first and second chips
        x1 = freq_bounds[2,1]
        x2 = freq_bounds[1,2]
        scale_dist = (freq-x1)/(x2-x1)
        weights = [1-scale_dist,scale_dist,0]
    else
        #then between second and third chips
        x1 = freq_bounds[2,2]
        x2 = freq_bounds[1,3]
        scale_dist = (freq-x1)/(x2-x1)
        weights = [0,1-scale_dist,scale_dist]
    end
    return weights
end
