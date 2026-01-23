using HDF5
using Statistics
using CairoMakie
using Printf

"""
    z2v(z; c=299792.458)

Convert redshift to velocity.
"""
function z2v(z; c=299792.458)
    return ((z + 1)^2 - 1) / ((z + 1)^2 + 1) * c
end

"""
    pix2v(x; delLog=6e-6)

Convert pixel offset to velocity.
"""
function pix2v(x; delLog=6e-6)
    z = 10^(x * delLog) - 1
    return z2v(z)
end

"""
    average_rvs(rv_pixoff_file::HDF5.File, full_list_file::HDF5.File)

Average out RVs across sdss_id + adjfiberindx pairs.

Returns: avg_rv, std_rv, nexp, fiber_ids, sdss_ids
"""
function average_rvs(rv_pixoff_file::HDF5.File, full_list_file::HDF5.File)
    sdss_id_all = read(full_list_file, "sdss_id")
    sdss_id_unq = unique(sdss_id_all)
    
    avg_rv = Float64[]
    std_rv = Float64[]
    nexp = Int[]
    fiber_ids = Int[]
    sdss_ids = typeof(sdss_id_all[1])[]
    
    rv_pixoff = read(rv_pixoff_file, "RV_pixoff_final")
    rv = pix2v.(rv_pixoff)
    adjfiberindx_all = read(full_list_file, "adjfiberindx")
    
    for sdss_id in sdss_id_unq
        idx = findall(sdss_id_all .== sdss_id)
        fids = unique(adjfiberindx_all[idx])
        
        if length(fids) == 1
            push!(avg_rv, mean(filter(!isnan, rv[idx])))
            push!(std_rv, std(filter(!isnan, rv[idx])))
            push!(nexp, length(idx))
            push!(fiber_ids, fids[1])
            push!(sdss_ids, sdss_id)
        else
            for f in fids
                idx_sub = findall((sdss_id_all .== sdss_id) .& (adjfiberindx_all .== f))
                push!(avg_rv, mean(filter(!isnan, rv[idx_sub])))
                push!(std_rv, std(filter(!isnan, rv[idx_sub])))
                push!(nexp, length(idx_sub))
                push!(fiber_ids, f)
                push!(sdss_ids, sdss_id)
            end
        end
    end
    
    return avg_rv, std_rv, nexp, fiber_ids, sdss_ids
end

"""
    plot_hist(avg_rv::Vector{Float64}, save_path::String)

Plot histogram of average RVs.
"""
function plot_hist(avg_rv::Vector{Float64}, save_path::String)
    fig = Figure(size=(1400, 1400))
    ax = Axis(fig[1, 1],
              xlabel="Radial Velocity (Avg. by sdss_id + adjfiberindx) [km/s]",
              ylabel="N",
              xlabelsize=28,
              ylabelsize=28,
              xticklabelsize=22,
              yticklabelsize=22)
    
    hist!(ax, avg_rv, bins=100)
    
    save(joinpath(save_path, "daily_rv_avg_hist.png"), fig)
    println("Saved: daily_rv_avg_hist.png")
end

"""
    plot_rv_vs_fiber(avg_rv::Vector{Float64}, fiber_ids::Vector{Int}, save_path::String)

Plot 2D histogram of RV vs fiber ID using heatmap.
"""
function plot_rv_vs_fiber(avg_rv::Vector{Float64}, fiber_ids::Vector{Int}, save_path::String)
    # Filter out NaN values
    valid_idx = .!isnan.(avg_rv)
    avg_rv_clean = avg_rv[valid_idx]
    fiber_ids_clean = fiber_ids[valid_idx]
    
    println("Creating 2D histogram with $(length(avg_rv_clean)) points...")
    
    # Create bins manually
    nbins = 30
    rv_min, rv_max = extrema(avg_rv_clean)
    fib_min, fib_max = extrema(fiber_ids_clean)
    
    rv_edges = range(rv_min, rv_max, length=nbins+1)
    fib_edges = range(fib_min, fib_max, length=nbins+1)
    
    # Compute 2D histogram manually
    H = zeros(Int, nbins, nbins)
    
    for (rv, fib) in zip(avg_rv_clean, fiber_ids_clean)
        rv_idx = searchsortedfirst(rv_edges, rv) - 1
        fib_idx = searchsortedfirst(fib_edges, fib) - 1
        
        rv_idx = clamp(rv_idx, 1, nbins)
        fib_idx = clamp(fib_idx, 1, nbins)
        
        H[rv_idx, fib_idx] += 1
    end
    
    # Replace zeros with NaN for log scale
    H_log = log10.(Float64.(H))
    H_log[H .== 0] .= NaN
    
    # Create bin centers
    rv_centers = (rv_edges[1:end-1] .+ rv_edges[2:end]) ./ 2
    fib_centers = (fib_edges[1:end-1] .+ fib_edges[2:end]) ./ 2
    
    fig = Figure(size=(2200, 1800))
    ax = Axis(fig[1, 1],
              xlabel="Radial Velocity (Avg. by sdss_id)",
              ylabel="adjfiberindx",
              xlabelsize=28,
              ylabelsize=28,
              xticklabelsize=22,
              yticklabelsize=22)
    
    # Note: No transpose here - H_log dimensions match (rv, fib)
    hm = heatmap!(ax, rv_centers, fib_centers, H_log, colormap=:viridis)
    Colorbar(fig[1, 2], hm, label="log10(N)", labelsize=28, ticklabelsize=22)
    
    save(joinpath(save_path, "daily_rv_avg_vs_fiber.png"), fig)
    println("Saved: daily_rv_avg_vs_fiber.png")
end

"""
    plot_rv_vs_fiber_violin(avg_rv::Vector{Float64}, fiber_ids::Vector{Int}, 
                           tele::String, save_path::String; std_thresh::Float64=15.0)

Plot range plot of average RVs per fiber with median line, 16th-84th percentile shaded bar, and extrema whiskers. 
Creates 6 subplots in one PNG file.
"""
function plot_rv_vs_fiber_violin(avg_rv::Vector{Float64}, 
                                fiber_ids::Vector{Int}, 
                                tele::String, 
                                save_path::String; 
                                std_thresh::Float64=15.0)
    if tele == "apo"
        fiber_plot = collect(1:300)
    else
        fiber_plot = collect(301:600)
    end
    
    chunk = 50
    nrows = length(fiber_plot) ÷ chunk
    
    println("Creating range plot with $nrows subplots using CairoMakie...")
    
    # Create figure with subplots - wider and shorter
    fig = Figure(size=(10000, 800 * nrows))
    
    # Matplotlib colors
    blue_color = RGBf(0.12, 0.47, 0.71)
    red_color = RGBf(1.0, 0.0, 0.0)
    
    for i in 1:nrows
        println("  Processing panel $i/$nrows...")
        fids = fiber_plot[(i-1)*chunk+1 : i*chunk]
        
        ax = Axis(fig[i, 1],
                  xlabel="adjfiberindx",
                  ylabel="Radial Velocity (Avg. by sdss_id)",
                  limits=((0.5, length(fids) + 0.5), (-315, 315)),
                  xticks=(1:length(fids), string.(fids)),
                  xticklabelrotation=0.0,
                  xlabelsize=32,
                  ylabelsize=32,
                  xticklabelsize=24,
                  yticklabelsize=24)
        
        for (j, f) in enumerate(fids)
            rvi = avg_rv[fiber_ids .== f]
            if length(rvi) > 0
                # Calculate statistics
                median_val = median(rvi)
                p16 = quantile(rvi, 0.16)  # 16th percentile
                p84 = quantile(rvi, 0.84)  # 84th percentile
                min_val = minimum(rvi)
                max_val = maximum(rvi)
                stdval = std(rvi)  # Still use std for coloring
                
                # Determine color based on std
                color = (!isnan(stdval) && stdval < std_thresh) ? red_color : blue_color
                
                # Draw shaded bar for 16th-84th percentile range
                poly!(ax, 
                     Point2f[(j-0.4, p16), (j+0.4, p16), 
                            (j+0.4, p84), (j-0.4, p84)],
                     color=(color, 0.6),
                     strokecolor=color,
                     strokewidth=1.5)
                
                # Draw thin lines (whiskers) from percentile edges to extrema
                linesegments!(ax, [Point2f(j, min_val), Point2f(j, p16)],
                             color=:black,
                             linewidth=1.5)
                linesegments!(ax, [Point2f(j, p84), Point2f(j, max_val)],
                             color=:black,
                             linewidth=1.5)
                
                # Mark extrema with small horizontal caps
                linesegments!(ax, [Point2f(j-0.2, min_val), Point2f(j+0.2, min_val)],
                             color=:black,
                             linewidth=2)
                linesegments!(ax, [Point2f(j-0.2, max_val), Point2f(j+0.2, max_val)],
                             color=:black,
                             linewidth=2)
                
                # Draw horizontal line for median (on top)
                linesegments!(ax, [Point2f(j-0.4, median_val), Point2f(j+0.4, median_val)],
                             color=:black,
                             linewidth=3)
                
                # Add count annotation with larger font
                text!(ax, j, -305, 
                     text="n=$(length(rvi))", 
                     align=(:center, :bottom),
                     fontsize=20)
            end
        end
        
        # Add grid
        vlines!(ax, 0.5:1:length(fids)+0.5, color=(:gray, 0.2), linewidth=0.5)
        hlines!(ax, -300:100:300, color=(:gray, 0.3), linewidth=0.5)
    end
    
    println("  Saving combined plot...")
    save(joinpath(save_path, "daily_rv_range_plot_avg_rv_fiber.png"), fig)
    
    if isfile(joinpath(save_path, "daily_rv_range_plot_avg_rv_fiber.png"))
        file_size = stat(joinpath(save_path, "daily_rv_range_plot_avg_rv_fiber.png")).size
        println("✓ Saved: daily_rv_range_plot_avg_rv_fiber.png ($(round(file_size/1024/1024, digits=2)) MB)")
    end
    
    # Print summary of low-std fibers
    low_std_count = 0
    println("\nFibers with std < $std_thresh km/s (shown in red):")
    for f in fiber_plot
        rvi = avg_rv[fiber_ids .== f]
        if length(rvi) > 0
            stdval = std(filter(!isnan, rvi))
            if !isnan(stdval) && stdval < std_thresh
                println("  Fiber $f: n=$(length(rvi)), std=$(round(stdval, digits=2)) km/s")
                low_std_count += 1
            end
        end
    end
    println("Total: $low_std_count fibers with low std")
end

# Main execution
function main()
    # Grab the files
    path = "."
    pix_off_file = "arMADGICS_out_RV_pixoff_final.h5"
    rv_pixoff = h5open(joinpath(path, pix_off_file), "r")
    
    full_list_file = "full_list_info.h5"
    full_list = h5open(joinpath(path, full_list_file), "r")
    
    # Get average RVs
    avg_rv, std_rv, nexp, fiber_ids, sdss_ids = average_rvs(rv_pixoff, full_list)
    
    # Create summary plots
    plot_path = "plots"
    save_path = joinpath(path, plot_path)
    
    # Create directory if it doesn't exist
    if !isdir(save_path)
        mkdir(save_path)
    end
    
    plot_hist(avg_rv, save_path)
    plot_rv_vs_fiber(avg_rv, fiber_ids, save_path)
    
    # Get telescope info
    tele_data = read(full_list, "tele")
    tele = String(tele_data[1])
    
    plot_rv_vs_fiber_violin(avg_rv, fiber_ids, tele, save_path)
    
    # Close files
    close(rv_pixoff)
    close(full_list)
    
    println("\nAll plots created successfully!")
end

# Run main if script is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
