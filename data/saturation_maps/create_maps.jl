# invoke with:
# julia --project="../../" create_maps.jl
# from this directory

using FITSIO, ProgressMeter
using Statistics: mean
include("../../src/utils.jl") # for safe_jldsave
include("../../src/makie_plotutils.jl") # for plotting the saturation map

"""
Get the final counts from an apz file. We could use apz2cube in ar3D.jl, but this is faster and uses
less memory.
"""
function last_read_from_apz(fname)
    counts = FITS(fname) do f
        avg_dcounts = read(f[2])

        total_counts = read(f[3])
        @showprogress desc="    Reading $(splitpath(fname)[end])" for i in 4:(length(f))
            total_counts .+= read(f[i]) .+ avg_dcounts
        end
        total_counts
    end

    # make sure it's safe, then convert to UInt16
    @assert all(counts .>= 0)
    @assert all(counts .<= 2^16)
    UInt16.(counts)
end

function save_plot(saturation_values, tel, j)
    fig = Figure(size = (840, 800), fontsize = 24)
    ax = Axis(fig[1, 1], title = "Saturation Map - $tel chip $("RGB"[j])")
    hm = heatmap!(ax, saturation_values, colorrange = (0, 2^16))
    Colorbar(fig[1, 2], hm)
    save("saturation_map_$(tel)_chip$("RGB"[j]).png", fig)
end

# LCO

metadata = [
    ("lco", 60755, 51930001, 51930005) # APO is handled below because the apz files are missing frames
]
for (tel, sjd, first_exposure, last_exposure) in metadata
    # there are 5 samples of each saturation map (5 extra long exposures).
    # take the mean of the last read for each pixel to get the saturation map.
    saturation_maps = map(1:3) do j
        println("Processing $tel chip $("RGB"[j])")

        saturation_values_samples = map(first_exposure:last_exposure) do i
            println("  saturation map for exposure $i (last one is $last_exposure)")

            file_slug = tel == "apo" ? "apR" : "asR"
            apz_file = "/uufs/chpc.utah.edu/common/home/sdss/sdsswork/data/apogee/$tel/$sjd/$file_slug-$("abc"[j])-$i.apz"
            last_read_from_apz(apz_file)
        end

        saturation_values = UInt16.(floor.(mean(saturation_values_samples)))

        println("    plotting saturation map for $tel chip $("RGB"[j])")
        save_plot(saturation_values, tel, j)

        println("    saving saturation map for $tel chip $("RGB"[j])")
        safe_jldsave("saturation_map_$(tel)_chip$("RGB"[j]).h5";
            no_metadata = true,
            saturation_values = saturation_values)
    end
    break
end

# APO

datadir = "/uufs/chpc.utah.edu/common/home/sdss42/sdsswork/users/u0914351/60823/ics"
# there are more exposures
# https://data.sdss5.org/sas/sdsswork/mwm/sandbox/airflow-ApogeeReduction.jl/daily/ApogeeReduction.jl/metadata/observing_log_viewer/?sjd=60823&site=apo
# but let's keep it simple
exposures = [52610013]
n_reads = [450]

println("reading in final APO counts")
fused_counts = map(exposures, n_reads) do exposure, n_read
    println("processing exposure APO exposure $exposure")

    path = joinpath(datadir, "apRaw-$exposure-$(lpad(n_read, 3, '0')).fits")
    if !isfile(path)
        throw("missing $path")
    end

    total_counts = FITS(path) do f
        read(f[1])
    end

    total_counts
end
#mean_counts = round.(UInt16, mean(fused_counts))
mean_counts = fused_counts[1]

for (j, xcoords) in enumerate([1:2048, 2049:4096, 4097:6144])
    save_plot(mean_counts[xcoords, :], "apo", j)

    safe_jldsave(
        "saturation_map_apo_chip$("RGB"[j]).h5";
        saturation_values = mean_counts[xcoords, :],
        no_metadata = true)
end
