# invoke with:
# julia --project="../../" create_maps.jl
# from this directory

using Statistics: mean
include("../../src/ar3D.jl") # for apz2cube
include("../../src/utils.jl") # for safe_jldsave
include("../../src/makie_plotutils.jl") # make a plot of the saturation map

metadata = [("lco", 60755, 51930001, 51930005), ("apo", 60823, 52610012, 52610015)]

for (tel, sjd, first_exposure, last_exposure) in metadata
    # there are 5 samples of each saturation map (5 extra long exposures).
    # take the mean of the last read for each pixel to get the saturation map.
    saturation_maps = map(1:3) do j
        println("Processing $tel chip $("RGB"[j])")

        saturation_values_samples = map(first_exposure:last_exposure) do i
            println("  saturation map for exposure $i (last one is $last_exposure)")

            file_slug = tel == "apo" ? "apR" : "asR"
            apz_file = "/uufs/chpc.utah.edu/common/home/sdss/sdsswork/data/apogee/$tel/$sjd/$file_slug-$("abc"[j])-$i.apz"
            println("    reading $apz_file")
            data = apz2cube(apz_file)[1]
            data[:, :, end]
        end

        saturation_values = mean(saturation_values_samples)

        println("    plotting saturation map for $tel chip $("RGB"[j])")

        # Plot the saturation map for this chip
        fig = Figure(size = (840, 800), fontsize = 24)
        ax = Axis(fig[1, 1], title = "Saturation Map - $tel chip $("RGB"[j])")
        hm = heatmap!(ax, saturation_values,
            colorrange = (0, 2^16))
        Colorbar(fig[1, 2], hm)
        save("saturation_map_$(tel)_chip$("RGB"[j]).png", fig)

        println("    saving saturation map for $tel chip $("RGB"[j])")
        safe_jldsave("saturation_map_$(tel)_chip$("RGB"[j]).jld2",
            saturation_values = saturation_values)
    end
end
