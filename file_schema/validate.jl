using YAML, HDF5

"""
Checks that a file matches a spec.
"""
function validate(path, spec_path)
    spec = YAML.load_file(spec_path)
    layout = h5open(h5structure, path, "r") 

    @show spec
    println()
    @show layout

    # TODO validate filename?
    # TODO make sure everything is in the glossary

    same_up_to_order(spec["structure"], layout)
end

function h5structure(h::HDF5.File)
    map(keys(h)) do k
        if h[k] isa HDF5.Dataset
            # TODO handle checking structure of info inside datasets
            k
        else
            #TODO handle groups
        end
    end
end

function same_up_to_order(spec, actual)
    # TODO recursively handle groups, etc.
    sort(spec) == sort(actual)
end

println(validate("/uufs/chpc.utah.edu/common/home/sdss42/sdsswork/users/u6039752-1/working/2024_09_24/outdir/apred/59142/ap2D_apo_59142_a_35800020_INTERNALFLAT.jld2", 
"files/ap2D.yaml"))