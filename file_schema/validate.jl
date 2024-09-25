using YAML, HDF5

function validate(path, spec_path)
    spec = YAML.load_file(spec_path)
    layout = h5open(h5structure, path, "r") 

    # TODO make sure these strctures are the same
    # TODO make sure everything is in the glossary
end

function h5structure(h::Union{HDF5.File, HDF5.Group})
    Dict([k=>h5structure(h[k]) for k in keys(h)])
end
function h5structure(h::HDF5.Dataset)
    #TODO
end

