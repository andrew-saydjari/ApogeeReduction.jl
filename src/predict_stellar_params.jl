#!/usr/bin/env julia
# predict_stellar_params_pycall.jl
# Use Python's XGBoost from Julia via PyCall to avoid the Julia XGBoost bug

using ArgParse
using HDF5
using DataFrames
using PyCall

# Import Python's xgboost
const xgb = pyimport("xgboost")
const np = pyimport("numpy")

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--input"
            help = "Path to arMADGICS starLineCof file"
            arg_type = String
            required = true
        "--models"
            help = "Directory containing XGBoost JSON models"
            arg_type = String
            default = "data/stellar_params/xgboost_models"
        "--outdir"
            help = "Output directory for predictions"
            arg_type = String
            required = true
        "--features_dset"
            help = "The feature key in the HDF5 file for the star coefficients"
            arg_type = String
            default = "x_starLineCof_v0"
        "--params"
            help = "Comma-separated list of stellar parameters (e.g. teff,logg,fe_h,raw_alpha_m_atm)"
            arg_type = String
            default = "teff,logg,fe_h,raw_alpha_m_atm"
    end
    return parse_args(s)
end


# load features matrix from HDF5
function load_features_from_h5(h5path::String, dset_path::String="x_starLineCof_v0")
    X_raw = h5open(h5path,"r") do f
            read(f[dset_path])
    end
    # ensure matrix shape: rows = samples, cols = features
    X_raw = Array{Float64}(X_raw)
    if ndims(X_raw) == 1
        X_raw = reshape(X_raw, :, 1)
    end
    # Transpose to get (n_samples, n_features)
    X = permutedims(Array{Float64,2}(X_raw))
    return X
end


function load_all_models(param_names::Vector{String}, models_dir::String)
    models = Dict{String, PyObject}()
    for param in param_names
        model_file = joinpath(models_dir, "model_$(param).json")
        if !isfile(model_file)
            error("Model file not found for parameter '$param': $model_file")
        end
        # Load using Python's XGBoost
        println("Loading model: $param")
        bst = xgb.Booster()
        bst.load_model(model_file)
        models[param] = bst
    end
    return models
end


function predict_all(models::Dict{String, PyObject}, X::AbstractMatrix)
    n = size(X, 1)
    df = DataFrame(Index = 1:n)
    
    # Convert Julia array to numpy array
    X_np = np.array(X)
    
    # Create DMatrix
    dmat = xgb.DMatrix(X_np)
    
    for (name, booster) in models
        println("Predicting for: $name")
        # Predict using Python's XGBoost
        preds_np = booster.predict(dmat)
        # Convert back to Julia array
        preds = convert(Array, preds_np)
        df[!, Symbol(name)] = preds
    end
    return df
end


function save_predictions_h5(h5path::AbstractString, df::DataFrame, param_names::Vector{String};
                             mode::String="w")
    @assert mode in ("w", "a")

    # open file (w = truncate/create; a = read/write append)
    open_mode = mode == "w" ? "w" : "r+"
    h5open(h5path, open_mode) do f
        # write each parameter as its own dataset
        for pname in param_names
            col = pname
            if !(col in names(df))
                @warn "Parameter '$pname' not present in DataFrame; skipping HDF5 dataset creation"
                continue
            end
            arr = convert(Vector{Float64}, df[:, col])

            # overwrite dataset if exists
            if haskey(f, pname)
                delete!(f, pname)
            end
            write(f, pname, arr)
        end
    end
end


# main function
function main()
    args = parse_commandline()
    input_h5 = args["input"]
    models_dir = args["models"]
    outdir = args["outdir"]
    features_dset = args["features_dset"]
    params = split(args["params"], ",")
    params = strip.(params)
    params = String.(params)

    println("="^60)
    println("Using Python XGBoost via PyCall")
    println("Python XGBoost version: ", xgb.__version__)
    println("="^60)

    # Load features
    X = load_features_from_h5(input_h5, features_dset)
    println("Loaded features: size=$(size(X))")

    # Load models
    println("\nLoading models from: $models_dir")
    models = load_all_models(params, models_dir)
    println("Found models: ", join(keys(models), ", "))

    # Predict
    println("\nRunning predictions for ", join(params, ", "))
    df_preds = predict_all(models, X)

    # Create output directory
    mkpath(outdir)

    # write to HDF5 file
    h5_outfn = joinpath(outdir, "xgboost_stellar_params.h5")
    save_predictions_h5(h5_outfn, df_preds, params; mode="w")
    println("\nWrote predictions to: $h5_outfn")
end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end