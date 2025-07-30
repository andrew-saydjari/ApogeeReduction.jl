using Pkg
Pkg.activate(dirname(@__DIR__))
using ApogeeReduction
Pkg.activate(@__DIR__)
using Documenter

DocMeta.setdocmeta!(ApogeeReduction, :DocTestSetup, :(using ApogeeReduction); recursive=true)

makedocs(;
    modules=[ApogeeReduction],
    authors="Andrew Saydjari <saydjari@cfa.harvard.edu> and contributors",
    sitename="ApogeeReduction.jl",
    format=Documenter.HTML(;
        canonical="https://andrew-saydjari.github.io/ApogeeReduction.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Error Bits" => "error_bits.md",
    ],
)

deploydocs(;
    repo="github.com/andrew-saydjari/ApogeeReduction.jl",
    devbranch="main",
)
