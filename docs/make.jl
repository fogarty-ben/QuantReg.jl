if !haskey(ENV, "TRAVIS_CI")
    push!(LOAD_PATH, joinpath(@__DIR__, "../src"))
end
using Documenter, QuantReg, DataFrames, StatsBase, StatsModels

makedocs(sitename="QuantReg.jl",
         modules=[QuantReg, StatsBase],
         pages=[
             "Introduction" => "index.md",
             "Quickstart" => "quickstart.md",
             "Types" => "types.md",
             "Fitting models" => "fitting.md",
             "Inference" => "inference.md",
             "Example usage" => "example.md",
         ])

deploydocs(
    repo = "github.com/fogarty-ben/QuantReg.jl.git",
)