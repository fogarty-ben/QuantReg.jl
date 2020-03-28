if !haskey(ENV, "TRAVIS_CI")
    push!(LOAD_PATH, joinpath(@__DIR__, "../src"))
end
using Documenter, QuantReg

makedocs(sitename="QuantReg.jl",
         pages = [
             "index.md",
             "Quickstart Guide" => "quickstart.md",
         ])

deploydocs(
    repo = "github.com/fogarty-ben/QuantReg.jl.git",
)