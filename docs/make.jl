push!(LOAD_PATH, joinpath(@__DIR__, "../src"))
using Documenter, QuantReg

makedocs(sitename="QuantReg.jl Documentation")