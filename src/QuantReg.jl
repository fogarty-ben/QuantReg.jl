module QuantReg

using DataFrames, Distributions, Format, GLM, Gurobi, JuMP, LinearAlgebra, StatsBase, StatsModels,
       Statistics

export QuantRegModel, fit, compute_inf

include("qrmodel.jl")
include("fitqr.jl")
include("ci.jl")

end # module
