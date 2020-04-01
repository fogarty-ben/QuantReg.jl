module QuantReg

using DataFrames, Distributions, Format, GLM, JuMP, LinearAlgebra, QuantReg_jll, StatsBase,
      StatsModels, Statistics

if haskey(ENV, "GUROBI_HOME") # ensure Gurobi is installed before loading package
    using Gurobi
end

export coef, coefnames, coeftable, confint, dof, dof_residual, fitted, @formula, isfitted,
       islinear, nobs, modelmatrix, response, responsename, residuals, stderr

export
    # Model type
    QuantRegModel,

    # functions
    rq,
    fit,
    fit!,
    compute_inf,
    compute_inf!

include("qrmodel.jl")
include("fitqr.jl")
include("inference.jl")

end # module
