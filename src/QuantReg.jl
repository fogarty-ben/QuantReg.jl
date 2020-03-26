module QuantReg

using DataFrames, Distributions, Format, GLM, Gurobi, JuMP, LinearAlgebra, StatsBase, StatsModels,
       Statistics

export coef, coefnames, coeftable, dof, dof_residual, fitted, @formula, isfitted, islinear,
       nobs, modelmatrix, response, responsename, residuals

export
    # Model type
    QuantRegModel,

    # functions
    rq,
    compute_inf,
    compute_inf!


include("qrmodel.jl")
include("fitqr.jl")
include("ci.jl")

end # module
