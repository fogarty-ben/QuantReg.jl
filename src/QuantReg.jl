module QuantReg

using DataFrames, Distributions, Format, GLM, JuMP, LinearAlgebra, StatsBase, StatsModels,
       Statistics

if haskey(ENV, "GUROBI_HOME") # ensure Gurobi is installed before loading package
    using Gurobi
end

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
