module QuantReg

using DataFrames, StatsBase, StatsModels, SparseArrays

"""
    QuantRegModel(formula, data; τ=0.5, method="br")

Quantile Regression model at the τth percentile fitting `data` according to `formula` using
`method`.

The default method is the Barrodale-Roberts simplex (`"br"`) written by Roger Kronker. Other
available methods include:
- 'gur': Fit using Gurobi
"""
struct QuantRegModel <: StatisticalModel
    formula::FormulaTerm
    data::DataFrame
    mf::ModelFrame
    τ::Number
    method::String
    y::Vector
    X::Array
    fitted::Bool
    coefficients::Union{Vector, Nothing}
    residuals::Union{Vector, Nothing}
    dual::Union{Vector, Nothing}
    fittedvalues::Union{Vector, Nothing}
end

function QuantRegModel(formula, data; τ=0.5, method="br")
    mf = ModelFrame(formula, data)
    QuantRegModel(formula, data, mf, τ, method, response(mf), ModelMatrix(mf).m, false, 
                  nothing, nothing, nothing, nothing)
end

end # module
