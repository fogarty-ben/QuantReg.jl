using DataFrames, StatsBase, StatsModels, SparseArrays

"""
    QuantRegModel(formula, data; τ=0.5, method="br")

Quantile Regression model at the τth percentile fitting `data` according to `formula` using
`method`.

The default method is the Barrodale-Roberts simplex (`"br"`) written by Roger Kronker. Other
available methods include:
- `"br"``: Fit using Barrodlae-Roberts simplex (FORTRAN code written by Roger Kronker)
- `"gurobi"``: Fit using Gurobi
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

# should this be convenience constructor or function rq
function QuantRegModel(formula, data; τ::Number=0.5, method::String="br")
    mf = ModelFrame(formula, data)
    QuantRegModel(formula, data, mf, τ, method, response(mf), ModelMatrix(mf).m, false, 
                  nothing, nothing, nothing, nothing)
end