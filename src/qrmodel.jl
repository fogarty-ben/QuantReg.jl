using DataFrames, StatsBase, StatsModels, SparseArrays

"""
    QuantRegResp

Contains response from fitting a quantile regression model.
"""
mutable struct QuantRegResp
    fitted::Bool
    coefficients::Union{Vector, Nothing}
    residuals::Union{Vector, Nothing}
    dual::Union{Array, Nothing}
    fittedvalues::Union{Vector, Nothing}
end

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
    mm::ModelMatrix
    τ::Number
    method::String
    response::QuantRegResp
end

function QuantRegModel(formula, data; τ::Number=0.5, method::String="br")
    mf = ModelFrame(formula, data)
    mm = ModelMatrix(mf)
    mr = QuantRegResp(false, nothing, nothing, nothing, nothing)
    QuantRegModel(formula, data, mf,mm, τ, method, mr)
end

function QuantRegModel(model::QuantRegModel, resp::QuantRegResp)
    QuantRegModel(model.formula, model.data, model.mf, model.mm,
                  model.τ, model.method, resp)
end