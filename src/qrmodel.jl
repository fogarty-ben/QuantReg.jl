using DataFrames, StatsBase, StatsModels, SparseArrays

"""
    QuantRegResp

Contains response from fitting a quantile regression model.
"""
mutable struct QuantRegFit
    computed::Bool
    coef::Union{Vector, Nothing}
    resid::Union{Vector, Nothing}
    dual::Union{Array, Nothing}
    yhat::Union{Vector, Nothing}
end

mutable struct QuantRegInf
    computed::Bool
    rankscore::Bool
    α::Number
    hs::Union{Nothing, Bool}
    iid::Bool
    interpolate::Bool
    tcrit::Bool
    lowerci::Union{Nothing, Array{Number}}
    upperci::Union{Nothing, Array{Number}}
    σ::Union{Nothing, Number}
    t::Union{Nothing, Number}
    p::Union{Nothing, Number}
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
    fit::QuantRegFit
    inf::QuantRegInf
end

# update default constructor
function QuantRegModel(formula, data; τ::Number=0.5, method::String="br", rankscore=false,
                       α=0.05, hs=true, iid=true, interp=true, tcrit=true)
    formula = apply_schema(formula, schema(formula, data), QuantRegModel)
    mf = ModelFrame(formula, data)
    mm = ModelMatrix(mf)
    mfit = QuantRegFit(false, nothing, nothing, nothing, nothing)
    minf = QuantRegInf(false, rankscore, α, hs, iid, interp, tcrit,
                       nothing, nothing, nothing, nothing, nothing)
    QuantRegModel(formula, data, mf, mm, τ, method, mfit, minf)
end

# update handling inference
function QuantRegModel(model::QuantRegModel, τ::Number)
    mfit = QuantRegFit(false, nothing, nothing, nothing, nothing)
    QuantRegModel(model.formula, model.data, model.mf, model.mm,
                  τ, model.method, mfit, model.inf)
end

function copy(model::QuantRegModel)
    mfit = QuantRegFit(model.fit.computed, model.fit.coef, model.fit.resid, model.fit.dual,
                       model.fit.yhat)
    minf = QuantRegInf(model.inf.computed, model.inf.rankscore, model.inf.α, model.inf.hs,
                       model.inf.iid, model.inf.interpolate, model.inf.tcrit,
                       model.inf.lowerci, model.inf.upperci, model.inf.σ, model.inf.t,
                       model.inf.p)
    mframe = ModelFrame(model.formula, model.data, QuantRegModel)
    mmatrix = ModelMatrix(mframe)
    QuantRegModel(model.formula, model.data, mframe, mmatrix, model.τ, model.method, mfit,
                  minf)
end

implicit_intercept(::QuantRegModel) = true

#coef(model::QuantRegModel) = model.fit.coef
#coefnames(model)