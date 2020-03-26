"""
    QuantRegFit

Contains response from fitting a quantile regression model.
"""
mutable struct QuantRegFit
    computed::Bool
    coef::Union{Vector, Nothing}
    resid::Union{Vector, Nothing}
    dual::Union{Array, Nothing}
    yhat::Union{Vector, Nothing}
end

"""
    QuantRegFit

Contains specs and results of computing inference for a quantile regression model.
"""
mutable struct QuantRegInf
    computed::Bool
    exact::Bool
    α::Number
    hs::Union{Nothing, Bool}
    iid::Bool
    interpolate::Bool
    tcrit::Bool
    lowerci::Union{Nothing, Array{Number}}
    upperci::Union{Nothing, Array{Number}}
    σ::Union{Nothing, Array{Number}}
    t::Union{Nothing, Array{Number}}
    p::Union{Nothing, Array{Number}}
end

"""
    QuantRegModel(formula, data; τ=0.5, method="br")

Quantile Regression model at the τth percentile fitting `data` according to `formula` using
`method`.

The default method is the Barrodale-Roberts simplex (`"br"`) written by Roger Kronker. Other
available methods include:
- `"br"`: Fit using Barrodlae-Roberts simplex (FORTRAN code originally pary of R quantreg
package)
- `"fn"``: Fit using Frisch-Newton algorithm (FORTRAN code originally part of R quantreg
package)
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

function QuantRegModel(formula::FormulaTerm, data::DataFrame, τ::Number;
                       method::String="br", exact::Union{Nothing, Bool}=nothing,
                       α::Number=0.05, hs::Bool=true, iid::Union{Nothing, Bool}=nothing,
                       interpolate::Bool=true, tcrit::Bool=true)
    formula = apply_schema(formula, schema(formula, data),
                                       QuantRegModel)                               
    mf = ModelFrame(formula, data)
    mm = ModelMatrix(mf)
    if exact == nothing
        if size(mm.m)[1] <= 1000
            exact = true
        else
            exact = false
        end
    end
    if iid == nothing
        if exact
            iid = true
        else
            iid = false
        end
    end
    mfit = QuantRegFit(false, nothing, nothing, nothing, nothing)
    minf = QuantRegInf(false, exact, α, hs, iid, interpolate, tcrit,
                       nothing, nothing, nothing, nothing, nothing)
    QuantRegModel(formula, data, mf, mm, τ, method, mfit, minf)
end

function QuantRegModel(model::QuantRegModel; τ::Union{Nothing, Number}=nothing,
                       method::Union{Nothing, String}=nothing,
                       exact::Union{Nothing, Bool}=nothing,
                       α::Union{Nothing, Number}=nothing, hs::Union{Nothing, Bool}=nothing,
                       iid::Union{Nothing, Bool}=nothing,
                       interpolate::Union{Nothing, Bool}=nothing,
                       tcrit::Union{Nothing, Bool}=nothing)
    fitchanged = any(map(x -> x != nothing, [τ, method]))
    infchanged = any(map(x -> x != nothing, [exact, α, hs, iid, interpolate, tcrit]))
    if fitchanged
        mfit = QuantRegFit(false, nothing, nothing, nothing, nothing)
        minf = QuantRegInf(false,
                           exact == nothing ? model.inf.exact : exact,
                           α == nothing ? model.inf.α : α,
                           hs == nothing ? model.inf.hs : hs,
                           iid == nothing ? model.inf.iid : iid,
                           interpolate == nothing ? model.inf.interpolate : interpolate,
                           tcrit == nothing ? model.inf.tcrit : tcrit,
                           nothing, nothing, nothing, nothing, nothing)
        newmodel = QuantRegModel(model.formula, model.data, model.mf, model.mm,
                                 τ == nothing ? model.τ : τ,
                                 method == nothing ? model.method : method,
                                 mfit, minf)
    elseif infchanged
        minf = QuantRegInf(false,
                           exact == nothing ? model.inf.exact : exact,
                           α == nothing ? model.inf.α : α,
                           hs == nothing ? model.inf.hs : hs,
                           iid == nothing ? model.inf.iid : iid,
                           interpolate == nothing ? model.inf.interpolate : interpolate,
                           tcrit == nothing ? model.inf.tcrit : tcrit,
                           nothing, nothing, nothing, nothing, nothing)
        newmodel = QuantRegModel(model.formula, model.data, model.mf, model.mm,
                                 τ == nothing ? model.τ : τ,
                                 method == nothing ? model.method : method,
                                 model.fit, minf)
    else
        newmodel = copy(model)
    end
    
    newmodel
end

function copy(model::QuantRegModel)
    mfit = QuantRegFit(model.fit.computed, model.fit.coef, model.fit.resid, model.fit.dual,
                       model.fit.yhat)
    minf = QuantRegInf(model.inf.computed, model.inf.exact, model.inf.α, model.inf.hs,
                       model.inf.iid, model.inf.interpolate, model.inf.tcrit,
                       model.inf.lowerci, model.inf.upperci, model.inf.σ, model.inf.t,
                       model.inf.p)
    mframe = ModelFrame(model.formula, model.data)
    mmatrix = ModelMatrix(mframe)
    QuantRegModel(model.formula, model.data, mframe, mmatrix, model.τ, model.method, mfit,
                  minf)
end

implicit_intercept(::QuantRegModel) = true

StatsBase.coef(model::QuantRegModel) = model.fit.computed ? model.fit.coef :
                             error("Model hasn't been fit.")
StatsBase.coefnames(model::QuantRegModel) = coefnames(model.mf)
StatsBase.dof(model::QuantRegModel) = size(model.mm.m)[2]
StatsBase.dof_residual(model::QuantRegModel) = size(model.mm.m)[1] - size(model.mm.m)[2]
StatsBase.isfitted(model::QuantRegModel) = model.fit.computed
StatsBase.islinear(::QuantRegModel) = true
StatsBase.nobs(model::QuantRegModel) = size(model.mm.m)[1]
StatsBase.stderr(model::QuantRegModel) = model.inf.computed ? (!model.inf.exact ? model.inf.σ :
                               ["NA" for i=1:size(model.mm.m)[2]]) :
                               error("Inference hasn't been computed")                            
StatsBase.fitted(model::QuantRegModel) = model.fit.computed ? model.fit.yhat :
                               error("Model hasn't been fit.")
StatsBase.modelmatrix(model::QuantRegModel) = model.mm
StatsBase.response(model::QuantRegModel) = response(model.mf)
StatsBase.responsename(model::QuantRegModel) = terms(model.formula)[1]
StatsBase.residuals(model::QuantRegModel) = model.fit.computed ? model.fit.resid :
                                  error("Model hasn't been fit.")

"""
    coeftable(model::QuantRegModel)

Generate a coefficient table from a QuantRegModel.
"""
function coeftable(model::QuantRegModel)
    if !model.fit.computed
        error("Model hasn't been fit.")
    else
        vals = hcat(model.fit.coef)
        header = ["Estimate"]
        if model.inf.computed & model.inf.exact
            vals = hcat(vals, transpose(model.inf.lowerci), transpose(model.inf.upperci))
            cistring = pyfmt("3.1d", (1 - model.inf.α) * 100)
            if model.inf.interpolate
                header = vcat(header, [cistring * "% CI Lower CI", cistring * "% CI Upper"])
            else
                header = vcat(header, [cistring * "% CI Lower", cistring* "% CI lower",
                                      cistring * "% CI upper", cistring * "% Upper",])
            end
        elseif model.inf.computed & !model.inf.exact
            vals = hcat(vals, model.inf.σ, model.inf.t, model.inf.p)
            if model.inf.tcrit
                header = vcat(header, ["Std. Error", "t", "P(>|t|)"])
            else
                header = vcat(header, ["Std. Error", "z", "P(>|z|)"])
            end
        end 
    end
    CoefTable(vals, header, coefnames(model))
end

"""
    show(io::IO, model:QuantRegModel:: multiple::Bool=false)

Display quantreg model.
"""
function Base.show(io::IO, model::QuantRegModel; multiple=false)
    if !multiple
        print(string(model.formula) * ", " * "τ=" * string(model.τ))
    else
        print("τ=" * string(model.τ))
        if model.inf.computed & model.inf.exact
            print(", α=" * string(model.inf.α))
        end
    end
    print("\n")
    if model.fit.computed
        Base.show(io, coeftable(model))
        if !multiple
            dof_total = dof(model) + dof_residual(model)
            println("\n\nDegrees of freedom: " * string(dof_total) * " total; " * 
                    string(dof_residual(model)) * " residual")
        end
    else
        print("Unfitted.")     
    end
end

"""
QuantRegModels(models::Dict)

Wrapper for a dictionary containing multiple QuantRegModel at different quantiles.
"""
struct QuantRegModels
    models::Dict
end

"""
    show(io::IO, model:QuantRegModel:: multiple::Bool=false)

Display quantreg model.
"""
function Base.show(io::IO, models::QuantRegModels)
    headerprinted = false
    for (τ, model) in models.models
        println()
        if !headerprinted
            println(string(model.formula) * ", α=" * string(model.inf.α))
            dof_total = dof(model) + dof_residual(model)
            println("Degrees of freedom: " * string(dof_total) * " total; " * 
            string(dof_residual(model)) * " residual")
            headerprinted = true
        end
        println()
        show(io, model, multiple=true)
    end

end

function QuantRegModels()
    QuantRegModels(Dict())
end

Base.getindex(X::QuantRegModels, i) = X.models[i]
Base.append!(X::QuantRegModels, model::QuantRegModel) = setindex!(X.models, model, model.τ)
taus(X::QuantRegModels) = keys(X.models)

"""
    rq(formula::FormulaTerm, data::DataFrame; kargs)

Generate, fit, and compute inference for the specified quantile regression model.

# Keyword Arguments
"""
function rq(formula::FormulaTerm, data::DataFrame; kargs...)
    kwargs = Dict(kargs)
    τ = pop!(kwargs, :τ, 0.5)
    if typeof(τ) <: Number
        model = QuantRegModel(formula, data, τ; kwargs...)
        fit!(model)
        compute_inf!(model)
        model
    elseif typeof(τ) <: Array
        if length(τ) == 0
            error("No τ values specified.")
        end
        models = QuantRegModels()
        for tau in τ
            model = QuantRegModel(formula, data, tau; kwargs...)
            fit!(model)
            compute_inf!(model)
            append!(models, model)
        end
        models
    else
        error("Invalid τ specification.")
    end
end