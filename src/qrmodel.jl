"""
    QuantReg.QuantRegFit(computed::Bool, coef::Union{Vector{<:Number}, Nothing},
                         resid::Union{Vector{<:Number}, Nothing}, 
                         dual::Union{Vector{<:Number}, Nothing},
                         yhat::Union{Vector{<:Number}, Nothing})

Stores results of fitting a quantile regression model.

`coef`, `resid`, `dual`, and `yhat` should be set to `nothing` until a model is fit.
""" 
mutable struct QuantRegFit
    computed::Bool
    coef::Union{Vector{<:Number}, Nothing}
    resid::Union{Vector{<:Number}, Nothing}
    dual::Union{Vector{<:Number}, Nothing}
    yhat::Union{Vector{<:Number}, Nothing}
end

"""
    deepcopy(fit::QuantRegFit)

Create a deep copy of `fit`.
"""
function Base.deepcopy(fit::QuantRegFit)
    newfit = QuantRegFit(fit.computed,
                         fit.coef == nothing ? nothing : copy(fit.coef),
                         fit.resid == nothing ? nothing : copy(fit.resid),
                         fit.dual == nothing ? nothing : copy(fit.dual),
                         fit.yhat == nothing ? nothing : copy(fit.yhat))

    newfit
end

"""
    QuantReg.QuantRegInf(computed::Bool,
                         lowerci::Union{Nothing, Vector{<:Number}, Matrix{<:Number}}
                         upperci::Union{Nothing, Vector{<:Number}, Matrix{<:Number}}
                         σ::Union{Nothing, Vector{<:Number}}
                         teststat::Union{Nothing, Vector{<:Number}}
                         p::Union{Nothing, Vector{<:Number}})

Stores results of computing inference for a quantile regression
model.

`lowerci`, `upperci`, `σ`, `teststat` should be set to `nothing` until inference is
computed.
"""
mutable struct QuantRegInf
    computed::Bool
    lowerci::Union{Nothing, Vector{<:Number}, Matrix{<:Number}}
    upperci::Union{Nothing, Vector{<:Number}, Matrix{<:Number}}
    σ::Union{Nothing, Vector{<:Number}}
    teststat::Union{Nothing, Vector{<:Number}}
    p::Union{Nothing, Vector{<:Number}}
end


"""
    deepcopy(inf::QuantRegInf)

Create a deep copy of `inf`.
"""
function Base.deepcopy(inf::QuantRegInf)
    newinf = QuantRegInf(inf.computed,
                         inf.lowerci == nothing ? nothing : copy(inf.lowerci),
                         inf.upperci == nothing ? nothing : copy(inf.upperci),
                         inf.σ  == nothing ? nothing : copy(inf.σ),
                         inf.teststat  == nothing ? nothing : copy(inf.teststat),
                         inf.p  == nothing ? nothing : copy(inf.p))
                         
    newinf
end

"""
    QuantRegModel(formula::FormulaTerm, data::DataFrame, mf::ModelFrame, mm::ModelMatrix,
                  τ::Number, fitmethod::String, invers::Bool, α::Number, hs::Bool,
                  iid::Bool, interpolate::Bool, tcrit::Bool, fit::QuantRegFit,
                  inf::QuantRegInfe)

Contains a quantile regression model at the `τ`th quantile fitting `data` according to
`formula`.

Use of this default constructor is not recommended. If rank test inversion is used to
compute inference, then the Hall-Sheather bandwidths flag will always be set to true,
overriding user choices.
"""
struct QuantRegModel <: StatisticalModel
    formula::FormulaTerm
    data::DataFrame
    mf::ModelFrame
    mm::ModelMatrix
    τ::Number
    fitmethod::String
    invers::Bool
    α::Number
    hs::Bool
    iid::Bool
    interpolate::Bool
    tcrit::Bool
    fit::QuantRegFit
    inf::QuantRegInf

    function QuantRegModel(formula::FormulaTerm, data::DataFrame, mf::ModelFrame,
                           mm::ModelMatrix, τ::Number, fitmethod::String, invers::Bool,
                           α::Number, hs::Bool, iid::Bool, interpolate::Bool, tcrit::Bool,
                           fit::QuantRegFit, inf::QuantRegInf)
        if invers
            hs = true
        end
        new(formula, data, mf, mm, τ, fitmethod, invers, α, hs, iid, interpolate, tcrit,
            fit, inf)
    end
end

"""
    QuantRegModel(formula::FormulaTerm, data::DataFrame; τ::Number=0.5,
                  fitmethod::String="br",invers::Union{Nothing, Bool}=nothing,
                  α::Number=0.05, hs::Bool=true, iid::Union{Nothing, Bool}=nothing,
                  interpolate::Bool=true, tcrit::Bool=true)

Construct a quantile regression model at the `τ`th quantile fitting `data` according to
`formula`.

In any call, `formula` and `data` must be set. The logic for selecting values for
unspecified parameters is as follows:
- `τ`: default is 0.5, i.e. median regression
- `fitmethod`: default method is the Barrodale-Roberts simplex
- `invers`: default setting is true for datasets with 1000 or fewer observations and false
    for larger datasets
- `α`: default setting is 0.05
- `hs`: default is true; always set to true when computing inference via rank test inversion
- `iid`: default is true if computing inference via rank test inversion and false otherwise
- `interpolate`: default is true
- `tcrit`: default is true
"""
function QuantRegModel(formula::FormulaTerm, data::DataFrame; τ::Number=0.5,
                        fitmethod::String="br", invers::Union{Nothing, Bool}=nothing,
                        α::Number=0.05, hs::Union{Nothing, Bool}=nothing, iid::Union{Nothing, Bool}=nothing,
                        interpolate::Bool=true, tcrit::Bool=true)
    formula = apply_schema(formula, schema(formula, data), QuantRegModel)                               
    mf = ModelFrame(formula, data)
    mm = ModelMatrix(mf)
    if invers == nothing
        if size(mm.m)[1] <= 1000
            invers = true
        else
            invers = false
        end
    end
    if iid == nothing
        if invers
            iid = true
        else
            iid = false
        end
    end
    if (invers) | (hs == nothing)
        hs = true
    end
    mfit = QuantRegFit(false, nothing, nothing, nothing, nothing)
    minf = QuantRegInf(false, nothing, nothing, nothing, nothing, nothing)
    QuantRegModel(formula, data, mf, mm, τ, fitmethod,  invers, α, hs, iid, interpolate,
                  tcrit, mfit, minf)
end

"""
    QuantRegModel(model::QuantRegModel; τ::Union{Nothing, Number}=nothing,
                  fitmethod::Union{Nothing, String}=nothing,
                  invers::Union{Nothing, Bool}=nothing, α::Union{Nothing, Number}=nothing,
                  hs::Union{Nothing, Bool}=nothing, iid::Union{Nothing, Bool}=nothing,
                  interpolate::Union{Nothing, Bool}=nothing,
                  tcrit::Union{Nothing, Bool}=nothing)

Construct a new quantile regression model by changing one or more parameters in an
existing model.

If `τ` or `method` are unchanged, then the model fit is retained from the passed model but
the model inference type is not.
"""
function QuantRegModel(model::QuantRegModel; τ::Union{Nothing, Number}=nothing,
                       fitmethod::Union{Nothing, String}=nothing,
                       invers::Union{Nothing, Bool}=nothing,
                       α::Union{Nothing, Number}=nothing, hs::Union{Nothing, Bool}=nothing,
                       iid::Union{Nothing, Bool}=nothing,
                       interpolate::Union{Nothing, Bool}=nothing,
                       tcrit::Union{Nothing, Bool}=nothing)
    fitchanged = any(map(x -> x != nothing, [τ, fitmethod]))
    infchanged = any(map(x -> x != nothing, [invers, α, hs, iid, interpolate, tcrit]))
    if fitchanged
        mfit = QuantRegFit(false, nothing, nothing, nothing, nothing)
        minf = QuantRegInf(false, nothing, nothing, nothing, nothing, nothing)
        newmodel = QuantRegModel(model.formula, model.data, model.mf, model.mm,
                                 τ == nothing ? model.τ : τ,
                                 fitmethod == nothing ? model.fitmethod : fitmethod,
                                 invers == nothing ? model.invers : invers,
                                 α == nothing ? model.α : α,
                                 hs == nothing ? model.hs : hs,
                                 iid == nothing ? model.iid : iid,
                                 interpolate == nothing ? model.interpolate : interpolate,
                                 tcrit == nothing ? model.tcrit : tcrit,
                                 mfit, minf)
    elseif infchanged
        minf = QuantRegInf(false, nothing, nothing, nothing, nothing, nothing)
        newmodel = QuantRegModel(model.formula, model.data, model.mf, model.mm,
                                 τ == nothing ? model.τ : τ, model.fitmethod,
                                 invers == nothing ? model.invers : invers,
                                 α == nothing ? model.α : α,
                                 hs == nothing ? model.hs : hs,
                                 iid == nothing ? model.iid : iid,
                                 interpolate == nothing ? model.interpolate : interpolate,
                                 tcrit == nothing ? model.tcrit : tcrit,model.fit, minf)
    else
        error("No parameters specified to be updated.")
    end
    
    newmodel
end

"""
    deepcopy(model::QuantRegModel)

Create a deep copy of `model` (excluding data and formula).
"""
function Base.deepcopy(model::QuantRegModel)
    mfit = deepcopy(model.fit)
    minf = deepcopy(model.inf)
    mframe = ModelFrame(model.formula, model.data)
    mmatrix = ModelMatrix(mframe)
    
    QuantRegModel(model.formula, model.data, mframe, mmatrix, model.τ, model.fitmethod,
                  model.invers, model.α, model.hs, model.iid, model.interpolate, model.tcrit,
                  mfit, minf)
end

# Informs how schema apply a function to a model
implicit_intercept(::QuantRegModel) = true

# Implementing StatBase abstraction functions
StatsBase.coef(model::QuantRegModel) = model.fit.computed ? model.fit.coef :
                                       error("Model hasn't been fit.")
StatsBase.coefnames(model::QuantRegModel) = coefnames(model.mf)
StatsBase.confint(model::QuantRegModel) = model.inf.computed ?
                                          hcat(model.inf.lowerci, model.inf.upperci) :
                                          error("Inference hasn't been computed.")
StatsBase.dof(model::QuantRegModel) = size(model.mm.m)[2]
StatsBase.dof_residual(model::QuantRegModel) = size(model.mm.m)[1] - size(model.mm.m)[2]
StatsBase.fitted(model::QuantRegModel) = model.fit.computed ? model.fit.yhat :
                                         error("Model hasn't been fit.")
StatsBase.isfitted(model::QuantRegModel) = model.fit.computed
StatsBase.islinear(::QuantRegModel) = true
StatsBase.nobs(model::QuantRegModel) = size(model.mm.m)[1]
StatsBase.stderr(model::QuantRegModel) = model.inf.computed ?
                                        (!model.inf.invers ? model.inf.σ :
                                         ["NA" for i=1:dof(model)]) :
                                        error("Inference hasn't been computed")                            
StatsBase.modelmatrix(model::QuantRegModel) = model.mm
StatsBase.response(model::QuantRegModel) = response(model.mf)
StatsBase.responsename(model::QuantRegModel) = terms(model.formula)[1]
StatsBase.residuals(model::QuantRegModel) = model.fit.computed ? model.fit.resid :
                                            error("Model hasn't been fit.")


function StatsBase.coeftable(model::QuantRegModel)
    if !model.fit.computed
        error("Model hasn't been fit.")
    else
        # Add coefficient estimates to table
        vals = hcat(model.fit.coef)
        header = ["Coefficient"]
        if model.inf.computed
            if !model.invers # Add asymptotic inf. to table
                vals = hcat(vals, model.inf.σ, model.inf.teststat, model.inf.p)
                if model.tcrit
                    header = vcat(header, ["Std. Error", "t", "P(>|t|)"])
                else
                    header = vcat(header, ["Std. Error", "z", "P(>|z|)"])
                end
            end 

            vals = hcat(vals, model.inf.lowerci, model.inf.upperci)
            cistring = pyfmt("3.1d", (1 - model.α) * 100)
            if model.invers & !model.interpolate
                header = vcat(header, [cistring * "% CI Lower", cistring* "% CI lower",
                                       cistring * "% CI upper", cistring * "% Upper",])
            else
                header = vcat(header, [cistring * "% CI Lower", cistring * "% CI Upper"])
            end
        end
    end
    
    CoefTable(vals, header, coefnames(model))
end

"""
    show(io::IO, model::QuantRegModel)

Display quantreg model.
"""
function Base.show(io::IO, model::QuantRegModel)
    println()
    println(string(model.formula) * ", " * "τ=" * string(model.τ))
    if model.fit.computed
        Base.show(io, coeftable(model))
            dof_total = dof(model) + dof_residual(model)
            println("\n\nDegrees of freedom: " * string(dof_total) * " total; " * 
                    string(dof_residual(model)) * " residual")
    else
        print("Unfitted.")     
    end
end

"""
QuantRegModels()

Wrapper containing multiple QuantRegModel at different quantiles indexed by quantile. 

This method partically implements some behaviors of a dictionary (for easy indexing) and
some behaviors of an array (for easy, consistent appends). This type is not intended to be
directly created by most end users.

For example, if a models::QuantRegModels contained models with `τ`=0.25, 0.5 and 0.75, these
models could be accessed as `models[0.25], models[0.5],` and `models[0.75]` respectively.
"""
struct QuantRegModels
    models::Dict
    function QuantRegModels()
        new(Dict())
    end
end

"""
    show(io::IO, models:QuantRegModels)

Display each model in `models`.
"""
function Base.show(io::IO, models::QuantRegModels)
    for τ in sort(collect(keys(models.models)))
        show(io, models[τ])
    end

end

"""
    getindex(X::QuantRegModels, τ::Number)

Returns the model in `X` fit at the τth quantile.
"""
Base.getindex(X::QuantRegModels, τ::Number) = X.models[τ]

"""
    getindex(X::QuantRegModels, τ::Number)

Check if `X` contains a model at the τth percentile.
"""
hastau(X::QuantRegModels, τ::Number) = haskey(X.models, τ)

"""
    append!(X::QuantRegModels, model::QuantRegModel)

Add `model` to `X` in-place; throws an error if `X` already contains a model with the same
τ value as `model`.
"""   
Base.append!(X::QuantRegModels, model::QuantRegModel) = !hastau(X, model.τ) ? 
                                                        setindex!(X.models, model, model.τ) :
                                                        error("Object already contains a" *
                                                              "model with τ=" *
                                                              model.τ)
                                                            

"""
    taus(X::QuantRegModels)

Retrive all the `τ` values for models stored in `X`.
"""
taus(X::QuantRegModels) = keys(X.models)

"""
    rq(formula::FormulaTerm, data::DataFrame; kwargs)

Generate, fit, and compute inference for the specified quantile regression model.

Acceptable kwargs are the same as those accepted by [QuantRegModel](@ref).
"""
function rq(formula::FormulaTerm, data::DataFrame; τ=0.5, kwargs...)
    kwargs = Dict(kwargs)
    if typeof(τ) <: Number
        model = QuantRegModel(formula, data; τ=τ, kwargs...)
        fit!(model)
        compute_inf!(model)
        model
    elseif typeof(τ) <: Array
        if length(τ) == 0
            error("No τ values specified.")
        end
        models = QuantRegModels()
        for tau in τ
            model = QuantRegModel(formula, data; τ=tau, kwargs...)
            fit!(model)
            compute_inf!(model)
            append!(models, model)
        end
        models
    else
        error("Invalid τ specification.")
    end
end