"""
    QuantRegFit(computed::Bool, method::String, coef::Union{Vector, Nothing},
                resid::Union{Vector, Nothing}, dual::Union{Array, Nothing},
                yhat::Union{Vector, Nothing})

Contains specifications for results of fitting a quantile regression model.

`computed` and `method` store specifications that will alwyas be set to their final values.
`coef`, `resid`, `dual`, and `yhat` store results of fitting a model and are set to
`nothing` until the model is fit.
"""
mutable struct QuantRegFit
    computed::Bool
    method::String
    coef::Union{Vector, Nothing}
    resid::Union{Vector, Nothing}
    dual::Union{Array, Nothing}
    yhat::Union{Vector, Nothing}
end

"""
    copy(fit::QuantRegFit)

Create a deep copy of `fit`.
"""
function Base.copy(fit::QuantRegFit)
    newfit = QuantRegFit(fit.computed, fit.method,
                         fit.coef == nothing ? nothing : copy(fit.coef),
                         fit.resid == nothing ? nothing : copy(fit.resid),
                         fit.dual == nothing ? nothing : copy(fit.dual),
                         fit.yhat == nothing ? nothing : copy(fit.yhat))

    newfit
end

"""
    QuantRegInf(computed::Bool, invers::Bool, α::Number, hs::Union{Nothing, Bool}, iid::Bool,
    interpolate::Bool, tcrit::Bool, lowerci::Union{Nothing, Array{Number}},
    upperci::Union{Nothing, Array{Number}}, σ::Union{Nothing, Array{Number}},
    teststatistic::Union{Nothing, Array{Number}}, p::Union{Nothing, Array{Number}})

Contains specifcations for and results of computing inference for a quantile regression
model.

`computed` `invers`, `α`, `hs`, `iid`, `interpolate`, and `trcit` store specifications that
will alwyas be set to their final values. `lowerci`, `upperci`, `σ`, `teststatistic`, and
`p` store results of fitting a model and are set to `nothing` until the model is fit.
"""
mutable struct QuantRegInf
    computed::Bool
    invers::Bool
    α::Number
    hs::Union{Nothing, Bool}
    iid::Bool
    interpolate::Bool
    tcrit::Bool
    lowerci::Union{Nothing, Array{Number}}
    upperci::Union{Nothing, Array{Number}}
    σ::Union{Nothing, Array{Number}}
    teststatistic::Union{Nothing, Array{Number}}
    p::Union{Nothing, Array{Number}}
end


"""
    copy(inf::QuantRegInf)

Create a deep copy of `inf`.
"""
function Base.copy(inf::QuantRegInf)
    newinf = QuantRegInf(inf.computed, inf.invers, inf.α, inf.hs, inf.iid, inf.interpolate,
                         inf.tcrit,
                         inf.lowerci == nothing ? nothing : copy(inf.lowerci),
                         inf.upperci == nothing ? nothing : copy(inf.upperci),
                         inf.σ  == nothing ? nothing : copy(inf.σ),
                         inf.teststatistic  == nothing ? nothing : copy(inf.teststatistic),
                         inf.p  == nothing ? nothing : copy(inf.p))
                         
    newinf
end

"""
    QuantRegModel(formula::FormulaTerm, data::DataFrame; τ::Number=0.5, method::String="br",
                  invers::Union{Nothing, Bool}=nothing, α::Number=0.05, hs::Bool=true,
                  iid::Union{Nothing, Bool}=nothing, interpolate::Bool=true,
                  tcrit::Bool=true)

Quantile Regression model at the `τ`th quantile fitting `data` according to `formula`.

In any call, `formula` and `data` must be set. Additional parameters and their defaults are
as specified below:
- `τ`: the quantile to fit the model at; default is 0.5, i.e. median regression
- `method`: the method to fit the model with; avaliable options are:
   - `"br"`: fit using Barrodlae-Roberts simplex (default method); see [`fitbr!`] for
   details
   - `"fn"``: Fit using Frisch-Newton algorithm; see [`fitfn!`] for details
   - `"gurobi"``: Fit using Gurobi (must have Gurobi installed); see [`fitgurobi!``]
- `invers`: if true, compute confidence intervals by inverting a rank test (otherwise use an
asymptotic esimtate of the covariance matrix); default setting is datasets with 1000 or
fewer observations and false for larger datasets
- `α`: size of test for computing inference; default setting is 0.05
- `hs`: if true, use Hall Sheather bandwidth when computing sparsity esimtates
(otherwise, use Bofinger bandwidth); default is true
- `iid`: if true, assume model errors are iid (otherwise, assume that the conditional
quantile function is locally (in tau) linear (in x)); default is true if using rank test
inversion and false if using an asymptotic estimate of the covariance matrix
- `interpolate`: if true, interpolate the confidence intervals produced by rank test
inversion inference (otherwise, print values just above and below); default is true
- `tcrit`: if true, use a Student's t distribution for calculating critical points
(otherwise use a normal distribution); default is true
"""
struct QuantRegModel <: StatisticalModel
    formula::FormulaTerm
    data::DataFrame
    mf::ModelFrame
    mm::ModelMatrix
    τ::Number
    fit::QuantRegFit
    inf::QuantRegInf

end

function QuantRegModel(formula::FormulaTerm, data::DataFrame; τ::Number=0.5,
                        method::String="br", invers::Union{Nothing, Bool}=nothing,
                        α::Number=0.05, hs::Bool=true, iid::Union{Nothing, Bool}=nothing,
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
    mfit = QuantRegFit(false, method, nothing, nothing, nothing, nothing)
    minf = QuantRegInf(false, invers, α, hs, iid, interpolate, tcrit,
            nothing, nothing, nothing, nothing, nothing)
    QuantRegModel(formula, data, mf, mm, τ, mfit, minf)
end

"""
    QuantRegModel(model::QuantRegModel; τ::Union{Nothing, Number}=nothing,
                  method::Union{Nothing, String}=nothing,
                  invers::Union{Nothing, Bool}=nothing, α::Union{Nothing, Number}=nothing,
                  hs::Union{Nothing, Bool}=nothing, iid::Union{Nothing, Bool}=nothing,
                  interpolate::Union{Nothing, Bool}=nothing,
                  tcrit::Union{Nothing, Bool}=nothing)

Construct a new QuantileRegression model by changing a parameter in an existing model.

If `τ` or `method` are unchanged, then the model fit is retained from the passed model but
the model inference type is not.
"""
function QuantRegModel(model::QuantRegModel; τ::Union{Nothing, Number}=nothing,
                       method::Union{Nothing, String}=nothing,
                       invers::Union{Nothing, Bool}=nothing,
                       α::Union{Nothing, Number}=nothing, hs::Union{Nothing, Bool}=nothing,
                       iid::Union{Nothing, Bool}=nothing,
                       interpolate::Union{Nothing, Bool}=nothing,
                       tcrit::Union{Nothing, Bool}=nothing)
    fitchanged = any(map(x -> x != nothing, [τ, method]))
    infchanged = any(map(x -> x != nothing, [invers, α, hs, iid, interpolate, tcrit]))
    if fitchanged
        mfit = QuantRegFit(false, method == nothing ? model.fit.method : method,
                           nothing, nothing, nothing, nothing)
        minf = QuantRegInf(false,
                           invers == nothing ? model.inf.invers : invers,
                           α == nothing ? model.inf.α : α,
                           hs == nothing ? model.inf.hs : hs,
                           iid == nothing ? model.inf.iid : iid,
                           interpolate == nothing ? model.inf.interpolate : interpolate,
                           tcrit == nothing ? model.inf.tcrit : tcrit,
                           nothing, nothing, nothing, nothing, nothing)
        newmodel = QuantRegModel(model.formula, model.data, model.mf, model.mm,
                                 τ == nothing ? model.τ : τ, mfit, minf)
    elseif infchanged
        minf = QuantRegInf(false,
                           invers == nothing ? model.inf.invers : invers,
                           α == nothing ? model.inf.α : α,
                           hs == nothing ? model.inf.hs : hs,
                           iid == nothing ? model.inf.iid : iid,
                           interpolate == nothing ? model.inf.interpolate : interpolate,
                           tcrit == nothing ? model.inf.tcrit : tcrit,
                           nothing, nothing, nothing, nothing, nothing)
        newmodel = QuantRegModel(model.formula, model.data, model.mf, model.mm,
                                 τ == nothing ? model.τ : τ, model.fit, minf)
    else
        error("No parameters specified to be updated.")
    end
    
    newmodel
end

"""
    copy(model::QuantRegModel)

Create a deep copy of `model` (excluding data and formula).
"""
function Base.copy(model::QuantRegModel)
    mfit = copy(model.fit)
    minf = copy(model.inf)
    mframe = ModelFrame(model.formula, model.data)
    mmatrix = ModelMatrix(mframe)
    
    QuantRegModel(model.formula, model.data, mframe, mmatrix, model.τ, mfit,
                  minf)
end

# Informs how schema apply a function to a model
implicit_intercept(::QuantRegModel) = true

# Implementing StatBase absetraction functions
StatsBase.coef(model::QuantRegModel) = model.fit.computed ? model.fit.coef :
                             error("Model hasn't been fit.")
StatsBase.coefnames(model::QuantRegModel) = coefnames(model.mf)
StatsBase.dof(model::QuantRegModel) = size(model.mm.m)[2]
StatsBase.dof_residual(model::QuantRegModel) = size(model.mm.m)[1] - size(model.mm.m)[2]
StatsBase.fitted(model::QuantRegModel) = model.fit.computed ? model.fit.yhat :
                               error("Model hasn't been fit.")
StatsBase.isfitted(model::QuantRegModel) = model.fit.computed
StatsBase.islinear(::QuantRegModel) = true
StatsBase.nobs(model::QuantRegModel) = size(model.mm.m)[1]
StatsBase.stderr(model::QuantRegModel) = model.inf.computed ? (!model.inf.invers ? model.inf.σ :
                               ["NA" for i=1:size(model.mm.m)[2]]) :
                               error("Inference hasn't been computed")                            

StatsBase.modelmatrix(model::QuantRegModel) = model.mm
StatsBase.response(model::QuantRegModel) = response(model.mf)
StatsBase.responsename(model::QuantRegModel) = terms(model.formula)[1]
StatsBase.residuals(model::QuantRegModel) = model.fit.computed ? model.fit.resid :
                                  error("Model hasn't been fit.")

"""
    coeftable(model::QuantRegModel)

Generate a coefficient table from a QuantRegModel.
"""
function StatsBase.coeftable(model::QuantRegModel)
    if !model.fit.computed
        error("Model hasn't been fit.")
    else
        # Add coefficient estimates to table
        vals = hcat(model.fit.coef)
        header = ["Coefficient"]
        if model.inf.computed & model.inf.invers # Add rank test inversion inf. to table
            vals = hcat(vals, transpose(model.inf.lowerci), transpose(model.inf.upperci))
            cistring = pyfmt("3.1d", (1 - model.inf.α) * 100)
            if model.inf.interpolate
                header = vcat(header, [cistring * "% CI Lower CI", cistring * "% CI Upper"])
            else
                header = vcat(header, [cistring * "% CI Lower", cistring* "% CI lower",
                                      cistring * "% CI upper", cistring * "% Upper",])
            end
        elseif model.inf.computed & !model.inf.invers # Add asymptotic inf. to table
            vals = hcat(vals, model.inf.σ, model.inf.teststatistic, model.inf.p)
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
    show(io::IO, model:QuantRegModel)

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
end

function QuantRegModels()
    QuantRegModels(Dict())
end

"""
    show(io::IO, models:QuantRegModels)

Display each model in `models`.
"""
function Base.show(io::IO, models::QuantRegModels)
    for (τ, model) in models.models
        show(io, model)
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
    rq(formula::FormulaTerm, data::DataFrame; kargs)

Generate, fit, and compute inference for the specified quantile regression model.

# Keyword Arguments
"""
function rq(formula::FormulaTerm, data::DataFrame; τ=0.5, kargs...)
    kwargs = Dict(kargs)
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