using LinearAlgebra, Distributions, Statistics, GLM
const rqbrlib = joinpath(@__DIR__, "FORTRAN/rqbr.dylib")

"""
    fit(model)

Fit a quantile regression model.
"""
function fit(model::QuantRegModel)
    if !model.fit.computed & 0 <= model.τ <= 1
        fitfxn = nothing
        if model.method == "br"
            fittedmodel = fitbr(model)
        elseif model.method == "gurobi"
            fitfxn = fitgurobi
        end
        fittedmodel 
    elseif model.τ < 0 | model.τ > 1
        error("Error: τ must be in [0, 1]")
    else
        return model
    end
end

"""
    fitbr(model, α=0.5, ci=false, iid=true, interp=true, tcrit=true), 

Fit quantile regresion model using the Barrodale-Roberts method.

# Arguments
- `model`: QuantRegModel to fit
- `α`: test size,
- `ci`: calculate confidence intervals with rank inversion
- `iid`: whether to pase rank inversion on an i.i.d. error model
- `interp`: iterpolate discrete rank order ci
- `tcrit`: use student's t dist for crit values (as opposed to normal dist)
"""
function fitbr(model::QuantRegModel; α=0.1, ci=false, iid=true, interp=true, tcrit=true)
    big = prevfloat(Inf)
    X = model.mm.m
    y = response(model.mf)
    n, k = size(X)
    if rank(X) != k
        error("Singular design matrix: cannot use Barrodale-Roberts method")
    end
    ny = size(y)[1]
    nsol = 2
    ndsol = 2
    if k == 1
        error("Cannot compute rankscore inference for model with one predictor")
    end
    if ci
        lci1 = true
        if tcrit
            dist = TDist(n - k)
        else
            dist = Normal()
        end
        cutoff = quantile(dist, 1 - α/2)
        if iid
            qn = 1 ./ diag(inv(transpose(X) * X))
        else
            band = compute_bandwidth(model.τ, n, hs=model.inf.hs)
            if model.τ + band > 1 || model.τ - band < 0
                error("Cannot compute CI: one of bandwidth bounds outside [0, 1].")
            end
            ϵ = eps()^(1/2) # question as to whether need this
            ubmodel = QuantRegModel(model, model.τ + band)
            ub = fitbr(ubmodel, ci=false).fit.coef
            lbmodel = QuantRegModel(model, model.τ - band)
            lb = fitbr(lbmodel, ci=false).fit.coef
            δypred = model.mm.m * (ub - lb)
            if any(δypred .<= 0)
                @warn sum(δypred <= 0) * "non-positive fis. See" *
                      "http://www.econ.uiuc.edu/~roger/research/rq/FAQ #7."
            end
            f = (max.(0, (2 * band)./(δypred .- ϵ)))
            regressors = filter(x -> (x != ConstantTerm(-1)) & (x != ConstantTerm(0)),
                                terms(model.formula)[2:end])
            enum_regressors = convert(Array{Integer}, [1:1:length(regressors);])
            qn = map(i -> compute_nid_qn(i, model.data, regressors, f), enum_regressors) 
        end
    else
        lci1 = false
        qn = zeros(k)
        cutoff = 0
    end

    β, μ, d, ci = fitbrfortran(n, k, X, y, model.τ, nsol, ndsol, qn, cutoff, lci1)
    mfit = QuantRegFit(true, β, μ, d, response(model.mf) - μ)
    QuantRegModel(model, mfit)
end

function compute_nid_qn(i, data, regressors, weights)
    model = lm(regressors[i] ~ foldl(+, vcat([ConstantTerm(0)], regressors[Not(i)])), data,
               wts=weights)
    resid = residuals(model)
    sum(resid .^ 2)
end

"""
    fitbrfortran()

Call rqbr RATFOR code.

# Arguments
- `n`: number of observations
- `k`: number of parameters
- `X`: x matrix
- `y`: y vector
- `τ`: the desired quantile
- `nsol`: an estimated (row) dimension of the primal solution array
- `ndsol1`: an estimated (row) dimension of the dual solution array
- `qn`: residuals variances from the projection of each column of X on remaining columns
- `cutoff`: the critical point for testing
- `lci1`: whether to calculate CI
"""
function fitbrfortran(n, k, X, y, τ, nsol, ndsol, qn, cutoff, lci1)
    tol = eps()^(2/3) # floating point tolerance
    ift = 1
    β = zeros(k)
    μ = zeros(n)
    s = zeros(n)
    wa = zeros(n + 5, k + 4)
    wb = zeros(n)
    sol = zeros(k + 3, nsol)
    dsol = zeros(n, ndsol)
    lsol = 0
    h = zeros(k, nsol)
    ci = zeros(4, k)
    tnmat = zeros(4, k)
    big = prevfloat(Inf)
    ccall(("rqbr_", rqbrlib), Cvoid,
          (Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Int32},
           Ptr{Float64}, Ptr{Float64}, Ref{Float64}, Ref{Float64}, Ref{Int32}, Ptr{Float64},
           Ptr{Float64}, Ptr{Int32}, Ptr{Float64}, Ptr{Float64}, Ref{Int32}, Ref{Int32},
           Ptr{Float64}, Ptr{Float64}, Ref{Int32}, Ptr{Int32}, Ptr{Float64}, Ref{Float64},
           Ptr{Float64}, Ptr{Float64}, Ref{Float64}, Ref{Int32}),
          n, k, n + 5, k + 3, k + 4, X, y, τ, tol, ift, β, μ, s, wa, wb, nsol, ndsol, sol,
          dsol, lsol, h, qn, cutoff, ci, tnmat, big, lci1)
    return β, μ, dsol, ci
end

"""
    fitgurobi(model) 

Fit quantile regresion model using Gurobi.
"""
function fitgurobi
end
