using LinearAlgebra, Distributions, Statistics, Gurobi, JuMP
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
function fitbr(model::QuantRegModel; ci=false)
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
        lci1, qn, cutoff = compute_inf_rs(model)
    else
        lci1 = false
        qn = zeros(k)
        cutoff = 0
    end
    β, μ, d, ci = fitbrfortran(n, k, X, y, model.τ, nsol, ndsol, qn, cutoff, lci1)
    mfit = QuantRegFit(true, β, μ, d, response(model.mf) - μ)
    QuantRegModel(model, mfit)
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
function fitgurobi(model::QuantRegModel)
    optimizer = Gurobi.Optimizer()
    lp = direct_model(optimizer)
    
    X = model.mm.m
    y = response(model.mf)
    n, k = size(X)

    @variable(lp, β[1:k])
    @variable(lp, u[1:n] >= 0)
    @variable(lp, v[1:n] >= 0)
    exp = AffExpr(0)
    for i in 1:n
        add_to_expression!(exp, model.τ, u[i])
        add_to_expression!(exp, (1 - model.τ), v[i])
    end
    @objective(lp, Min, exp)
    @constraint(lp, feas[i=1:n], sum(X[i, j] * β[j] for j in 1:k) - u[i] + v[i] == y[i])
    optimize!(lp)
    β = value.(β)
    μ = value.(v) - value.(u)
    d = dual.(feas) .+ 0.5
    print(d)
    mfit = QuantRegFit(true, β, μ, d, y - μ)
    QuantRegModel(model, mfit)
end
