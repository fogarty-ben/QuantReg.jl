const rqbrlib = joinpath(@__DIR__, "FORTRAN/rqbr.dylib")
const rqfnblib = joinpath(@__DIR__, "FORTRAN/rqfnb.dylib")

"""
    fit!(model::QuantRegModel)

In-place version of `[fit(model)]`.
"""
function StatsBase.fit!(model::QuantRegModel)
    if !model.fit.computed & 0 <= model.τ <= 1
        fitfxn = nothing
        if model.method == "br"
            fitbr!(model)
        elseif model.method == "gurobi"
            fitgurobi!(model)
        elseif model.method == "fn"
            fitfn!(model)
        else
            @error("Fitting method " * model.method * " unsupported.")
        end
        model.fit.computed = true
        model
    elseif model.τ < 0 || model.τ > 1
        error("Error: τ must be in [0, 1]")
    else
        @warn("Model already fitted.")
    end
    model
end

"""
    fit(model::QuantRegModel)

Fit `model` according to `model.method`.
"""
fit(model::QuantRegModel) = fit!(copy(model))

"""
    fitbr!(model::QuantRegModel; ci::Bool=false), 

Fit `model` using the Barrodale-Roberts method.

If ci is false, `model.fit` is updated in place to reflect the fit produced by running the
Barrodale-Roberts simplex, and confidence intervals are not computed. Otherwise, confidence
intervals are computed, and `model.inf` is updated in place to reflect the confidence
invervals produced by this method, but `model.fit` is not updated.

This fitting method leverages public domain FORTRAN code written by Roger Koenker for the R
`quantreg` package.
"""
function fitbr!(model::QuantRegModel; ci=false)
    big = prevfloat(Inf) 
    n, k = size(model.mm.m)
    if rank(model.mm.m) != k
        error("Fitting error: singular design matrix.")
    end
    nsol = 2 # guess for primal solution row dimension
    ndsol = 2 # guess for dual solution row dimension
    if ci
        # Initialize values needed to calculate confidence intervals
        lci1, qn, cutoff = init_ci_invers(model)
    else
        # Dummy values if not calculating confidence intervals
        lci1 = false
        qn = zeros(k)
        cutoff = 0
    end
    β, μ, d, confint, tnmat, flag = fitbrfortran(n, k, model.mm.m, response(model.mf),
                                                 model.τ, nsol, ndsol, qn, cutoff, lci1)
    if flag[1] != 0
        @warn("Solution may be non-unique. See " *
              "http://www.econ.uiuc.edu/~roger/research/rq/FAQ #1/2.")
    end
    
    if !ci # update model.fit
        model.fit.computed = true
        model.fit.coef = β
        model.fit.resid = μ
        model.fit.dual = d
        model.fit.yhat = response(model.mf) - μ
    else # update model.inf
        write_ci!(model, confint, tnmat, cutoff)
    end

    model
end

"""
Need to update from here down.

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
function fitbrfortran(n::Integer, k::Integer, X::Matrix{Number}, y::Vector{Number},
                      τ::Number, nsol::Integer, ndsol::Integer, qn::Vector{Number},
                      cutoff::Number, lci1::Bool)
    tol = eps()^(2/3) # floating point tolerance
    ift = [1]
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
           Ptr{Float64}, Ptr{Float64}, Ref{Float64}, Ref{Float64}, Ptr{Int32}, Ptr{Float64},
           Ptr{Float64}, Ptr{Int32}, Ptr{Float64}, Ptr{Float64}, Ref{Int32}, Ref{Int32},
           Ptr{Float64}, Ptr{Float64}, Ref{Int32}, Ptr{Int32}, Ptr{Float64}, Ref{Float64},
           Ptr{Float64}, Ptr{Float64}, Ref{Float64}, Ref{Int32}),
          n, k, n + 5, k + 3, k + 4, X, y, τ, tol, ift, β, μ, s, wa, wb, nsol, ndsol, sol,
          dsol, lsol, h, qn, cutoff, ci, tnmat, big, lci1)
    
          β, μ, dsol, ci, tnmat, ift
end

"""
    fitgurobi(model) 

Fit quantile regresion model using Gurobi.
"""
function fitgurobi!(model::QuantRegModel)
    optimizer = Gurobi.Optimizer(OutputFlag=0)
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
    @constraint(lp, feas[i=1:n], sum(X[i, j] * β[j] for j in 1:k) + u[i] - v[i] == y[i])
    optimize!(lp)
    β = value.(β)
    μ = value.(u) - value.(v)
    d = dual.(feas) .+ 0.5
    
    model.fit.computed = true
    model.fit.coef = β
    model.fit.resid = μ
    model.fit.dual = d
    model.fit.yhat = y - μ

    model
end

"""
    fitfn(model::QuantRegModel)

Fit model using Frish-Newton algorithm.
"""
function fitfn!(model::QuantRegModel)
    ϵ = eps()^(1/2)
    if model.τ < ϵ || model.τ > 1- ϵ
        error("Ccannot use Barrodale-Roberts method for τ extremely close to 0 or 1.")
    end
    n, k = size(model.mm.m)
    a = convert(Array, transpose(model.mm.m))
    y = -1 .* response(model.mf)
    rhs = (1 - model.τ) .* mapslices(sum, model.mm.m, dims=1)
    dsol = ones(n)
    μ = ones(n)
    β = 0.99995
    wn = zeros(10 * n)
    wn[1 : n] .= (1 - model.τ)
    wp = zeros(k, k + 3)
    nit = zeros(3)
    info = [1]
    ccall(("rqfnb_", rqfnblib), Cvoid,
          (Ref{Int32}, Ref{Int32}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
           Ptr{Float64}, Ref{Float64}, Ref{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Int32},
           Ptr{Int32}),
          n, k, a, y, rhs, dsol, μ, β, ϵ, wn, wp, nit, info)
    if info[1] != 0
        error("Fitting error: singular design matrix in stepy.")
    end
    model.fit.computed = true
    model.fit.coef = -1 .* wp[:, 1]
    model.fit.yhat = model.mm.m * model.fit.coef
    model.fit.resid = response(model.mf) - model.fit.yhat

    model
end
