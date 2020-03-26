const rqbrlib = joinpath(@__DIR__, "FORTRAN/rqbr.dylib")
const rqfnblib = joinpath(@__DIR__, "FORTRAN/rqfnb.dylib")

"""
    fit(model)

Fit a quantile regression model.
"""
function fit!(model::QuantRegModel)
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

fit(model::QuantRegModel) = fit!(copy(model))

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
function fitbr!(model::QuantRegModel; ci=false)
    big = prevfloat(Inf)
    X = model.mm.m
    y = response(model.mf)
    n, k = size(X)
    if rank(X) != k
        error("Fitting error: singular design matrix.")
    end
    ny = size(y)[1]
    nsol = 2
    ndsol = 2
    if ci
        lci1, qn, cutoff = compute_inf_exact(model)
    else
        lci1 = false
        qn = zeros(k)
        cutoff = 0
    end
    β, μ, d, confint, tnmat, flag = fitbrfortran(n, k, X, y, model.τ, nsol, ndsol, qn,
                                                 cutoff, lci1)
    
    if !ci
        if flag[1] != 0
            @warn("Solution may be non-unique.")
        end
        model.fit.computed = true
        model.fit.coef = β
        model.fit.resid = μ
        model.fit.dual = d
        model.fit.yhat = y - μ
    else
        write_ci!(model, confint, tnmat, cutoff)
    end

    model
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
    @constraint(lp, feas[i=1:n], sum(X[i, j] * β[j] for j in 1:k) - u[i] + v[i] == y[i])
    optimize!(lp)
    β = value.(β)
    μ = value.(v) - value.(u)
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
    model.fit.resid = μ
    model.fit.dual = dsol
    model.fit.yhat = response(model.mf) - μ

    model
end
