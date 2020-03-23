using LinearAlgebra, Distributions, Statistics, Libdl

"""
    fit(model)

Fit a quantile regression model.
"""
function fit(model::QuantRegModel)
    if !model.response.fitted & 0 <= model.τ <= 1
        fitfxn = nothing
        if model.method == "br"
            print("here")
            fittedmodel = fitbr(model)
        elseif model.method == "gurobi"
            fitfxn = fitgurobi
        end
    elseif model.τ < 0 | model.τ > 1
        error("Error: τ must be in [0, 1]")
    fittedmodel
    end
end

"""
    fitbr(model, α=0.5, ci=false, iid=true, interp=true, tcrit=true), 

Fit quantile regresion model using the Barrodale-Roberts method.

# Arguments
- `model`: QuantRegModel to fit
- `α`: test size,
- `ci`: calculate confidence intervals with rank inversion,
- `iid`: whether to pase rank inversion on an i.i.d. error model
- `interp`: iterpolate discrete rank order ci
- `tcrit`: use student's t dist for crit values (as opposed to normal dist)
"""
function fitbr(model::QuantRegModel; α=0.5, ci=false, iid=true, interp=true, tcrit=true)
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
    # This is where Kronker checks (& corrects) for out of bounds tau
    if k == 1
        ci = false
    end
    if ci
        lci = true
        if (tcrit)
            dist = TDist(n-k)
        else
            dist = Normal()
        end
        cutoff = quantile(dist, 1 - α/2)
        if (!iid) # Implement different form of testing w/ density estimation
        else
            qn = 1 / diag(transpose(X) * X)
        end
    else
        lci1 = false
        qn = zeros(k)
        cutoff = 0
    end
    β, μ, d= fitbrfortran(n, k, X, y, model.τ, nsol, ndsol, cutoff, lci1)

    resp = QuantRegResp(true, β, μ, d, response(model.mf) - μ)
    QuantRegModel(model, resp)
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
- `cutoff`: the critical point for testing
- `lci1`: whether to calculate CI
"""
function fitbrfortran(n, k, X, y, τ, nsol, ndsol, cutoff, lci1)
    tol = eps()^(2/3) # floating point tolerance
    ift = 0
    β = zeros(k)
    μ = zeros(n)
    s = zeros(n)
    wa = zeros(n + 5, k + 4)
    wb = zeros(n)
    sol = zeros(k + 3, nsol)
    dsol = zeros(n, ndsol)
    lsol = 0
    h = zeros(k, nsol)
    qn = zeros(k)
    ci = zeros(4, k)
    tnmat = zeros(4, k)
    big = prevfloat(Inf)
    ccall(("rqbr_", "rqbr"), Cvoid, (Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Int32},
                          Ptr{Float64}, Ptr{Float64}, Ref{Float64}, Ref{Float64},
                          Ref{Int32}, Ptr{Float64}, Ptr{Float64}, Ptr{Int32}, Ptr{Float64},
                          Ptr{Float64}, Ref{Int32}, Ref{Int32}, Ptr{Float64}, Ptr{Float64},
                          Ref{Int32}, Ptr{Int32}, Ptr{Float64}, Ref{Float64}, Ptr{Float64},
                          Ptr{Float64}, Ref{Float64}, Ref{Int32}),
          n, k, n + 5, k + 3, k + 4, X, y, τ, tol, ift, β, μ, s, wa, wb, nsol, ndsol, sol,
          dsol, lsol, h, qn, cutoff, ci, tnmat, big, lci1)
    return β, μ, dsol
end

# Fit using Gurobi
function fitgurobi
end