using Distributions, StatsModels, LinearAlgebra

"""
    calcbandwidth(τ, n, hs=false)

Calculate the bandwidth for nid inference.

# Arguments
- `τ`: quantile
- `n`: sample size
- `hs`: if true, use Hall Sheather bandwidth
"""
function compute_bandwidth(τ, n; hs=true)
    α = 0.05
    dist = Normal()
    q = quantile(dist, τ)
    d = pdf(dist, q)
    if hs
        qsize = quantile(dist, 1 - α/2)
        print(n^(-1/3) * qsize^(2/3) * ((1.5 * d^2)/(2 * q^2 + 1))^(1/3))
        n^(-1/3) * qsize^(2/3) * ((1.5 * d^2)/(2 * q^2 + 1))^(1/3)
    else
        n^(-1/5) * ((4.5 * d^4)/(2 * q^2 + 1)^2)^(1/5)
    end
end

"""
    compute_inf(model, se; hs=true)

Do inference for quantile regression model.

# Arguments
- `model`: quantile regression model
- `method`: means of computing inference; available options:
  - `iid`: asymptotic covariance matrix with assumption errors are iid
  - `nid`: asymptotic covariance matrix with assumption errors are nid
- `hs`: if true, use Hall Sheather bandwidth (only applicable for iid and nid methods)
"""
function compute_inf(model::QuantRegModel; se="nid", hs=true)
    if !model.fit.computed
        error("Model must be fitted before calculating confidence interval.")
    elseif se == "iid"
        σ = compute_iid_inf(model, hs)
    elseif se == "nid"
        σ = compute_nid_inf(model, hs)
    end   
end  

"""
    compute_iid_inf(model, hs)

Compute inference for a quantile regression model under iid assumption.
"""
function compute_iid_inf(model::QuantRegModel, hs)
    n, k = size(model.mm.m)
    resid = model.fit.resid
    ϵ = eps()^(1/2)
    fr = qr(model.mm.m).R
    fnorminv = fr \ I
    fnorminv = fnorminv * transpose(fnorminv)
    pz = sum(abs.(resid) .< ϵ)
    band = max(k + 1, ceil(n * compute_bandwidth(model.τ, n; hs=hs)))
    ir = convert(Array{Integer}, [(pz + 1):1:(band + pz + 1);])
    residorder = sort(resid[sortperm(abs.(resid))][ir])
    xt = ir/(n-k)
    df = DataFrame(residorder=residorder, xt=xt)
    auxmodel = QuantRegModel(@formula(residorder ~ xt), df)
    auxmodel = fit(auxmodel)
    sparsity = auxmodel.fit.coef[2]
    cov = sparsity^2 * fnorminv * model.τ * (1 - model.τ)
    scale = 1/sparsity
    sqrt.(diag(cov))
end

"""
    compute_nid_inf(model, hs)

Calculate a confidence interval for a quantile regression model under nid assumption.
"""
function compute_nid_inf(model::QuantRegModel, hs)
    n = length(response(model.mf))
    band = compute_bandwidth(model.τ, n, hs=hs)
    if model.τ + band > 1 || model.τ - band < 0
        error("Cannot compute CI: one of bandwidth bounds outside [0, 1].")
    end
    ϵ = eps()^(1/2)
    ubmodel = QuantRegModel(model, model.τ + band)
    ub = fit(ubmodel).fit.coef
    lbmodel = QuantRegModel(model, model.τ - band)
    lb = fit(lbmodel).fit.coef
    δypred = model.mm.m * (ub - lb)
    if any(δypred .<= 0)
        @warn sum(δypred <= 0) * "non-positive fis. See" *
              "http://www.econ.uiuc.edu/~roger/research/rq/FAQ #7."
    end
    f = (max.(0, (2 * band)./(δypred .- ϵ)))
    fr = qr(sqrt.(f) .* model.mm.m).R
    fnorminv = fr \ I
    fnorminv = fnorminv * transpose(fnorminv)
    cov = model.τ * (1 - model.τ) .* fnorminv * transpose(model.mm.m) * model.mm.m * fnorminv
    scale = mean(f)
    sqrt.(diag(cov))
end