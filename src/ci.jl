using Distributions, StatsModels, LinearAlgebra, GLM

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
        n^(-1/3) * qsize^(2/3) * ((1.5 * d^2)/(2 * q^2 + 1))^(1/3)
    else
        n^(-1/5) * ((4.5 * d^4)/(2 * q^2 + 1)^2)^(1/5)
    end
end

"""
    write_ci(model::QuantRegModel, ci::Array{Number, 2})

Write CI results to a model.
"""
function write_ci!(model::QuantRegModel, ci::Array{Float64, 2}, tnmat::Array{Float64, 2},
                   cutoff::Float64)
    if model.inf.interpolate
        model.inf.lowerci = ci[2,:] .- (abs.(ci[1,:] .- ci[2,:]) .* (cutoff .- abs.(tnmat[2,:]))) ./ abs.(tnmat[1,:] .- tnmat[2,:])
        model.inf.upperci = ci[3,:] .+ (abs.(ci[4,:] .- ci[3,:]) .* (cutoff .- abs.(tnmat[3,:]))) ./ abs.(tnmat[4,:] .- tnmat[3,:])
    else
        model.inf.lowerci = ci[1:2,:]
        model.inf.upperci = ci[3:4,:]
    end
end

"""
    compute_inf(model)

Compute inference for quantile regression model.
"""
function compute_inf!(model::QuantRegModel)
    if !model.fit.computed
        error("Model must be fitted before calculating confidence interval.")
    elseif model.inf.computed
        return model
    elseif model.inf.rankscore
        model = fitbr!(model; ci=true)
    else
        if model.inf.iid
            σ = compute_inf_asy_iid(model)
        else
            σ = compute_inf_asy_nid(model)
        end
        model.inf.σ = σ
        model.inf.t = model.fit.coef ./ model.inf.σ
        n, k = size(model.mm)
        dist = TDist(n - k)
        model.inf.p = 2 * (1 - cdf(dist, model.inf.t))
    end
    model.inf.computed = true
end  

compute_inf(model::QuantRegModel) = compute_inf!(copy(model))

"""
    compute_rs_nid_qn(i::Integer, data::DataFrame, regressors::Array{Term},
                      weights::Array{Float64})

Computes residuals variances from the projection of each column of X on remaining columns
for rankscore inference under the n.i.d. assumption.
"""
function compute_rs_nid_qn(i, data, regressors, weights)
    model = lm(regressors[i] ~ foldl(+, vcat([ConstantTerm(0)], regressors[Not(i)])), data,
               wts=weights)
    resid = residuals(model)
    sum(resid .^ 2)
end

"""
    compute_inf_rankscore(model, hs)

Compute inference for a quantile regression model under iid assumption.
"""
function compute_inf_rs(model::QuantRegModel)
    X = model.mm.m
    n, k = size(X)
    lci1 = true
    if model.inf.tcrit
        dist = TDist(n - k)
    else
        dist = Normal()
    end
    cutoff = quantile(dist, 1 - model.inf.α/2)
    if model.inf.iid
        qn = 1 ./ diag(inv(transpose(X) * X))
    else
        band = compute_bandwidth(model.τ, n, hs=model.inf.hs)
        if model.τ + band > 1 || model.τ - band < 0
            error("Cannot compute CI: one of bandwidth bounds outside [0, 1].")
        end
        ϵ = eps()^(1/2) # question as to whether need this
        ubmodel = QuantRegModel(model, model.τ + band)
        ub = fitbr!(ubmodel, ci=false).fit.coef
        lbmodel = QuantRegModel(model, model.τ - band)
        lb = fitbr!(lbmodel, ci=false).fit.coef
        δypred = model.mm.m * (ub - lb)
        if any(δypred .<= 0)
            @warn sum(δypred <= 0) * "non-positive fis. See" *
                    "http://www.econ.uiuc.edu/~roger/research/rq/FAQ #7."
        end
        f = (max.(0, (2 * band)./(δypred .- ϵ)))
        regressors = filter(x -> (x != ConstantTerm(-1)) & (x != ConstantTerm(0)),
                            terms(model.formula)[2:end])
        enum_regressors = convert(Array{Integer}, [1:1:length(regressors);])
        qn = map(i -> compute_rs_nid_qn(i, model.data, regressors, f), enum_regressors) 
    end
    return lci1, qn, cutoff
end

"""
    compute_inf_asy_iod(model, hs)

Compute asymptotic inference for a quantile regression model under iid assumption.
"""
function compute_inf_asy_iid(model::QuantRegModel)
    n, k = size(model.mm.m)
    resid = model.fit.resid
    ϵ = eps()^(1/2)
    fr = qr(model.mm.m).R
    fnorminv = fr \ I
    fnorminv = fnorminv * transpose(fnorminv)
    pz = sum(abs.(resid) .< ϵ)
    band = max(k + 1, ceil(n * compute_bandwidth(model.τ, n; hs=model.inf.hs)))
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
    compute_inf_asy_nid(model, hs)

Compute asymptotic inference for a quantile regression model under nid assumption.
"""
function compute_inf_asy_nid(model::QuantRegModel)
    n = length(response(model.mf))
    band = compute_bandwidth(model.τ, n, hs=model.inf.hs)
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