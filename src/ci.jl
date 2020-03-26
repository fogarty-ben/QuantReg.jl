"""
    calcbandwidth(τ::Number, n::Integer, α::Number; hs=false)

Calculate the τ bandwidth for sparsity estimation when calculating inference under the
assumption that the conditional quantile function is locally (in tau) linear (in x).  

# Arguments
- `τ`: quantile
- `n`: sample size
- `α`: alpha level for intended confidence interval
- `hs`: if true, use Hall Sheather bandwidth; otherwise, use Bofinger bandwidth
"""
function compute_bandwidth(τ, n, α; hs=true)
    dist = Normal()
    q = quantile(dist, τ)
    d = pdf(dist, q)
    if hs
        qsize = quantile(dist, 1 - α/2)
        n^(-1/3) * qsize^(2/3) * ((1.5 * d^2)/(2 * q^2 + 1))^(1/3) # Hall Sheather
    else
        n^(-1/5) * ((4.5 * d^4)/(2 * q^2 + 1)^2)^(1/5) # Bofinter
    end
end

"""
    write_ci(model::QuantRegModel, ci::Array{Number, 2})

Write confidence intervals to `model.inf`.

`ci` should be a set of  confidence interval matrix produced by a call to
[`fitbr(model; ci=true)`].
"""
function write_ci!(model::QuantRegModel, ci::Array{Float64, 2}, tnmat::Array{Float64, 2},
                   cutoff::Float64)
    if model.inf.interpolate
        model.inf.lowerci = ci[2:2,:] .- (abs.(ci[1:1,:] .- ci[2:2,:]) .* 
                            (cutoff .- abs.(tnmat[2:2,:]))) ./ abs.(tnmat[1:1,:] .-
                             tnmat[2:2,:])
        model.inf.upperci = ci[3:3,:] .+ (abs.(ci[4:4,:] .- ci[3:3,:]) .* 
                            (cutoff .- abs.(tnmat[3:3,:]))) ./ abs.(tnmat[4:4,:] .- 
                             tnmat[3:3,:])
    else
        model.inf.lowerci = ci[1:2, :]
        model.inf.upperci = ci[3:4, :]
    end
end

"""
    compute_inf!(model)

In-place version of [`compute_inf!(model)`]
"""
function compute_inf!(model::QuantRegModel)
    if !model.fit.computed
        error("Model must be fitted before calculating confidence interval.")
    elseif model.inf.computed
        @info("Inference already computed, skipping recomputation.")
    elseif model.inf.exact
        fitbr!(model; ci=true)
    else
        if model.inf.iid
            σ = compute_σ_iid_asy(model)
        else
            σ = compute_σ_nid_asy(model)
        end
        # Write standard errors, test statistic, and p-values to to model.inf
        model.inf.σ = σ
        model.inf.teststat = model.fit.coef ./ model.inf.σ
        n, k = size(model.mm)
        dist = TDist(n - k)
        model.inf.p = 2 .* (1 .- cdf.(dist, abs.(model.inf.t)))
    end
    model.inf.computed = true
    
    model
end  

"""
    compute_inf(model)

Compute inference for `model` as specified in `model.inf`.
"""
compute_inf(model::QuantRegModel) = compute_inf!(copy(model))

"""
    compute_exact_nid_qn(i::Integer, data::DataFrame, regressors::Array{Term},
                      weights::Array{Float64})

Computes residuals variances from the projection of each column of X on remaining columns
for exact inference under the n.i.d. assumption.

This function should be of little interest to end users as it is primarily a helper function
for computing inference with a rank test inversion and with n.i.d. errors.

# Arguments
- `i`: enumeration of column to project
- `data`: data used to fit the model
- `regressors`: list of regressors from the model
- `weights`: weight to use in calculating projection as a linear regression
"""
function compute_exact_nid_qn(i::Integer, data::DataFrame, regressors::Array{Term},
                              weights::Array{Float64, 1})
    formula = regressors[i] ~ foldl(+, vcat([ConstantTerm(0)], regressors[Not(i)]))
    model = lm(formula, data, wts=weights)
    resid = residuals(model)
    rv = sum(resid .^ 2)

    rv
end

"""
    init_inf_invers(model::QuantRegModel)

Initialize necessary values for calculating confidence intervals for model with a rank test
inversion. 

This function should be of little interest to end users as it is solely used as a subroutine
when calling [`fitbr(model; ci=true)`].
"""
function init_ci_invers(model::QuantRegModel)
    n, k = size(model.mm.m)
    if k == 1
        error("CI calculation: cannot compute inference with a rank test inversion for a " *
              "model with one predictor.")
    end
    
    # Set confidience interval flag
    lci1 = true

    # Get critical point for the confidence interval
    if model.inf.tcrit
        dist = TDist(n - k)
    else
        dist = Normal()
    end
    cutoff = quantile(dist, 1 - model.inf.α/2)
    
    # Calculate resid. var. from projecting each predictor column on others
    if model.inf.iid # assuming errors are iid
        qn = 1 ./ diag(inv(transpose(model.mm.m) * model.mm.m))
    else # assuming errors are nid
        band = compute_bandwidth(model.τ, n, model.inf.α; hs=model.inf.hs)
        if model.τ + band > 1 || model.τ - band < 0
            error("Cannot compute CI: one of bandwidth bounds outside [0, 1].")
        end
        ubmodel = QuantRegModel(model, model.τ + band)
        ub = fitbr!(ubmodel, ci=false).fit.coef
        lbmodel = QuantRegModel(model, model.τ - band)
        lb = fitbr!(lbmodel, ci=false).fit.coef
        dypred = model.mm.m * (ub - lb)
        if any(dypred .<= 0)
            @warn sum(dypred <= 0) * "non-positive fis. See" *
                    "http://www.econ.uiuc.edu/~roger/research/rq/FAQ #7."
        end
        
        ϵ = eps()^(1/2)
        f = (max.(0, (2 * band)./(dypred .- ϵ)))
        regressors = filter(x -> (x != ConstantTerm(-1)) & (x != ConstantTerm(0)),
                            terms(model.formula)[2:end])
        enum_regressors = convert(Array{Integer}, [1:1:length(regressors);])
        qn = map(i -> compute_exact_nid_qn(i, model.data, regressors, f), enum_regressors) 
    end

    lci1, qn, cutoff
end

"""
    compute_σ_iid_asy(model::QuantRegModel)

Compute standard errors for a `model` using an estimate of the asymptotic covariance matrix
under the assumtion that errors are iid.
"""
function compute_σ_iid_asy(model::QuantRegModel)
    # Prep for auxillary model
    n, k = size(model.mm.m)
    ϵ = eps()^(1/2)
    fr = qr(model.mm.m).R
    fnorminv = fr \ I
    fnorminv = fnorminv * transpose(fnorminv)
    pz = sum(abs.(model.fit.resid) .< ϵ)
    band = max(k + 1, ceil(n * compute_bandwidth(model.τ, n, model.inf.α; hs=model.inf.hs)))
    ir = convert(Array{Integer}, [(pz + 1):1:(band + pz + 1);])
    residorder = sort(model.fit.resid[sortperm(abs.(model.fit.resid))][ir])
    xt = ir/(n-k)
    
    # Generate and fit auxillary model
    df = DataFrame(residorder=residorder, xt=xt)
    auxmodel = QuantRegModel(@formula(residorder ~ xt), df, 0.5)
    fit!(auxmodel)
    
    # Use auxillary model to calculate standard errors
    sparsity = auxmodel.fit.coef[2]
    cov = sparsity^2 * fnorminv * model.τ * (1 - model.τ)
    scale = 1/sparsity
    σ = sqrt.(diag(cov))

    σ
end

"""
    compute_σ_nid_asy(model::QuantRegModel)

Compute standard errors for a `model` using an estimate of the asymptotic covariance matrix
under the assumption that the conditional quantile function is locally (in tau) linear
(in x).  
"""
function compute_σ_nid_asy(model::QuantRegModel)
    # Generate and fit auxillary models
    band = compute_bandwidth(model.τ, length(response(model.mf)), model.inf.α;
                             hs=model.inf.hs)
    if model.τ + band > 1 || model.τ - band < 0
        error("Cannot compute CI: one of bandwidth bounds outside [0, 1].")
    end
    ϵ = eps()^(1/2)
    ubmodel = QuantRegModel(model, model.τ + band)
    ub = fit(ubmodel).fit.coef
    lbmodel = QuantRegModel(model, model.τ - band)
    lb = fit(lbmodel).fit.coef
    
    # Use auxillary models to calculate standard errors
    δypred = model.mm.m * (ub - lb)
    if any(δypred .<= 0)
        @warn sum(δypred <= 0) * "non-positive fis. See" *
              "http://www.econ.uiuc.edu/~roger/research/rq/FAQ #7."
    end
    f = (max.(0, (2 * band)./(δypred .- ϵ)))
    fr = qr(sqrt.(f) .* model.mm.m).R
    fnorminv = fr \ I
    fnorminv = fnorminv * transpose(fnorminv)
    cov = model.τ * (1 - model.τ) .* fnorminv * transpose(model.mm.m) * model.mm.m *
          fnorminv
    scale = mean(f)
    σ = sqrt.(diag(cov))

    σ
end