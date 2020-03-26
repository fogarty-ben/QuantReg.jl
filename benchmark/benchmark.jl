using DataFrames, Distributions, QuantReg

"""
    generatedata(n:;Integer, k::Integer; scale<:Number=10)

Generate a random, heteroskedastic n x k dataset.

# Arguments
- `n`: sample size
- `k`: number of regressors
- `scale` for all generated variables; must be positive
"""
function generatedata(n::Integer, k::Integer; scale::Number=10)
    @assert scale > 0
    
    # generate arbitrary covariance matrix
    n_entries = convert(Int64, k * (k + 1) / 2)
    unif = Uniform(-scale/2, scale/2)
    entries = rand(unif, n_entries)
    σ = zeros(k, k)
    for j in 1:k
        for i in 1:k
            if  i >= j
                σ[i, j] = pop!(entries)
            else
                σ[i, j] = σ[j, i]
            end
        end
    end
    σ = transpose(σ) * σ

    # generate design matrix
    norm = MvNormal([0 for i in 1:k], σ) 
    X = transpose(rand(norm, n))

    # generate perfect outcomes
    β = rand(unif, k)
    Y = X * β

    # add heteroskedastic errors
    βhet = rand(unif, k)
    σresid = X * βhet
    σresid = σresid .- minimum(σresid) .+ 1 # ensure all positive
    μnorm = MvNormal([0 for i in 1:n], sqrt.(σresid))
    Y = Y + rand(μnorm)

    df = hcat(DataFrame(y=Y), DataFrame(X))
    df
end

"""
    simstep(n:;Integer, k::Integer; scale<:Number=10)

Execute one step in the simulation for an n x k simulated dataset.
"""
function simstep(n::Integer, k::Integer; scale::Number=10)
    dataset = generatedata(n, k, scale)
    
end