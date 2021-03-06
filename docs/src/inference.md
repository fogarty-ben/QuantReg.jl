# Computing Inference for Quantile Regression Models

QuantReg.jl currently contains two primary methods for computing inference for a quantile
regression model: (1) via a rank test inversion and (2) via an asymptotic estimate of the
covariance matrix. The choice between these methods is determined by the `invers` flag in
a `QuantRegModel` object:
- `invers`: if true, compute confidence intervals by inverting a rank test (only recommended
    for datasets with <= 1000 observations); otherwise, use an asymptotic esimtate of the
    covariance matrix

Within each of these primary methods, there are a number of sub-options:
- `α`: size of the test for computing inference; be aware that this can affect results
    beyond the confidence intervals if the Hall Sheather bandwidth is being used
- `hs`: if true, use the Hall Sheather bandwidth when computing sparsity esimtates (the
    Hall Sheather bandwidth is always used for inference via rank test inversion)
- `iid`:  if true, assume regression errors are iid; otherwise, assume that the conditional
    quantile function is locally (in tau) linear (in x)
- `interpolate`: if true, interpolate the confidence intervals produced by rank test
    inversion inference; otherwise report values just above and just below each cutoff (only
    applicable if inference is calculated via rank test inversion)
- `tcrit`: if true, use a Student's t distribution for calculating critical points; 
    otherwise use a standard normal distribution

## Generic functions

Much like the `fit!` and `fit` functions, inference for a `QuantRegModel` object can be
computed using the generic `compute_inf!` and `compute_inf` functions. The `compute_inf!`
function computes inference in place, modifying the `QuantRegInf` object stored at
`model.inf`. The `compute_inf()` method creates a deep copy of the passed model, computes
inference for the copied model, and returns the new model. For example:

```
julia> rq(@formula(Y ~ X1 + X2 + X3), df; τ=0.80, fitmethod="fn", α=0.20, invers=true,
          iid=false)
```

creates an unfitted a quantile regression model at the 0.80th quantile and specifies for it
to be fit via the Frisch-Newton. Inference is to be conducted via a rank test inversion
under the assumpton that the conditional quantile function is locally (in tau) linear (in x)
to produce 80% Confidence Intervals. Then:

```
julia> fit!(model)
```

fits the model in place via the specified method. And finally:

```
julia> compute_inf!(model)
```

computes inference in place as specified. Note that a model *must* be fit before inference
can be computed.

```@docs
QuantReg.compute_inf!
QuantReg.compute_inf
```

Ultimately, these two functions simply serve as a one-stop wrappers for computing inference
and call different subroutines depending on the specifications of the model.

## Rank Test Inversion

To compute confidence intervals by inverting a rank test, this algorithm leverages the same
function and FORTRAN routine for fitting a model with the Barrodale-Roberts simplex 
algorithm. In this case, however, the keyword argument, `ci` is set to true. One can see the
documentation for this function in the [Barrodale-Roberts Simplex Method](@ref) section 
of the [Fitting Quantile Regression Models](@ref) section.

There are three subroutines for `fitbr!` specific to computing inference via rank test
inversion. End users should not need to call these directly, but they are included here
for completeness.

```@docs
QuantReg.compute_nid_qn_invers
QuantReg.init_ci_invers
QuantReg.write_ci!
```

## Asymptotic Estimate of Covariance Matrix

The package can compute inference via an asymptotic estimate of the covariance matrix.
All code for this method is written directly in Julia using one of two functions,
`compute_σ_iid_asy`  if regression errors are assumed to be iid and `compute_σ_nid_asy`
otherwise. These functions are closely modeled on analogous code in the R quantreg pakage.

```@docs
QuantReg.compute_σ_iid_asy
QuantReg.compute_σ_nid_asy
```

## Generic Subroutines

In addition, there is a subroutine, `compute_bandwidth` that is used to calculate the τ
bandwidth for sparsity estimation in both primary methods for calculating inference.

```@docs
QuantReg.compute_bandwidth
```