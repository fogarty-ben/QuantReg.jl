# Quickstart Guide

## Installation

QuantReg.jl is available via the Julia general repository. To use it on your machine, you
will need to run the following commands:

```
julia> using Pkg
julia> Pkg.add(QuantReg)
```

The Gurobi.jl dependency package may fail to build if Gurobi is not installed on the
machine; this is okay. QuantReg.jl can now be loaded in the REPL with the following command:

```
julia> using QuantReg
```

Note that pre-compilation may take some time the first time the package is loaded; this is
also to be expected.

## Running your first model

The easiest way to generate, fit, and compute inference for a quantile regression model is
via the `rq` command.

```
julia> using CSV, QuantReg
julia> df = CSV.read("/path/to/data.csv") # load dataset to fit
julia> model = rq(@formula(Y ~ X1 + X2 + X3), df; τ=0.50) # fit dataset
```

The resulting object `model` contains a median regression fitted and with inference computed
according to default settings. For more information on default settings, see the
[QuantRegModel](@ref) setting.

Within the result object, `model.fit` contains the results of fitting the model, including
the coefficients, residuals, fitted values, and, for some fitting methods, the solution to
the dual problem. For more information on how to access these values, see the
[QuantRegFit](@ref) section.

`model.inf` contains the results of computing inference for the model, which encompasses
confidence interval bounds if a rank test inversion is used to compute inference and
encompasses standard errors, test statistics, p-values, and confidence interval bounds if
inference is computed using asymptotic estimates of the covariance matrix. For more
information on how to access these values, see the [QuantRegInf](@ref) section.

## Running models at multiple quantiles

It is also possible to fit more than one model simultaneously by passing an array of τ 
values. 

```
julia> models = rq(@formula(Y ~ X1 + X2 + X3), df; τ=[0.25, 0.50, 0.75])
```

This call yields a QuantRegModels object containing three quantile regression models, fit at
the 0.25th, 0.50th, and 0.75th quantiles. The individual models can be accessed by using the
τ values as indices. For example,

```
julia> models[0.25]
```

will return the quantile regression model fit at the 0.25th quantile. For more information
on QuantRegModels objects, see the [QuantRegModels](@ref) section.

## Configuring models

Users can provide more detailed model specifications as keyword arguments to the quantile
regression command. Any specification field in the [QuantRegModel](@ref) type can
be accepted as a keyword argument. The following command fits a quantile regression model at
the 0.80th quantile using the Frisch-Newton and computes 80% confidence intervals under the
assumption that the conditional quantile function is locally (in tau) linear (in x) (i.e.
regression errors are not i.i.d.).

```
julia> rq(@formula(Y ~ X1 + X2 + X3), df; τ=0.80, fitmethod="fn", α=0.20, iid=false)
```

```@docs
rq
```

## Under the hood of rq

In code, the `rq` function is really a wrapper for three functions that constitute a common
workflow. It first constructs a new `QuantRegModel` object according to the specifications
provided. Then, it fits the model in place according to the specifications in the model with
`fit!`. This function call updates the `QuantRegFit` object stored in the model. Finally, it
computes inference for the model according to the specifications in the model with
`compute_inf!`. This function call updates the `QuantRegInf` object stored in the model.

For example, the call:

```
julia> model = rq(@formula(Y ~ X1 + X2 + X3), df; τ=0.80, fitmethod="fn", α=0.20, iid=false)
```

is equivalent to:

```
julia> model = QuantRegModel(@formula(Y ~ X1 + X2 + X3), df; τ=0.80, fitmethod="fn", α=0.20, 
                             iid=false)
julia> fit!(model)
julia> compute_inf!(model)
```