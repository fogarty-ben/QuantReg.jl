# Quickstart

## Installation

QuantReg.jl is not currently available via the Julia general repository and must be
installed from GitHub.

From the Pkg REPL, run the following command:

```
pkg> add https://github.com/fogarty-ben/QuantReg.jl
```
Alternatively, run the following commands from the Julia REPL:

```
julia> using Pkg;  
julia> Pkg.add(PackageSpec(url="https://github.com/fogarty-ben/QuantReg.jl", rev="master"))
```

## Running your first model

The easiest way to generate, fit, and compute inference for a quantile regression model is
via the `rq` command.

```
julia> using CSV, QuantReg
julia> df = CSV.read("/path/to/data.csv") # load dataset to fit
julia> model = rq(@formula(Y ~ X1 + X2 + X3), df; τ=0.50)
```

The resulting object `model` contains a median regression fitted and with inference computed
according to default settings. For more information on default settings, see the
[type(QuantRegModel)](@ref) reference.

Withing the result object, `model.fit` contains the results of fitting the model, including
the coefficients, residuals, fitted values, and, for some fitting methods, the solution to
the dual problem. For more information on how to access these values, see the
[type(QuantRegFit)](@ref) reference.

`model.inf` containts the results of computing inference for the model, encompasses
confidence interval bounds if a rank test inversion is used to compute inference and
encompasses standard errors, test statistics, p-values, and confidence interval bounds if
inference is computed using asymptotic estimates of the covariance matrix. For more
information on how to access these values, see the [type(QuantRegInf)](@ref) reference.

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

will return the quantile regression model fit at the 0.25th quantile.

## Configuring models

Users can provide more detailed model specifications as keyword arguments to the quantile
regression command. Any specification fields in the [QuantRegModel][@ref] type can be
accepted as a keyword argument. The following command fits a quantile regression model at
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

In code, the `rq` functions is really a wrapper for three functions that constitute a common
work flow. It first constructs a `QuantRegModel` object according to the specifications
provided. Then, it fits the model in place according to the specifications in the generated
model object with [func(fit!)](@ref). This function call updates the `QuantRegFit` object
stored in the model object. Finally, it computes inference for the model according to the
specifications in the generated model object with [func(compute_inf!)](@ref). This function
call updates the `QuantRegInf` object stored in the model object.