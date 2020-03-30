# Quickstart

## Installation

QuantReg.jl is not currently available via the Julia general repository and is pre-release.
To use it on your machine, you will need to run the following commands.

Using git, clone the QuantReg.jl repository:

```
% git clone https://github.com/fogarty-ben/QuantReg.jl
```

After cloning the repository, `cd` into `QuantReg.jl/src/FORTRAN`:

```
% cd QuantReg.jl/src/FORTRAN
```

Once inside this folder, users need to build the two FORTRAN dependences, `rqbr.f` and
`rqfnb.f`. I recommend using the gfortran compiler as in the following lines, replacing 
`${dlext}` with `dll` on Windows, `dylib` on OSX, and `so` on Linux:

```
% gfortran -o rqbr.${dlext} -fPIC -shared -std=legacy rqbr.f 
% gfortran -o rqfnb.${dlext} -fPIC -shared -std=legacy rqfnb.f
```

This should produce two new files, `rqbr.${dlext}` and `rqfnb.${dlext}` in the
`QuantReg/src/FORTRAN` directory.

Now, the package is ready do be loaded into Julia; first, `cd` into the root directory of
the project. If you are still in the `QuantReg/src/FORTRAN` directory, run the following
command:

```
% cd ../..
```

Then, open the Julia REPL in project mode using the following command:

```
% julia --project
```

Once the Julia REPL opens, install dependencies for the project as follows:

```
julia> using Pkg; 
julia> Pkg.instantiate()
```

The Gurobi.jl package may fail to build if Gurobi is not installed on machine; this is okay.
QuantReg.jl can now be loaded in the REPL with the following command:

```
using QuantReg
```

Note that pre-compilation may take some time the first time the package is loaded; this is
also to be expected. After performing these steps once, users should be able to load the
project simply by starting the Julia REPL in project mode from the root of the QuantReg.jl
project directory and loading the package. There should be no need to rebuild the FORTRAN
libraries or redownload the package dependencies.

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
[QuantRegModel](@ref) reference.

Withing the result object, `model.fit` contains the results of fitting the model, including
the coefficients, residuals, fitted values, and, for some fitting methods, the solution to
the dual problem. For more information on how to access these values, see the
[QuantRegFit](@ref) reference.

`model.inf` containts the results of computing inference for the model, encompasses
confidence interval bounds if a rank test inversion is used to compute inference and
encompasses standard errors, test statistics, p-values, and confidence interval bounds if
inference is computed using asymptotic estimates of the covariance matrix. For more
information on how to access these values, see the [QuantRegInf](@ref) reference.

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
regression command. Any specification fields in the [QuantRegModel][@ref] type can
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

In code, the `rq` functions is really a wrapper for three functions that constitute a common
work flow. It first constructs a `QuantRegModel` object according to the specifications
provided. Then, it fits the model in place according to the specifications in the generated
model object with `fit!`. This function call updates the `QuantRegFit` object
stored in the model object. Finally, it computes inference for the model according to the
specifications in the generated model object with `compute_inf!`(@ref). This
function call updates the `QuantRegInf` object stored in the model object.

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