# QuantReg.jl Types

The QuantReg.jl package introduces four types: [QuantRegFit](@ref), [QuantRegInf](@ref),
[QuantRegModel](@ref), and [QuantRegModels](@ref). For most users, a deep understanding of
these types is unnecessary, but the structure of these types will be helpful for advanced
users, students investigating the properties of quantile regression, or developers wishing
to extend the package.

These four types work together as follows. A `QuantRegModel` is an immutable object and
contains all of the specification necessary to fit and compute inference for a model plus
a `QuantRegFit` object and a `QuantRegInf` object that respectively store the results of
fitting the model and computing inference. A `QuantRegModel` object stores multiple 
`QuantRegModel` objects, indexed by their τ values; note that a `QuantRegModels` object
cannot contain multple `QuantRegModel` objects with the same τ value.

## QuantRegFit

`QuantRegFit` objects store the results of fitting a quantile regression model and are
mutable. Each instance contains the following fields:

- `computed`: denotes whether or not the parent model has been fit
- `coef`: stores the coefficients for a fitted model or nothing if
    model hasn't been fit
- `resid`: stores the residuals for each observation or nothing if
    the model hasn't been fit
- `dual`: stores the solution to the dual problem if the fitting
    method produces the dual solution or nothing if the model hasn't been fit or the fitting
    method does not produce the dual solution
- `yhat`: stores the fitted values for the dataset used to build the
    model or nothing if the model hasn't been fit

`QuantRegFit` has one construtor, the default constructor in which all field values are
specified directly. `QuantRegFit` objects are primarily constructed in the creation of
`QuantRegModel` objects.

```@docs
QuantReg.QuantRegFit
```

Additionally, there is a method to create deep copies of `QuantRegFit` objects.

```@docs
deepcopy(::QuantReg.QuantRegFit)
```

## QuantRegInf

Similarly, `QuantRegInf` objects store the results of computing inference for a quantile
regression model and are mutable. Each instance contains the following fields:

- `computed`: denotes whether or not inference has been computed for the parent model
- `lowerci`: lower bounds for the specified confidence intervals or nothing if inference
    hasn't been computed for the model; a (k x 1) matrix unless inference is computed via 
    rank test inversion and interpolate is set to false in which case, a (k x 2) matrix
- `upperci`: upper bounds for the specified confidence intervals or nothing if the inference
    hasn't been computed for the model; a (k x 1) matrix unless inference is computed via
    rank test inversion and interpolate is set to false in which case, a (k x 2) matrix
- `σ`: standard errors for the coefficients or nothing if the inference hasn't been computed
    for the model or is calculated via a rank test inversion
- `teststat`: t or z statistics for the coefficients or nothing if inference hasn't been
    computed for the model or is calculated via a rank test inversion
- `p`: p values for the coefficients or nothing if inference hasn't been computed for the
    model or is calculated via a rank test inversion

`QuantRegInf` also has one construtor, the default constructor in which all field values are
specified directly. `QuantRegFit` objects are primarily constructed in the creation of
`QuantRegModel` objects.

```@docs
QuantReg.QuantRegInf
```

Additionally, there is a method to create deep copies of `QuantRegInf` objects.

```@docs
deepcopy(::QuantReg.QuantRegInf)
```

# QuantRegModel

`QuantRegModel` is the type with which end users are most likely to interact.
`QuantRegModel` objects are immutable and contain the following fields:

- `formula`: the formula to fit the data
- `data`: DataFrame object containing data to fit the model with
- `mf`: ModelFrame object created by apply `formula` to `data`
- `mm`: ModelMatrix object created from `mf`
- `τ`: quantile to fit the model at
- `fitmethod`: the algorithm/method to fit the model; available options are:
  - `"br"`: fit using Barrodale-Roberts simplex (default method); see [`fitbr!`] for 
   details
   - `"fn"`: Fit using Frisch-Newton algorithm; see [`fitfn!`] for details
   - `"gurobi"`: Fit using Gurobi (must have Gurobi installed); see [`fitgurobi!``]
- `invers`: if true, compute confidence intervals by inverting a rank test (only recommended
    for datasets with <= 1000 observations); otherwise, use an asymptotic esimtate of the
    covariance matrix
- `α`: size of the test for computing inference
- `hs`: if true, use Hall Sheather bandwidth when computing sparsity esimtates (not
    applicable if inference is calculated via rank test inversion with iid regression errors)
- `iid`:  if true, assume regression errors are iid; otherwise, assume that the conditional
    quantile function is locally (in tau) linear (in x)
- `interpolate`: if true, interpolate the confidence intervals produced by rank test
    inversion inference; otherwise report values just above and just below each cutoff (only
    applicable if inference is calculated via rank test inversion)
- `tcrit`: if true, use a Student's t distribution for calculating critical points; 
    otherwise use a normal distribution
- `fit`: `QuantRegFit` object storing the results of fitting the model
- `inf`: `QuantRegInf` object storing the results of computing inference for the model

### Constructors

There are multiple ways to construct a `QuantRegModel` object. The preferred for end users
constructor accepts a formula and dataset as positional arguments and model specifications
as keyword arguments. Depending on the passed arguments, the constructor intelligently
choses defaults for unspecified fields as described below and constructs empty `QuantRegFit`
and `QuantRegInf` objects.

```@docs
QuantReg.QuantRegModel(::StatsModels.FormulaTerm, ::DataFrames.DataFrame; ::Number,
                       ::String, ::Union{Nothing, Bool}, ::Number, ::Bool,
                       ::Union{Nothing, Bool}, ::Bool, ::Bool)
```

Additionally, one can construct a new `QuantRegModel` object from an existing
`QuantRegModel` object by specifying parameters to change.

```@docs
QuantReg.QuantRegModel(::QuantRegModel; ::Union{Nothing, Number}, ::Union{Nothing, String},
              ::Union{Nothing, Bool}, ::Union{Nothing, Number}, hs::Union{Nothing, Bool},
              ::Union{Nothing, Bool}, ::Union{Nothing, Bool}, ::Union{Nothing, Bool})
```

One particularly convient use caase for this constructor would be changing the test size or
inference type for an existing model. For example, consider generating the following model:

```
julia> model = rq(@formula(Y ~ X1 + X2 + X3), df; τ=0.80, fitmethod="fn", α=0.20, iid=false)
```

Suppose after fitting this model, the use wanted to calculate 90% confidence intervals
without recomputing the model fit. Such a change is nontrivial relative to standard linear
regression because the Hall-Sheather bandwidth depends on the the confidence intervlals (as
does the rank test inversion, if calcuating in that manner). To compute this change, the
user could make employ the above constructor as follows:

```
julia> newmodel = QuantRegModel(model; α=0.10) # Generate model with new test size
julia> compute_inf!(newmodel) # Compute inference for the modified model in-place
```

Finally, the `QuantRegModel` object also has a default constructor. Use of the default
constructor is *strongly* discouraged.

```@docs
QuantReg.QuantRegModel
```

### Methods

There is a method to create a deep copy of a `QuantRegModel` object.

```@docs
deepcopy(::QuantReg.QuantRegModel)
```

There are also methods for the following StatsBase abstractions:

```@docs
StatsBase.coef
StatsBase.coefnames
StatsBase.coeftable
StatsBase.confint
StatsBase.dof
StatsBase.dof_residual
StatsBase.fitted
StatsBase.isfitted
StatsBase.islinear
StatsBase.nobs
StatsBase.stderr
StatsBase.modelmatrix
StatsBase.response
StatsBase.responsename
StatsBase.residuals
```

As well as a method to nicely print the results of a model:

```@docs
Base.show(::IO, ::QuantRegModel)
```

## QuantRegModels

QuantReg.jl also contains a `QuantRegModels` type that enables users to store multiple
models at different τ values. This type is a hybrid between a dictionary and a list, so that
models are always indexed by the τ value for which they are fit. A single `QuantRegModels`
object cannot contain multiple models with the same τ value. 

!!! note
    The primary intent of these objects is to store quantile regression models with the
    exact same specifications aside from τ, though this is not strictly enforced. The
    internal implementation may change in the future, however, so other uses are *strongly*
    discouraged.

Currently, a `QuantRegModel` has a single field, models, which stores all the
`QuantRegModel` objects in a dictionary indexed by τ value. This dictionary should not be
modified directly.

### Constructor 

This type has one constructor, which always creates an empty `QuantRegModels` object.

```@docs
QuantReg.QuantRegModels
```

### Methods

As a hybrid of an array and a dictionary, some methods intended for arrays and some intended
for dictionaries are implemented for `QuantRegModels`.

For example, to add a new model to a `QuantRegModels`, one should use `append!`. This method
ensures that the underlying dictionary exclusively contains `QuantRegModel` objects,
that the models are indexed by their τ value, and that each `QuantRegModels` object contains
no more than one `QuantRegModel` for a given τ value. 

```@docs
Base.append!(::QuantReg.QuantRegModels, ::QuantReg.QuantRegModel)
```

A `getindex` method is implemented for this type, so users can access a particular model via
the following syntax

```
julia> models[τ]
```

where models is a `QuantRegModels` object and τ  is a τ value for a `QuantRegModel` stored
in models.

```@docs
Base.getindex(::QuantReg.QuantRegModels, ::Number)
```

There are also two methods, `taus` and `hastau`, which respectively display all the τ values
for which a given `QuantRegModels` object contains a model and checks whther a given
`QuantRegModels` object contains a model for a specific value of τ.

```@docs
QuantReg.taus
QuantReg.hastau
```

Lastly, there is also an implementation of `show` which nicely displays all the models in a
`QuantRegModels` object.

```@docs
Base.show(::IO, ::QuantReg.QuantRegModels)
```
