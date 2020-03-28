# Fitting Quantile Regression Models

QuantReg.jl currently contains three methods for fitting quantile regression models: (1) the
Barrodale-Roberts simple algorithm, (2) the Frisch-Newton algorithm, and (3) via the
industrial solver Gurobi. The first two methods leverage FORTRAN code originally
written for the R quantreg package. The third method was developed specifically for this
package and only works if Gurobi is already installed on the user's machine.

## Generic functions

A `QuantRegModel` object, named `model`, can be fit using the generic `fit!()` function, and
indeed this is the function that most end users will leverage. This functions fits the
object in place, modifying the `QuantRegFit` object stored at `model.fit`. For example,

```
julia> model = QuantRegModel(@formula(Y ~ X1 + X2 + X3), df; τ=0.80, fitmethod="fn")
```

creates an unfitted a quantile regression model at the 0.80th quantile and specifies for it
to be fit via the Frisch-Newton. Then,

```
julia> fit!(model)
```

fits the model in-place via the specified method.

Additionaly, there is a `fit()` method that creates a deep copy of the passed model, fits
that copy, and returns the new, fitted model.

```@docs
QuantReg.fit!
QuantReg.fit
```

Ultimately, these two functions simply serve as a one-stop wrapper for fitting a model and
call different subroutines depending on the value of `model.fitmethod`.

## Choosing a fitting method

* TODO: Write about simulation results *

## Barrodale-Roberts Simplex Method

The Barrodale-Roberts simplex algorithm implemented can be used for both fitting the model
and for computing inference via a rank test inversion. When fitting, the following method is
called with the keyword argument, `ci`, set to its default value of false.

To select this method when specifying a quantile regression model, set the keyword argument
`fitmethod="br"` in either the `rq` function of the `QuantRegModel` constructor.

```@docs
QuantReg.fitbr!
```

This method relies on a subroutine, `fitbrfortran`, that handles the FORTRAN call.

```@docs
QuantReg.fitbr!
```

## Frisch-Newton Algorithm

The Frisch-Newton algorithm implemented is an interior point method that can be used to
solve quantile regression problems. Unlike the other two methods in this package, this
method does not produce solutions to the dual problem.

To select this method when specifying a quantile regression model, set the keyword argument
`fitmethod="fn"` in either the `rq` function of the `QuantRegModel` constructor.

```@docs
QuantReg.fitfn!
```

## Gurobi

The only fitting method available via this package that is not available via the R quantreg
package is fitting via Gurobi. This fitting method is only available to users who already
have a licensed version of Gurobi installed on their machine and the Gurobi.jl package
added in Julia. Note that Gurobi.jl is not installed as a dependency of QuantReg.jl so that
the package is accessible to users without a licensed version of Gurobi (the package
Gurobi.jl cannot be properly installed without a licensed version of Gurobi).

To select this method when specifying a quantile regression model, set the keyword argument
`fitmethod="fn"` in either the `gurobi` function of the `QuantRegModel` constructor.

```@docs
QuantReg.fitgurobi!
```