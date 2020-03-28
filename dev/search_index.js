var documenterSearchIndex = {"docs":
[{"location":"quickstart/#Quickstart-1","page":"Quickstart Guide","title":"Quickstart","text":"","category":"section"},{"location":"quickstart/#Installation-1","page":"Quickstart Guide","title":"Installation","text":"","category":"section"},{"location":"quickstart/#","page":"Quickstart Guide","title":"Quickstart Guide","text":"QuantReg.jl is not currently available via the Julia general repository and must be installed from GitHub.","category":"page"},{"location":"quickstart/#","page":"Quickstart Guide","title":"Quickstart Guide","text":"From the Pkg REPL, run the following command:","category":"page"},{"location":"quickstart/#","page":"Quickstart Guide","title":"Quickstart Guide","text":"pkg> add https://github.com/fogarty-ben/QuantReg.jl","category":"page"},{"location":"quickstart/#","page":"Quickstart Guide","title":"Quickstart Guide","text":"Alternatively, run the following commands from the Julia REPL:","category":"page"},{"location":"quickstart/#","page":"Quickstart Guide","title":"Quickstart Guide","text":"julia> using Pkg;  \njulia> Pkg.add(PackageSpec(url=\"https://github.com/fogarty-ben/QuantReg.jl\", rev=\"master\"))","category":"page"},{"location":"quickstart/#Running-your-first-model-1","page":"Quickstart Guide","title":"Running your first model","text":"","category":"section"},{"location":"quickstart/#","page":"Quickstart Guide","title":"Quickstart Guide","text":"The easiest way to generate, fit, and compute inference for a quantile regression model is via the rq command.","category":"page"},{"location":"quickstart/#","page":"Quickstart Guide","title":"Quickstart Guide","text":"julia> using CSV, QuantReg\njulia> df = CSV.read(\"/path/to/data.csv\") # load dataset to fit\njulia> model = rq(@formula(Y ~ X1 + X2 + X3), df; τ=0.50)","category":"page"},{"location":"quickstart/#","page":"Quickstart Guide","title":"Quickstart Guide","text":"The resulting object model contains a median regression fitted and with inference computed according to default settings. For more information on default settings, see the type(QuantRegModel) reference.","category":"page"},{"location":"quickstart/#","page":"Quickstart Guide","title":"Quickstart Guide","text":"Withing the result object, model.fit contains the results of fitting the model, including the coefficients, residuals, fitted values, and, for some fitting methods, the solution to the dual problem. For more information on how to access these values, see the type(QuantRegFit) reference.","category":"page"},{"location":"quickstart/#","page":"Quickstart Guide","title":"Quickstart Guide","text":"model.inf containts the results of computing inference for the model, encompasses confidence interval bounds if a rank test inversion is used to compute inference and encompasses standard errors, test statistics, p-values, and confidence interval bounds if inference is computed using asymptotic estimates of the covariance matrix. For more information on how to access these values, see the type(QuantRegInf) reference.","category":"page"},{"location":"quickstart/#Running-models-at-multiple-quantiles-1","page":"Quickstart Guide","title":"Running models at multiple quantiles","text":"","category":"section"},{"location":"quickstart/#","page":"Quickstart Guide","title":"Quickstart Guide","text":"It is also possible to fit more than one model simultaneously by passing an array of τ  values. ","category":"page"},{"location":"quickstart/#","page":"Quickstart Guide","title":"Quickstart Guide","text":"julia> models = rq(@formula(Y ~ X1 + X2 + X3), df; τ=[0.25, 0.50, 0.75])","category":"page"},{"location":"quickstart/#","page":"Quickstart Guide","title":"Quickstart Guide","text":"This call yields a QuantRegModels object containing three quantile regression models, fit at the 0.25th, 0.50th, and 0.75th quantiles. The individual models can be accessed by using the τ values as indices. For example,","category":"page"},{"location":"quickstart/#","page":"Quickstart Guide","title":"Quickstart Guide","text":"julia> models[0.25]","category":"page"},{"location":"quickstart/#","page":"Quickstart Guide","title":"Quickstart Guide","text":"will return the quantile regression model fit at the 0.25th quantile.","category":"page"},{"location":"quickstart/#Configuring-models-1","page":"Quickstart Guide","title":"Configuring models","text":"","category":"section"},{"location":"quickstart/#","page":"Quickstart Guide","title":"Quickstart Guide","text":"Users can provide more detailed model specifications as keyword arguments to the quantile regression command. Any specification fields in the QuantRegModel type can be accepted as a keyword argument. The following command fits a quantile regression model at the 0.80th quantile using the Frisch-Newton and computes 80% confidence intervals under the assumption that the conditional quantile function is locally (in tau) linear (in x) (i.e. regression errors are not i.i.d.).","category":"page"},{"location":"quickstart/#","page":"Quickstart Guide","title":"Quickstart Guide","text":"julia> rq(@formula(Y ~ X1 + X2 + X3), df; τ=0.80, fitmethod=\"fn\", α=0.20, iid=false)","category":"page"},{"location":"quickstart/#","page":"Quickstart Guide","title":"Quickstart Guide","text":"A full description of the available configurations and their defaults is available in the type(QuantRegModel) reference.","category":"page"},{"location":"quickstart/#Under-the-hood-1","page":"Quickstart Guide","title":"Under the hood","text":"","category":"section"},{"location":"quickstart/#","page":"Quickstart Guide","title":"Quickstart Guide","text":"In code, the rq functions is really a wrapper for three functions that constitute a common work flow. It first constructs a QuantRegModel object according to the specifications provided. Then, it fits the model in place according to the specifications in the generated model object with func(fit!). This function call updates the QuantRegFit object stored in the model object. Finally, it computes inference for the model according to the specifications in the generated model object with func(compute_inf!). This function call updates the QuantRegInf object stored in the model object.","category":"page"},{"location":"#QuantReg.jl-Documentation-1","page":"QuantReg.jl Documentation","title":"QuantReg.jl Documentation","text":"","category":"section"},{"location":"#","page":"QuantReg.jl Documentation","title":"QuantReg.jl Documentation","text":"This package provides types and functions for specifying, fitting, and computing inference for quantile regression models. This package was originally modeled after the R package quantreg, authored by Roger Koenker. Currently, the packages implements many of the abstractions for statistical models from StatsBase.jl and requires tabular data in the form of a DataFrame.","category":"page"},{"location":"#","page":"QuantReg.jl Documentation","title":"QuantReg.jl Documentation","text":"Contributions to this packages that incorporate more features from the R quantreg package, add new features, or refine existing features are welcomed and encouraged on the project's GitHub. Please feel free to open an issue or pull request!","category":"page"}]
}
