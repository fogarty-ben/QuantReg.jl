# QuantReg.jl
[![Build Status](https://travis-ci.org/fogarty-ben/QuantReg.jl.png?branch=master)](https://travis-ci.org/{ORG-or-USERNAME}/{REPO-NAME})
[![Stable Docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://fogarty-ben.github.io/QuantReg.jl/stable)


This package provides types and functions for specifying, fitting, and computing inference
for quantile regression models. This package was originally modeled after the [R package
quantreg, authored by Roger Koenker.](https://cran.r-project.org/web/packages/quantreg/index.html)
Currently, the packages implements many of the abstractions for statistical models from StatsBase.jl
and requires tabular data in the form of a `DataFrame`.

Contributions to this packages that incorporate more features from the R `quantreg` package,
add new features, or refine existing features are welcomed and encouraged
[on the project's GitHub.](https://github.com/fogarty-ben/QuantReg.jl) Please feel free to
open an issue or pull request!
