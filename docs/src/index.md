# QuantReg.jl Documentation

This package provides types and functions for specifying, fitting, and computing inference
for quantile regression models. This package was modeled after and employs FORTRAN libraries
from the [R package quantreg, authored by Roger Koenker.](https://cran.r-project.org/web/packages/quantreg/index.html) 
Currently, the package implements many of the abstractions for statistical models from 
StatsBase.jl and requires tabular data in the form of a DataFrame.

Contributions to this package that incorporate more features from the R quantreg package,
add new features, or refine existing features are welcomed and encouraged
[on the project's GitHub.](https://github.com/fogarty-ben/QuantReg.jl) Please feel free to
open an issue or pull request!

!!! note
    This project is currently pre-release pending a solution to reliably cross-compiling
    the underlying FORTRAN libraries from quantreg. This issue should be resolved with a
    few small changes pending the merge of
    [this pull request to Yggdrasil](https://github.com/JuliaPackaging/Yggdrasil/pull/683)
    into master.
    
    For the time being, however, this repository contains the FORTRAN source needed to
    execute the project. Users should build this source into dynamic libraries using the
    commands in the [Quickstart Guide](@ref) to run the program.
