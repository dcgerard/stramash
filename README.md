
<!-- README.md is generated from README.Rmd. Please edit that file -->
STRAMASH: Specifying The Right Amount of Mixture-likelihood for Adaptive SHrinkage
==================================================================================

This package is similar to and depends on the `ashr` package (<https://github.com/stephens999/ashr>) but also allows the error distribution to be any mixture of normals or uniforms. See that link for details on the model.

Right now, I would only suggest using a mixture of normals as your error distribution. The mixture of uniforms seems to work OK if you have a huge grid, but not well otherwise.

Installation
============

To install, run the following code in R

``` r
install.packages("devtools")
devtools::install_github("dcgerard/stramash")
```
