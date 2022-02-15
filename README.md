
<!-- README.md is generated from README.Rmd. Please edit that file -->

# remotePARTS

<!-- badges: start -->

[![License: GPL
v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

[![Lifecycle:
maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)

<!-- [![Travis build status](https://travis-ci.com/morrowcj/remotePARTS.svg?branch=master)](https://travis-ci.com/morrowcj/remotePARTS) -->
<!-- [![Travis build status](https://travis-ci.com/morrowcj/remotePARTS.svg?branch=master)](https://travis-ci.org/github/morrowcj/remotePARTS) -->

[![R-CMD-check](https://github.com/morrowcj/remotePARTS/workflows/R-CMD-check/badge.svg)](https://github.com/morrowcj/remotePARTS/actions)
<!-- badges: end -->

`remotePARTS` is an `R` package that contains tools for running
spatio-temporal auto regression analyses on large remotely-sensed data
sets by partitioning data into manageable chunks.

## Description

This package is based on the PARTS method for analyzing spatially
autocorrelated time series (Ives et al., 2021).

## Instalation

To install the package and it’s dependencies, use the following R code:

``` r
install.packages("devtools") # ensure you have the latest devtools
devtools::install_github("morrowcj/remotePARTS")
```

Then, upon successful installation, load the package with
`library(remotePARTS)`.

The latest version of
[Rtools](https://cran.r-project.org/bin/windows/Rtools/) is required for
Windows and C++11 is required for other systems.

## Example usage

For examples on how to use `remotePARTS` in it’s current state, see the
`Alaska` vignette by using the following R code:

``` r
vignette("Alaska")
```

Note that the vignette needs to be build when installing and may require
the `build_vignettes = TRUE` argument when istalling with
`install_github()`.

# References

Ives, Anthony R., et al. “Statistical inference for trends in
spatiotemporal data.” Remote Sensing of Environment 266 (2021): 112678.
