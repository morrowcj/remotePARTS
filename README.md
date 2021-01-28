
<!-- README.md is generated from README.Rmd. Please edit that file -->

# remotePARTS

<!-- badges: start -->

[![Lifecycle:
maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
<!-- badges: end -->

`remotePARTS` is an `R` package that contains tools for running
spatio-temporal auto regression analyses on large remotely-sensed data
sets by partitioning data into managable chunks.

## Description

This package is based on the PARTS method for analyzing spatially
autocorrelated time series (Ives et al., in prep).

## Instalation

To install the package and it’s dependencies, use the following R code:

``` r
install.packages("remotes") # or install.packages("devtools")
remotes::install_github("morrowcj/remotePARTS", dependencies = "Depends", build_vignettes = TRUE)
```

Then, upon successful installation, load the package with
`library(remotePARTS)`.

The latest version of
[Rtools](https://cran.r-project.org/bin/windows/Rtools/) is required for
Windows and C++11 is required for other systems.

### Windows installation notes/troubleshooting:

On my Windows 10 PC, I had to change the permission settings for R in
order for `install_github()` to work:

1)  right click on “C:\\Program Files\\R\\R-4.0.2\\library\\base”

2)  click properties

3)  Select “Security” Tab

4)  find and select select “Users” in the “Group or user names” scroll
    menu

5)  tick “Full control”

## Example usage

For examples on how to use `remotePARTS` in it’s current state, see the
`Alaska` vignette by using the following R code:

``` r
vignette("Alaska")
```

## Testing

`remotePARTS` is currently in early development. Stability and
efficiency tests are ongoing and improvements occur incrementally.
Automated tests have exists for some, but not all, of the functions and,
at present, tests have only been conducted in limited environments.

    Test Environments:
    
    1. Windows 10.0.19041 x64

Please report any bugs or suggested features as git issues.

## Data formats

In general, this package requires data to be contained in R matrices.
Data stored in other formats such as image files or rasters will need to
be converted. Furthermore, it is highly recommended that map pixels
without data values (e.g. water in analyses of land patterns) be removed
entirely. This is especially true when using the partitioned form of the
analysis because the matrix multiplication can’t handle missing data or
unequal dimensions.

The following resource demonstrates how to manipulate raster files in R:
[Geospatial raster with R data
carpentry](http://datacarpentry.org/r-raster-vector-geospatial/)

## Planned Features / To-do

Since this package in developmental stages, there are many features that
are currently unimplemented. This section will keep track of the
features and design implementations that I plan to include or change in
the next version as well.

  - [ ] allow users to, optionally, input parameters (e.g. `r` and `a`
    in the exponential-power function) instead of fitting ML parameters.

  - [ ] add example for testing “is there an overall time trend” to the
    vignette

  - [ ] make providing distance matrix **optional** instead of required
    for `fitGLS.partition_rcpp()` and the partitioned method as a whole.
    (**note** this has been demonstrated in the vignette but has not yet
    been turned into a function)

  - [ ] include parallelization and distributed comptuting options. If
    these are not natively implemented (i.e. using openMP in C++), then
    examples of how to make it work with other parallelization tools
    should be provided.

  - [ ] more explicit handling of missing data: “How should a constant
    time series be treated?”; “What happends if there is a missing data
    point within a single time series?”

  - [ ] possibly change the CLS function so that it reads 1 line of data
    at a time to save memory. Also, using `RcppEigen::fastLM()` may be
    better than `lm()` in terms of speed.

  - [x] **Break up C++ functions into more than 1 file**

  - [ ] **WRITE TESTS FOR EVERY FUNCTION\!\!\!** - most functions have
    at least some automated tests but more are needed.

  - [ ] update documentation for every function to include output format
    and example usage. - most function have decent documentation but
    some of the documentation is inconsistent across the package and
    some functions are under-documented.

  - [ ] replace `fitGLS` with `fitGLS2` code and change the way all
    other C++ functions handle lists. The C++ code in `fitGLS2` modifies
    an existing list made in R rather than building one within the C++
    code.

  - [ ] create S3 constructors and methods (i.e. `print()`, `summary()`,
    etc.) for all relevant functions.
    
      - [x] `remoteGLS()` constructs a GLS object and
        `print.remoteGLS()` prints a compact and summarized display.
    
      - [x] `remoteCLS()` constructor and methods for CLS objects. Both
        the pixel-level `fitCLS` and the map-level `fitCLS.map` have S3
        methods.
    
      - [x] `remoteAR()` constructor and methods for AR\_REML objects.
        Both pixel and map-level `fitAR()` have S3 methods.
    
      - [ ] `PARTmat()`: partition matrix for the partitioned GLS
        method.

If there are any additional features that you would like to see
implemented, or bugs/issues that you ran into, please submit an issue on
github.

## Installation (OLD)

Currently to install this package, the best way is to install with the
`remotePARTS_[version].tar.gz` file created with `R CMD check`.

Once a user has the tar.gz file they can install it with

    install.packages("remotePARTS_[version].tar.gz", repos = NULL, type = "source")

and then from the R console with load it with

    library(remotePARTS)

<!-- Eventually, the following lines should replace the above installation info: -->

<!-- You can install the released version of remotePARTS from 
[CRAN](https://CRAN.R-project.org) with: -->

<!-- ``` r -->

<!-- install.packages("remotePARTS") -->

<!-- ``` -->

<!-- And the development version from [GitHub](https://github.com/) with: -->

<!-- ``` r -->

<!-- # install.packages("devtools") -->

<!-- devtools::install_github("morrowcj/remotePARTS") -->

<!-- ``` -->

<!-- ## Example -->

<!-- This is a basic example which shows you how to solve a common problem: -->

<!-- ```{r example} -->

<!-- library(remotePARTS) -->

<!-- ## basic example code -->

<!-- ``` -->

<!-- What is special about using `README.Rmd` instead of just `README.md`? You can include R chunks like so: -->

<!-- ```{r cars} -->

<!-- summary(cars) -->

<!-- ``` -->

<!-- You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. -->

<!-- You can also embed plots, for example: -->

<!-- ```{r pressure, echo = FALSE} -->

<!-- plot(pressure) -->

<!-- ``` -->

<!-- In that case, don't forget to commit and push the resulting figure files, so they display on GitHub! -->
