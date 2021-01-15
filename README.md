
<!-- README.md is generated from README.Rmd. Please edit that file -->

# remoteSTAR

<!-- badges: start -->

<!-- badges: end -->

The goal of remoteSTAR is to …

## Description

`remoteSTAR` is currently in early development. Organization is less
than ideal, official unit tests are not present, and C++ code has
occasionally exhibited bugs.

This package is not stable, consider it a beta. Please report any
comments or bugs either directly to me `morrow5@wisc.edu` or through
github: <https://github.com/morrowcj/remoteSTAR> (it is a private repo
so please send me a github username and I’ll add you to the repo).

**Note** I plan to eventually show examples of importing rasters. Here
is a lovely resource for manipulating rasters in R, including how to
convert them to data frames: [Geospatial raster with R data
carpentry](http://datacarpentry.org/r-raster-vector-geospatial/)

## Installation

Currently to install this package, the best way is to install with the
`remoteSTAR_[version].tar.gz` file created with `R CMD check`.

Once a user has the tar.gz file they can install it with

    install.packages("remoteSTAR_[version].tar.gz", repos = NULL, type = "source")

and then from the R console with load it with

    library(remoteSTAR)

<!-- Eventually, the following lines should replace the above installation info: -->

<!-- You can install the released version of remoteSTAR from 
[CRAN](https://CRAN.R-project.org) with: -->

<!-- ``` r -->

<!-- install.packages("remoteSTAR") -->

<!-- ``` -->

<!-- And the development version from [GitHub](https://github.com/) with: -->

<!-- ``` r -->

<!-- # install.packages("devtools") -->

<!-- devtools::install_github("morrowcj/remoteSTAR") -->

<!-- ``` -->

## Example usage

For examples on how to use `remoteSTAR` in it’s current state, see the
`Alaska` vignette by using the following R code:

    vignette("Alaska")

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

  - [ ] **Break up C++ functions into more than 1 file**

  - [ ] **WRITE TESTS FOR EVERY FUNCTION\!\!\!**

  - [ ] update documentation for every function to include output format
    and example usage.

  - [ ] replace `fitGLS` with `fitGLS2` code and change the way all
    other C++ functions handle lists. The C++ code in `fitGLS2` modifies
    an existing list made in R rather than building one within the C++
    code.

  - [ ] create S3 constructors and methods (i.e. `print()`, `summary()`,
    etc.) for all functions
    
      - [x] `remoteGLS()` constructs a GLS object and
        `print.remoteGLS()` prints a compact and summarized display.
    
      - [ ] `remoteCLS()` constructor and methods for CLS objects
    
      - [ ] `remoteAR()` constructor and methods for AR\_REML objects
    
      - [ ] `PARTmat()`: partition matrix for the partitioned GLS
        method.

If there are any additional features that you would like to see
implemented, or bugs/issues that you ran into, please submit an issue on
github.

<!-- ## Example -->

<!-- This is a basic example which shows you how to solve a common problem: -->

<!-- ```{r example} -->

<!-- library(remoteSTAR) -->

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

# <!-- In that case, don't forget to commit and push the resulting figure files, so they display on GitHub! -->