# v1.0.4.2

* updated the parallel defaults for all functions. Now, `ncores=NA` is equivalent to `ncores=1` instead of using
the `C++` compiler default for `Eigen`.

* Additionally, we removed parallel examples from `partGLS` due to uncertainty in our understanding of CRAN's
test setup using more cores than any of our test machines.


# v1.0.4

* fixed bug where R would crash if 0-intercept model given

* fixed bug where CPU usage would not respect `ncores` for some functions

* updated comments to cran to reflect new tests in preparation for submission.

# v1.0.3
* fixed bug where `fitGLS_opt()` would fail if one of the iterations in the 
optimization loop raises an error.
* added explicit call to "BFGS" method while using `fitGLS_opt()` in the 
vignette and lowered tolerance to improve speed.

# v1.0.2

* fixed a bug where `fitCor()` was calculating a full distance matrix instead
of one among only the sampled pixels.
* added option to suppress saving the `nls` model in `fitCor()`
* added $\rho$ to the output of `fitAR_map()`

# v1.0.1

fixed bug where `fitGLS_partition` would not work properly when `ncores > 1` and
`formula == formula0`. 

# v1.0.0

`remotePARTS` has been completely restructured from previous development builds. 
Please see the Alaska vignette with `vignette("Alaska")` for usage demonstrations. 
