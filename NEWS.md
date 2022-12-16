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
