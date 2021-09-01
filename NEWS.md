# v0.1.1.9002

* Added a `NEWS.md` file to track changes to the package.

* Fixed some installation bugs

* Fixed some build bugs (Travis CI)

* added `optimize_GLS()` to optimize the GLS over `r`, `a,` and `nugget`

# v0.1.1.9003

* added multicore function with `fitGLS.partition.mc()`, which uses 
functionality from the `foreach` and `doParallel` packages.

# v0.1.1.9004

* Fixed a bug with `fitGLS.partition()` functions whereby X0 could not be 
specified within the `part_f` function. `part_data()` and `part_csv()` have also
been corrected to output X0. 

# v0.1.1.9005

* updated fitGLS.partition functions to prevent them from hogging memory. This
means that the cholesky inversions cannot be reused and need to be recalculated,
leading to longer compute times. 

* fixed installation bug where vignette would not build.

* added multithreading capability, through Eigen, to the GLS suite of functions.
