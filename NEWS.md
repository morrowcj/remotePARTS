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

# v0.1.1.9006

* updated fitGLS.partition functions to prevent them from hogging memory. This
means that the cholesky inversions cannot be reused and need to be recalculated,
leading to longer compute times. 

* fixed installation bug where vignette would not build.

* added multithreading capability, through Eigen, to the GLS suite of functions, 
though this seems to have limited utility since eigen already uses all available
threads for matrix operations.

* added functionality that allows co-estimation of spatial autocorrelation
and nugget parameters from the data - i.e., fixed `optimize_GLS()`. 

* added lower and upper limits to optimize_GLS (and new partitioned version)

* added new vignette describing optimize_GLS's new funtionality

# v0.1.1.9008

* changed the way that `fitGLS.partiton.mc()` handles cross-paritioning to be
more reproducible

* fixed a typo in source code of `crosspart_worker()` that caused `rcoef` to be
calculated incorrectly. This change has major implications for p-values obtained
from `cor_t()` and `cor_t.test()` in previous versions. To update the GLS
object from previous versions without re-running the entire dataset, you can 
re-run a limmited dataset to calculate `rcoef` and then replace the relevant 
values in the old code: 

```
# Original partition matrix:
og.parts = sample_partitions(npix = nrow(data), partsize = 2000)
# Original GLS, using mincross = 6:
old.gls <- fitGLS.partition.mc(..., part_mat = og.parts, npart = ncol(og.parts),
                               partsize = nrow(og.parts), mincross = 6)
# New GLS to get rcoef, using npart = mincross:
new.gls = fitGLS.partition.mc(..., part_mat = og.parts[1:6], npart = 6,
                              partsize = nrow(og.parts), mincross = 6)
# Update relevant statistics in original gls:
old.gls$cross.stats <- new.gls$cross.stats
old.gls$overall.stats$rcoefmean <- new.gls$overall.stats$rcoefmean
old.gls$overall.stats$meanstats["rSSRmean"] = new.gls$overall.stats["rSSRmean"]
old.gls$overall.stats$meanstats["rSSEmean"] = new.gls$overall.stats["rSSEmean"]
# re-run t-test:
cor_t.test(old.gls)
```

# v0.1.1.9010

* Fixed a bug in `optimize_GLS` where `V.meth = 'exponential-power'` was not 
working properly.

# v0.1.1.9011

## Restructuring: 

* [x] update fitCLS
  - now handles formula
  - now outputs lm-style list object
* [ ] update fitCLS.map to handle new fitCLS output
* [ ] update fitAR to have similar output to fitCLS
* [ ] update fitAR.map to 
  - handle new fitAR output and 
  - have similar output to fitCLS.map

