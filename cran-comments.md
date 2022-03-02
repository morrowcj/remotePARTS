## R CMD check results
There were no ERRORs or WARNINGs.

There was 1 NOTE: 

  * "Maintainer: 'Clay Morrow <morrowcj@outlook.com>'"

  This is a first-time submission.
  
### Test Environments
This version of the package was tested, using github actions, with the code

```
rcmdcheck::rcmdcheck(args = c("--no-manual", "--as-cran"), error_on = "warning", check_dir = "check")
```

Checks were conducted on the following systems:

  * windows-latest (release)
  
  * macOS-latest (release)
  
  * ubuntu-20.04 (release)
  
Additionally, the package was tested on the development version of R for Windows 
(on 2022-03-01) with 

```
devtools::check_win_devel()
```

*R-devel results here*

## Downstream dependencies
There are currently no downstream dependencies for this package.

## Other Comments

* Some users have had trouble building the vignette on personal machines 
(using `devtools::install_github()`), but I have not had any trouble on any 
of the systems I have personally tested on. 

* Some examples in `partGLS.Rd` and `fitGLS_opt.Rd` can take a bit longer to run
(10-20 seconds total) since the functions documented in these files need to 
invert multiple large matrices. 
