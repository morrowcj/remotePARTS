## R CMD check results

There were no ERRORs or WARNINGs.

There was 1 NOTE associated with the first-time submission: 

```
  Maintainer: 'Clay Morrow <morrowcj@outlook.com>'

  New submission
```

### Test Environments

This version of the package was tested, using github actions, with the following call to `rcmdcheck` on 2022-06-30:

```
rcmdcheck::rcmdcheck(args = c("--no-manual", "--as-cran"), error_on = "warning", check_dir = "check")
```

Checks were conducted on the following systems:

  * windows-latest (release)
  
  * macOS-latest (release)
  
  * ubuntu-20.04 (release)
  
Additionally, the package was checked on the development version of R for Windows 
(on 2022-06-30) with `devtools::check_win_devel()` and the same results as above
were returned. 

## Downstream dependencies

There are currently no downstream dependencies for this package.

## Other Comments

* On the ubuntu-20.04 (release) system, an additional NOTE occurred, indicating
that the libraries are much larger than on the other the systems:

```
  ❯ checking installed package size ... NOTE
    installed size is 32.7Mb
    sub-directories of 1Mb or more:
      libs  31.3Mb
```

* Some users have had trouble building the vignette on personal machines 
(using `devtools::install_github()`), but I have not had any trouble on any 
of the systems I have personally tested with. 

* Some examples in `partGLS.Rd` and `fitGLS_opt.Rd` can take a bit longer to run
(10-20 seconds total) since the functions documented in these files need to 
invert multiple large matrices. 
