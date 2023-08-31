## R CMD check results

There were no ERRORs or WARNINGs.

There was 1 NOTE associated with a URL that was flagged as invalid.

```
  Maintainer: 'Clay Morrow <morrowcj@outlook.com>'
    
  New submission

  Found the following (possibly) invalid URLs:
    URL: https://support.posit.co/hc/en-us/articles/200486498-Package-Development-Prerequisites
      From: README.md
      Status: 403
      Message: Forbidden
```

I checked the URL (on 2023-08-30), which is in the README file, and it works properly.

### Test Environments

This version of the package was tested on 2023-08-23, using `check-r-package` GitHub Action from the 
[r-lib repository](https://github.com/r-lib/actions), with the following call to `rcmdcheck`:

```
rcmdcheck::rcmdcheck(args = c("--no-manual", "--as-cran"), error_on = "error", check_dir = "check")
```

Checks were conducted on the following systems, using :

  * windows-latest (release)
  
  * macOS-latest (release)
  
  * ubuntu-latest (devel)
  
  * ubuntu-latest (oldrel-1)
  
  * ubuntu-latest (release)
  
----  
  
Additionally, the package was checked on the development version of R for Windows 
(on 2023-08-30) with `devtools::check_win_devel()` and the same results as above
were returned. 

## Downstream dependencies

There are currently no downstream dependencies for this package.

## Other Comments

* Some users have had trouble building the vignette on personal machines 
(using `devtools::install_github()`), but I have not had any trouble on any 
of the systems I have personally tested with. 

* Some examples in `partGLS.Rd` and `fitGLS_opt.Rd` can take a bit longer to run
(10-20 seconds total) since the functions documented in these files need to 
invert multiple large matrices. 

## Updates after intitial submission

After initial submission, I received correspondence from Benjamin Altmann at CRAN. In response to their requested 
changes I have:

* removed trailing white space from all lines in the DESCRIPTION file
* added `@return` values to all exported functions. 
* removed `:::` from all documentation
* removed all `@examples` from un-exported functions
* removed examples wrapped in `\dontrun{}` and replaced them with `\donttest{}` where necessary.
* removed commented out code in examples


