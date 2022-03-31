
#' @title Fit a PARTS GLS model, with maximum likelihood spatial parameters
#'
#' @details Estimate spatial parameters, via maximum likelihood, from data
#' rather than from time series residuals; Fit a GLS with these specifications.
#'
#' @param formula a model formula, passed to \code{fitGLS}
#' @param data an optional data frame environment in which to search for
#' variables given by \code{formula}; passed to \code{fitGLS}
#' @param coords a numeric coordinate matrix or data frame, with two columns and
#' rows corresponding to each pixel
#' @param distm_FUN a function to calculate a distance matrix from \code{coords}
#' @param covar_FUN a function to estimate distance-based covariances
#' @param start a named vector of starting values for each parameter to be estimated;
#' names must match the names of arguments in \code{covar_FUN} or "nugget"
#' @param fixed an optional named vector of fixed parameter values; names
#' must match the names of arguments in \code{covar_FUN} or "nugget"
#' @param opt.only logical: if TRUE, execution will halt after estimating the parameters;
#' a final GLS will not be fit with the estimated parameters
#' @param formula0,save.xx,save.invchol,no.F arguments passed to \code{fitGLS}
#' for final GLS output
#' @param trans optional list of functions for transforming the values in
#' \code{start} or \code{fixed} in order to constrain the parameter space within
#' \code{optim}
#' @param backtrans optional list of functions for back-transforming parameters
#' to their correct scale (for use with \code{trans})
#' @param ... additional arguments passed to \code{stats::optim()}
#'
#' @details \code{fitGLS_opt} fits a GLS by estimating spatial parameters from
#' data. \code{\link{fitCor}}, combined with \code{\link{fitGLS}(nugget = NA)},
#' gives better estimates of spatial parameters, but time-series residuals may
#' not be available in all cases. In these cases, spatial parameters can be
#' estimated from distances among points and a response vector. Mathematical
#' optimization of the log likelihood of different GLS models are computed by
#' calling \code{optim()} on \code{fitGLS}.
#'
#' Distances are calculated with \code{distm_FUN} and a covariance matrix is
#' calculated from these distances with \code{covar_FUN}. Arguments to to
#' \code{covar_FUN}, except distances, are given by \code{start} and \code{fixed}.
#' Parameters specified in \code{start} will be be estimated while those given
#' by \code{fixed} will remain constant throughout fitting. Parameter names in
#' \code{start} and \code{fixed} should exactly match the names of arguments in
#' \code{covar_FUN} and should not overlap (though, \code{fixed} takes precedence).
#'
#' In addition to arguments of \code{covar_FUN} a "nugget" component can
#' also be occur in \code{start} or \code{fixed}. If "nugget" does not occur
#' in either vector, the GLS are fit with \code{nugget = 0}. A zero nugget also
#' allows much faster computation, through recycling the common
#' inverse cholesky matrix in each GLS computation. A non-zero nugget requires
#' inversion of a different matrix at each iteration, which can be
#' substantially slower.
#'
#' If \code{opt.only = FALSE}, the estimated parameters are used to fit the final
#' maximum likelihood GLS solution with \code{fitGLS()} and arguments
#' \code{formula0}, \code{save.xx}, \code{save.invchol}, and \code{no.F}.
#'
#' Some parameter combinations may not produce valid covariance matrices. During
#' the optimization step messages about non-positive definitive V may result on
#' some iterations. These warnings are produced by \code{fitGLS} and NA
#' log-likelihoods are returned in those cases.
#'
#' Note that \code{fitGLS_opt} fits multiple GLS models, which requires
#' inverting a large matrix for each one (unless a fixed 0 nugget is used).
#' This process is very computationally intensive and may take a long time to
#' finish depending upon your machine and the size of the data.
#'
#' Sometimes \code{optim} can have a difficult time finding a reasonable solution
#' and without any constraits on parameter space (with certain algorithms), results
#' may even be nonsensical. To combat this, \code{fitGLS_opt} has the arguments
#' \code{trans} and \code{backtrans} which allow you to transform
#' (and back-transform) parameters to a different scale. For example, you may
#' want to force the 'range' parameter between 0 and 1. The logit function can
#' do just that, as its limits are -Inf and Inf as x approaches 0 and 1,
#' respectively. So, we can set \code{trans} to the logit function:
#' \code{trans = list(range = function(x)log(x/(1-x)))}. Then we need to set
#' \code{backtrans} to the inverse logit function to return a parameter value
#' between 0 and 1: \code{backtrans = list(range = function(x)1/(1+exp(-x)))}.
#' This will force the optimizer to only search for the range parameter in the
#' space from 0 to 1. Any other constraint function can be used for \code{trans}
#' provided that there is a matching back-transformation.
#'
#'
#' @seealso \code{\link{fitCor}} for estimating spatial parameters from time
#' series residuals; \code{\link{fitGLS}} for fitting GLS and with the option
#' of estimating the maximum-likelihood nugget component only.
#'
#' @return If \code{opt.only = TRUE}, \code{fitGLS_opt} returns the
#' output from \code{stats::optim()}: see it's documentation for more details.
#'
#' Otherwise, a list with two elements is returned:
#'
#' \describe{
#'     \item{opt}{output from \code{optim}, as above}
#'     \item{GLS}{a "remoteGLS" object. See \code{\link{fitGLS}} for more details.}
#' }
#'
#' @examples
#' \donttest{
#' ## read data
#' data(ndvi_AK10000)
#' df = ndvi_AK10000[seq_len(200), ] # first 500 rows
#'
#' ## estimate nugget and range (very slow)
#' fitGLS_opt(formula = CLS_coef ~ 0 + land, data = df,
#'             coords = df[, c("lng", "lat")], start = c(range = .1, nugget = 0),
#'             opt.only = TRUE)
#'
#' ## estimate range only, fixed nugget at 0, and fit full GLS (slow)
#' fitGLS_opt(formula = CLS_coef ~ 0 + land, data = df,
#'              coords = df[, c("lng", "lat")],
#'              start = c(range = .1), fixed = c("nugget" = 0),
#'              method = "Brent", lower = 0, upper = 1)
#'
#' ## constrain nugget to 0 and 1
#' logit <- function(p) {log(p / (1 - p))}
#' inv_logit <- function(l) {1 / (1 + exp(-l))}
#'
#' fitGLS_opt(formula = CLS_coef ~ 0 + land, data = df,
#'            coords = df[, c("lng", "lat")],
#'            start = c(range = .1, nugget = 1e-10),
#'            trans = list(nugget = logit), backtrans = list(nugget = inv_logit),
#'            opt.only = TRUE)
#' }
#' @export
fitGLS_opt <- function(formula, data = NULL, coords, distm_FUN = "distm_scaled",
                       covar_FUN = "covar_exp",
                       start = c(range = .01, nugget = 0),
                       fixed = c(), opt.only = FALSE,
                       formula0 = NULL, save.xx = FALSE, save.invchol = FALSE,
                       no.F = TRUE,
                       trans = list(), backtrans = list(),
                       ...){
  call = match.call()

  # transform variables to constrain them in optim, if needed
  is.trans = FALSE
  if(length(trans) > 0 & length(backtrans) > 0){
    stopifnot(length(trans) == length(backtrans))
    stopifnot(names(trans) == names(backtrans))
    stopifnot(all(names(trans) %in% c(names(start), names(fixed))))
    is.trans = TRUE
    for(n in names(trans)){
      if(n %in% names(start)){
        start[n] = do.call(trans[[n]], list(start[n]))
      }
      if(n %in% names(fixed)){
        fixed[n] = do.call(trans[[n]], list(fixed[n]))
      }
    }
  }

  # create a list of arguments to pass to do.call
  arg.list <- list(par = start[! names(start) %in% names(fixed)], #parameters to optimze
                   fn = fitGLS_opt_FUN, formula = formula, data = data,
                   coords = coords, covar_FUN = covar_FUN,
                   distm_FUN = distm_FUN,
                   is.trans = is.trans, backtrans = backtrans
  )

  # append fixed parameters, if they exist
  if(length(fixed) > 0){
    arg.list <- append(arg.list, list(fp = fixed))
  }

  # append arguments given by ... to the argument list
  arg.list = append(arg.list, list(...)) # add additional arguments to arg list

  # call optim, and pass arguments
  opt.out <- do.call(optim, args = arg.list)

  # back-transform the parameter values to their original scale
  if(is.trans){
    for(n in names(backtrans)){
      if(n %in% names(opt.out$par)){
        opt.out$par[n] = do.call(backtrans[[n]], list(opt.out$par[n]))
      }
      if(n %in% names(start)){
        start[n] = do.call(backtrans[[n]], list(start[n]))
      }
      if(n %in% names(fixed)){
        fixed[n] = do.call(backtrans[[n]], list(fixed[n]))
      }
    }
  }

  if(opt.only){
    return(opt.out)
  } else {
    names(opt.out$par) = names(start)
    spars = opt.out$par[!names(opt.out$par) %in% "nugget"]
    nug = ifelse(test = "nugget" %in% names(opt.out$par), yes = opt.out$par["nugget"],
                 no = ifelse(test = "nugget" %in% names(fixed), yes = fixed["nugget"],
                             no = 0))
    GLS.out = fitGLS(formula = formula, data = data, formula0 = formula0,
                     save.xx = save.xx, save.invchol = save.invchol,
                     logLik.only = FALSE, no.F = no.F, coords = coords,
                     distm_FUN = distm_FUN ,covar_FUN = covar_FUN,
                     nugget = nug, covar.pars = spars)
    GLS.out$call = call
    return(list(opt = opt.out, GLS = GLS.out))
  }
}

#' Function that fitGLS_opt optimizes over
#'
#' @param op a named vector of parameters to be optimized
#' @param fp a named vector of fixed parameters
#' @param coords a coordinate matrix
#' @param covar_FUN a covariance function
#' @param distm_FUN a distm function
#' @param formula GLS model formula
#' @param data data source
#' @param is.trans logical: are any of the values in \code{op} or \code{fp}
#' transformed, needing back-transformation?
#' @param backtrans optional: a named list of functions used to backtransform any element
#' of \code{op} or \code{fp}. Names must correspond to names in \code{op}
#' or \code{fp}.
#'
#' @return \code{fitGLS_opt_FUN} returns the negative log likelihood of a GLS,
#' given the parameters in \code{op} and \code{fp}
#'
#' @examples
#' \dontrun{
#' data(ndvi_AK10000)
#' df = ndvi_AK10000[seq_len(200), ] # first 500 rows
#' coords = df[, c("lng", "lat")]
#' remotePARTS:::fitGLS_opt_FUN(op = c(range = .1, nugget = .2),
#'                              formula = CLS_coef ~ 0 + land, data = df,
#'                              coords = coords)
#' remotePARTS:::fitGLS_opt_FUN(op = c(range = .1), fp = c(nugget = 0),
#'                              formula = CLS_coef ~ 0 + land, data = df,
#'                              coords = coords)
#'
#' logit <- function(p) {log(p / (1 - p))}
#' inv_logit <- function(l) {1 / (1 + exp(-l))}
#'
#' # input logit-transformed range parameters
#' remotePARTS:::fitGLS_opt_FUN(op = c(range = .1, nugget = logit(.2)),
#'                              formula = CLS_coef ~ 0 + land, data = df,
#'                              coords = coords, is.trans = TRUE,
#'                              backtrans = list(nugget = inv_logit))
#' # transformed range and nugget
#' remotePARTS:::fitGLS_opt_FUN(op = c(range = logit(.1), nugget = logit(.2)),
#'                              formula = CLS_coef ~ 0 + land, data = df,
#'                              coords = coords, is.trans = TRUE,
#'                              backtrans = list(nugget = inv_logit, range = inv_logit))
#' }
fitGLS_opt_FUN <- function(op, fp, formula, data = NULL, coords, covar_FUN = "covar_exp", distm_FUN = "distm_scaled",
                           is.trans = FALSE, backtrans = list()){
  ## combine the optimized and fixed parameters into one vector
  all.pars <- if(missing(fp)){op}else{c(op, fp)}

  if(is.trans == TRUE) {
    if (is.null(names(backtrans)) | !all(names(backtrans) %in% names(all.pars))) {
      stop("backtrans mut be a list with element names matching names in op or fp")
    }
    for (n in names(backtrans)){
      all.pars[n] = do.call(backtrans[[n]], list(all.pars[n]))
    }
  }

  ## extract the nugget
  nug = ifelse("nugget" %in% names(all.pars), all.pars["nugget"], 0)
  ## extract the non-nugget parameters
  sp.pars = as.list(all.pars[!names(all.pars) %in% "nugget"])
  ## calculate covariance
  cov.f = match.fun(covar_FUN)
  dist.f = match.fun(distm_FUN)
  V = dist.f(coords) # distance
  args = append(list(d = V), as.list(sp.pars))
  V = do.call(cov.f, args) # replace with covariance
  ## Calculate log-likelihood
  logLik = suppressWarnings(fitGLS(formula = formula, data = data, V = V, formula0 = NULL,
                  save.xx = FALSE, save.invchol = FALSE, logLik.only = TRUE, no.F = TRUE,
                  nugget = nug))
  return(-logLik)
}

