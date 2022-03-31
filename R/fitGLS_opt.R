
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
#' }
#' @export
fitGLS_opt <- function(formula, data = NULL, coords, distm_FUN = "distm_scaled",
                       covar_FUN = "covar_exp",
                       start = c(range = .01, nugget = 0),
                       fixed = c(), opt.only = FALSE,
                       formula0 = NULL, save.xx = FALSE, save.invchol = FALSE,
                       no.F = TRUE,
                       ...){
  call = match.call()


  # create a list of arguments to pass to do.call
  arg.list <- list(par = start[! names(start) %in% names(fixed)], #parameters to optimze
                   fn = fitGLS_opt_FUN, formula = formula, data = data,
                   coords = coords, covar_FUN = covar_FUN,
                   distm_FUN = distm_FUN
  )

  # append fixed parameters, if they exist
  if(length(fixed) > 0){
    arg.list <- append(arg.list, list(fp = fixed))
  }

  # append arguments given by ... to the argument list
  arg.list = append(arg.list, list(...)) # add additional arguments to arg list

  # call optim, and pass arguments
  opt.out <- do.call(optim, args = arg.list)

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
#' }
fitGLS_opt_FUN <- function(op, fp, formula, data = NULL, coords, covar_FUN = "covar_exp", distm_FUN = "distm_scaled"){
  ## combine the optimized and fixed parameters into one vector
  all.pars <- if(missing(fp)){op}else{c(op, fp)}
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
  logLik = fitGLS(formula = formula, data = data, V = V, formula0 = NULL,
                  save.xx = FALSE, save.invchol = FALSE, logLik.only = TRUE, no.F = TRUE,
                  nugget = nug)
  return(-logLik)
}

