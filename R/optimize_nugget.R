# optimize_nugget R function wrapper for .optimize_nugget_cpp

#' Find the maximum likelihood estimate of the nugget
#' @rdname optimize_nugget
#'
#' @param X numeric (double) nxp matrix
#' @param y numeric (double) nx1 column vector
#' @param V numeric (double) nxn matrix
#' @param lower lower boundary for nugget search
#' @param upper upper boundary for nugget search
#' @param tol desired accuracy of nugget search
#' @param debug logical: debug mode?
#' @param ncores an optional integer indicating how many CPU threads to use for
#' matrix calculations.
#'
#' @return maximum likelihood nugget estimate
#'
#' @details Finds the maximum likelihood nugget estimate via mathematical
#' optimization.
#'
#' To maximize efficiency, \code{optimize_nugget()} is implemented entirely
#' in C++. Optimization takes place via a C++ version of the \code{fmin} routine
#' (Forsythe et al 1977). Translated from http://www.netlib.org/fmm/fmin.f
#'
#' The function \code{LogLikGLS()} is optimized for \code{nugget}. Once the
#' \code{LogLikGLS()} functionality is absorbed by \code{fitGLS()}, it will
#' be used instead.
#'
#' @seealso \code{?stats::optimize()}
#'
#' @examples
#' \dontrun{
#' ## read data
#' data(ndvi_AK10000)
#' df = ndvi_AK10000[seq_len(200), ] # first 200 rows
#'
#' ## format data
#' X = stats::model.matrix(CLS_coef ~ 0 + land, data = df)
#'
#' ## fit covariance matrix
#' V = covar_exp(distm_scaled(cbind(df$lng, df$lat)), range = .01)
#'
#' ## find the ML nugget
#' remotePARTS:::optimize_nugget(X = X, V = V, y = df$CLS_coef, debug = TRUE)
#' }
optimize_nugget <- function(X, y, V, lower = 0.001, upper = 0.999,
                            tol = .Machine$double.eps^.25, debug = FALSE,
                            ncores = NA) {
  if(is.na(ncores)){
    ncores = 0L
  } else {
    ncores = as.integer(ncores)
  }

  # # coerce input to matrices
  X = as.matrix(X)
  X0 = diag(1)
  y = as.matrix(y)
  stopifnot(ncol(y) == 1)
  V = as.matrix(V)

  # checks
  ## check positive definitive
  if (!all(check_posdef(V))) {
    stop("V is not positive definitive")
  }
  ## check for correct dimensions
  if (!all.equal(length(y), nrow(X), nrow(V), ncol(V))) {
    stop("Input dimension mismatch")
  }
  ## check that all variables are numeric
  if (!all(is.double(y), is.double(X), is.double(V))) {
    stop("All inputs must be numeric (double precision)")
  }
  ## boundaries handling
  stopifnot(lower >= 0)
  stopifnot(upper <= 1)

  .Call(`_remotePARTS_optimize_nugget_cpp`, X, X0, V, y, lower, upper, tol,
        diag(1), FALSE, debug, ncores)
}
