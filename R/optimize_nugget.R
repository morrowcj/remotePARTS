# optimize_nugget R function wrapper for .optimize_nugget_cpp

#' Find the maximum likelihood estimate of the nugget
#' @rdname optimize_nugget
#'
#' @param X numeric (double) nxp matrix
#' @param V numeric (double) nxn matrix
#' @param y numeric (double) nx1 column vector
#' @param lower lower boundary for nugget search
#' @param upper upper boundary for nugget search
#' @param tol desired accuracy of nugget search
#' @param debug logical: debug mode?
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
#' @seealso [stats::optimize()]
#'
#' @examples
#'
#' @export
optimize_nugget <- function(X, V, y, lower = 0, upper = 1, tol = 1e-5,
                            debug = FALSE){
  # coerce input to matrices
  X = as.matrix(X)
  V = as.matrix(V)
  y = as.matrix(y)

  ## error handling
  stopifnot(all(is.double(X), is.double(V), is.double(y)))
  stopifnot(all.equal(nrow(X), nrow(V), nrow(y)))
  stopifnot(all(check_posdef(V)))
  stopifnot(lower >= 0)
  stopifnot(upper <= 1)

  ## execute the function
  return(.Call(`_remotePARTS_optimize_nugget_cpp`, X, V, y, lower, upper,
               tol, debug))
}

#' fitNugget - R version
#' @rdname optimize_nugget
#' @details \code{fitNugget} is the R-only version of \code{optimize_nugget()}.
#' It uses \code{fitGLS_R()} and \code{stats::optimize()} to obtain the ML
#' nugget.
#'
#' @export
fitNugget <-  function(X, V, y, lower = 0, upper = 1, tol = .00001){
  int = c(lower, upper)
  N.opt <- optimize(f = function(nug){return(fitGLS_R(X, V, y, nugget = nug)$logLik)},
                    interval = int, tol = tol, maximum = TRUE)
  if(N.opt$maximum < tol){
    N0.LL <- fitGLS_R(X, V, y, nugget = 0)$logLik
    if(N0.LL > N.opt$objective){
      N.opt$maximum <- 0
    }
  }
  return(N.opt$maximum)
}

#' fitNugget_Rcpp - Rcpp version
#' @rdname optimize_nugget
#'
#' @details \code{fitNugget_Rcpp} is a hybrid between \code{optimize_nugget()}
#' and \code{fitNugget()}. It optimizes the C++ function \code{LogLikGLS()} with
#' the R function \code{stats::optimize()}.
#'
#' After further testing, only the most efficeint of the 3 nugget optimizers will
#' remain in the package.
#'
#' @export
fitNugget_Rcpp <-  function(X, V, y, lower = 0, upper = 1, tol = .00001){
  int = c(0,1)
  N.opt <- optimize(f = function(nug){return(LogLikGLS(nugget = nug, X, V, y))},
                    interval = int, tol = tol, maximum = TRUE)
  if(N.opt$maximum < tol){
    N0.LL <- LogLikGLS(nugget = 0, X, V, y)
    if(N0.LL > N.opt$objective){
      N.opt$maximum <- 0
    }
  }
  return(N.opt$maximum)
}

