# optimize_nugget R function wrapper for .optimize_nugget_cpp

#' Find the maximum likelihood estimate of the nugget
#'
#' @details this is the C++ version of `optimize()` which is specific to
#' finding the nugget value that maximizes the log-likelihood of `fitGLS_cpp()`
#'
#' This function is a translation from the forchan algorithm fmin into C++:
#' http://www.netlib.org/fmm/fmin.f
#'
#' Note: this function actually uses `LogLikGLS_cpp()` which should be swapped
#' for `fitGLS_cpp()` once the correct functionality is added to the latter.
#'
#' @param X numeric (double) nxp matrix
#' @param V numeric (double) nxn matrix
#' @param y numeric (double) nx1 column vector
#' @param lower lower boundary for nugget search
#' @param upper upper boundary for nugget search
#' @param tol desired accuracy of nugget search
#' @param debug logical: debug mode?
#'
#' @examples #TBA
#'
#' @export
optimize_nugget <- function(X, V, y, lower = 0, upper = 1, tol = 1e-5,
                            debug = FALSE){
  # coerce input to matrices
  X = as(X, "matrix")
  V = as(V, "matrix")
  y = as(y, "matrix")

  ## error handling
  stopifnot(all(is.double(X), is.double(V), is.double(y)))
  stopifnot(all.equal(nrow(X), nrow(V), nrow(y)))
  stopifnot(all(check_posdef(V)))
  stopifnot(lower >= 0)
  stopifnot(upper <= 1)

  ## execute the function
  return(.Call(`_remoteSTAR_optimize_nugget_cpp`, X, V, y, lower, upper,
               tol, debug))
}




