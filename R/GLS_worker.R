# GLS worker wrapper

#' Worker function 1 for paritioned GLS
#'
#' @details this function is the first of 2 (maybe 3) worker functions that,
#' together, perform the partitioned GLS analysis.
#'
#' This function is simply a wrapper for fitGLS_cpp() that finds the MLE nugget
#' and adds it to the output.
#'
#' NOTE: eventually, the worker functions will perform the analysis using
#' multiple cores but that has not yet been implemented.
#'
#' @param y numeric vector
#' @param X numeric matrix
#' @param V numeric matrix
#' @param X0 numeric matrix
#' @param nug_l lower boundary for nugget optimization
#' @param nug_u upper boundary for nugget optimization
#' @param nug_tol tolerance of nugget optimization
#' @param save_xx logical: should xx, xx0, and invcholV be returned?
#'
#' @examples #TBA
#'
#' @export
GLS_worker <- function(y, X, V, X0, nug.l = 0, nug.u = 1, nug.tol = 1e-5,
                       save_xx = FALSE){
  ## coerce to matrices
  X = as(X, "matrix")
  V = as(V, "matrix")
  y = as(y, "matrix")
  X0 = as(X0, "matrix")

  ## error handling
  stopifnot(all(is.double(X), is.double(V), is.double(y), is.double(X0)))
  stopifnot(all.equal(nrow(X), nrow(V), nrow(X0), nrow(y)))
  stopifnot(all(check_posdef(V)))
  stopifnot(nug.l >= 0)
  stopifnot(nug.u <= 1)

  return(.Call(`_remoteSTAR_GLS_worker_cpp`, y, X, V, X0, nug.l, nug.u, nug.tol,
               save_xx))
}
