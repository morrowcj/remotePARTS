#' Fit GLS to remote sensing data
#'
#' @details see `fitGLS()`
#'
#' @rdname fitGLS
#'
#' @param X numeric (double) matrix
#' @param V numeric (double) matrix
#' @param y numeric (double) column vector
#' @param X0 numeric (double) matrix
#' @param nugget numeric nugget to add to V
#' @param save_xx logical: should xx, xx0, and invcholV be returned? This
#' functionality is meant for use with the partitioned GLS whereby these
#' values are used to calculate cross-partition statistics.
#' @param threads integer indicating the number of threads to use. This current
#' version does not have multi-thread functionality so this argument does
#' nothing yet.
#'
#' @examples #TBA
#'
#' @export
fitGLS <- function(X, V, y, X0, nugget = 0, save_xx = FALSE, threads = 1){


  ## coerce to matrices
  X = as(X, "matrix")
  V = as(V, "matrix")
  y = as(y, "matrix")
  X0 = as(X0, "matrix")

  ## error handling
  stopifnot(all(is.double(X), is.double(V), is.double(y), is.double(X0)))
  stopifnot(all.equal(nrow(X), nrow(V), nrow(X0), nrow(y)))
  stopifnot(all(check_posdef(V)))


  ## call the c++ function
  # return(.fitGLS_cpp) # contains additional function call
  out <- .Call(`_remoteSTAR_fitGLS_cpp`, X, V, y, X0, nugget, save_xx, threads)
  ## add in p values
  out$pval.t <- 2 * pt(abs(out$tstat), df = out$dft, lower.tail = F)
  out$pval.F <- pf(out$Fstat, df1 = out$df.F[1], df2 = out$df.F[2], lower.tail = F)
  return(out)
}

#' Caculate log-liklihood of GLS model
#'
#' @details this function is mostly meant to optimize the nugget for a paritular
#' set of data.
#'
#' Note: this function should be deprecated and simply added as functionality
#' to `.fitGLS_cpp()`.
#'
#' @rdname fitGLS
#'
#' @param nugget the nugget to add to V
#' @param X numeric matrix
#' @param V numeric matrix
#' @param y numeric vector
#'
#' @examples #TBA
LogLikGLS <- function(nugget, X, V, y){
  ## coerce to matrices
  X = as(X, "matrix")
  V = as(V, "matrix")
  y = as(y, "matrix")

  ## error handling
  stopifnot(all(is.double(X), is.double(V), is.double(y)))
  stopifnot(all.equal(nrow(X), nrow(V), nrow(y)))
  stopifnot(all(check_posdef(V)))

  return(.Call(`_remoteSTAR_LogLikGLS_cpp`, nugget, X, V, y))
}
