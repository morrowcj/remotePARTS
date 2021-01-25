# GLS worker function ----
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

  ## Execute C++ GLS worker function
  out <- .Call(`_remotePARTS_GLS_worker_cpp`, y, X, V, X0, nug.l, nug.u, nug.tol,
               save_xx)

  ## calculate p values outside of C++
  out$pval.t <- sapply(out$tstat, function(tval){ #t test
    2*pt(abs(tval), df = out$dft, lower.tail = FALSE)
  })
  out$pval.F <- pf(out$Fstat, out$df.F[1], out$df.F[2], lower.tail = FALSE) #F test

  return(out)
}

## cross-partition worker ----
#' Worker function 2 for partitioned GLS
#'
#' @details this is the second worker function for the partitioned GLS analysis.
#'
#' NOTE: currently, there is no parallel functionality and the partitioned
#' form of the GLS is not implemented entirely in C++. Instead, the R function
#' fitGLS.partition_rcpp() weaves between R and C++ on a single core. While
#' this method is still much faster than the purely R implementation, migration
#' to entirely C++ will greatly improve speed further. This migration requires
#' calculating geographic distances with C++ which I've not yet written.
#'
#' Additionally, there seems to be a memory-related issue with this code. I've
#' successfully used this function when partitions have 100 or fewer rows (too
#' small). However, larger partitions cause a fatal error that causes a crash.
#'
#' @param xxi numeric matrix xx from  partition i
#' @param xxj numeric matrix xx from  partition j
#' @param xxi0 numeric matrix xx0 from  partition i
#' @param xxj0 numeric matrix xx0 from  partition j
#' @param invChol_i numeric matrix invcholV from  partition i
#' @param invChol_j numeric matrix invcholV from  partition j
#' @param Vsub numeric variance matrix for Xij (upper block)
#' @param nug_i nugget from partition i
#' @param nug_j nugget from partition j
#' @param df1 first degree of freedom
#' @param df2 second degree of freedom
#'
#' @examples #TBA
#'
#' @export
crosspart_worker <- function(xxi, xxj, xxi0, xxj0, invChol_i, invChol_j, Vsub,
                             nug_i, nug_j, df1, df2){
  # coerce input to matrices
  xxi = as(xxi, "matrix")
  xxj = as(xxj, "matrix")
  xxi0 = as(xxi0, "matrix")
  xxj0 = as(xxj0, "matrix")
  invChol_i = as(invChol_i, "matrix")
  invChol_j = as(invChol_j, "matrix")
  Vsub = as(Vsub, "matrix")

  ## error handling
  stopifnot(all(is.double(xxi), is.double(xxj),
                is.double(xxi0), is.double(xxj0),
                is.double(invChol_i), is.double(invChol_j),
                is.double(Vsub)))
  stopifnot(all.equal(nrow(xxi), nrow(xxj), nrow(xxi0), nrow(xxj0)))
  stopifnot(all.equal(ncol(xxi), ncol(xxj)))
  stopifnot(all.equal(ncol(xxi0), ncol(xxj0)))
  # stopifnot(all(check_posdef(V)))

  return(.Call(`_remotePARTS_crosspart_worker_cpp`, xxi, xxj, xxi0, xxj0,
               invChol_i, invChol_j, Vsub, nug_i, nug_j,  df1, df2))
}
