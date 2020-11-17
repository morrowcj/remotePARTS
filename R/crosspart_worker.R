# crosspartition GLS worker

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
#' @param df1 first degree of freedom
#' @param df2 second degree of freedom
#'
#' @examples #TBA
#'
#' @export
crosspart_worker <- function(xxi, xxj, xxi0, xxj0, invChol_i, invChol_j, Vsub,
                             df1, df2){
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

  return(.Call(`_remoteSTAR_crosspart_worker_cpp`, xxi, xxj, xxi0, xxj0,
               invChol_i, invChol_j, Vsub,  df1, df2))
}
