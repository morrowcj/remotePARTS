# GLS worker function ----
#' Worker function 1 for paritioned GLS
#' @rdname GLS_worker
#' @family remoteGLS
#'
#' @param y numeric vector
#' @param X numeric matrix
#' @param V numeric matrix
#' @param X0 numeric matrix
#' @param nug_l lower boundary for nugget optimization
#' @param nug_u upper boundary for nugget optimization
#' @param nug_tol tolerance of nugget optimization
#' @param threads
#' @param save_xx logical: should xx, xx0, and invcholV be returned?
#'
#' @details \code{GLS_worker()} is meant to be called by other functions
#' and performs GLS similar to \code{fitGLS()}. However, unlike \code{fitGLS()},
#' the nugget is optimized each time the function is called. This additional
#' step is ideal for partitioned GLS.
#'
#' At present, \code{GLS_worker()} is only single-core natively but will have
#' multi-core functionality eventually.
#'
#' @return a list of remoteGLS output
#'
#' @seealso [fitGLS()] and [fitGLS.parition_rcpp()]
#'
#' @examples
#'
#' @export
GLS_worker <- function(y, X, V, X0, nug_l = 0, nug_u = 1, nug_tol = 1e-5,
                       save_xx = FALSE, threads = 1){
  ## coerce to matrices
  X = as.matrix(X)
  V = as.matrix(V)
  y = as.matrix(y)
  X0 = as.matrix(X0)

  ## error handling
  stopifnot(all(is.double(X), is.double(V), is.double(y), is.double(X0)))
  stopifnot(all.equal(nrow(X), nrow(V), nrow(X0), nrow(y)))
  stopifnot(all(check_posdef(V)))
  stopifnot(nug_l >= 0)
  stopifnot(nug_u <= 1)

  ## Execute C++ GLS worker function
  out <- .Call(`_remotePARTS_GLS_worker_cpp`, y, X, V, X0, nug_l, nug_u, nug_tol,
               save_xx, threads)

  ## calculate p values outside of C++
  out$pval.t <- sapply(out$tstat, function(tval){ #t test
    2*pt(abs(tval), df = out$dft, lower.tail = FALSE)
  })
  out$pval.F <- stats::pf(out$Fstat, out$df.F[1], out$df.F[2], lower.tail = FALSE) #F test
  class(out) <- append("remoteGLS", class(out))
  attr(out, "no_F") = FALSE

  return(out)
}

## cross-partition worker ----
#' Worker function 2 for partitioned GLS
#' @rdname crosspart_worker
#' @family remoteGLS
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
#' @return a list of cross-partition statistics
#'
#' @details Cross-partition statistics are calculated for a pair of partitions
#' i and j.
#'
#' @examples
#'
#' @export
crosspart_worker <- function(xxi, xxj, xxi0, xxj0, invChol_i, invChol_j, Vsub,
                             nug_i, nug_j, df1, df2){
  # coerce input to matrices
  xxi = as.matrix(xxi)
  xxj = as.matrix(xxj)
  xxi0 = as.matrix(xxi0)
  xxj0 = as.matrix(xxj0)
  invChol_i = as.matrix(invChol_i)
  invChol_j = as.matrix(invChol_j)
  Vsub = as.matrix(Vsub)

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


#' R version of cross-part worker
#' @rdname crosspart_worker
#'
#' @details \code{crosspart_worker_R} is not exported by default and is just
#' a reference function. is a purely R version and will be removed in future
#' implementations.
crosspart_worker_R <- function(xxi, xxj, xxi0, xxj0, invChol_i, invChol_j,
                               Vsub, df1, df2){
  np = nrow(xxi)
  # rescale nuggets
  # nug_i = ifelse(nug_i == 0, 0, (1 - nug_i)/nug_i)
  # nug_j = ifelse(nug_i == 0, 0, (1 - nug_i)/nug_i)

  # variance matrix with nuggets (ARE NEVER USED)
  # Vn <- diag(rep(c(nug_i, nug_j), each = np)) + Vij
  # Vn <- Vn[1:np, (np+1):(2*np)] # upper right block

  # calculate stats
  Rij <- crossprod(t(invChol_i), tcrossprod(Vsub, invChol_j))

  Hi <- xxi %*% solve(crossprod(xxi)) %*% t(xxi)
  Hj <- xxj %*% solve(crossprod(xxj)) %*% t(xxj)

  Hi0 <- xxi0 %*% solve(crossprod(xxi0)) %*% t(xxi0)
  Hj0 <- xxj0 %*% solve(crossprod(xxj0)) %*% t(xxj0)

  SiR <- Hi - Hi0
  SjR <- Hj - Hj0

  SiE <- diag(np) - Hi
  SjE <- diag(np) - Hj

  # rSSRij <- (SiR %*% (Rij %*% SjR %*% t(Rij)))/df1
  # rSSEij <- (SiE %*% (Rij %*% SjE %*% t(Rij)))/df2

  rSSRij <- matrix(SiR, nrow=1) %*%
    matrix(Rij %*% SjR %*% t(Rij), ncol=1)/df1

  rSSEij <- matrix(SiE, nrow=1) %*%
    matrix(Rij %*% SjE %*% t(Rij), ncol=1)/df2


  # output
  out_lst <- list("Rij" = Rij,
                  "Hi" = Hi,
                  "Hj" = Hj,
                  "Hi0" = Hi0,
                  "Hj0" = Hj0,
                  "SiR" = SiR,
                  "SjR" = SjR,
                  "rSSRij" = rSSRij,
                  "rSSEij" = rSSEij)
  return(out_lst)
}
