# crosspartition GLS worker

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
               invChol_i, invChol_j, Vsub, df1, df2))
}
