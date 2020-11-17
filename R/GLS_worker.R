# GLS worker wrapper

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
