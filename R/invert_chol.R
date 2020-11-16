#' Invert the cholesky decomposition of V
#'
#' @param M numeric (double), positive definite matrix
#' @param nugget numeric (double) nugget to add to M
#'
#' @return numeric matrix: inverse of the cholesky decomposition
#'
#' @examples
#' M <- crossprod(matrix(1:6, 3))
#' # without a nugget:
#' invert_chol(M)
#' # with a nugget:
#' invert_chol(M, nugget = 0.2)
#'
#' @export
invert_chol <- function(M, nugget = 0){
  # Coerce inputs to proper format ----

  if(!is.matrix(M)){M = as.matrix(M)}
  if(!is.double(M)){M = as.double(M)}

  # exception handling ----

  ## need to check if M is positive and square or C++ will crash.
  stopifnot(all(check_posdef(M)))

  # execute the C++ function ----
  # return(.invchol_cpp(M, nugget))
  return(.Call(`_remoteSTAR_invchol_cpp`, M, nugget))
}
