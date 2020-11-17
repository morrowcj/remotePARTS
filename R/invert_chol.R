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

  # exception handling ----
  if(!is.matrix(M)){stop("M is not of class 'matrix'")}
  if(!is.double(M)){stop("M is not of type 'double'")}
  if(!all(check_posdef(M))){stop("M is not positive definite")}

  # execute the C++ function ----
  # return(.invchol_cpp(M, nugget))
  return(.Call(`_remoteSTAR_invchol_cpp`, M, nugget))
}
