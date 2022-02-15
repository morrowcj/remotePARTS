
## C++ Version (preferred!) ----
#' Invert the cholesky decomposition of V
#' @rdname invert_chol
#'
#' @param M numeric (double), positive definite matrix
#' @param nugget numeric (double) nugget to add to M
#'
#' @return numeric matrix: inverse of the Cholesky decomposition (lower triangle)
#'
#' @details Calculates the inverse of the Cholesky decomposition of M which
#' should not be confused with the inverse of M *derived* from the
#' Cholesky decomposition (i.e. `chol2inv(M)`).
#'
#' @examples
#' M <- crossprod(matrix(1:6, 3))
#'
#' # without a nugget:
#' invert_chol(M)
#'
#' # with a nugget:
#' invert_chol(M, nugget = 0.2)
#' @export
invert_chol <- function(M, nugget = 0){

  # exception handling ----
  if(!is.matrix(M)){stop("M is not of class 'matrix'")}
  if(!is.double(M)){stop("M is not of type 'double'")}
  if(!all(check_posdef(M))){stop("M is not positive definite")}

  # execute the C++ function ----
  # return(.invchol_cpp(M, nugget))
  return(.Call(`_remotePARTS_invchol_cpp`, M, nugget))
}
