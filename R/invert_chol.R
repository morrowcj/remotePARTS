#' Invert the cholesky decomposition of V
#'
#' @param M numeric matrix
#' @param nugget numeric nugget to add to V
#' @param debug logical: enter debug mode?
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
invert_chol <- function(M, nugget = 0, debug = FALSE){
  # Coerce inputs to proper format ----

  if(!is.matrix(M)){M <- as.matrix(M)}
  if(!is.numeric(M)){M <- as.matrix(M)}

  # exception handling ----

  ## need to check if M is positive definitive or C++ will crash.
  if(nrow(M) != ncol(M)){
    stop("M not square")
  } # not square
  if(any(M[upper.tri(M)] != t(M[lower.tri(M)]))){
    stop("M not symmetric")
  } # not symetric
  tmp <- eigen(M, only.values = TRUE)$values
  if(any(eigen(M, only.values = TRUE)$values < 1e-8)){
    stop("M not positive definite")
  }

  # execute the C++ function ----
  return(.invchol_cpp(M, nugget))
}
