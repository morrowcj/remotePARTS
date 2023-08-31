#' @title Check if a matrix is positive definite
#'
#' @details check if a matrix is 1) square, 2) symmetric, and 3) positive
#' definite
#'
#' @param M numeric matrix
#'
#' @return returns a named logical vector with the following elements:
#'
#' \describe{
#'     \item{sqr}{logical: indicating whether \code{M} is square}
#'     \item{sym}{logical: indicating whether \code{M} is symmetric}
#'     \item{posdef}{logical: indicating whether \code{M} is positive-definitive}
#' }
#'
#' @examples
#'
#' # distance matrix
#' M = distm_scaled(expand.grid(x = 1:3, y = 1:3))
#'
#' # check if it is positive definitive
#' check_posdef(M)
#'
#' # check if the covariance matrix is positive definitive
#' check_posdef(covar_exp(M, .1))
#'
#' # non-symmetric matrix
#' check_posdef(matrix(1:9, 3, 3))
#'
#' # non-square matrix
#' check_posdef(matrix(1:6, 3, 2))
#'
#' @export
check_posdef <- function(M){
  res <- c(sqr = FALSE, sym = FALSE, posdef = FALSE)
  res["sqr"] = nrow(M) == ncol(M)
  if(res["sqr"]){
    res["sym"] = isSymmetric(M)
    if(res["sym"]){
      res["posdef"] = all(eigen(M, only.values = TRUE)$values > 1e-8)
    }
  }
  return(res)
}
