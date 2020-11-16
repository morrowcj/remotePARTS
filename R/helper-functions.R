#' check if matrix is positive definitive
#'
#' @param M numeric matrix
check_posdef <- function(M){

  res <- c(sqr = FALSE, sym = FALSE, posdef = FALSE)

  res["sqr"] = nrow(M) == ncol(M)
  # if(!res["sqr"]){
  #   stop("M not square")
  # } # not square

  if(res["sqr"]){

    res["sym"] = isSymmetric(M)
    # if(!res["sym"]){
    #   stop("M not symmetric")
    # } # not symetric

    if(res["sym"]){

      res["posdef"] = all(eigen(M, only.values = TRUE)$values > 1e-8)
      # if(!res["posdef"]){
      #   stop("M not positive definite")
      # }
    }
  }
  return(res)
}

