# partGLmethods

## Print method ----
#' @title S3 print method for "partGLS" objects
#'
#' @param x "partGLS" object
#' @param ... additional arguments passed to print
#'
#' @method print partGLS
#' @export
print.partGLS <- function(x, ...){
  cat("\nCoefficients:\n")
  print(x$overall$coefficients, ...)
  cat("\nCross coefficients:\n")
  print(x$overall$rcoefficients, ...)
  cat("\nCross-partition statistics:\n")
  print(c(rSSR = x$overall$rSSR, rSSE = x$overall$rSSE, Fstat = x$overall$Fstat), ...)
}

## Correlated chisqr test ----
#' @title Chisqr test for partitioned GLS
#'
#' @param Fmean mean value of F-statistic from correlated F-tests
#' @param rSSR
#' @param df1 first degree of freedom for F-tests
#' @param npart number of partitions
#'
#' @return a p-value for the correlated chisqr test
#'
#' @examples
#' remotePARTS:::part_chisqr(Fmean = 3.6, rSSR = .021, df1 = 2, npart = 5)
part_chisqr <- function(Fmean, rSSR, df1, npart){
  rZ <- rSSR^.5/df1
  v.MSR <- diag(df1) - rZ
  V.MSR <- kronecker(diag(npart),v.MSR) + rZ
  lambda <- eigen(V.MSR)$values
  pvalue <- suppressWarnings(CompQuadForm::imhof(q = npart * df1 * Fmean, lambda = lambda)$Qq)
  pvalue = ifelse(pvalue <= 1e-06, 1e-06, pvalue) # prevent from being negative/too low
  return(c("pval.chisqr" = pvalue))
}

#' Conduct a chi-squared test
#'
#' @description generic S3 method for a chi-squared test
#'
#' @param x object on which to conduct the test
#' @param ... additional arguments
#'
#' @export
chisqr <- function(x, ...){
  UseMethod("chisqr")
}

#' @title Conduct a chisqr test of "partGLS" object
#'
#' @description Conduct a correlated chi-square test on a partitioned GLS
#'
#' @param x "remoteGLS" object
#' @param ... additional arguments passed to print
#'
#' @return a p-value for the correlated chisqr test
#'
#' @method chisqr partGLS
#' @export
chisqr.partGLS <- function(x, ...){
 part_chisqr(Fmean = x$overall$Fstat, rSSR = x$overall$rSSR,
             df1 = x$overall$dfs[1],
             npart = x$overall$partdims["npart"])
}

## correlated t-test
#' @title Correlated t-test for paritioned GLS
#' @param coefs vector average GLS coefficients
#' @param part.SEs matrix of partition SEs for each coefficient (columns)
#' @param rcoef
#' @param df2 second degree of freedom from partitioned GLS
#' @param npart number of partitions
#'
#' @return coefficient table with estimates, standard errors, t-statistics, and p-values
part_ttest <- function(coefs, part.SEs, rcoef, df2, npart){

  secoef <- matrix(NA, length(coefs), 1)
  for(i in seq_len(length(coefs))){
    R <- (1 - rcoef[i]) * diag(npart) + rcoef[i] * matrix(1, npart, npart)
    secoef[i, ] <- (part.SEs[,i] %*% R %*% part.SEs[, i])^.5/npart
  }

  tscore <- coefs/secoef
  pvalue <- 2 * pt(abs(tscore), df=df2 * npart, lower.tail = FALSE)

  ttest <- cbind(coefs, secoef, tscore, pvalue)
  colnames(ttest) <- c("Est", "SE", "t.stat", "pval.t")

  return(p.t = ttest)
}

#' @title Conduct a t-test of "partGLS" object
#'
#' @description Conduct a correlated t-test of a partitioned GLS
#'
#' @param x "partGLS" object
#' @param ... additional arguments passed to print
#'
#' @return a coefficient table with estimates, standard errors, t-statistics, and p-values
#'
#' @method t.test partGLS
#' @export
t.test.partGLS <- function(x, ...){
  part_ttest(coefs = x$overall$coefficients, part.SEs = x$part$SEs,
             rcoef = x$overall$rcoefficients, df2 = x$overall$dfs[2],
             npart = x$overall$partdims["npart"])
}
