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
  cross_stats = c("rSSE" = x$overall$rSSE)
  if(!is.na(x$overall$rSSR) | !is.null(x$overall$rSSR)){
    cross_stats["rSSR"] = x$overall$rSSR
  }
  if(!is.na(x$overall$Fstat) | !is.null(x$overall$Fstat)){
    cross_stats["Fstat"] = x$overall$Fstat
  }
  print(c(rSSR = x$overall$rSSR, rSSE = x$overall$rSSE, Fstat = x$overall$Fstat), ...)

  if(!is.null(x$overall$t.test)){
    cat("\nT-test results:\n")
    print(x$overall$t.test)
  }
}

## Correlated chisqr test ----
#' @title Chisqr test for partitioned GLS
#'
#' @param Fmean mean value of F-statistic from correlated F-tests
#' @param rSSR correlation among partition regression sum of squares
#' @param df1 first degree of freedom for F-tests
#' @param npart number of partitions
#'
#' @return a p-value for the correlated chisqr test
#'
#' @examples
#' remotePARTS:::part_chisqr(Fmean = 3.6, rSSR = .021, df1 = 2, npart = 5)
part_chisqr <- function(Fmean, rSSR, df1, npart){
  checks = c(is.na(Fmean) | is.infinite(Fmean),
             is.na(rSSR) | is.infinite(rSSR),
             is.na(npart) | is.infinite(npart))
  if(any(checks)){
    stop("NA values supplied to part_chisqr")
  }
  rZ <- (rSSR/df1)^0.5
  v.MSR <- diag(df1) - rZ
  V.MSR <- kronecker(diag(npart),v.MSR) + rZ
  D.MSR <- t(chol(V.MSR, pivot = TRUE))
  if (attr(D.MSR, "rank") < npart*df1) {
    rank.deficient.MSR <- npart*df1-attr(D.MSR, "rank")
    v.MSR <- diag(df1) - (0.99/df1)
    V.MSR <- kronecker(diag(npart), v.MSR) + (0.99/df1)
  } else {
    rank.deficient.MSR <- 0
  }
  lambda <- eigen(V.MSR)$values
  pvalue <- suppressWarnings(CompQuadForm::imhof(q = npart * df1 * Fmean, lambda = lambda)$Qq)
  pvalue = ifelse(pvalue <= 1e-06, 1e-06, pvalue) # prevent from being negative/too low
  out = c("pval.chisqr" = pvalue)
  attr(out, 'rankdef.MSR') <- rank.deficient.MSR
  return(out)
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
#' @param part.covar_coef an array of covar_coef from each partition
#' @param rcoefficients an rcoefficeints array, one for each partition
#' @param df2 second degree of freedom from partitioned GLS
#' @param npart number of partitions
#'
#' @return a list whose first element is a coefficient table with estimates,
#' standard errors, t-statistics, and p-values and whose second element is a
#' matrix of correlations among coefficients.
part_ttest <- function(coefs, part.covar_coef, rcoefficients, df2, npart){
  secoef <- matrix(NA, length(coefs), 1)
  for(i in seq_len(length(coefs))){
    R <- (1 - rcoefficients[i,i]) * diag(npart) + rcoefficients[i,i] * matrix(1, npart, npart)
    part.SEs <- part.covar_coef[i,i,]^.5
    secoef[i, ] <- (part.SEs %*% R %*% part.SEs)^.5/npart
  }

  covar_coef <- matrix(NA, length(coefs), length(coefs))
  for(i in seq_len(length(coefs))) for(j in seq_len(length(coefs))){
    covar_coef[i, j] <- (sum(part.covar_coef[i, j, ]) +
            rcoefficients[i, j] * (npart - 1) * abs(sum(part.covar_coef[i, j, ])))/npart^2
  }
  rownames(covar_coef) <- names(coefs)
  colnames(covar_coef) <- names(coefs)

  # secoef <- diag(covar_coef)^.5
  tscore <- coefs/secoef
  pvalue <- 2 * pt(abs(tscore), df=df2 * npart, lower.tail = FALSE)

  ttest <- cbind(coefs, secoef, tscore, pvalue)
  colnames(ttest) <- c("Est", "SE", "t.stat", "pval.t")
  return(list(p.t = ttest, covar_coef = covar_coef))
}

#' @title Conduct a t-test of "partGLS" object
#'
#' @description Conduct a correlated t-test of a partitioned GLS
#'
#' @param x "partGLS" object
#' @param ... additional arguments passed to print
#'
#' @return a list whose first element is a coefficient table with estimates,
#' standard errors, t-statistics, and p-values and whose second element is a
#' matrix of correlations among coefficients.
#'
#' @method t.test partGLS
#' @export
t.test.partGLS <- function(x, ...){
  part_ttest(coefs = x$overall$coefficients,
                  part.covar_coef = x$part$covar_coef,
                  rcoefficients = x$overall$rcoefficients,
                  df2 = x$overall$dfs[2],
                  npart = x$overall$partdims["npart"])
}
