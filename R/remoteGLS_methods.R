## remoteGLS methods

## GLS Constructor ----
#' remoteGLS constructor (S3)
#'
#' @param formula optional argument specifying the GLS formula
#' @param formula0 optional argument specifying the null GLS formula
#' @param no.F optional argument specifying the no.F attribute
#' @return an empty S3 object of class "remoteGLS"
#'
#' @examples
#' # tmp <- remoteGLS() #empty remoteGLS object
#'
remoteGLS <- function(formula, formula0, no.F = FALSE){
  # create empty list
  elements = c("call", # model info
               "coefficients", "SSE", "MSE", "SE",
               "df_t", "logDetV",
               "tstat", "pval_t", "logLik", "nugget",
               "coefficients0", "SSE0", "MSE0", "SE0",
               "MSR", "df0", "LL0", "df_F",
               "Fstat","pval_F",
               "xx", "xx0", "invcholV", "formula", "formula0")
  GLS.obj = vector("list", length(elements))
  names(GLS.obj) = elements

  GLS.obj$call = match.call()

    if(!missing(formula)){
      GLS.obj$formula = deparse(as.formula(formula))
    }

    if(!missing(formula0)){
      GLS.obj$formula0 = deparse(as.formula(formula0))
    }

  class(GLS.obj) <- append("remoteGLS", class(GLS.obj))
  attr(GLS.obj, "no.F") = no.F

  return(GLS.obj) # return the empty list
}

## GLS print method ----
#' print method for remoteGLS
#'
#' @param x remoteGLS object
#' @param digits digits to print
#' @param ... additional arguments
#'
#' @return formatted output for remoteGLS object
#'
#' @method print remoteGLS
#' @export
print.remoteGLS <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  ## Coefficient table
  coefs = data.frame("Est" = x$coefficients, "t stat" = x$tstat, "pval t" = x$pval_t)

  ## Model stats
  if(!attr(x, "no.F")){
    mod.stats = data.frame("model" = c(x$formula, x$formula0),
                           "df_F" = x$df_F,
                           "SSE" = c(x$SSE, x$SSE0),
                           "MSE" = c(x$MSE, x$MSE0),
                           "logLik" = c(x$logLik, x$LL0),
                           "Fstat" = c(x$Fstat, NA),
                           "pval_F" = c(x$pval_F, NA))
  } else {
    mod.stats = data.frame("SSE" = x$SSE,
                           "MSE" = x$MSE,
                           "logLik" = x$logLik)
  }
  cat("t-tests:\n")
  print(coefs, digits = digits)

  if(!attr(x, "no.F")){
    cat("\nF-test:\n")
  } else {
    cat("\nModel statistics:\n")
  }
  print(mod.stats, digits = digits)
}
