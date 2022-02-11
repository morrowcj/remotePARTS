# This document contains the class methods for the remotePARTS classes
## remoteTS and mapTS

## print methods ----
#' @title S3 print method for remoteTS class
#'
#' @rdname TS.methods
#'
#' @param x remoteTS object
#' @param digits significant digits to show
#' @param signif.stars logical, passed to \code{stats::printCoefmat}
#' @param ..., additional parameters passed to \code{stats::printCoefmat}
#'
#'
#' @return returns formatted output
#'
#' @examples
#' # simulate dummy data
#'  time.points = 9 # time series length
#'  map.width = 5 # square map width
#'  coords = expand.grid(x = 1:map.width, y = 1:map.width) # coordinate matrix
#'  ## create empty spatiotemporal variables:
#'  X <- matrix(NA, nrow = nrow(coords), ncol = time.points) # response
#'  Z <- matrix(NA, nrow = nrow(coords), ncol = time.points) # predictor
#'  # setup first time point:
#'  Z[, 1] <- .05*coords[,"x"] + .2*coords[,"y"]
#'  X[, 1] <- .5*Z[, 1] + rnorm(nrow(coords), 0, .05) #x at time t
#'  ## project through time:
#'  for(t in 2:time.points){
#'    Z[, t] <- Z[, t-1] + rnorm(map.width^2)
#'    X[, t] <- .2*X[, t-1] + .1*Z[, t] + .05*t + rnorm(nrow(coords), 0 , .25)
#'  }
#'
#'  ## Pixel CLS
#'  tmp.df = data.frame(x = X[1, ], t = nrow(X), z = Z[1, ])
#'  CLS <- fitCLS(x ~ z, data = tmp.df)
#'  print(CLS)
#'  summary(CLS)
#'  residuals(CLS)
#'  coef(CLS)
#'  logLik(CLS)
#'  fitted(CLS)
#'  # plot(CLS) # doesn't work
#'
#'  ## Pixel AR
#'  AR <- fitAR(x ~ z, data = tmp.df)
#'  print(AR)
#'  summary(AR)
#'  coef(AR)
#'  residuals(AR)
#'  logLik(AR)
#'  fitted(AR)
#'  # plot(AR) # doesn't work
#'
#'  ## Map CLS
#'  CLS.map <- fitCLS_map(X, coords, y ~ Z, X.list = list(Z = Z), lag.x = 0, resids.only = TRUE)
#'  print(CLS.map)
#'  summary(CLS.map)
#'  residuals(CLS.map)
#'  # plot(CLS.map)# doesn't work
#'
#'  CLS.map <- fitCLS_map(X, coords, y ~ Z, X.list = list(Z = Z), lag.x = 0, resids.only = FALSE)
#'  print(CLS.map)
#'  summary(CLS.map)
#'  coef(CLS.map)
#'  residuals(CLS.map)
#'  # logLik(CLS.map) # doesn't work
#'  fitted(CLS.map)
#'  # plot(CLS.map) # doesn't work
#'
#'  ## Map AR
#'  AR.map <- fitAR_map(X, coords, y ~ Z, X.list = list(Z = Z), resids.only = TRUE)
#'  print(AR.map)
#'  summary(AR.map)
#'  residuals(AR.map)
#'  # plot(AR.map)# doesn't work
#'
#'  AR.map <- fitAR_map(X, coords, y ~ Z, X.list = list(Z = Z), resids.only = FALSE)
#'  print(AR.map)
#'  summary(AR.map)
#'  coef(AR.map)
#'  residuals(AR.map)
#'  # logLik(AR.map) # doesn't work
#'  fitted(AR.map)
#'  # plot(AR.map) # doesn't work
#'
#' @method print remoteTS
#' @export
print.remoteTS <- function(x, digits = max(3L, getOption("digits") - 3L),
                           signif.stars = getOption("show.signif.stars"), ...){
  # Function call
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  # Coefficient table
  cat("Coefficients:\n")
  coef.tab <- cbind("Estimate" = x$coefficients, "Std. Error" = x$SE,
                    "t value" = x$tstat, "Pr(>|t|)" = x$pval)
  stats::printCoefmat(coef.tab, digits = digits, signif.stars = signif.stars,
                      na.print = "NA", has.Pvalue = TRUE, P.values = TRUE, ...)
}

#' @title S3 summary method for remoteTS class
#'
#' @rdname TS.methods
#'
#' @param object remoteTS object
#' @param digits significant digits to show
#' @param signif.stars logical, passed to \code{stats::printCoefmat}
#' @param ..., additional parameters passed to \code{stats::printCoefmat}
#'
#' @return returns formatted output, including summary stats
#'
#' @method summary remoteTS
#' @export
summary.remoteTS <- function(object, digits = max(3L, getOption("digits") - 3L),
                             signif.stars = getOption("show.signif.stars"), ...){
  # Function Call
  cat("\nCall:\n", paste(deparse(object$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  # Residuals
  cat("Residuals:\n")
  if (object$df.residual > 5L) {
    names <- c("Min", "1Q", "Median", "3Q", "Max")
    resid.q <- if (length(dim(object$residuals)) == 2L)
      structure(apply(t(object$residuals), 1L, quantile),
                dimnames = list(names, dimnames(object$residuals)[[2L]]))
    else {
      zz <- zapsmall(quantile(object$residuals), digits + 1L)
      structure(zz, names = names)
    }
    print(resid.q, digits = digits, ...)
  }
  else if (object$df.residual > 0L) {
    print(object$residuals, digits = digits, ...)
  }
  else {
    cat("ALL", object$rank, "residuals are 0: no residual degrees of freedom!")
    cat("\n")
  }

  # Coefficient table
  cat("\nCoefficients:\n")
  coef.tab <- cbind("Estimate" = object$coefficients, "Std. Error" = object$SE,
                    "t value" = object$tstat, "Pr(>|t|)" = object$pval)
  stats::printCoefmat(coef.tab, digits = digits, signif.stars = signif.stars,
                      na.print = "NA", has.Pvalue = TRUE, P.values = TRUE, ...)

  # MSE
  cat("\nMean squared error:", round(object$MSE, digits))
  cat("\nLog-likelihood:", round(object$logLik, digits))
}

#' @title S3 print method for mapTS class
#'
#' @rdname TS.methods
#'
#' @param x mapTS object
#' @param digits significant digits to show
#' @param ..., additional parameters passed to further print methods
#'
#' @return returns formatted output
#'
#' @method print mapTS
#' @export
print.mapTS <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
  # Function Call
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  # Output
  if(attr(x, "resids.only")){
    cat("Time series residuals:\n")
    print(x$residuals, digits)
  } else {
    cat("Coefficients:\n")
    print(x$coefficients, digits)

    cat("\nTime series residuals:\n")
    print(x$residuals, digits)
  }
}

#' @title S3 summary method for mapTS class
#'
#' @rdname TS.methods
#'
#' @param object mapTS object
#' @param digits significant digits to show
#' @param CL confidence level (default = .95)
#' @param na.rm logical, should observations with NA be removed?
#' @param ..., additional parameters passed to further print methods
#'
#' @return returns formatted summary stats
#'
#' @method summary mapTS
#' @export
summary.mapTS <- function(object, digits = max(3L, getOption("digits") - 3L), CL = .95,
                          na.rm = TRUE, ...){

  # Function Call
  cat("\nCall:\n", paste(deparse(object$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  if(attr(object, "resids.only")){
    cat("Correlation among time series residuals:\n")
    tmp = summary(as.vector(cor(t(object$residuals))))
    names(tmp) <- c("Min", "1Q", "Median", "Mean", "3Q", "Max")
    print(tmp, digits = digits)
  } else {
    # Coefficients
    cat("Coefficients:\n")
    print(smry_funM(object$coefficients, CL = CL, na.rm = na.rm), digits = digits)

    cat("\nP-values:\n")
    print(smry_funM(object$pvals, CL = CL, na.rm = na.rm), digits = digits)

    cat("\nCorrelation among time series residuals:\n")
    tmp = summary(as.vector(cor(t(object$residuals))))
    names(tmp) <- c("Min", "1Q", "Median", "Mean", "3Q", "Max")
    print(tmp, digits = digits)

    cat("\nModel fit:\n")
    print(rbind(MSE = smry_funV(object$MSE, CL = CL, na.rm = na.rm),
                LogLik = smry_funV(object$logLik, CL = CL, na.rm = na.rm)),
          digits = digits)
  }
}

#' @title helper summary function (matrix)
#'
#' @rdname TS.methods
#'
#' @param x numeric matrix
#' @param CL confidence level (default = .95)
#' @param na.rm logical, should observations with NA be removed?
#'
#' @return summary statistics for each column including quartiles, mean, and
#' upper and lower confidence levels (given by CL)
smry_funM <- function(x, CL = .95, na.rm = TRUE){
  alph = 1 - CL
  mean <- apply(x, 2, mean, na.rm = na.rm)
  quarts <- apply(x, 2, quantile,
                  probs = c(0, .25, .5, .75, 1),
                  na.rm = na.rm)
  CI <- apply(x, 2, quantile, probs = c(alph/2, 1-(alph/2)),
              na.rm = na.rm)
  out <- rbind(quarts, mean, CI)
  row.names(out) <- c("Min", "1Q", "Median", "3Q", "Max", "Mean", "CL.lower", "CL.upper")
  return(t(out))
}

#' @title helper summary function (vector)
#'
#' @rdname TS.methods
#'
#' @param x numeric matrix
#' @param CL confidence level (default = .95)
#' @param na.rm logical, should observations with NA be removed?
#'
#' @return summary statistics including quartiles, mean, and upper and lower
#' confidence levels (given by CL)
smry_funV <- function(x, CL = .95, na.rm = TRUE){
  alph = 1 - CL
  mean = mean(x, na.rm = na.rm)
  quarts = quantile(x, probs = c(0, .25, .5, .75, 1))
  CI <- quantile(x, probs = c(alph/2, 1-(alph/2)))
  out <- c(quarts, mean, CI)
  names(out) <- c("Min", "1Q", "Median", "3Q", "Max", "Mean", "CL.lower", "CL.upper")
  return(out)
}

## summary ----
## summary.lm will work for fitCLS but NOT fitAR, so a new method should exist?

## coef ----
## coef() already works for fitAR, fitCLS, fitAR_map, and fitCLS_map

## residuals ----
## resid() already works for fitAR, fitCLS, fitAR_map, and fitCLS_map

## plot methods ----
## plot.lm will no work for either fitCLS or fitAR, so a new method should exist?
