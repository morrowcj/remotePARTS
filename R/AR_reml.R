## Pixel AR REML----
#' @title AR regressions by REML
#'
#' @description \code{fitAR} is used to fit AR(1) time series regression
#' analysis using restricted maximum likelihood
#'
#' @family remoteTS
#'
#' @param formula a model formula, as used by \code{stats::lm()}
#' @param data optional data environment to search for variables in \code{formula}.
#' As used by \code{lm()}
#'
#' @details
#' This function finds the restricted maximum likelihood (REML) to estimate
#' parameters for the regression model with AR(1) random error terms
#'
#' \deqn{y(t) =  X(t) \beta + \varepsilon(t)}{y(t) = X(t)*beta + e(t)}
#' \deqn{\varepsilon(t) =  \rho \varepsilon(t-1) + \delta(t)}{e(t) = rho*e(t-1) + delta(t)}
#'
#' where \eqn{y(t)} is the response at time \eqn{t};
#'
#' \eqn{X(t)} is a model matrix containing covariates;
#'
#' \eqn{\beta}{beta} is a vector of effects of \eqn{X(t)};
#' \eqn{\varepsilon(t)}{e(t)} is the autocorrelated random error;
#'
#' \eqn{\delta \sim N(0, \sigma)}{delta ~ N(0, sigma)} is a temporally independent
#' Gaussian random variable with mean zero and standard deviation
#' \eqn{\sigma}{sigma};
#'
#' and \eqn{\rho}{rho} is the AR(1) autoregression parameter
#'
#' \code{fitAR} estimates the parameter via mathematical optimization
#' of the restricted log-likelihood function calculated by the workhorse function
#' \code{remotePARTS:::AR_fun()}.
#'
#' @references
#'
#' Ives, A. R., K. C. Abbott, and N. L. Ziebarth. 2010. Analysis of ecological
#'
#'     time series with ARMA(p,q) models. Ecology 91:858-871.
#'
#' @return \code{fitAR} returns a list object of class "remoteTS", which contains
#' the following elements.
#'
#' \describe{
#'     \item{call}{the function call}
#'     \item{coefficients}{a named vector of coefficients}
#'     \item{SE}{the standard errors of parameter estimates}
#'     \item{tstat}{the t-statistics for coefficients}
#'     \item{pval}{the p-values corresponding to t-tests of coefficients}
#'     \item{MSE}{the model mean squared error}
#'     \item{logLik}{the log-likelihood of the model fit}
#'     \item{residuals}{the residuals: response minus fitted values}
#'     \item{fitted.values}{the fitted mean values}
#'     \item{rho}{The AR parameter, determined via REML}
#'     \item{rank}{the numeric rank of the fitted model}
#'     \item{df.residual}{the residual degrees of freedom}
#'     \item{terms}{the \code{stats::terms} object used}
#' }
#'
#' Output is structured similarly to an "lm" object.
#'
#' @seealso \code{\link{fitAR_map}} to easily apply \code{fit_AR} to many pixels;
#' \code{\link{fitCLS}} and \code{\link{fitCLS_map}} for conditional least squares
#' time series analyses.
#'
#' @examples
#'
#' # simulate dummy data
#' t = 1:30 # times series
#' Z = rnorm(30) # random independent variable
#' x = .2*Z + (.05*t) # generate dependent effects
#' x[2:30] = x[2:30] + .2*x[1:29] # add autocorrelation
#'
#' # fit the AR model, using Z as a covariate
#' (AR = fitAR(x ~ Z))
#'
#' # get specific components
#' AR$residuals
#' AR$coefficients
#' AR$pval
#'
#' # now using time as a covariate
#' (AR.time <- fitAR(x ~ t))
#'
#' # source variable from a dataframe
#' df = data.frame(y = x, t.scaled = t/30, Z = Z)
#' fitAR(y ~ t.scaled + Z, data = df)
#'
#' ## Methods
#' summary(AR)
#' residuals(AR)
#' coefficients(AR)
#'
#' @export
fitAR <- function(formula, data = NULL){
  # structure data from input
  call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- quote(stats::model.frame) # rename the function call
  mf$drop.unused.levels <- TRUE
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  resp <- stats::model.response(mf, "numeric")
  if (is.matrix(resp)){stop("response is a matrix: must be a vector")} # throw dimension error
  X <- stats::model.matrix(mt, mf, contrasts = NULL)

  # return(list(call, mf = mf, mt = mt, resp = resp, X = X))

  # Optimize over AR_fun() for par
  opt <- stats::optim(fn = AR_fun, par = 0.2, y = resp, X = X, logLik.only = TRUE,
                      method = "Brent", upper = 1, lower = -1,
                      control = list(maxit = 10^4))

  b <- opt$par # optimized parameter

  # Perform the AR regression with the optimized parameter and return results
  AR.out = AR_fun(par = b, y = resp, X = X, logLik.only = FALSE)
  AR.out$call = call

  AR.out$rho = b
  AR.out$rank = ncol(X)
  AR.out$df.residual = nrow(X) - ncol(X)
  AR.out$terms = terms(as.formula(formula))

  class(AR.out) <- append("lm", class(AR.out))
  class(AR.out) <- append("remoteTS", class(AR.out))

  return(AR.out)
}

## AR workhorse function ----
#' @rdname fitAR
#' @family remoteTS
#'
#' @param par AR parameter value
#' @param y vector of time series (response)
#' @param X model matrix (predictors)
#' @param logLik.only logical: should only the partial log-likelihood be computed
#'
#' @details \code{AR_fun} is the work horse behind \code{fitAR} that is called
#' by \code{optim} to estimate the autoregression parameter \eqn{\rho}{rho}.
# ' \code{AR_fun} calculates the restricted log likelihood function.
#'
#' @return
#' When \code{logLik.only == F}, \code{AR_fun} returns the output described in
#' \code{?fitAR}. When \code{logLik.only == T}, it returns a quantity that is
#' linearly and negatively related to the restricted log likelihood
#' (i.e., partial log-likelihood).
#'
#' @examples
#'
#' # simulate dummy data
#' # x = rnorm(31)
#' # time = 1:30
#' # x = x[2:31] + x[1:30] + 0.3*time #AR(1) process + time trend
#' # U = stats::model.matrix(formula(x ~ time))
#'
#' # fit an AR
#' # remotePARTS:::AR_fun(par = .2, y = x, X = U, logLik.only = FALSE)
#'
#' # get the partial logLik of the AR parameter, given the data.
#' # remotePARTS:::AR_fun(par = .2, y = x, X = U, logLik.only = FALSE)
#'
#' # # show that minimizing the partial logLik maximizes the true logLik (NOT RUN)
#' # n = 100
#' # out.mat = matrix(NA, nrow = n, ncol = 3,
#' #                  dimnames = list(NULL, c("par", "logLik", "partialLL")))
#' # out.mat[, "par"] = seq(-10, 10, length.out = n)
#' # for (i in seq_len(n) ) {
#' #    p = out.mat[i, "par"]
#' #    out.mat[i, "logLik"] = remotePARTS:::AR_fun(par = p, y = x, X = U, logLik.only = TRUE)
#' #    out.mat[i, "partialLL"] = remotePARTS:::AR_fun(par = p, y = x, X = U,
#' #                                                   logLik.only = FALSE)$logLik
#' # }
#' # plot(x = out.mat[, "partialLL"], y = out.mat[, "logLik"])
AR_fun <- function(par, y, X, logLik.only = TRUE) {
  call = match.call()

  # setup variables
  b <- par # parameter of interest
  n.obs <- length(y) # number of time points
  q <- ncol(X) # number of covariates in model matrix
  B <- diag(n.obs) # n length identity matrix
  diag(B[-1, ]) <- -b # set sub-diagonal to -b
  iS <- diag(n.obs) # another identity matrix
  iS[1, 1] <- (1 - b^2) # set first element to 1-b^2
  iV <- t(B) %*% iS %*% B

  # remove NA elements
  if(any(is.na(y))){
    iV <- iV[!is.na(y), !is.na(y)]
    X <- X[!is.na(y),]
    y <- y[!is.na(y)]
  }

  # log determinant of iV
  logdetV <- -determinant(iV)$modulus[1]

  # Regression calculations
  ## solve Am + B for m where A = (X'*iV*X) and B = (X'*iV*y)
  ### beta is the effect of the covariates
  beta <- solve(t(X) %*% iV %*% X, t(X) %*% iV %*% y)
  # remove the effect of the covariates from y
  H <- y - X %*% beta
  # estimate the variance
  s2 <- (t(H) %*% iV %*% H)/(n.obs - q)
  # calculate the partial log-likelihood of b given y and X
  logLik <- 0.5 * ((n.obs - q) * log(s2) + logdetV +
                 determinant(t(X) %*% iV %*% X)$modulus[1] +
                 (n.obs - q))

  if(logLik.only){
    # return log-likelihood
    return(as.vector(logLik))

  } else {

    # important stats
    MSE <- as.numeric(s2) #MSE
    s2beta <- MSE * solve(t(X) %*% iV %*% X) #SE
    t.stat = (abs(beta) / diag(s2beta)^0.5)
    pval = 2 * stats::pt(q = t.stat, df = n.obs - q,
                         lower.tail = FALSE )

    yhat = X %*% beta

    # log likelihood without constants (i.e. s2) - no parameter dependancy
    logLik <- 0.5 * (n.obs - q) * log(2 * pi) +
      determinant(t(X) %*% X)$modulus[1] - logLik

    # collect output
    out.list <- list(call = call,
                     coefficients = beta[, 1],
                     SE = diag(s2beta),
                     tstat = t.stat[, 1],
                     pval = pval[, 1],
                     MSE = MSE,
                     logLik = logLik[, 1],
                     residuals = as.vector(H),
                     fitted.values = as.vector(yhat))

    return(out.list)
  }
}


## Full Map AR ----
#' @title Map-level AR REML
#'
#' @description \code{fitAR_map} is used to fit AR REML regression to each spatial
#' location (pixel) within spatiotemporal data.
#'
#' @family remoteTS
#'
#' @param Y a spatiotemporal response variable: a numeric matrix or data frame
#' where columns correspond to time points and rows correspond to pixels.
#' @param coords a numeric coordinate matrix or data frame, with two columns and
#' rows corresponding to each pixel
#' @param formula a model formula, passed to \code{fitAR()}: the left side
#' of the formula should always be "y" and the right hand side should refer to
#' variables in \code{X.list}
#' @param X.list a named list of temporal or spatiotemporal predictor variables:
#' elements must be either numeric vectors with one element for each time point
#' or a matrix/data frame with rows corresponding to pixels and columns
#' corresponding to time point. These elements must be named and referred to
#' in \code{formula}
#' @param resids.only logical: should output beyond coordinates and residuals be
#' withheld? Useful when passing output to \code{fitCor()}
#'
#' @details \code{fitAR_map} is a wrapper function that applies \code{fitAR} to
#' many pixels.
#'
#' The function loops through the rows of \code{Y}, matched with rows of
#' spatiotemporal predictor matrices. Purely temporal predictors, given by
#' vectors, are used for all pixels. These predictor variables, given by the
#' right side of \code{formula} are sourced from named elements in \code{X.list}.
#'
#' @seealso \code{\link{fitAR}} for fitting AR REML to individual time series and \code{\link{fitCLS}}
#' & \code{\link{fitCLS_map}} for time series analysis based on conditional least squares.
#'
#' @return \code{fitCLS_map} returns a list object of class "mapTS".
#'
#' The output will always contain at least these elements:
#'
#' \describe{
#'    \item{call}{the function call}
#'    \item{coords}{the coordinate matrix or dataframe}
#'    \item{residuals}{time series residuals: rows correspond to pixels
#'    (\code{coords})}
#' }
#'
#' When \code{resids.only = FALSE}, the output will also contain the following
#' components. Matrices have rows that correspond to pixels and columns that
#' correspond to time points and vector elements correspond to pixels.
#'
#' \describe{
#'    \item{coefficients}{a numeric matrix of coefficeints}
#'    \item{SEs}{a numeric matrix of coefficient standard errors}
#'    \item{tstats}{a numeric matrix of t-statistics for coefficients}
#'    \item{pvals}{a numeric matrix of p-values for coefficients t-tests}
#'    \item{MSEs}{a numeric vector of MSEs}
#'    \item{logLiks}{a numeric vector of log-likelihoods}
#'    \item{fitted.values}{a numeric matrix of fitted values}
#' }
#'
#' An attribute called "resids.only" is also set to match the value of
#' \code{resids.only}
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
#' # # visualize dummy data (NOT RUN)
#' # library(ggplot2);library(dplyr)
#' # data.frame(coords, X) %>%
#' #   reshape2::melt(id.vars = c("x", "y")) %>%
#' #   ggplot(aes(x = x, y = y, fill = value)) +
#' #   geom_tile() +
#' #   facet_wrap(~variable)
#'
#' # fit AR, showing all output
#' fitAR_map(X, coords, formula = y ~ t, resids.only = TRUE)
#'
#' # fit AR with temporal and spatiotemporal predictors
#' (AR.map <- fitAR_map(X, coords, formula = y ~ t + Z, X.list = list(t = 1:ncol(X),
#'                      Z = Z), resids.only = FALSE))
#' ## extract some values
#' AR.map$coefficients # coefficients
#' AR.map$logLik # log-likelihoods
#'
#' ## Methods
#' summary(AR.map)
#' residuals(AR.map)
#' coefficients(AR.map)
#'
#' @export
fitAR_map <- function(Y, coords, formula = "y ~ t",
                      X.list = list(t = 1:ncol(Y)),
                      resids.only = FALSE){
  call = match.call()

  ## input checks
  if (!is.matrix(Y) & !is.data.frame(Y)) { # Y
    stop(paste("Y must be a matrix or dataframe with ncol(Y) time points",
               "and nrow(Y) pixels"))
  }
  if (!is.matrix(coords) & !is.data.frame(coords) |
      nrow(coords) != nrow(Y) | ncol(coords) != 2){
    stop(paste("coords must be a matrix or dataframe with 2 columns",
               "and rows equal to nrow(Y)"))
  }
  if (is.null(names(X.list))) {
    stop(paste("X.list must be a named list, with element names",
               "corresponding to the right hand side of formula"))
  }

  ## data dimensions
  n.pix = nrow(Y) # pixels
  n.time = ncol(Y) # time points

  ## setup output
  residuals = matrix(NA, nrow = n.pix, ncol = n.time)

  ## apply CLS to map
  for (i in 1:n.pix) {

    ## extract one pixel of data from Y into a dataframe
    df = data.frame(y = Y[i, ]) # times series at pixel i

    ## Parse covariate list, filling in df
    for(j in names(X.list)) {
      if (length(X.list[[j]]) == ncol(Y)) { ## purely temporal predictor
        df[, j] = X.list[[j]]

      } else if (all(dim(X.list[[j]]) == c(nrow(Y), ncol(Y)))) { ## spatio-temporal predictor
        df[, j] = X.list[[j]][i, ]

      } else { ## dimension mismatch
        stop(paste0("all elements j of X.list must satisfy one of the following conditions:",
                    "\n  1) dim(X.list[[j]]) == c(nrow(Y), ncol(Y))",
                    "\n  2) length(X.list[[j]]) == ncol(Y)",
                    "\n  culprit: X.list[['", j,"']]"))
      }
    }

    ## fit CLS
    ar = fitAR(formula = as.formula(formula), data = df)

    ## build output elements
    residuals[i, ] <- ar$residuals
    if (!resids.only){ # additional output
      if(i == 1) {
        ## output setup
        coefficients <- matrix(NA, nrow = n.pix, ncol = length(ar$coefficients),
                               dimnames = list(NULL ,names(ar$coefficients)))
        SEs <- matrix(NA, nrow = n.pix, ncol = length(ar$SE),
                      dimnames = list(NULL, names(ar$SE)))
        tstats <- matrix(NA, nrow = n.pix, ncol = length(ar$tstat),
                         dimnames = list(NULL, names(ar$tstat)))
        pvals <- matrix(NA, nrow = n.pix, ncol = length(ar$pval),
                        dimnames = list(NULL, names(ar$pval)))
        MSEs <- vector("numeric", n.pix)
        logLiks <- vector("numeric", n.pix)
        fitted.values <- matrix(NA, nrow = n.pix, ncol = length(ar$fitted.values))
      }
      ## add to output
      coefficients[i, ] <- ar$coefficients
      SEs[i, ] <- ar$SE
      tstats[i, ] <- ar$tstat
      pvals[i, ] <- ar$pval
      MSEs[i] <- ar$MSE
      logLiks[i] <- ar$logLik
      fitted.values[i, ] <- ar$fitted.values
    }
  }

  ## conditional return
  if(resids.only){
    out <- list(call = call, coords = coords, residuals = residuals)
    attr(out, "resids.only") = TRUE
  } else {
    out <- list(call = call, coords = coords, coefficients = coefficients, SEs = SEs, tstats = tstats,
                pvals = pvals, MSEs = MSEs, logLiks = logLiks, fitted.values = fitted.values,
                residuals = residuals)
    attr(out, "resids.only") = FALSE
  }

  class(out) <- append("mapTS", class(out))

  return(out)
}
