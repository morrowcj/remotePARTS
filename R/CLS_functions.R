## fitCLS ----
#' @title CLS for time series
#'
#' @description \code{fitCLS} is used to fit conditional least squares regression
#' to time series data.
#'
#' @family remoteTS
#'
#' @param formula a model formula, as used by \code{stats::lm()}
#' @param data optional data environment to search for variables in \code{formula}.
#' As used by \code{lm()}
#' @param lag.y an integer indicating the lag (in time steps) between y and y.0
#' @param lag.x an integer indicating the lag (in time steps) between y and the
#' independent variables (except y.0).
#' @param model logical, should the used model matrix be returned? As used by
#' \code{lm()}
#' @param y logical, should the used response variable be returned? As used by
#' \code{lm()}
#' @param debug logical debug mode
#'
#' @details
#' This function regresses the response variable (y) at time t, conditional on the
#' response at time t-\code{lag.y } and the specified dependent variables (X) at
#' time t-\code{lag.x}:
#'
#' \deqn{y(t) = y(t - lag.y) + X(t - lag.x) + \varepsilon}{y(t) = y(t - lag.y) + X(t - lag.x) + e}
#'
#' where \eqn{y(t)} is the response at time \eqn{t};
#'
#' \eqn{X(t)} is a model matrix containing covariates;
#'
#' \eqn{\beta}{beta} is a vector of effects of \eqn{X(t)};
#'
#' and \eqn{\varepsilon(t)}{t(t)} is a temporally independent Gaussian random
#' variable with mean zero and standard deviation \eqn{\sigma}{sigma}
#'
#' \code{stats::lm()} is then called, using the above equation.
#'
#' @return \code{fitCLS} returns a list object of class "remoteTS", which
#' inherits from  "lm". In addition to the default "lm" components, the output
#' contains these additional list elements:
#'
#' \describe{
#'    \item{tstat}{the t-statistics for coefficients}
#'    \item{pval}{the p-values corresponding to t-tests of coefficients}
#'    \item{MSE}{the model mean squared error}
#'    \item{logLik}{the log-likelihood of the model fit}
#' }
#'
#' @seealso \code{\link{fitCLS_map}} to easily apply \code{fitCLS} to many pixels;
#' \code{\link{fitAR}} and \code{\link{fitAR_map}} for AR time series analyses.
#'
#' @examples
#'
#' # simulate dummy data
#' t = 1:30 # times series
#' Z = rnorm(30) # random independent variable
#' x = .2*Z + (.05*t) # generate dependent effects
#' x[2:30] = x[2:30] + .2*x[1:29] # add autocorrelation
#' x = x + rnorm(30, 0, .01)
#' df = data.frame(x, t, Z) # collect in data frame
#'
#' # fit a CLS model with previous x, t, and Z as predictors
#' ## note, this model does not follow the underlying process.
#' ### See below for a better fit.
#' (CLS <- fitCLS(x ~ t + Z, data = df))
#'
#' # extract other values
#' CLS$MSE #MSE
#' CLS$logLik #log-likelihood
#'
#' # fit with no lag in independent variables (as simulated):
#' (CLS2 <- fitCLS(x ~ t + Z, df, lag.x = 0))
#' summary(CLS2)
#'
#' # no lag in x
#' fitCLS(x ~ t + Z, df, lag.y = 0)
#'
#' # visualize the lag
#' ## large lag in x
#' fitCLS(x ~ t + Z, df, lag.y = 2, lag.x = 0, debug = TRUE)$lag
#' ## large lag in Z
#' fitCLS(x ~ t + Z, df, lag.y = 0, lag.x = 2, debug = TRUE)$lag
#'
#' # # throws errors (NOT RUN)
#' # fitCLS(x ~ t + Z, df, lag.y = 28) # longer lag than time
#' # fitCLS(cbind(x, rnorm(30)) ~ t + Z, df) # matrix response
#'
#' ## Methods
#' summary(CLS)
#' residuals(CLS)
#'
#' @export
fitCLS <- function(formula, data = NULL, lag.y = 1, lag.x = 1, debug = FALSE,
                   model = FALSE, y = FALSE){
  # structure data from input
  call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  mf$formula <- update(formula, . ~ . -1)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- quote(stats::model.frame) # rename the function call
  mf$drop.unused.levels <- TRUE
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  resp <- stats::model.response(mf, "numeric")
  if (is.matrix(resp)){stop("response is a matrix: must be a vector")} # throw dimension error
  stopifnot(max(lag.x, lag.y) + 2 < length(resp)) # throw lag error
  X <- stats::model.matrix(mt, mf, contrasts = NULL)

  # variable names
  resp.name = all.vars(update(formula, . ~ 1))
  # temporal indices
  current <- (max(lag.x, lag.y) + 2):length(resp)
  previous <- current - lag.y
  lagged <- current - lag.x
  # shift data
  y.i <- resp[current] # response at current time
  y.0 <- resp[previous] # response at previous time
  X.lag <- X[lagged, ] # dependent variables (lagged)
  # stack shifted data into dataframe
  tmp <- data.frame(y.i, y.0, X.lag)
  prev.name <- paste0(resp.name, ".0")
  names(tmp) <- c(resp.name, prev.name, colnames(X))
  # fit the model
  update.form <- paste0(". ~ ", prev.name, " + .")
  fm <- lm(update(formula, update.form), tmp, model = model, y = y)

  smry <- suppressWarnings(summary(fm))
  ## update the output
  fm$call <- call # update call
  fm$SE <- smry$coefficients[, "Std. Error"] #standard errors
  fm$tstat <- smry$coefficients[, "t value"] #t-statistic
  fm$pval <- smry$coefficients[, "Pr(>|t|)"] #p-value
  fm$MSE <- suppressWarnings(anova(fm))["Residuals", "Mean Sq"]#MSE
  fm$logLik <- logLik(fm) # include log-likelihood

  class(fm) <- append("remoteTS", class(fm))

  # return statement
  if(debug) { # debug mode
    mf.new <- stats::model.frame(update(formula, update.form), tmp)
    return(list(mf = mf, y = y, X = X, resp.name = resp.name, lag = cbind(current, previous, lagged),
                tmp = tmp, mf.new = mf.new, fm = fm))
  } else { # default mode

    return(fm)
  }
}

## fitCLS_map ----
#' @title Map-level CLS for time series
#'
#' @description \code{fitCLS_map} is used to fit conditional least squares
#' regression to each spatial location (pixel) within spatiotemporal data.
#'
#' @family remoteTS
#'
#' @param Y a spatiotemporal response variable: a numeric matrix or data frame
#' where columns correspond to time points and rows correspond to pixels.
#' @param coords a numeric coordinate matrix or data frame, with two columns and
#' rows corresponding to each pixel
#' @param formula a model formula, passed to \code{fitCLS()}: the left side
#' of the formula should always be "y" and the right hand side should refer to
#' variables in \code{X.list}
#' @param X.list a named list of temporal or spatiotemporal predictor variables:
#' elements must be either numeric vectors with one element for each time point
#' or a matrix/data frame with rows corresponding to pixels and columns
#' corresponding to time point. These elements must be named and referred to
#' in \code{formula}
#' @param lag.y the lag between y and y.0, passed to \code{fitCLS()}
#' @param lag.x the lag between y and predictor variables, passed to
#' \code{fitCLS()}
#' @param resids.only logical: should output beyond coordinates and residuals be
#' withheld? Useful when passing output to \code{fitCor()}
#'
#' @details \code{fitCLS_map} is a wrapper function that applies
#' \code{fitCLS()} to many pixels.
#'
#' The function loops through the rows of \code{Y}, matched with rows of
#' spatiotemporal predictor matrices. Purely temporal predictors, given by
#' vectors, are used for all pixels. These predictor variables, given by the
#' right side of \code{formula} are sourced from named elements in \code{X.list}.
#'
#' @seealso \code{\link{fitCLS}} for fitting CLS on individual time series and
#' \code{\link{fitAR}} and \code{\link{fitAR_map}} for AR REML time series analysis.
#'
#' @return \code{fitCLS_map} returns a list object of class "mapTS".
#'
#' The output will always contain at least these elements:
#'
#' \describe{
#'    \item{call}{the function call}
#'    \item{coords}{the coordinate matrix or dataframe}
#'    \item{residuals}{time series residuals: rows correspond to pixels (\code{coords})}
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
#'
#' # simulate dummy data
#' time.points = 9 # time series length
#' map.width = 5 # square map width
#' coords = expand.grid(x = 1:map.width, y = 1:map.width) # coordinate matrix
#' ## create empty spatiotemporal variables:
#' X <- matrix(NA, nrow = nrow(coords), ncol = time.points) # response
#' Z <- matrix(NA, nrow = nrow(coords), ncol = time.points) # predictor
#' # setup first time point:
#' Z[, 1] <- .05*coords[,"x"] + .2*coords[,"y"]
#' X[, 1] <- .5*Z[, 1] + rnorm(nrow(coords), 0, .05) #x at time t
#' ## project through time:
#' for(t in 2:time.points){
#'   Z[, t] <- Z[, t-1] + rnorm(map.width^2)
#'   X[, t] <- .2*X[, t-1] + .1*Z[, t] + .05*t + rnorm(nrow(coords), 0 , .25)
#' }
#'
#' # # visualize dummy data (NOT RUN)
#' # library(ggplot2);library(dplyr)
#' # data.frame(coords, X) %>%
#' #   reshape2::melt(id.vars = c("x", "y")) %>%
#' #   ggplot(aes(x = x, y = y, fill = value)) +
#' #   geom_tile() +
#' #   facet_wrap(~variable)
#'
#' # fit CLS, showing all output
#' fitCLS_map(X, coords, formula = y ~ t, resids.only = TRUE)
#'
#' # fit CLS with temporal and spatiotemporal predictors
#' (CLS.map <- fitCLS_map(X, coords, formula = y ~ t + Z,
#'                        X.list = list(t = 1:ncol(X), Z = Z),
#'                        resids.only = FALSE))
#' ## extract some values
#' CLS.map$coefficients # coefficients
#' CLS.map$logLik # log-likelihoods
#'
#' ## Methods
#' summary(CLS.map)
#' residuals(CLS.map)
#' coefficients(CLS.map)
#'
#' @export
fitCLS_map <- function(Y, coords, formula = "y ~ t",
                       X.list = list(t = 1:ncol(Y)),
                       lag.y = 1, lag.x = 0, resids.only = FALSE){
  call = match.call()

  ## input checks
  if (!is.matrix(Y) & !is.data.frame(Y)) { # Y
    stop(paste("Y must be a matrix or dataframe with ncol(Y) time points",
               "and nrow(Y) pixels"))
  }
  Y = as.matrix(Y)
  if (!is.matrix(coords) & !is.data.frame(coords) |
      nrow(coords) != nrow(Y) | ncol(coords) != 2){
    stop(paste("coords must be a matrix or dataframe with 2 columns",
               "and rows equal to nrow(Y)"))
  }
  coords = as.matrix(coords)
  if (is.null(names(X.list))) {
    stop(paste("X.list must be a named list, with element names",
               "corresponding to the right hand side of formula"))
  }

  ## data dimensions
  n.pix = nrow(Y) # pixels
  n.time = ncol(Y) # time points

  ## setup output
  residuals = matrix(NA, nrow = n.pix, ncol = length((max(lag.x, lag.y) + 2):n.time) )

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
    cls = fitCLS(formula = as.formula(formula), data = df, lag.y = lag.y, lag.x = lag.x,
                 debug = FALSE, model = FALSE)

    ## build output elements
    residuals[i, ] <- cls$residuals
    if (!resids.only){ # additional output
      if(i == 1) {
        ## output setup
        coefficients <- matrix(NA, nrow = n.pix, ncol = length(cls$coefficients),
                               dimnames = list(NULL ,names(cls$coefficients)))
        SEs <- matrix(NA, nrow = n.pix, ncol = length(cls$SE),
                      dimnames = list(NULL, names(cls$SE)))
        tstats <- matrix(NA, nrow = n.pix, ncol = length(cls$tstat),
                         dimnames = list(NULL, names(cls$tstat)))
        pvals <- matrix(NA, nrow = n.pix, ncol = length(cls$pval),
                        dimnames = list(NULL, names(cls$pval)))
        MSEs <- vector("numeric", n.pix)
        logLiks <- vector("numeric", n.pix)
        fitted.values <- matrix(NA, nrow = n.pix, ncol = length(cls$fitted.values))
      }
      ## add to output
      coefficients[i, ] <- cls$coefficients
      SEs[i, ] <- cls$SE
      tstats[i, ] <- cls$tstat
      pvals[i, ] <- cls$pval
      MSEs[i] <- cls$MSE
      logLiks[i] <- cls$logLik
      fitted.values[i, ] <- cls$fitted.values
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
