# pixel level CLS ----
#' @title Pixel CLS
#' @description Fit a CLS regression to a time series
#' @family remoteCLS
#'
#' @param x numeric vector of length p containing time series data for one
#' location (i.e. a pixel)
#' @param t numeric vector of length p containing the values for time.
#' @param Z (optional) a numeric vector/matrix with p rows and q columns. Each
#' column contains a different covariate measured at the same site as x at the
#' same p time points.
#' @param save_AR.df should the auto-regression data frame be returned?
#' default: FALSE
#'
#' @return \code{remoteCLS.pixel} object. A list with the following elements:
#' \describe{
#'     \item{\code{$call}}{the matched call to this function}
#'     \item{\code{$fm}}{the model object fit using \code{stats::lm()}}
#'     \item{\code{$AR.df}}{an optionally returned AR data frame built with
#'     \code{AR_df()}}
#' }
#'
#' @details CLS is fit to a time series \code{x} by first building an AR
#' data frame with \code{AR_df()}. Then, \code{x} at time j (\code{x.tj})
#' is regressed on \code{x} at time i (\code{x.ti}) and \code{t} at time j
#' (\code{t.j}).
#'
#' By default the print.remoteCLS() method does not *show* all output,
#' even if the ret_* flags are set to TRUE. Use the S3 \code{$} operator to
#' access elements by name or use one of the extractor functions
#' (shown in examples below).
#'
#' @seealso [AR_df()] for AR data frames, [fitCLS.map()] for fitting full-map
#' time series CLS, [fitAR()] & [fitAR.map()] for using AR REML instead of CLS.
#'
#' @export
#'
#' @examples
#' data(ndvi_AK3000)
#' x = unlist(ndvi_AK3000[1, -c(1:6)]) # time series for first pixel only
#' t = 1:length(x) # time points
#'
#' ## CLS without covariates
#' fitCLS(x = x, t = t)
#'
#' ## now with 2 covariates
#' Z = cbind(rnorm(length(x)), rnorm(length(x)))
#' fitCLS(x = x, t = t, Z = Z)
#'
#' ## get the AR data frame with save_AR.df = TRUE
#' fitCLS(x = x, t = t, save_AR.df = TRUE)$AR.df
#'
#' ## print a summary
#' summary(fitCLS(x = x, t = t))
#'
#' ## use get_fm() to extract $fm, which is a stats::lm() object
#' fm = get_fm(fitCLS(x = x, t = t))
#' ## then treat it just like lm()
#' summary(fm)
#' resid(fm)
fitCLS <- function(x, t, Z = NULL, save_AR.df = FALSE) {
  stopifnot(length(x) == length(t)) # check
  # build AR data frame
  AR.df <- AR_df(x, t, Z)
  # fit the model
  fm <- if(is.null(Z)) {
     lm(x.tj ~ x.ti + t.j, data = AR.df)
  } else {
     lm(x.tj ~ x.ti + t.j + z, data = AR.df)
  }
  # build output list
  out <- list(call = match.call(),
              fm = fm)
  # add AR.df to output list
  if(save_AR.df){
    out$AR.df = AR.df
  } else {out$AR.df = NULL}
  # set remoteCLS class
  class(out) <- c("remoteCLS", "pixel")

  return(out)
}

## map level CLS ----
#' @title Map CLS
#' @description Fit a CLS model to an entire map or map subset
#'
#' @param X \eqn{n x p} matrix of data with 1-n rows of pixels and 1-p columns of time
#' points
#' @param t time vector of length p
#' @param ret_xi.coef logical: return/save coefficient table for xi?
#' @param ret_int.coef logical: return/save coefficient table for intercept?
#' @param ret_MSE logical: return/save model MSE?
#' @param ret_resid logical: return/save model residuals?
#'
#' @return \code{remoteCLS.map} object. A list with the following elements:
#' \describe{
#'    \item{\code{$call}}{the matched call to this function}
#'    \item{\code{$time.coef}}{coefficient table for the effect of time
#'    (1 row per pixel)}
#'    \item{\code{$mean}}{mean of X at each pixel, averaged across time}
#'    \item{\code{$xi.coef}}{optional coefficient table for the effect of x at
#'    time i
#'    on x at time j (one row per pixel)}
#'    \item{\code{$int.coef}}{optional coefficient table for the intercept}
#'    \item{\code{$MSE}}{model MSE for each pixel}
#'    \item{\code{$residuals}}{matrix of model residuals (1 row per pixel)}
#' }
#'
#' @details \code{fitCLS.map()} is a vectorized version of \code{fitCLS()},
#' which is called internally for each row of \code{X}.
#'
#' By default the print.remoteCLS() method does not *show* all output,
#' even if the ret_* flags are set to TRUE. Use the S3 \code{$} operator to
#' access elements by name or use one of the extractor functions
#' (shown in examples below).
#'
#' @seealso [fitCLS()] for fitting CLS on individual time series and [fitAR()]
#' & [fitAR.map()] for AR REML time series analysis.
#'
#' @export
#'
#' @examples
#' data(ndvi_AK3000)
#' X = (ndvi_AK3000[1:20, -c(1:6)]) # first 20 pixels of NDVI data
#' t = 1:ncol(X) # time points
#'
#' # fit CLS (only save time coefficient)
#' CLS <- fitCLS.map(X, t)
#' CLS # print.remoteCLS() only prints model info and time coefficients by default
#' coef(CLS) # coefficient table as data frame
#'
#' # fit again, but save all possible output
#' CLS.full <- fitCLS.map(X, t, ret_xi.coef = TRUE, ret_int.coef = TRUE,
#'                        ret_MSE = TRUE, ret_resid = TRUE)
#'
#' # printing CLS.full still only prints model info and time coefficients by default:
#' CLS.full
#'
#' # But all coefficeints are still contained in the output and can be accessed:
#' coef(CLS.full) # time coefficient table; also coef(CLS.full, "time") or CLS.full$time.coef
#' coef(CLS.full, "xi") # extract coefficient table for xi; also CLS.full$xi.coef
#' coef(CLS.full, "intercept") # extract the coefficient table for intercept; also CLS.full$int.coef
#'
#' resid(CLS.full) # residual matrix; also CLS.full$residuals
#'
#' CLS.full$MSE # vector of MSEs
#'
#' summary(CLS.full) # summaries
fitCLS.map <- function(X, t, ret_xi.coef = FALSE, ret_int.coef = FALSE,
                       ret_MSE = TRUE, ret_resid = TRUE){
  stopifnot(ncol(X) == length(t))
  n.pixels = nrow(X)
  n.time = length(t)

  ## run CLS
  cls.list <- lapply(1:n.pixels, function(x){
    fitCLS(unlist(X[x,]), t, save_AR.df = TRUE)
  })

  ## get coefficients from cls
  coef.list <- lapply(cls.list, function(x){
    as.data.frame(summary(get_fm(x))$coefficient)
  })

  ## setup coefficient tables
  out.list <- list(call = match.call(),
                   time.coef = as.data.frame(matrix(NA, ncol = 4, nrow = n.pixels)))
  colnames(out.list$time.coef) <- c("Est", "SE", "t","p.t")
  out.list$mean = rowMeans(X)
  if (ret_xi.coef){
    out.list$xi.coef <- as.data.frame(matrix(NA, ncol = 4, nrow = n.pixels))
    colnames(out.list$xi.coef) <- c("Est", "SE", "t","p.t")
  }
  if (ret_int.coef){
    out.list$int.coef <- as.data.frame(matrix(NA, ncol = 4, nrow = n.pixels))
    colnames(out.list$int.coef) <- c("Est", "SE", "t","p.t")
  }
  if (ret_MSE) {out.list$MSE = vector("numeric", n.pixels)}
  if (ret_resid) {out.list$residuals = matrix(NA, ncol = n.time - 1, nrow = n.pixels)}

  # build the tables
  for (i in 1:n.pixels){
    out.list$time.coef[i, ] <- unlist(coef.list[[i]]["t.j", ])
    if (ret_xi.coef) {out.list$xi.coef[i, ] <- unlist(coef.list[[i]]["x.ti", ])}
    if (ret_int.coef) {out.list$int.coef[i, ] <- unlist(coef.list[[i]]["(Intercept)", ])}
    if (ret_MSE) {out.list$MSE[i] <- sigma(get_fm(cls.list[[i]]))^2}
    if (ret_resid) {out.list$residuals[i, ] <- residuals(get_fm(cls.list[[i]]))}
  }


  class(out.list) <- c("remoteCLS", "map")

  return(out.list)
}

# # cls_star DEPRICATED ----
# #' CLS regression of remote sensing data
# #'
# #' @description
# #' Fit a constrained linear regression model to remote sensing data `X`, testing
# #' the effects of time `t`, `x_{i, t-1}`, and, optionally `z_{i, t-1}` on
# #' `x_{i, t}`.
# #'
# #' @details
# #' This function is a wrapper that calls the single-site function `fitCLS()`
# #' for each site. At present, matrices `X` and `Z` must not contain any `NA`s.
# #'
# #' @param X n x p numeric matrix (usually of of remote sensing observations)
# #'  taken from n sites (rows) and p time points (columns).
# #' @param t numeric vector of length p containing the values for time. Recommended:
# #' `scale(1:ncol(X))`.
# #' @param Z.list Optional list with q elements. Each element of this list
# #' should be an n x p numeric matrix (`z_{i, t}`) with values corresponding to a
# #' covariate measured at each site:time combination.
# #'
# #' @return an n x 8 matrix with the following columns:
# #' `site` the site number from 1-n
# #' `mean` the site average of `x_{i}` (i.e. `rowMeans(X)`)
# #' `Est` the coefficient estimate for the effect of time
# #' `SE` standard error of the coefficient estimate
# #' `t` student's t-test statistic for the effect of time
# #' `p` p-value of the two-tailed t test for the effect of time
# #' `MSE` site-specific mean squared error
# #' `x_t0.EST` the coefficient estimate for the effect of x_{i, t-1} on x_{i, t}
# #'
# #' @export
# #'
# #' @examples #TBA
# #' # starting parameters
# #' n = 25; p = 20; beta = .2
# #'
# #' # simulate data
# #' time = scale(1:p)
# #' X = matrix(rnorm(n*p), ncol = p) # X without time effect
# #' # effect of time from 1 - p
# #' betaT = matrix(rep(time, each = n) * beta, ncol = p)
# #'
# #' # fit the model
# #' cls_star(X*betaT, time)
# cls_star <- function(X, t, Z.list = NULL){
#   # data to fill
#   tmp <- data.frame(pixel = 1:nrow(X), mean = rowMeans(X),
#                     Est = NA, SE = NA, t = NA, p = NA, MSE = NA,
#                     x_t0.EST = NA)
#
#   for (i in 1:nrow(X)) {
#     if (is.null(Z.list) | missing(Z.list)) {
#       fit <- fitCLS(X[i, ], t, Z = NULL)
#     } else {
#       stopifnot(sapply(Z.list, function(Y)dim(Y)==dim(X)))
#       z <-  as.matrix(sapply(Z.list, function(Y){Y[i, ]}))
#       fit <- fitCLS(X[i, ], t, z)
#     }
#     # d[i, "mean"] <- mean(X[i, ])
#     tmp[i, 3:6] <-  fit$coef[3, ]
#     tmp[i, "MSE"] <- fit$MSE
#     tmp[i, "x_t0.EST"] <- fit$coef[2, 1]
#   }
#   return(tmp)
# }
