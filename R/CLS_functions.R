# fitCLS ----
#' Fit a CLS regression to a time series
#'
#' @description see `?cls_star()`
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
#' @return list of 3 elements: \code{$call} the function call used to produce
#' the output, \code{$fm} the model object fit using \code{stats::lm()}, and, if
#' if \code{save_AR.df = TRUE}, \code{$AR.df} an AR data frame built with
#' \code{AR_df()} .
#'
#' @seealso [AR_df()] for AR data frames and [fitCLS.map()] for fitting full-map
#' time series
#'
#' @export
#'
#' @examples
#' data(ndvi_AK3000)
#' x = unlist(ndvi_AK3000[1, -c(1:6)]) # time series for first pixel only
#' t = 1:length(x) # time points
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

## Function to replace cls_star ----
#' Title
#'
#' @param X
#' @param t
#' @param ret_xi.coef
#' @param ret_int.coef
#' @param ret_MSE
#' @param ret_resid
#'
#' @return
#' @export
#'
#' @examples
fitCLS.map <- function(X, t, ret_xi.coef = FALSE, ret_int.coef = FALSE,
                       ret_MSE = TRUE, ret_resid = TRUE){
  n.pixels = nrow(X)
  n.time = length(t)
  stopifnot(ncol(X) == n.time)

  ## run CLS
  cls.list <- lapply(1:n.pixels, function(x){
    fitCLS(X[x,], t, save_AR.df = TRUE)
  })

  ## get coefficients from cls
  coef.list <- lapply(cls.list, function(x){
    as.data.frame(summary(get_fm(x))$coefficient)
  })

  ## setup coefficient tables
  out.list <- list(call = match.call(),
                   time.coef = as.data.frame(matrix(NA, ncol = 4, nrow = n.pixels)))
  names(out.list$time.coef) <- c("Est", "SE", "t","p.t")
  out.list$mean = rowMeans(X)
  if (ret_xi.coef){
    out.list$xi.coef <- as.data.frame(matrix(NA, ncol = 4, nrow = n.pixels))
    names(out.list$xi.coef) <- c("Est", "SE", "t","p.t")
  }
  if (ret_int.coef){
    out.list$int.coef <- as.data.frame(matrix(NA, ncol = 4, nrow = n.pixels))
    names(out.list$int.coef) <- c("Est", "SE", "t","p.t")
  }
  if (ret_MSE) {out.list$MSE = vector("numeric", n.pixels)}
  if (ret_resid) {out.list$residuals = matrix(NA, ncol = n.time - 1, nrow = n.pixels)}

  # build the tables
  for (i in 1:n.pixels){
    out.list$time.coef[i, ] <- coef.list[[i]]["t.j", ]
    if (ret_int.coef) {out.list$int.coef[i, ] <- coef.list[[i]]["(Intercept)", ]}
    if (ret_xi.coef) {out.list$xi.coef[i, ] <- coef.list[[i]]["x.ti", ]}
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
