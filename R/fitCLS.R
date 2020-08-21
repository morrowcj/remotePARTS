# fitCLS ----
#' Fit a CLS regression to a time series
#'
#' @description see `?cls_star()`
#'
#' @param x numeric vector of length p containing remote sensing data
#' @param t numeric vector of length p containing the values for time. Recommended:
#' `scale(1:ncol(X))`
#' @param Z (optional) a numeric vector/matrix with p columns and q rows. Each
#' row contains a different covariate measured at the same site as x at the same
#' p time points.
#'
#' @return a list of 3 or 4 elements:
#' `coef` the model coefficients from `summary(lm(...))$coef`
#' `MSE` means squared error of the regression model
#' `resids` the residual errors of the model
#' `Z` the coefficients for optional covariates
#'
#' @export
#'
#' @examples

fitCLS <- function(x, t, Z = NULL) {
  # variables
  t_n = length(t)

  # covar handling
  if (is.null(Z) | missing(Z)) {
    Zs = NA
  } else {
    Zs = as.matrix(Z)
    stopifnot(dim(Z)[1] == t_n)
    Zs = Zs[2:t_n, ]
  }

  # data frame for the model
  tmp <- data.frame(
    x_t = x[2:t_n], #X_{t}
    x_t0 = x[1:(t_n - 1)], #X_{t-1}
    time = t[2:t_n], #T_{t}
    Z = Zs)

  # fit the model
  if(is.null(Z) | missing(Z)){
    fm <- lm(x_t ~ x_t0 + time, data = tmp)
  } else {
    fm <- lm(x_t ~ x_t0 + time + Z, data = tmp)
  }

  out <- list(coef = summary(fm)$coef, MSE = summary(fm)$sigma^2,
              resids = resid(fm))

  # add Z to output
  if (!is.null(Z) & !missing(Z)){
    out$Z <- summary(fm)$coef[-c(1:3), ]
  } else {out$Z <-  NULL}

  return(out)
}

# cls_star ----
#' CLS regression of remote sensing data
#'
#' @description
#' Fit a constrained linear regression model to remote sensing data `X`, testing
#' the effects of time `t`, `x_{i, t-1}`, and, optionally `z_{i, t-1}` on
#' `x_{i, t}`.
#'
#' @details
#' This function is a wrapper that calls the single-site function `fitCLS()`
#' for each site. At present, matrices `X` and `Z` must not contain any `NA`s.
#'
#' @param X n x p numeric matrix (usually of of remote sensing observations)
#'  taken from n sites (rows) and p time points (columns).
#' @param t numeric vector of length p containing the values for time. Recommended:
#' `scale(1:ncol(X))`.
#' @param Z.list Optional list with q elements. Each element of this list
#' should be an n x p numeric matrix (`z_{i, t}`) with values corresponding to a
#' covariate measured at each site:time combination.
#'
#' @return an n x 8 matrix with the following columns:
#' `site` the site number from 1-n
#' `mean` the site average of `x_{i}` (i.e. `rowMeans(X)`)
#' `Est` the coefficient estimate for the effect of time
#' `SE` standard error of the coefficient estimate
#' `t` student's t-test statistic for the effect of time
#' `p` p-value of the two-tailed t test for the effect of time
#' `MSE` site-specific mean squared error
#' `x_t0.EST` the coefficient estimate for the effect of x_{i, t-1} on x_{i, t}
#'
#' @export
#'
#' @examples
#' # starting parameters
#' n = 25; p = 20; beta = .2
#'
#' # simulate data
#' time = scale(1:p)
#' X = matrix(rnorm(n*p), ncol = p) # X without time effect
#' # effect of time from 1 - p
#' betaT = matrix(rep(time, each = n) * beta, ncol = p)
#'
#' # fit the model
#' cls_star(X*betaT, time)

cls_star <- function(X, t, Z.list = NULL){
  # data to fill
  tmp <- data.frame(site = 1:nrow(X), mean = rowMeans(X),
                  Est = NA, SE = NA, t = NA, p = NA, MSE = NA,
                  x_t0.EST = NA)

  for (i in 1:nrow(X)) {
    if (is.null(Z.list) | missing(Z.list)) {
      fit <- fitCLS(X[i, ], t, Z = NULL)
    } else {
      stopifnot(sapply(Z.list, function(Y)dim(Y)==dim(X)))
      z <-  as.matrix(dplyr::bind_cols(lapply(Z.list, function(Y){Y[i, ]})))
      fit <- fitCLS(X[i, ], t, z)
    }
    # d[i, "mean"] <- mean(X[i, ])
    tmp[i, 3:6] <-  fit$coef[3, ]
    tmp[i, "MSE"] <- fit$MSE
    tmp[i, "x_t0.EST"] <- fit$coef[2, 1]
  }
  return(tmp)
}
