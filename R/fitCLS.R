# fit_cls ----
## fit a cls to a timeseries
#' Title
#'
#' @param x
#' @param t
#' @param covars
#'
#' @return
#' @export
#'
#' @examples
fit_cls <- function(x, t, covars = NULL) {
  # variables
  t_n = length(t)

  # covar handling
  if (is.null(covars) | missing(covars)) {
    Z = NA
  } else {
    Z = as.matrix(covars)
    stopifnot(dim(Z)[1] == t_n)
    Z = Z[2:t_n, ]
  }

  # data frame for the model
  tmp <- data.frame(
    x_t = x[2:t_n], #X_{t}
    x_t0 = x[1:(t_n - 1)], #X_{t-1}
    time = t[2:t_n], #T_{t}
    Z = Z)

  # fit the model
  if(is.null(covars) | missing(covars)){
    fm <- lm(x_t ~ x_t0 + time, data = tmp)
  } else {
    fm <- lm(x_t ~ x_t0 + time + Z, data = tmp)
  }

  out <- list(summary = summary(fm)$coef, MSE = summary(fm)$sigma^2,
              resids = resid(fm))

  # add Z to output
  if (!is.null(covars) & !missing(covars)){
    out$Z <- summary(fm)$coef[-c(1:3), ]
  }

  return(out)
}

# cls_star ----
## this is a multi-site wrapper for fit_cls
#' Title
#'
#' @param X
#' @param t
#' @param covar.list
#'
#' @return
#' @export
#'
#' @examples
cls_star <- function(X, t, covar.list = NULL){
  # data to fill
  tmp <- data.frame(site = 1:nrow(X), mean = rowMeans(X),
                  Est = NA, SE = NA, t = NA, p = NA, MSE = NA,
                  x_t0.EST = NA)

  for (i in 1:nrow(X)) {
    if (is.null(covar.list) | missing(covar.list)) {
      fit <- fit_cls(X[i, ], t, covars = NULL)
    } else {
      stopifnot(sapply(covar.list, function(Y)dim(Y)==dim(X)))
      z <-  as.matrix(dplyr::bind_cols(lapply(covar.list, function(Y){Y[i, ]})))
      fit <- fit_cls(X[i, ], t, z)
    }
    # d[i, "mean"] <- mean(X[i, ])
    tmp[i, 3:6] <-  fit$summary[3, ]
    tmp[i, "MSE"] <- fit$MSE
    tmp[i, "x_t0.EST"] <- fit$summary[2, 1]
  }
  return(tmp)
}
