## AR function ----
#' restricted maximum likelihood of an AR model
#'
#' @param par AR parameter value
#' @param x vector of time series (response)
#' @param U model matrix (predictors)
#'
#' @return
#'
#' if \code{LL.only = TRUE}: returns log-likelihood of \code{par}
#' given \code{x} and \code{u}.
#'
#' otherwise, returns a "remoteAR.pixel" object which is a list of a
#' regression coefficeint table (\code{$coef}), the REML AR parameter (\code{$b}),
#' model MSE (\code{$MSE}),
#' estimated coefficeint covariance matrix (\code{$s2beta}),
#' and the model log-likelihood (\code{$logLik}).
#'
#' @details used by \code{fitAR()}
#'
#' @export
#' @examples
#' time = 1:30
#' x = rnorm(31)
#' x = x[2:31] + x[1:30] + 0.3*time #AR(1) process + time trend
#' U = model.matrix(formula(x ~ time))
#' AR_funct(par = .2, x, U, LL.only = TRUE)
#' AR_funct(par = .2, x, U, LL.only = FALSE)
AR_funct <- function(par, x, U, LL.only = TRUE) {
  b <- par # parameter of interest
  n.obs <- length(x) # number of time points
  q <- ncol(U) # number of covariates in model matrix
  B <- diag(n.obs) # n length identity matrix
  diag(B[-1, ]) <- -b # set sub-diagonal to negative b
  iS <- diag(n.obs) # another identity matrix
  iS[1, 1] <- (1 - b^2) # set first element to 1-b^2
  iV <- t(B) %*% iS %*% B # matrix multiply B'*iS*B

  # handle missing/NA values by removing all NA elements from each object
  if(any(is.na(x))){
    iV <- iV[!is.na(x), !is.na(x)]
    U <- U[!is.na(x),]
    x <- x[!is.na(x)]
  }

  # log determinant of iV
  logdetV <- -determinant(iV)$modulus[1]

  # solve Am + B for m where A = (U'*iV*U) and B = (U'*iV*x)
  ## beta is the effect of the covariates
  beta <- solve(t(U) %*% iV %*% U, t(U) %*% iV %*% x)
  # remove the effect of the covariates from x
  H <- x - U %*% beta
  # estimate the variance
  s2 <- (t(H) %*% iV %*% H)/(n.obs - q)
  # calculate the log-likelihood of b given x and U
  LL <- 0.5 * ((n.obs - q) * log(s2) + logdetV +
                 determinant(t(U) %*% iV %*% U)$modulus[1] +
                 (n.obs - q))
  #show(c(LL,b))

  if(LL.only){
  # return log-likelihood
  return(as.vector(LL))
  } else {
    MSE <- as.numeric(s2) #MSE
    s2beta <- MSE * solve(t(U) %*% iV %*% U) #SE
    t.stat = (abs(beta) / diag(s2beta)^0.5)
    pval = 2 * pt(q = t.stat, df = n.obs - q,
                  lower.tail = FALSE )

    ## log likelihood without constants (i.e. s2) - no parameter dependancy
    logLik <- 0.5 * (n.obs - q) * log(2 * pi) +
      determinant(t(U) %*% U)$modulus[1] - LL

    coef.tab <- data.frame("Est" = beta,
                      "SE" = diag(s2beta),
                      "t.stat" = t.stat,
                      "p.val" = pval)

    out.list = list(call = match.call(),
                    coef = coef.tab,
                    b = b,
                    MSE = MSE,
                    s2beta = s2beta,
                    resids = as.vector(H),
                    logLik = as.vector(logLik))

    class(out.list) <- c("remoteAR", "pixel")

    return(out.list)
  }
}

## AR Wrapper----
#' Fit an Auto-regressive time series analysis using restricted maximum
#' likelihood
#'
#' @param formula model formula
#' @param data object in which the model will first look for data
#'
#' @return a "remoteAR.pixel" object. See [AR_funct()] for more details.
#'
#' @details [fitAR()] is a wrapper function for [AR_funct()].
#'
#' By default, the print.remoteAR() method does not show all output.
#' to access individual components, use \code{names()} to see element names
#' and the S3 \code{$} operator to access them.
#'
#' @export
#' @examples
#' time = 1:30
#' x = rnorm(31)
#' x = x[2:31] + x[1:30] + 0.3*time #AR(1) process + time trend
#' fitAR(x ~ time)
fitAR <- function(formula, data){
  ## Produce model matrix and model frame from call ----
  call <- match.call() # function call
  mf <- match.call(expand.dots = FALSE) # don't expand ...
  m <- match(c("formula", "data", "subset",
               "weights", "na.action", "offset"),
             names(mf), 0L) # match arguments provided by call
  mf <- mf[c(1L, m)] #function name, plus arguments matched
  mf$drop.unused.levels <- TRUE # show that we dropped levels
  mf[[1L]] <- quote(stats::model.frame) # rename the function call
  mf <- eval(mf, parent.frame()) # evaluate the model frame with the data
  mt <- attr(mf, "terms") # model terms
  y <- model.response(mf, "numeric") # response (vector)
  # w <- as.vector(model.weights(mf)) # model.weights
  # offset <- model.offset(mf) # model offset
  if (is.matrix(y)){stop("response is a matrix: must be a vector")}
  ny <- length(y)
  U <- model.matrix(mt, mf, contrasts = NULL) # create model matrix

  ## Optimize AR_funct() for par ----
  opt <- optim(fn = AR_funct, par = 0.2, x = y, U = U, LL.only = TRUE,
               method = "Brent", upper = 1, lower = -1,
               control = list(maxit = 10^4))

  b <- opt$par # optimized parameter

  ## Perform the AR regression with the optimized parameter ----
  AR.out = AR_funct(par = b, x = y, U = U, LL.only = FALSE)
  AR.out$call = call

  return(AR.out)
}

#' Fit AR REML models to a time series matrix
#'
#' @param X nxp time series response matrix with p columns corresponding to time
#'  points and n columns corresponding to the number of pixels
#' @param t p length temporal response vector
#' @param Z
#' @param ret_int.coef should the intercept coeffients be returned? logical
#' @param ret_AR.par should the AR parameter estimates be returned? logical
#' @param ret_MSE should the model MSEs be returned? logical
#' @param ret_resid should the model residuals be returned? logical
#' @param ret_logLik should the model log-likelihoods be returned? logical
#'
#' @return a list with the following elements: the initial function call
#' (\code{$call}), a coefficient matrix for the temporal variable
#' (\code{$time.coef}), an optional vector of the AR parameters (\code{$AR.par}),
#' an optional vector of MSEs (\code{$MSE}), an optional vector of
#' log-likelihoods (\code{$logLik}), and an optional nxp matrix of model
#' residuals (\code{$resids}).
#'
#' @details by default the print.remoteAR() method does not show all output.
#' to access individual components, use \code{names()} to see element names
#' and the S3 \code{$} operator to access them.
#'
#' @export
#'
#' @examples
#' t = 1:30; n.pix = 10
#' X = matrix(rnorm(length(t)*n.pix), ncol = length(t))
#' fitAR.map(X, t) # only $call and $time.coef are printed by print.remoteAR()
#' summary(fitAR.map(X, t))
#'
#' coef(fitAR.map(X, t)) # data frame of time coefficeints (alternatively fitAR.map(X, t)$time.coef)
#' fitAR.map(X, t)$AR.par # AR parameters
#' fitAR.map(X, t)$MSE # model MSEs
#' fitAR.map(X, t)$logLik # model log-likelihoods
#' resid(fitAR.map(X, t)) # matrix of model residuals (alternatively fitAR.map(X, t)$resids)
fitAR.map <- function(X, t, Z = NULL,
                      ret_int.coef = FALSE, ret_AR.par = TRUE,
                      ret_MSE = TRUE, ret_resid = TRUE, ret_logLik = TRUE){
  stopifnot(ncol(X) == length(t))
  n.pixels = nrow(X)
  n.time = length(t)

  ## Run fitAR
  if (!is.null(Z)){
    message("handling of Z not yet implemented")
  }
  AR.list = lapply(1:n.pixels, function(x){
    y = X[x, ]
    return(fitAR(y ~ t))
  })

  ## Extract coefficients
  coef.list <- lapply(AR.list, coef)

  ## Initialize output
  out.list <- list(call = match.call(),
                   time.coef = matrix(NA, ncol = 4, nrow = n.pixels))
  colnames(out.list$time.coef) <- c("Est", "SE", "t.stat", "p.val")
  if (ret_int.coef) {
    out.list$int.coef = matrix(NA, ncol = 4, nrow = n.pixels)
    colnames(out.list$int.coef) <- c("Est", "SE", "t.stat", "p.val")
  }
  if (ret_AR.par) {out.list$AR.par = numeric(n.pixels)}
  if (ret_MSE) {out.list$MSE = numeric(n.pixels)}
  if (ret_logLik) {out.list$logLik = numeric(n.pixels)}
  if (ret_resid) {out.list$resids = matrix(NA, ncol = n.time, nrow = n.pixels)}

  ## Fill output
  for (i in 1:n.pixels) {
    out.list$time.coef[i, ] <- unlist(coef.list[[i]]["t", ])
    if (ret_int.coef) {
      out.list$int.coef[i, ] <- unlist(coef.list[[i]]["(Intercept)", ])
    }
    if (ret_AR.par) {out.list$AR.par[i] = AR.list[[i]]$b}
    if (ret_MSE) {out.list$MSE[i] = AR.list[[i]]$MSE}
    if (ret_logLik) {out.list$logLik[i] = AR.list[[i]]$logLik}
    if (ret_resid) {out.list$resids[i, ] = AR.list[[i]]$resids}
  }

  class(out.list) <- c("remoteAR", "map")

  return(out.list)
}
