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
#'
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
#' @export
#'
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


fitAR.map <- function(X, t, Z = NULL){
  print("TBA")
}
