
## C++ Versions ----
#' Fit GLS to remote sensing data
#' @rdname fitGLS
#'
#' @param X model matrix (double)
#' @param V varcov matrix (double)
#' @param y response vector (double)
#' @param X0 null model matrix (double)
#' @param nugget nugget added to V (double)
#' @param save_xx logical: should xx, xx0, and invcholV be returned? This
#' functionality is meant for use with the partitioned GLS whereby these
#' values are used to calculate cross-partition statistics.
#' @param threads integer indicating the number of threads to use. Currently
#' this parameter does nothin but multi-core functionality will be added soon.
#'
#' @return remoteGLS object
#'
#' @details
#'
#' @examples
#'
#' @export
fitGLS <- function(X, V, y, X0, nugget = 0, save_xx = FALSE, threads = 1){


  ## coerce to matrices
  X = as(X, "matrix")
  V = as(V, "matrix")
  y = as(y, "matrix")
  X0 = as(X0, "matrix")

  ## error handling
  stopifnot(all(is.double(X), is.double(V), is.double(y), is.double(X0)))
  stopifnot(all.equal(nrow(X), nrow(V), nrow(X0), nrow(y)))
  stopifnot(all(check_posdef(V)))


  ## call the c++ function
  # return(.fitGLS_cpp) # contains additional function call
  out <- .Call(`_remotePARTS_fitGLS_cpp`, X, V, y, X0, nugget, save_xx, threads)
  ## add in p values
  out$pval.t <- 2 * pt(abs(out$tstat), df = out$dft, lower.tail = F)
  out$pval.F <- pf(out$Fstat, df1 = out$df.F[1], df2 = out$df.F[2], lower.tail = F)
  class(out) <- "remoteGLS"
  GLS$model.info$call <- match.call()
  GLS$nugget = nugget

  return(out)
}


#' Alternative fitGLS function
#' @rdname fitGLS
#'
#' @param formula formula to build the model with
#' @param data object containing the data
#' @param form.0 null model formula (default: "y ~ 1")
#' @param contrasts optional linear contrasts to use
#' @param ... additional arguments passed to \code{\link{optimize_nugget}}
#'
#' @details \code{fitGLS2()} first creates an empty remoteGLS object
#' (with \code{remoteGLS()}) and then fills in the elements by modifying
#' the remoteGLS object in the C++ function.
#'
#' After further testing to determine which function is most efficient,
#' either \code{fitGLS()} or \code{fitGLS2()} will be deprecated in favor
#' of the other. Though, the formula notation of \code{fitGLS2()} will likely
#' be kept in either case.
#'
#' \code{fitGLS2()} uses a slightly different C++ function from \code{fitGLS()}
#' to achieve the different behavior.
#'
#' @export
#'
#' @examples
fitGLS2 <- function(formula, data, V, nugget = 0, form.0 = NULL,save_xx = FALSE,
                    threads = 1, contrasts = NULL, ...){

  ## Parse formula arguments to make model matrix ----
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
  X <- model.matrix(mt, mf, contrasts) # create model matrix
  rm(mf) # delete the large model frame from memory

  ## Handle missing nugget (NULL or NA) ----
  if (is.null(nugget) || is.na(nugget)){
    nugget = optimize_nugget(X, V, y, ...)
  }

  ## Build null model----
  if (is.null(form.0)){
    form.0 = formula(y ~ 1)
  } else {
    form.0 = formula(form.0)
  }
  # conditionally assign X0
  X0 <- if (missing(data) || is.null(data))
    model.matrix(form.0)
  else
    model.matrix(form.0, data)

  ## coerce to matrices ----
  # X = as(X, "matrix")
  # V = as(V, "matrix")
  # # y = as(y, "matrix")
  # X0 = as(X0, "matrix")

  ## error handling ----
  stopifnot(all(is.double(X), is.double(V), is.double(y), is.double(X0)))
  stopifnot(all.equal(nrow(X), nrow(V), nrow(X0), ny))
  stopifnot(all(check_posdef(V)))

  GLS <- remoteGLS(form = formula)
  GLS$model.info$call <- call
  GLS$nugget = nugget

  ## Run GLS ----
  .Call(`_remotePARTS_fitGLS2_cpp`, GLS, X, V, y, X0, nugget, save_xx, threads)

  # add in p values
  GLS$pval.t <- 2 * pt(abs(GLS$tstat), df = GLS$dft, lower.tail = F)
  GLS$pval.F <- pf(GLS$Fstat, df1 = GLS$df.F[1], df2 = GLS$df.F[2], lower.tail = F)

  ## Return ----
  return(GLS)


  # return(list(test.out = list(X = X, y = y, X0 = X0),
  #             out = out))

}


#' Caculate log-liklihood of GLS model
#' @rdname fitGLS
#'
#' @details \code{LogLikGLS()} returns only the log-likelihood and is used to
#' by other functions for optimization but should be deprecated and simply
#' added as functionality to \code{fitGLS()} and/or \code{fitGLS2()} in future
#' implementations.
#'
#' @examples #TBA
LogLikGLS <- function(nugget, X, V, y){
  ## coerce to matrices
  X = as(X, "matrix")
  V = as(V, "matrix")
  y = as(y, "matrix")

  ## error handling
  stopifnot(all(is.double(X), is.double(V), is.double(y)))
  stopifnot(all.equal(nrow(X), nrow(V), nrow(y)))
  stopifnot(all(check_posdef(V)))

  return(.Call(`_remotePARTS_LogLikGLS_cpp`, nugget, X, V, y))
}

## R Versions ----
#' @rdname fitGLS
#'
#' @details \code{fitGLS()} is a C++ implementation of \code{fitGLS_R()}.
#'
#' @export
fitGLS_R <- function(X, V, y, X0 = NULL, nugget = 0){
  stopifnot(all.equal(nrow(X), ncol(V), nrow(V), length(y)))
  n <- nrow(X)
  invcholV <- invert_cholR(V, nugget)
  xx <- invcholV %*% X
  yy <- invcholV %*% y
  varX <- crossprod(xx)
  beta <- as.numeric(solve(varX, crossprod(xx, yy)))
  SSE <- as.numeric(crossprod(yy - xx %*% beta))
  MSE <- SSE/(n - ncol(xx))
  varcov <- MSE * solve(varX)
  se <- diag(varcov)^0.5
  t <- beta/se
  df.t <- n - ncol(xx)
  p.t <- 2 * pt(abs(t), df = df.t, lower.tail = F) # compute in R not C++

  logdetV <- -2 * sum(log(diag(invcholV)))
  logLik <- -0.5 * (n * log(2 * pi) + n * log((n-ncol(xx))*MSE/n) + logdetV + n)

  if(is.null(X0) | missing(X0)){
    X0 <- matrix(1, nrow = n, ncol = 1)
  }
  xx0 <- invcholV %*% X0
  varX0 <- crossprod(xx0)
  betahat0 <- as.numeric(solve(varX0, crossprod(xx0, yy)))
  SSE0 <- as.numeric(crossprod(yy - xx0 %*% betahat0))
  df0 <- length(betahat0)
  MSE0 <- SSE0/(n - ncol(xx))
  MSR <- (SSE0 - SSE)/(ncol(xx) - ncol(xx0))
  logLik0 <- -0.5 * (n * log(2 * pi) + n * log((n-df0)*MSE0/n) + logdetV + n)

  varX0 <- t(xx0) %*% xx0
  varcov0 <- MSE0 * solve(varX0)
  se0 <- diag(varcov0)^0.5

  if (ncol(xx) > 1) { # when/why would ncol(xx) == 1 ever ?
    FF <- (n - ncol(xx))/(ncol(xx) - ncol(xx0)) * (SSE0 - SSE)/SSE
    df1.F <- ncol(xx) - ncol(xx0)
    df2.F <- n - ncol(xx)
    p.F <- pf(FF, df1 = df1.F, df2 = df2.F, lower.tail = F)
    df.F <- c(df1.F, df2.F)
  } else {
    FF <- (n - 1) * (SSE0 - SSE)/SSE
    df1.F <- 1
    df2.F <- n - 1
    p.F <- pf(FF, df1 = df1.F, df2 = df2.F, lower.tail = F)
    df.F <- c(df1.F, df2.F)
  }

  return(list(betahat = beta, VarX = varX, SSE = SSE, MSE = MSE, varcov = varcov,
              SE = se, tstat = t, pval.t = p.t, dft = df.t, logDetV = logdetV,
              logLik = logLik, betahat0 = betahat0, SE0 = se0, SSE0 = SSE0,
              SSR = SSE0 - SSE, MSE0 = MSE0, MSR = MSR, df0 = df0,
              logLik0 = logLik0, Fstat = FF, pval.F = p.F, df.F = df.F
  ))
}
