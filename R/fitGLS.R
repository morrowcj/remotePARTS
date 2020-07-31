## invert_choldec ----

invert_cholR <- function(M, nugget = NULL){
  stopifnot(nrow(M) == ncol(M))

  n = nrow(M)

  # handle nugget
  if(!missing(nugget) & !is.null(nugget)){
    print("using nugget")
    M <- (1 - nugget) * M + nugget * diag(n)
  }

  # return the result
  return(t(backsolve(chol(M), diag(n))))
}


## fitGLS ----
## This function calls invert_choldec()

## fitDistVar ----
fitDistVar_R <- function(Dist, spatialcor, fun = "exponential"){
  switch(fun,
         # exponential (with alias)
         "exponential" = exp(-Dist/spatialcor), ## This version yeilds non positive-definitive matrix!!!
         "exp" = exp(-Dist/spatialcor),
         # exponential power (with alias)
         "exponential-power" = exp(-(Dist/spatialcor[1])^spatialcor[2]),
         "exp-pwr" = exp(-(Dist/spatialcor[1])^spatialcor[2]),
         # taper-spherical (aliases)
         "taper-spherical" = taper_sphere(Dist, spatialcor),
         "taper" = taper_sphere(Dist, spatialcor),
         "sphr" = taper_sphere(Dist, spatialcor))
}

## fitNugget ----
## This function calls fitGLS() and invert_choldec()

## Original Code ----

V.fit <- function(Dist, spatialcor, FUN = "exponential") {

  if (FUN == "exponential")
    return(exp(-Dist/spatialcor))

  if (FUN == "exponential-power")
    return(exp(-(Dist/spatialcor[1])^spatialcor[2]))

  if (FUN == "taper-spherical")
    return(taper.spherical(Dist, spatialcor))

}

nugget.fit.funct <- function(nugget, formula, data, V, verbose = FALSE) {
  n <- ncol(V)
  invcholV <- t(backsolve(chol((1 - nugget) * V + nugget * diag(n)), diag(n)))
  z <- GLS.fit(formula, data = data, invcholV = invcholV)
  if(verbose == TRUE) show(c(z$logLik, nugget))
  return(z$logLik)
}

nugget.fit <- function(formula, data, V, nugget.tol = 0.00001, interval = c(0, 1), verbose = FALSE) {
  opt.nugget <- optimize(nugget.fit.funct, formula, data = data, V = V, interval = interval, maximum = T, tol = nugget.tol, verbose = verbose)
  # check at the zero boundary
  if(opt.nugget$maximum < nugget.tol){
    nugget0.fit <- nugget.fit.funct(0, formula, data, V)
    if(nugget0.fit > opt.nugget$objective) opt.nugget$maximum <- 0
  }
  return(opt.nugget$maximum)
}

GLS.fit <- function(formula, formula0 = NULL, data, V = NULL, invcholV = NULL, save.invcholV = F) {

  mf <- model.frame(formula = formula, data = data)
  x <- model.matrix(attr(mf, "terms"), data = mf)
  y <- model.response(mf)
  n <- length(y)

  if (is.null(invcholV)) {
    if (is.null(V)) {
      invcholV <- diag(n)
    } else {
      invcholV <- t(backsolve(chol(V), diag(n)))
    }
  }

  xx <- invcholV %*% x
  yy <- invcholV %*% y

  coef <- as.numeric(solve(crossprod(xx), crossprod(xx,yy))) # B(XX' %*% XX)' = (XX' %*% YY)
  names(coef) <- colnames(x)
  varX <- t(xx) %*% xx
  SSE <- as.numeric(crossprod(yy - xx %*% coef))
  MSE <- SSE/(n - ncol(xx))

  varcov <- MSE * solve(varX)
  se <- diag(varcov)^0.5
  t <- coef/se
  df.t <- n - ncol(xx)
  p.t <- 2 * pt(abs(t), df = df.t, lower.tail = F)

  logdetV <- -2 * sum(log(diag(invcholV)))
  logLik <- -0.5 * (n * log(2 * pi) + n * log((n-ncol(xx))*MSE/n) + logdetV + n)

  if (is.null(formula0)) {
    x0 <- matrix(1, nrow = n, ncol = 1)
    xx0 <- invcholV %*% x0
    coef0 <- solve(crossprod(xx0), crossprod(xx0, yy))
    SSE0 <- as.numeric(crossprod(yy - xx0 %*% coef0))
    df0 <- ncol(coef0)
  } else {
    mf0 <- model.frame(formula = formula0, data = data)
    x0 <- model.matrix(attr(mf0, "terms"), data = mf0)
    xx0 <- invcholV %*% x0

    if(any(xx0 != 0)){
      coef0 <- solve(crossprod(xx0), crossprod(xx0, yy))
      SSE0 <- as.numeric(crossprod(yy - xx0 %*% coef0))
      df0 <- ncol(coef0)
    }else{
      SSE0 <- as.numeric(crossprod(yy))
      df0 <- 1
      coef0 <- NA
    }
  }
  MSE0 <- SSE0/(n - ncol(xx))
  MSR <- (SSE0 - SSE)/(ncol(xx) - ncol(xx0))
  logLik0 <- -0.5 * (n * log(2 * pi) + n * log((n-df0)*MSE0/n) + logdetV + n)

  if(any(xx0 != 0)){
    varX0 <- t(xx0) %*% xx0
    varcov0 <- MSE0 * solve(varX0)
    se0 <- diag(varcov0)^0.5
  }else{
    varcov0 <- NULL
    se0 <- NULL
  }

  if (ncol(xx) > 1) {
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

  if(!save.invcholV) invcholV <- NULL

  return(list(coef = coef, se = se, t = t, df.t = df.t, p.t = p.t, F = FF, df1.F = df1.F, df2.F = df2.F, p.F = p.F, logLik = logLik, logLik0 = logLik0, MSE = MSE, MSE0 = MSE0, MSR = MSR, SSE = SSE, SSE0 = SSE0, SSR = SSE0 - SSE, coef0 = coef0, se0 = se0, varX = varX, varcov = varcov, varcov0 = varcov0, invcholV = invcholV, xx=xx, xx0=xx0, yy=yy))
}
