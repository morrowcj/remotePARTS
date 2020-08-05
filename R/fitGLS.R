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

fitGLS <- function(X, V, y, X0 = NULL){
  stopifnot(all.equal(nrow(X), ncol(V), nrow(V), length(y)))
  n <- nrow(X)
  invcholV <- invert_cholR(V)
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
              logLik = logLik, betahat0 = betahat0, SSE0 = SSE0, MSE0 = MSE0,
              MSR = MSR, df0 = df0, logLik0 = logLik0, Fstat = FF, pval.F = p.F,
              df.F = df.F
              ))
}

# Rcpp::sourceCpp("Cpp/GLS_Chol.cpp")
## This function can be depracated if I implement t-test and F-test in C++
wrap_glscpp <- function(X, V, y, X0 = NULL){
  if(is.null(X0) | missing(X0)){
    X0 <- matrix(as.double(1), nrow = nrow(X), ncol = 1)
  }
  out <- fitGLS_cpp(X = X, V = V, y = y, X0 = X0)

  out$pval.t <- 2 * pt(abs(out$tstat), df = out$dft, lower.tail = F)
  out$pval.F <- pf(out$Fstat, df1 = out$df.F[1], df2 = out$df.F[2], lower.tail = F)
}

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

nugget.fit <- function(formula, data, V, nugget.tol = 0.00001,
                       interval = c(0, 1), verbose = FALSE) {
  opt.nugget <- optimize(nugget.fit.funct, formula, data = data, V = V, interval = interval, maximum = T, tol = nugget.tol, verbose = verbose)
  # check at the zero boundary
  if(opt.nugget$maximum < nugget.tol){
    nugget0.fit <- nugget.fit.funct(0, formula, data, V)
    if(nugget0.fit > opt.nugget$objective) opt.nugget$maximum <- 0
  }
  return(opt.nugget$maximum)
}

GLS.fit <- function(formula, formula0 = NULL, data, V = NULL, invcholV = NULL,
                    save.invcholV = F) {

  mf <- model.frame(formula = formula, data = data) # model frame from data
  x <- model.matrix(attr(mf, "terms"), data = mf) # model matrix X from data
  y <- model.response(mf) # model response
  n <- length(y) # n observations

  # compute the inverse of the chol decomp: t(Uinv)
  if (is.null(invcholV)) {
    if (is.null(V)) {
      invcholV <- diag(n)
    } else {
      invcholV <- t(backsolve(chol(V), diag(n))) #t(Uinv)
    }
  }

  xx <- invcholV %*% x #t(Uinv) %*% X
  yy <- invcholV %*% y #t(Uinv) %*% y

  # betahat
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

GLS.partition.data <- function(formula, formula0 = NULL, data,
                               spatial.autocor.FUN = "exponential-power",
                               spatialcor = spatialcor, est.nugget = T,
                               npart = 10, partition = NULL,
                               nugget.interval = c(0,1), fixed.nugget = NULL,
                               nugget.tol = 0.00001, min.num.cross.part = 5,
                               verbose = F, rm.spatial.autocorrelation = F) {

  n <- nrow(data)
  if (!is.null(partition)) {
    npart <- nrow(partition)
    nn <- n - (n%%npart)
    n.p <- nn/npart
    pick <- partition
  } else {
    nn <- n - (n%%npart)
    n.p <- nn/npart
    pick <- matrix(sample(n)[1:nn], nrow = npart)
  }

  mf <- model.frame(formula = formula, data = data)
  df2 <- n.p - (ncol(model.matrix(attr(mf, "terms"), data = mf)) - 1)
  mf0 <- model.frame(formula = formula0, data = data)
  df0 <- n.p - (ncol(model.matrix(attr(mf0, "terms"), data = mf0)) - 1)
  df1 <- df0 - df2

  if(min.num.cross.part > npart) min.num.cross.part <- npart
  if(ncol(mf) <= 1){
    min.num.cross.part <- NA
  }

  SSR.part <- NULL
  SSE.part <- NULL
  SSE0.part <- NULL
  coef.part <- NULL
  coef0.part <- NULL
  se.part <- NULL
  se0.part <- NULL
  F.part <- NULL
  p.F.part <- NULL
  logLik.part <- NULL
  logLik0.part <- NULL
  nugget.part <- NULL
  invcholV.part <- list(NULL)
  xx.part <- list(NULL)
  xx0.part <- list(NULL)
  for (i in 1:npart) {

    data.part <- data[pick[i,],]

    if(rm.spatial.autocorrelation == F){
      # create distance matrix in kilometers
      location <- data.part[,c('lng','lat')]
      Dist.part <- geosphere::distm(location, fun=distGeo)/1000

      Vp <- V.fit(Dist.part, spatialcor = spatialcor, FUN = spatial.autocor.FUN)
      if (is.null(fixed.nugget) & est.nugget) {
        nugget <- nugget.fit(formula, data.part, Vp, interval = nugget.interval,
                             verbose = verbose)
        nugget.interval <- c(0, max(1000*nugget.tol, min(100*nugget,1)))
        Vp <- (1 - nugget) * Vp + nugget * diag(n.p)
      } else {
        if (is.null(fixed.nugget)) {
          nugget <- 0
          Vp <- Vp
        } else {
          nugget <- fixed.nugget[i]
          Vp <- (1 - nugget) * Vp + nugget * diag(n.p)
        }
      }
    }else{
      Vp <- diag(n.p)
    }
    invcholV <- t(backsolve(chol(Vp), diag(n.p)))
    z.part <- GLS.fit(formula, formula0, data = data.part, invcholV = invcholV,
                      save.invcholV = F)

    SSR.part <- c(SSR.part, z.part$SSR)
    SSE.part <- c(SSE.part, z.part$SSE)
    SSE0.part <- c(SSE0.part, z.part$SSE0)
    coef.part <- cbind(coef.part, z.part$coef)
    coef0.part <- cbind(coef0.part, z.part$coef0)
    se.part <- cbind(se0.part, z.part$se)
    se0.part <- cbind(se.part, z.part$se0)
    F.part <-  c(F.part, z.part$F)
    p.F.part <-  c(p.F.part, z.part$p.F)
    logLik.part <- c(logLik.part, z.part$logLik)
    logLik0.part <- c(logLik0.part, z.part$logLik0)
    nugget.part <- c(nugget.part, nugget)

    if(!is.na(min.num.cross.part) && i <= min.num.cross.part){
      invcholV.part[[i]] <- invcholV
      xx.part[[i]] <- z.part$xx
      xx0.part[[i]] <- z.part$xx0
    }
    if(verbose) {
      show(paste0("partition ",i," of ", npart))
      show(z.part$coef)
      show(z.part$p.F)
    }
  }

  if(!is.na(min.num.cross.part)){
    rSSE.part <- matrix(NA, nrow=npart, ncol=npart)
    rSSR.part <- matrix(NA, nrow=npart, ncol=npart)
    for (i in 1:(min.num.cross.part-1)) for (j in (i+1):min.num.cross.part) {
      data.part <- data[c(pick[i,], pick[j,]),]

      # create distance matrix in kilometers
      location <- data.part[,c('lng','lat')]
      Dist.part <- geosphere::distm(location, fun=distGeo)/1000
      Vpick <- V.fit(Dist.part, spatialcor = spatialcor,
                     FUN = spatial.autocor.FUN)
      Vnugget <- diag(c(rep((1-nugget.part[i])/nugget.part[i], n.p),
                        rep((1-nugget.part[j])/nugget.part[j], n.p)))
      Vnugget[is.infinite(Vnugget)] <- 0
      Vpick <- Vpick + Vnugget

      xx1 <- xx.part[[i]]
      xx2 <- xx.part[[j]]
      xx10 <- xx0.part[[i]]
      xx20 <- xx0.part[[j]]

      Rij <- crossprod(t(invcholV.part[[i]]),
                       tcrossprod(Vpick[1:n.p, (n.p+1):(2*n.p)],
                                  invcholV.part[[j]]))
      H1 <- xx1 %*% solve(t(xx1) %*% xx1) %*% t(xx1)
      H2 <- xx2 %*% solve(t(xx2) %*% xx2) %*% t(xx2)

      if(!is.na(xx10[1])){
        H10 <- xx10 %*% solve(t(xx10) %*% xx10) %*% t(xx10)
        H20 <- xx20 %*% solve(t(xx20) %*% xx20) %*% t(xx20)
      }else{
        H10 <- 0
        H20 <- 0
      }

      S1R <- H1 - H10
      S2R <- H2 - H20

      S1E <- diag(n.p) - H1
      S2E <- diag(n.p) - H2

      rSSR.part[i,j] <- matrix(S1R, nrow=1) %*% matrix(Rij %*% S2R %*%
                                                         t(Rij), ncol=1)/df1
      rSSE.part[i,j] <- matrix(S1E, nrow=1) %*% matrix(Rij %*% S2E %*%
                                                         t(Rij), ncol=1)/df2
    }

    rSSR <- mean(rSSR.part, na.rm=T)
    rSSE <- mean(rSSE.part, na.rm=T)

    Fmean <- mean(F.part)
  }else{
    rSSR.part <- NA
    rSSE.part <- NA

    rSSR <- NA
    rSSE <- NA

    Fmean <- NA
  }

  coef <- rowMeans(coef.part)
  coef0 <- rowMeans(coef0.part)

  return(list(coef = coef, Fmean = Fmean, df1 = df1, df2 = df2,
              SSR.part = SSR.part, SSE.part = SSE.part, SSE0.part = SSE0.part,
              logLik.part = logLik.part, logLik0.part = logLik0.part,
              nugget = mean(nugget.part), nugget.part = nugget.part,
              F.part = F.part, p.F.part = p.F.part, coef.part=coef.part,
              se.part=se.part, coef0.part=coef0.part, se0.part=se0.part,
              rSSR = rSSR, rSSE = rSSE, rSSR.part = rSSR.part,
              rSSE.part = rSSE.part, npart = npart, partition = pick,
              spatial.autocor.FUN = "exponential-power",
              spatialcor = spatialcor))
}

