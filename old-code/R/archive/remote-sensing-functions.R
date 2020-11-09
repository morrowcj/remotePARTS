# Remote sensing functions to use with this package
## Authors: Anthony Ives, Clay Morrow

## AR.reml ----
#' @title Perform AR(1) REML on a model with covariates u
#' @description
#' This function performs AR(1) REML using ... [Description]
#' last modified on 01-April-2020
#'
#' @param data
#' @param par
#' @param x
#' @param u
#'
#' @return
#' @export
#'
#' @examples
AR.reml <- function(formula,
                    data = list(),
                    par, # strength of AR(1) correlation [0,1]
                    x, # vector of responses (time-series)
                    u){ # matrix of covariates

  AR.reml.funct <- function(par = par,
                            x = x, #n x 1 response variables
                            u = u){#n x q covariates
    b <- par #b parameter?? Correlation of x_i and x_i+1 (and x_i-1)
    n.obs <- length(x) #n number of observations
    q <- dim(u)[2] #number of covariate paramters
    B <- diag(n.obs) #n x n identity matrix
    diag(B[-1, ]) <- -b #make -b the 1st subdiagonal

    iS <- diag(n.obs) #n x n identity matrix
    iS[1, 1] <- (1 - b^2) #set the upper left value to 1-b^2
    iV <- t(B) %*% iS %*% B #create block-diagonal covariance matrix
    logdetV <- -determinant(iV)$modulus[1] #calculate the negative log determinant modulus

    A <- t(u) %*% iV %*% u # Variance component for u
    B <- t(u) %*% iV %*% x # Variance component for x
    beta <- solve(A, B) #solve 'B = Ax' for x (slope coefficients: effect of U_i on X)
    H <- x - u %*% beta #standardized response: x - effects of u

    s2 <- (t(H) %*% iV %*% H)/(n.obs - q) # Total Variance Component
    LL <- 0.5 * ((n.obs - q) * log(s2) + logdetV + # Log-liklihood? p x p matrix...
                   determinant(t(u) %*% iV %*% u)$modulus[1] + (n.obs - q))
    #show(c(LL,b))
    return(LL)
  }

  mf <- model.frame(formula = formula, data = data)
  u <- model.matrix(attr(mf, "terms"), data = mf) #variables
  x <- model.response(mf) #response variables

  q <- dim(u)[2]
  # define optimizer function using AR.reml.funct
  opt <- optim(fn = AR.reml.funct, par = 0.2, method = "Brent", upper = 1,
               lower = -1, control = list(maxit = 10^4), x = x, u = u)
  b <- opt$par #maximum likelihood (reml) parameter estimate for b.

  # repeat AR.reml.funct steps with this b
  n.obs <- length(x)
  q <- dim(u)[2]
  B <- diag(n.obs)
  diag(B[-1, ]) <- -b

  iS <- diag(n.obs)
  iS[1, 1] <- (1 - b^2)
  iV <- t(B) %*% iS %*% B
  logdetV <- -determinant(iV)$modulus[1]

  beta <- solve(t(u) %*% iV %*% u, t(u) %*% iV %*% x)
  H <- x - u %*% beta

  MSE <- as.numeric((t(H) %*% iV %*% H)/(n.obs - q)) #mean square error
  s2beta <- MSE * solve(t(u) %*% iV %*% u)

  # calculate p values and log likelihoods
  Pr <- 1:q
  for (i in 1:q) Pr[i] <- 2 * pt(abs(beta[i])/s2beta[i, i]^0.5, df = n.obs - q,
                                 lower.tail = F)

  logLik <- 0.5 * (n.obs - q) * log(2 * pi) +
    determinant(t(u) %*% u)$modulus[1] - opt$value

  return(list(beta = beta, #linear regression slope coefficient estimates
              b = b, # AR(1) correlation paramter
              MSE = MSE, #mean square error
              s2beta = s2beta, #variance component of beta
              Pr = Pr, #p values for beta
              logLik = logLik)) #log likelihoods of beta
}

## simX ----
#' @title Simulate a collection of correlated time series
#' @description
#' simX will simulate ... [Description]
#'
#' @param formula
#' @param data
#' @param coef
#' @param b
#' @param s
#' @param Dr
#' @param t.scale
#' @param n
#' @param n.obs
#' @param n.burn
#' @param seed
#'
#' @return
#' @export
#'
#' @examples
simX <- function(formula,
                 data = data.frame(rep(1, n)),
                 coef,
                 b,
                 s,
                 Dr = NULL,
                 t.scale, #time scale
                 n, #number of independent X's to simulate
                 n.obs, #number of observations
                 n.burn, #burn-in iterations
                 seed = 0){
  set.seed(seed=seed)
  mf <- model.frame(formula = formula, data = data)
  u <- model.matrix(attr(mf, "terms"), data = mf)
  if (!is.matrix(coef))
    coef <- matrix(coef, ncol = 1)

  if (nrow(coef) != ncol(u)){
    stop(paste("Length of coef must equal the number of independent variables",
               "(including the intercept)."))
  }

  XX <- matrix(0, nrow = n, ncol = n.obs) #'empty' X design matrix
  x <- matrix(0, nrow = n, ncol = 1)
  d <- 0
  for (t in 1:(n.burn + n.obs)){
    if (is.null(Dr)){
      e <- rnorm(n, sd = s)
    } else{
      e <- Dr %*% rnorm(n, sd = s)
    }
    if (t <= n.burn){
      d <- b * d + e
      x <- d
    } else{
      d <- b * d + e
      x <- t.scale[t - n.burn] * as.numeric(u %*% coef) + d
    }
    if (t > n.burn)
      XX[, t - n.burn] <- as.matrix(x)
  }
  XX <- XX - rowMeans(XX)
  return(XX)
}

## CLS.fit ----
#' @title Fit a CLS model to a matrix of time series
#' @description
#' CLS.fit fits a CLS model using lm() to a matrix of time series located in
#' each row
#'
#' @param X
#' @param t.scale
#'
#' @return
#' @export
#'
#' @examples
CLS.fit <- function(X, t.scale){

  n <- dim(X)[1] #rows (observations) in X

  # CLS for entire map
  d <- data.frame(site = 1:n)
  for (i in 1:dim(X)[1]){
    x <- X[i, ]

    d$mean[i] <- mean(x)

    z.CLS <- lm(x[2:length(x)] ~ x[1:(length(x) - 1)] + t.scale[2:length(x)])
    d$c[i] <- summary(z.CLS)$coef[3, 1] #estimate of time effect
    d$t[i] <- summary(z.CLS)$coef[3, 3] #t-value
    d$p[i] <- summary(z.CLS)$coef[3, 4] #p-value
    d$b[i] <- summary(z.CLS)$coef[2, 1] #estimate of x_t-1 effect
    d$MSE[i] <- summary(z.CLS)$sigma^2
  }
  return(d)
}

## LS.fit ----
#' @title Fit a LS model to a matrix of time series
#' @description
#' LS.fit fits a LS model using lm() to a matrix of time series located in
#' each row
#'
#' @param X
#' @param t.scale
#'
#' @return
#' @export
#'
#' @examples
LS.fit <- function(X, t.scale){

  n <- dim(X)[1]

  # LS for entire map
  d <- data.frame(site = 1:n)
  for (i in 1:dim(X)[1]){
    x <- X[i, ]

    d$mean[i] <- mean(x)

    z.LS <- lm(x ~ t.scale)
    d$c.LS[i] <- summary(z.LS)$coef[2, 1] #effect of time
    d$p.LS[i] <- summary(z.LS)$coef[2, 4] #p value
  }
  return(d)
}

## taper.spherical ----
#' @name taper.spherical
#' @rdname taper.spherical
#' @title Create tapered distance matrices
#' @description
#' taper.spherical and taper.spherical.dif are different formulations for
#' tapered distance matrices in which cells more distant than beta are set to
#' zero.
#'
#' @param d
#' @param beta
#' @param b
#' @param core
#' @return
NULL
#' @rdname taper.spherical
#' @export
#'
#' @examples
#' taper.spherical(d, beta)
taper.spherical <- function(d, beta){
  x <- d
  x[d > beta] <- 0
  x[d <= beta] <- ((1 - d[d <= beta]/beta)^2) * (1 + d[d <= beta]/(2 * beta))
  x
}
#' @rdname taper.spherical
#' @export
#'
#' @examples
#' taper.spherical.dif(d, cor, b)
taper.spherical.dif <- function(d, cor, b){
  beta <- exp(-b)
  x <- d
  x[d > beta] <- 0
  x[d <= beta] <- ((1 - d[d <= beta]/beta)^2) * (1 + d[d <= beta]/(2 * beta))
  cor - x
}

## spatialcor.fit ----
#' @name spatialcor.fit
#' @rdname spatialcor.fit
#' @title Fit a distance matrix to a set of time series using CLS
#'
#' @description
#' spatialcor.fit fits an exponential or taper-spherical distance matrix to a
#' spatial set of time series by performing CLS, calculating the correlation
#' matrix of residuals, and then fitting the correlations using nls.
#'
#' spatialcor.fit.data fits an exponential or taper-spherical distance matrix
#' to a spatial set of time series by performing CLS, calculating the
#' correlation matrix of residuals, and then fitting the correlations using nls.
#'
#' NOTE: In contrast to spatialcor.fit, spatialcor.fit.data returns the spatial
#' correlation scaled to the input distance matrix in km (from distm)
#'
#' @param X
#' @param t.scale
#' @param Dist
#' @param r.start
#' @param fit.n.sample
#' @param FUN
#' @param plot.fig
#' @param col.plot
#'
# #' @param X
# #' @param t.scale
#' @param data
# #' @param r.start
# #' @param fit.n.sample
# #' @param FUN
# #' @param plot.fig
# #' @param col.plot
#'
#' @return
NULL
#' @rdname spatialcor.fit
#' @export
#'
#' @examples
spatialcor.fit <- function(X, t.scale, Dist, r.start = 0.1, fit.n.sample,
                           FUN = "exponential", plot.fig = F, col.plot = NULL){

  n <- nrow(X)

  # subsample for r.fit
  fit.pick <- sample.int(n = n, size = fit.n.sample)

  resid <- matrix(0, nrow = fit.n.sample, ncol = n.obs - 1)
  for (i in 1:fit.n.sample){
    x <- X[fit.pick[i], ]

    z.CLS <- lm(x[2:length(x)] ~ x[1:(length(x) - 1)] + t.scale[2:length(x)])
    resid[i, ] <- z.CLS$resid
  }
  cor.resid <- cor(t(resid))
  dist <- Dist[fit.pick, fit.pick]/max(Dist)

  # colors for plotting
  if (is.null(col.plot)){
    col.plot <- "black"
  } else{
    col.plot <- col.plot[fit.pick]
  }

  cor.resid[lower.tri(cor.resid)] <- NA
  dist[lower.tri(dist)] <- NA

  v.cor.resid <- matrix(cor.resid, ncol = 1)
  v.dist <- matrix(dist, ncol = 1)
  v.cor.resid <- v.cor.resid[!is.na(v.cor.resid)]
  v.dist <- v.dist[!is.na(v.dist)]

  w <- as.data.frame(cbind(v.dist, v.cor.resid))
  names(w) <- c("dist", "cor")

  if (FUN == "exponential"){
    fit <- nls(cor ~ exp(-dist/r), data = w, start = list(r = r.start),
               nls.control(maxiter = 500))
    spatialcor <- coef(fit) * max(Dist)
  }
  if (FUN == "taper-spherical"){
    fit <- nls(~taper.spherical.dif(d = dist, cor = cor, b = b), data = w,
               start = list(b = 0.5), nls.control(maxiter = 500))
    spatialcor <- exp(-coef(fit)) * max(Dist)
  }
  if (plot.fig){
    plot(dist * max(Dist), cor.resid, pch = 20, cex = 0.5, col = col.plot)
    x.dist <- (1:fit.n.sample)/fit.n.sample * max(Dist)
    if (FUN == "exponential")
      lines(x.dist, exp(-x.dist/spatialcor), col = "red", lty = 2)
    if (FUN == "taper-spherical")
      lines(x.dist, taper.spherical(d = x.dist, beta = spatialcor), col = "red",
            lty = 2)
  }
  return(list(spatialcor = spatialcor, spatialcor.sigma = summary(fit)$sigma))
}
#' @rdname spatialcor.fit
#' @export
#'
#' @examples
spatialcor.fit.data <- function(X, t.scale, data, r.start = 0.1, fit.n.sample,
                                FUN = "exponential", plot.fig = F,
                                col.plot = NULL){

  n <- nrow(X)

  # subsample for r.fit
  fit.pick <- sample.int(n = n, size = fit.n.sample)

  resid <- matrix(0, nrow = fit.n.sample, ncol = n.obs - 1)
  for (i in 1:fit.n.sample){
    x <- X[fit.pick[i], ]

    z.CLS <- lm(x[2:length(x)] ~ x[1:(length(x) - 1)] + t.scale[2:length(x)])
    resid[i, ] <- z.CLS$resid
  }
  cor.resid <- cor(t(resid))

  # create distance matrix in kilometers
  location <- data[fit.pick,c('lng','lat')]
  Dist <- geosphere::distm(location, fun=distGeo)/1000
  dist <- Dist/max(Dist)

  # colors for plotting
  if (is.null(col.plot)){
    col.plot <- "black"
  } else{
    col.plot <- col.plot[fit.pick]
  }

  cor.resid[lower.tri(cor.resid)] <- NA
  dist[lower.tri(dist)] <- NA

  v.cor.resid <- matrix(cor.resid, ncol = 1)
  v.dist <- matrix(dist, ncol = 1)
  v.cor.resid <- v.cor.resid[!is.na(v.cor.resid)]
  v.dist <- v.dist[!is.na(v.dist)]

  w <- as.data.frame(cbind(v.dist, v.cor.resid))
  names(w) <- c("dist", "cor")

  if (FUN == "exponential"){
    fit <- nls(cor ~ exp(-dist/r), data = w, start = list(r = r.start))
    spatialcor <- coef(fit) * max(Dist)
  }
  if (FUN == "taper-spherical"){
    fit <- nls(~taper.spherical.dif(d = dist, cor = cor, b = b), data = w,
               start = list(b = 0.5))
    spatialcor <- exp(-coef(fit)) * max(Dist)
  }
  if (plot.fig){
    plot(dist * max(Dist), cor.resid, pch = 20, cex = 0.5, col = col.plot)
    x.dist <- (1:fit.n.sample)/fit.n.sample * max(Dist)
    if (FUN == "exponential")
      lines(x.dist, exp(-x.dist/spatialcor), col = "red", lty = 2)
    if (FUN == "taper-spherical")
      lines(x.dist, taper.spherical(d = x.dist, beta = spatialcor),
            col = "red", lty = 2)
  }
  return(list(spatialcor = spatialcor, spatialcor.sigma = summary(fit)$sigma))
}

## V.fit ----
#' @title Construct covariance matrix for remoteSTAR structures
#' @description
#' V.fit constructs the covariance matrix from a distance matrix for exponential
#' and taper-spherical structures
#'
#' @param Dist
#' @param spatialcor
#' @param FUN
#'
#' @return
#' @export
#'
#' @examples
V.fit <- function(Dist, spatialcor, FUN = "exponential"){

  if (FUN == "exponential")
    return(exp(-Dist/spatialcor))

  if (FUN == "taper-spherical")
    return(taper.spherical(Dist, spatialcor))

}

## GLS.fit ----
#' @title Fit a GLS to remoteSTAR data
#' @description
#' GLS.fit fits a GLS to the data given a specified V or invcholV
#'
#' @param formula
#' @param formula0
#' @param data
#' @param V
#' @param invcholV
#'
#' @return
#' @export
#'
#' @examples
GLS.fit <- function(formula, formula0 = NULL, data, V = NULL, invcholV = NULL){

  mf <- model.frame(formula = formula, data = data)
  x <- model.matrix(attr(mf, "terms"), data = mf)
  y <- model.response(mf)
  n <- length(y)

  if (is.null(invcholV)){
    if (is.null(V)){
      invcholV <- diag(n)
    } else{
      invcholV <- t(backsolve(chol(V), diag(n)))
    }
  }

  xx <- invcholV %*% x
  yy <- invcholV %*% y

  coef <- as.numeric(solve(crossprod(xx), crossprod(xx,yy)))
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

  if (is.null(formula0)){
    x0 <- matrix(1, nrow = n, ncol = 1)
    xx0 <- invcholV %*% x0
    coef0 <- solve(crossprod(xx0), crossprod(xx0, yy))
    SSE0 <- as.numeric(crossprod(yy - xx0 %*% coef0))
    df0 <- ncol(coef0)
  } else{
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

  if (ncol(xx) > 1){
    FF <- (n - ncol(xx))/(ncol(xx) - ncol(xx0)) * (SSE0 - SSE)/SSE
    df1.F <- ncol(xx) - ncol(xx0)
    df2.F <- n - ncol(xx)
    p.F <- pf(FF, df1 = df1.F, df2 = df2.F, lower.tail = F)
    df.F <- c(df1.F, df2.F)
  } else{
    FF <- (n - 1) * (SSE0 - SSE)/SSE
    df1.F <- 1
    df2.F <- n - 1
    p.F <- pf(FF, df1 = df1.F, df2 = df2.F, lower.tail = F)
    df.F <- c(df1.F, df2.F)
  }

  return(list(coef = coef, se = se, t = t, df.t = df.t, p.t = p.t, F = FF,
              df1.F = df1.F, df2.F = df2.F, p.F = p.F, logLik = logLik,
              logLik0 = logLik0, MSE = MSE, MSE0 = MSE0, MSR = MSR, SSE = SSE,
              SSE0 = SSE0, SSR = SSE0 - SSE, coef0 = coef0, se0 = se0,
              varX = varX, varcov = varcov, varcov0 = varcov0,
              invcholV = invcholV, xx=xx, xx0=xx0, yy=yy))
}

## nugget.fit ----
#' @name nugget.fit
#' @rdname nugget.fit
#' @title fit a nugget to the covariance matrix V
#' @description
#' nugget.fit fits a nugget to the covariance matrix V. Note that this involves
#' refitting the data using GLS.fit. nugget.fit uses the function nuget.fit.funct.
#' @param nugget
#' @param formula
#' @param data
#' @param V
#' @param verbose
#'
# #' @param formula
# #' @param data
# #' @param V
#' @param nugget.tol
#' @param interval
# #' @param verbose
#'
#' @return
NULL
#' @rdname nugget.fit
#' @export
#'
#' @examples
nugget.fit.funct <- function(nugget, formula, data, V, verbose = FALSE){
  n <- ncol(V)
  invcholV <- t(backsolve(chol((1 - nugget) * V + nugget * diag(n)), diag(n)))
  z <- GLS.fit(formula, data = data, invcholV = invcholV)
  if(verbose == TRUE) show(c(z$logLik, nugget))
  return(z$logLik)
}
#' @rdname nugget.fit
#' @export
#'
#' @examples
nugget.fit <- function(formula, data, V, nugget.tol = 0.00001,
                       interval = c(0, 1),
                       verbose = FALSE){
  opt.nugget <- optimize(nugget.fit.funct, formula, data = data, V = V,
                         interval = interval, maximum = T, tol = nugget.tol,
                         verbose = verbose)
  # check at the zero boundary
  if(opt.nugget$maximum < nugget.tol){
    nugget0.fit <- nugget.fit.funct(0, formula, data, V)
    if(nugget0.fit > opt.nugget$objective) opt.nugget$maximum <- 0
  }
  return(opt.nugget$maximum)
}

## GLS.partition.data ----
#' @title Fit a GLS to specified partitions
#' @description
#' GLS.partition.data fits a GLS to a specified number of random partitions of
#' the data after fitting the nugget. It calls GLS.fit().
#' @param formula
#' @param formula0
#' @param data
#' @param r
#' @param est.nugget
#' @param npart
#' @param partition
#' @param fixed.nugget
#' @param nugget.tol
#' @param min.num.rSS
#'
#' @return
#' @export
#'
#' @examples
GLS.partition.data <- function(formula, formula0 = NULL, data, r,
                               est.nugget = T, npart = 10, partition = NULL,
                               fixed.nugget = NULL, nugget.tol = 0.00001,
                               min.num.rSS = 12){

  max.offdiag.matrices <- ceiling((1 + (1 + 8*min.num.rSS)^.5)/2)

  n <- nrow(data)
  if (!is.null(partition)){
    npart <- nrow(partition)
    nn <- n - (n%%npart)
    n.p <- nn/npart
    pick <- partition
  } else{
    nn <- n - (n%%npart)
    n.p <- nn/npart
    pick <- matrix(sample(n)[1:nn], nrow = npart)
  }

  mf <- model.frame(formula = formula, data = data)
  df2 <- n.p - (ncol(model.matrix(attr(mf, "terms"), data = mf)) - 1)
  mf0 <- model.frame(formula = formula0, data = data)
  df0 <- n.p - (ncol(model.matrix(attr(mf0, "terms"), data = mf0)) - 1)
  df1 <- df0 - df2

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
  interval <- c(0,1)
  for (i in 1:npart){

    data.part <- data[pick[i,],]

    # create distance matrix in kilometers
    location <- data.part[,c('lng','lat')]
    Dist.part <- geosphere::distm(location, fun=distGeo)/1000

    Vp <- V.fit(Dist.part, spatialcor = r, FUN = "exponential")
    if (is.null(fixed.nugget) & est.nugget){
      nugget <- nugget.fit(formula, data.part, Vp, interval = interval)
      interval <- c(0, max(1000*nugget.tol, min(100*nugget,1)))
      Vp <- (1 - nugget) * Vp + nugget * diag(n.p)
    } else{
      if (is.null(fixed.nugget)){
        nugget <- 0
        Vp <- Vp
      } else{
        nugget <- fixed.nugget[i]
        Vp <- (1 - nugget) * Vp + nugget * diag(n.p)
      }
    }
    z.part <- GLS.fit(formula, formula0, data = data.part, V = Vp)

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
    if(i <= max.offdiag.matrices){
      invcholV.part[[i]] <- z.part$invcholV
      xx.part[[i]] <- z.part$xx
      xx0.part[[i]] <- z.part$xx0
    }
  }

  rSSE.part <- matrix(NA, nrow=npart, ncol=npart)
  rSSR.part <- matrix(NA, nrow=npart, ncol=npart)
  for (i in 1:min(max.offdiag.matrices-1,(npart-1))){
    for (j in (i+1):min(max.offdiag.matrices,npart)){
      data.part <- data[c(pick[i,], pick[j,]),]

      # create distance matrix in kilometers
      location <- data.part[,c('lng','lat')]
      Dist.part <- geosphere::distm(location, fun=distGeo)/1000
      Vpick <- V.fit(Dist.part, spatialcor = r, FUN = "exponential")
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
  }

  coef <- rowMeans(coef.part)
  coef0 <- rowMeans(coef0.part)
  rSSR <- mean(rSSR.part, na.rm=T)
  rSSE <- mean(rSSE.part, na.rm=T)

  Fmean <- mean(F.part)

  return(list(coef = coef, Fmean = Fmean, df1 = df1, df2 = df2,
              SSR.part = SSR.part, SSE.part = SSE.part, SSE0.part = SSE0.part,
              logLik.part = logLik.part, logLik0.part = logLik0.part,
              nugget = mean(nugget.part), nugget.part = nugget.part,
              F.part = F.part, p.F.part = p.F.part, coef.part=coef.part,
              se.part=se.part, coef0.part=coef0.part, se0.part=se0.part,
              rSSR = rSSR, rSSE = rSSE, rSSR.part = rSSR.part,
              rSSE.part = rSSE.part, npart = npart, partition = pick))
}

## correlated.F.bootstrap ----
#' @title bootstrap multiple-comparison corrected p-values
#' @description
#' correlated.F.bootstrap gives the multiple-comparisons-corrected p-values
#' and also performs a bootstrap to get the distribution of a correlated F test.
#' It might be possible to work this out algebraically.
#'
#' @param Fmean.obs
#' @param rSSR
#' @param rSSE
#' @param df1
#' @param df2
#' @param npart
#' @param nboot
#'
#' @return
#' @export
#'
#' @examples
correlated.F.bootstrap <- function(Fmean.obs, rSSR, rSSE, df1, df2, npart,
                                   nboot = 2000){
  part <- rep(1:npart, each=df1)

  rZ <- rSSR^.5/df1
  v.MSR <- diag(df1) - rZ
  v.MSR <- kronecker(diag(npart),v.MSR) + rZ
  D.MSR <- t(chol(v.MSR, pivot=T))
  if(attr(D.MSR, "rank") < npart){
    rank.MSR <- attr(D.MSR, "rank")
    v.MSR <- diag(df1) - .99/df1
    v.MSR <- kronecker(diag(npart),v.MSR) + rZ
    D.MSR <- t(chol(v.MSR, pivot=T))
  }else{
    rank.MSR <- NA
  }

  v.MSE <- (1-rSSE) * diag(npart) + rSSE
  D.MSE <- t(chol(v.MSE))

  count <- 0
  for(boot in 1:nboot){
    Z1 <- D.MSR %*% rnorm(npart*df1)
    MSR.boot <- aggregate(Z1^2, by=list(part), FUN=sum)[,2]/df1
    MSE.boot <- 1 + D.MSE %*% rnorm(npart, mean=0, sd=(2*df2)^.5/df2)
    if(Fmean.obs < mean(MSR.boot/MSE.boot)) count <- count + 1
  }
  return(list(pvalue = count/nboot, nboot = nboot, rank.MSR = rank.MSR))
}

## GLS.partition.pvalue ----
#' @title Wrapper for correlated.F.bootstrap()
#' @description
#' GLS.partition.pvalue is a wrapper for correlated.F.bootstrap()
#'
#' @param z
#' @param nboot
#'
#' @return
#' @export
#'
#' @examples
GLS.partition.pvalue <- function(z, nboot = 2000){
  if(is.finite(z$rSSR)){
    p.Fmean <- correlated.F.bootstrap(Fmean.obs = z$Fmean, rSSR = z$rSSR,
                                      rSSE = z$rSSE, df1 = z$df1, df2 = z$df2,
                                      npart = z$npart, nboot = nboot)
  }else{
    p.Fmean <- list(NA,NA,NA)
  }

  p.Fhochberg <- min(p.adjust(z$p.F.part, "hochberg"))
  p.Fhommel <- min(p.adjust(z$p.F.part, "hommel"))
  p.Ffdr <- min(p.adjust(z$p.F.part, "fdr"))

  return(list(p.Fmean = p.Fmean, p.Fhochberg = p.Fhochberg,
              p.Fhommel = p.Fhommel, p.Ffdr = p.Ffdr))
}
