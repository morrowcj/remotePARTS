# taper_sphere() ----
## performs a taper-spherical transofrmation of d
## if a correlation 'cor' is given, the transformation is subtracted from cor
taper_sphere <- function(d, beta, cor = NULL){

  ## conditional beta
  if (!missing(cor) && !is.null(cor)){
    beta <- exp(-beta)
  }

  # taper d
  x <- ifelse(test = d > beta,
              yes = 0,
              no = ((1 - (d/beta))^2) * (1 + (d/(2 * beta)))
  )

  ## conditional return
  if (!missing(cor) && !is.null(cor)){
    return(cor - x) # taper.spherical.diff
  } else {
    return(x) # taper.spherical
  }
}

# scale_dist() ----
## calculates a scaled distance matrix
scale_dist <- function(X, location, scl = 1000){
  D <- geosphere::distm(location)/scl
  m = max(D)
  d = D/m
  attr(d, "max.dist") <- m
  return(d)
}

# fit_spatialcor() ----
X = Xmat
location = dat[, c("lng", "lat")]
t = t.scale
data = data
r.start = .1
fit.n = 100
dist.scl = 1000
covars = NULL
fun = "exp"
Dist = geosphere::distm(location)
scale.dist = TRUE


fit_spatialcor <- function(X, t, r.start = 0.1, a.start = 1,
                           fit.n = 100, fun = "exp",
                           dist, location,
                           scale.dist = TRUE, dist.scl = 1000, U = NULL,
                           covars = NULL, plot.fig = FALSE, cols.plot = NULL,
                           ... ## additional arguments to pass to nls()
                           ){

  sub.inx <- sample.int(nrow(X), fit.n) # index for subsetting the data
  X.sub <- X[sub.inx, ]
  # calculate the residuals
  resids <- sapply(seq_len(fit.n), function(i){
    x <- X.sub[i, ]
    lm(x[2:length(x)] ~ x[1:(length(x) - 1)] + t[2:length(x)])$resid
  })

  # resids <- matrix(0, ncol = fit.n, nrow = ncol(X.sub) - 1)
  # for (i in 1:fit.n) {
  #   x <- X.sub[i, ]
  #   resids[, i] <- lm(x[2:length(x)] ~ x[1:(length(x) - 1)] + t[2:length(x)])$resid
  #   # resids[, i] <- z.CLS$resid
  # }
  # cor.resids <- cor(resids)

  cor.resids <- cor(resids) # correlation of residuals

  # scale the distance matrix, if asked
  if (scale.dist) {
    dist <- scale_dist(X.sub, location[sub.inx, ], scl = dist.scl)
    max.d = attr(dist, "max.dist")
  } else {
    dist <- dist[sub.inx, sub.inx]
    max.d = max(dist)
    dist <- dist/max.d
  }

  x.dist <- (1:fit.n)/fit.n * max.d

  # convert matrices to vectors and combine to data frame
  v.cor.resids <- cor.resids[upper.tri(cor.resids, diag = TRUE)]
  v.dist <- dist[upper.tri(dist, diag = TRUE)]
  w = data.frame(dist = v.dist, cor = v.cor.resids)


  # perform nls regression to get the spatial correlation
  if (fun == "exp" || fun == "exponential") {
    fit <- nls(cor ~ exp(-dist/r), data = w, start = list(r = r.start))
    spatialcor <- coef(fit) * max.d
    if (plot.fig){y <- exp(-x.dist/spatialcor)}
  }
  if (fun ==  "exp-pwr" || fun == "exponential-power") {
    fit <- nls(cor ~ exp(-(dist/r)^a), data = w,
               start = list(r = r.start, a = a.start), nls.control(maxiter = 500))
    spatialcor <- coef(fit)
    spatialcor[1] <- spatialcor[1] * max.d
    if (plot.fig){y <- exp(-(x.dist/spatialcor[1])^spatialcor[2])}
  }
  if (fun == "sphr" || fun ==  "taper" || fun == "taper-spherical") {
    # fit <- nls(~taper.spherical.dif(d = dist, cor = cor, b = b), data = w,
    #            start = list(b = 0.5))
    fit <- nls(~taper_sphere(d = dist, beta = b, cor = cor), data = w,
               start = list(b = 0.5))
    spatialcor <- exp(-coef(fit)) * max.d
    if (plot.fig){y <- taper_sphere(d = x.dist, beta = spatialcor)}
  }


  if (plot.fig) {

    # colors for plotting
    if (is.null(cols.plot) | missing(cols.plot)) {
      cols.plot <- "black"
    } else {
      cols.plot <- cosl.plot[sub.inx]
    }

    plot(dist * max.d, cor.resids, pch = 20, cex = 0.5, col = cols.plot)
    lines(x.dist, y, col = "red", lty = 1, lwd = 2)
  }

  return(list(mod = fit, spatialcor = spatialcor, logLik = logLik(fit)))
}

## Original code ----

taper.spherical <- function(d, beta) {
  x <- d
  x[d > beta] <- 0
  x[d <= beta] <- ((1 - d[d <= beta]/beta)^2) * (1 + d[d <= beta]/(2 * beta))
  x
}

taper.spherical.dif <- function(d, cor, b) {
  beta <- exp(-b)
  x <- d
  x[d > beta] <- 0
  x[d <= beta] <- ((1 - d[d <= beta]/beta)^2) * (1 + d[d <= beta]/(2 * beta))
  cor - x
}
spatialcor.fit <- function(X, t.scale, Dist, r.start = 0.1, fit.n.sample,
                           FUN = "exponential", plot.fig = F,
                           col.plot = NULL) {

  n <- nrow(X)

  # subsample for r.fit
  fit.pick <- sample.int(n = n, size = fit.n.sample)

  resid <- matrix(0, nrow = fit.n.sample, ncol = n.obs - 1)
  for (i in 1:fit.n.sample) {
    x <- X[fit.pick[i], ]

    z.CLS <- lm(x[2:length(x)] ~ x[1:(length(x) - 1)] + t.scale[2:length(x)])
    resid[i, ] <- z.CLS$resid
  }
  cor.resid <- cor(t(resid))
  dist <- Dist[fit.pick, fit.pick]/max(Dist)

  # colors for plotting
  if (is.null(col.plot)) {
    col.plot <- "black"
  } else {
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

  if (FUN == "exponential") {
    fit <- nls(cor ~ exp(-dist/r), data = w, start = list(r = r.start), nls.control(maxiter = 500))
    spatialcor <- coef(fit) * max(Dist)
  }
  if (FUN == "exponential-power") {
    fit <- nls(cor ~ exp(-(dist/r)^a), data = w, start = list(r = r.start, a = 1), nls.control(maxiter = 500))
    spatialcor <- coef(fit) * max(Dist)
  }
  if (FUN == "taper-spherical") {
    fit <- nls(~taper.spherical.dif(d = dist, cor = cor, b = b), data = w, start = list(b = 0.5), nls.control(maxiter = 500))
    spatialcor <- exp(-coef(fit)) * max(Dist)
  }
  if (plot.fig) {
    plot(dist * max(Dist), cor.resid, pch = 20, cex = 0.5, col = col.plot)
    x.dist <- (1:fit.n.sample)/fit.n.sample * max(Dist)
    if (FUN == "exponential")
      lines(x.dist, exp(-x.dist/spatialcor), col = "red", lty = 1, lwd = 2)
    if (FUN == "exponential-power")
      lines(x.dist, exp(-(x.dist/spatialcor[1])^spatialcor[2]), col = "red", lty = 1, lwd = 2)
    if (FUN == "taper-spherical")
      lines(x.dist, taper.spherical(d = x.dist, beta = spatialcor), col = "red", lty = 1, lwd = 2)
  }
  return(list(spatialcor = spatialcor, spatialcor.sigma = summary(fit)$sigma))
}

spatialcor.fit.U <- function(X, U=NULL, Dist, r.start = 0.1, fit.n.sample,
                             FUN = "exponential", plot.fig = F,
                             col.plot = NULL) {

  n <- nrow(X)

  # subsample for r.fit
  fit.pick <- sample.int(n = n, size = fit.n.sample)

  resid <- matrix(0, nrow = fit.n.sample, ncol = n.obs - 1)
  for (i in 1:fit.n.sample) {
    x <- X[fit.pick[i], ]

    if(is.null(U)){
      z.CLS <- lm(x[-1] ~ x[-length(x)])
    }else{
      u <- U[i, ]
      z.CLS <- lm(x[-1] ~ x[-length(x)] + u[-length(x)])
    }
    resid[i, ] <- z.CLS$resid
  }
  cor.resid <- cor(t(resid))
  dist <- Dist[fit.pick, fit.pick]/max(Dist)

  # colors for plotting
  if (is.null(col.plot)) {
    col.plot <- "black"
  } else {
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

  if (FUN == "exponential") {
    fit <- nls(cor ~ exp(-dist/r), data = w, start = list(r = r.start), nls.control(maxiter = 500))
    spatialcor <- coef(fit) * max(Dist)
  }
  if (FUN == "taper-spherical") {
    fit <- nls(~taper.spherical.dif(d = dist, cor = cor, b = b), data = w, start = list(b = 0.5), nls.control(maxiter = 500))
    spatialcor <- exp(-coef(fit)) * max(Dist)
  }
  if (plot.fig) {
    plot(dist * max(Dist), cor.resid, pch = 20, cex = 0.5, col = col.plot)
    x.dist <- (1:fit.n.sample)/fit.n.sample * max(Dist)
    if (FUN == "exponential")
      lines(x.dist, exp(-x.dist/spatialcor), col = "red", lty = 2)
    if (FUN == "taper-spherical")
      lines(x.dist, taper.spherical(d = x.dist, beta = spatialcor), col = "red", lty = 2)
  }
  return(list(spatialcor = spatialcor, spatialcor.sigma = summary(fit)$sigma))
}


# NOTE: This returns the spatial correlation scaled to the input distance matrix in km (from distm)
spatialcor.fit.data <- function(X, t.scale, data, r.start = 0.1, a.start = 1,
                                fit.n.sample, FUN = "exponential", plot.fig = F,
                                col.plot = NULL) {

  n <- nrow(X)
  n.obs <- ncol(X)

  # subsample for r.fit
  fit.pick <- sample.int(n = n, size = fit.n.sample)

  resid <- matrix(0, nrow = fit.n.sample, ncol = n.obs - 1)
  for (i in 1:fit.n.sample) {
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
  if (is.null(col.plot)) {
    col.plot <- "black"
  } else {
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

  if (FUN == "exponential") {
    fit <- nls(cor ~ exp(-dist/r), data = w, start = list(r = r.start))
    spatialcor <- coef(fit) * max(Dist)
  }
  if (FUN == "exponential-power") {
    #fit <- nls(cor ~ exp(-(dist/r)^a), data = w, start = list(r = r.start, a = a.start), nls.control(maxiter = 500))
    #browser()
    # f <- function(x) {
    # r <- x[1]
    # a <- x[2]
    # return(sum(exp((-w$dist/r)^a) - w$cor)^2)
    # }
    # fit <- optim(f, par = c(r.start, a.start))

    fit <- nls(cor ~ exp(-(dist/r)^a), data = w, start = list(r = r.start, a = a.start), nls.control(maxiter = 500))
    spatialcor <- coef(fit)
    spatialcor[1] <- spatialcor[1] * max(Dist)
  }
  if (FUN == "taper-spherical") {
    fit <- nls(~taper.spherical.dif(d = dist, cor = cor, b = b), data = w, start = list(b = 0.5))
    spatialcor <- exp(-coef(fit)) * max(Dist)
  }
  if (plot.fig) {
    plot(dist * max(Dist), cor.resid, pch = 20, cex = 0.5, col = col.plot)
    x.dist <- (1:fit.n.sample)/fit.n.sample * max(Dist)
    if (FUN == "exponential")
      lines(x.dist, exp(-x.dist/spatialcor), col = "red", lty = 1, lwd = 2)
    if (FUN == "exponential-power")
      lines(x.dist, exp(-(x.dist/spatialcor[1])^spatialcor[2]), col = "red", lty = 1, lwd = 2)
    if (FUN == "taper-spherical")
      lines(x.dist, taper.spherical(d = x.dist, beta = spatialcor), col = "red", lty = 1, lwd = 2)
  }
  return(list(spatialcor = spatialcor, spatialcor.sigma = summary(fit)$sigma, logLik = logLik(fit)))
}

