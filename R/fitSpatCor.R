# taper_sphere() ----
## performs a taper-spherical transofrmation of d
## if a correlation 'cor' is given, the transformation is subtracted from cor

#' Title
#'
#' @param d
#' @param beta
#' @param cor
#'
#' @return
#' @export
#'
#' @examples
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

#' Title
#'
#' @param location
#' @param scl
#'
#' @return
#' @export
#'
#' @examples
scale_dist <- function(location, scl = 1000){
  D <- geosphere::distm(location)/scl
  m = max(D)
  d = D/m
  attr(d, "max.dist") <- m
  return(d)
}

# fit_spatialcor() ----
# X = Xmat
# location = dat[, c("lng", "lat")]
# t = t.scale
# data = data
# r.start = .1
# fit.n = 100
# dist.scl = 1000
# covars = NULL
# fun = "exp"
# Dist = geosphere::distm(location)
# scale.dist = TRUE


#' Title
#'
#' @param X
#' @param t
#' @param r.start
#' @param a.start
#' @param fit.n
#' @param fun
#' @param dist
#' @param location
#' @param scale.dist
#' @param dist.scl
#' @param U
#' @param covars
#' @param plot.fig
#' @param cols.plot
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
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
    dist <- scale_dist(location[sub.inx, ], scl = dist.scl)
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
