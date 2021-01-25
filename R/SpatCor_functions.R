# taper_sphere() ----
## performs a taper-spherical transofrmation of d
## if a correlation 'cor' is given, the transformation is subtracted from cor

#' Taper-spherical transformation of distance matrix
#'
#' @description
#'
#' Note: this documentation is incomplete - arguments need better documenting
#'
#' @param d a numeric distance matrix or vector
#' @param beta beta parameter
#' @param cor optional correlation parameter to be used to calculate Delta d
#'
#' @return taper-spherical transformation of d
#' @export
#'
#' @examples #TBA
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

#' Scale a distance matrix by its maximum value
#'
#' @details a distance matrix is fit using `geosphere::distm()` and then divided
#' by itKs maximum distance.
#'
#' @param location n x 2 numeric matrix with latitude and longitude
#' (respectively) for sites (rows)
#'
#' @param scl scalar by which to divide the initial distance matrix. This
#' determines the units of the `max.dist` attribute. `scl = 1000` (default)
#' corresponds to distances in km.
#'
#' @return a scaled distance matrix with a `max.dist` attribute containing
#' the maximum distance by which the initial values were divided.
#'
#' @export
#'
#' @examples #TBA
scale_dist <- function(location, scl = 1000){
  D <- geosphere::distm(location)/scl
  m = max(D)
  d = D/m
  attr(d, "max.dist") <- m
  return(d)
}

# fit_spatialcor() ----

#' Title
#'
#' @param X n x p numeric matrix (usually of of remote sensing observations)
#'  taken from n sites (rows) and p time points (columns).
#' @param t numeric vector of length p containing the values for time. Recommended:
#' `scale(1:ncol(X))`.
#' @param r.start numeric starting point for r parameter which is mathematically optimized
#' @param a.start numeric starting point for a parameter
#' @param fit.n integer signifying the size of the random subset of X rows to be
#' used in the estimation.
#' @param fun character string indicating which function to use. Currently the available options are
#' "exponential", "exponential-power", and "taper-spherical".
#'
#' Note: Currently, only "exponential-power" is working properly. I'm not sure why. - morrowcj
#'
#' @param dist n x n numeric distance matrix
#' @param location n x 2 numeric matrix containing latitude and longitude respectively
#' @param scale.dist logical: should the distance matrix be scaled?
#' @param dist.scl if `scale.dist == TRUE`, `scl` see `?scale.dist()`.
#' @param U CURRENTLY UNUSED ARGUMENT
#' @param covars CURRENTLY UNUSED ARGUMENT
#' @param plot.fig logical: should a figure be plotted displaying results?
#'
#' Note: this functionality should be replaced by defining a `plot.method()`
#' for an appropriate S3 class.
#'
#' @param cols.plot optional vector of colors to use in the plot (size n)
#' @param ... additional arguments passed to `nls()`
#'
#' @return a list containing
#' mod: the nls() fit object
#' spatialcor: the maximum likelihood spatial coefficient parameter(s)
#' logLik: the log-likelihood of the fit
#' @export
#'
#' @examples #TBA
fit_spatialcor <- function(X, t, r.start = 0.1, a.start = 1,
                           fit.n = 1000, method = "exp",
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
  if (method == "exp" || method == "exponential") {
    fit <- nls(cor ~ exp(-dist/r), data = w, start = list(r = r.start))
    spatialcor <- coef(fit) * max.d
    if (plot.fig){y <- exp(-x.dist/spatialcor)}
  }
  if (method ==  "exp-pwr" || method == "exponential-power") {
    fit <- nls(cor ~ exp(-(dist/r)^a), data = w,
               start = list(r = r.start, a = a.start), nls.control(maxiter = 500))
    spatialcor <- coef(fit)
    spatialcor[1] <- spatialcor[1] * max.d
    if (plot.fig){y <- exp(-(x.dist/spatialcor[1])^spatialcor[2])}
  }
  if (method == "sphr" || method ==  "taper" || method == "taper-spherical") {
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
