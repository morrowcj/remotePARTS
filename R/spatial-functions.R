
## distance in km ----
#' @title Calculate distance matrix in kilometers
#' @family distm
#'
#' @param coords coordinate matrix
#' @param coords2 optional coordinate matrix
#'
#' @return distance matrix
#'
#' @details this function is simply a wrapper for \code{geosphere::distm()}
#' @seealso \code{?geosphere::distm()}
#'
#' @export
distm_km <- function(coords, coords2 = NULL){
  if(is.null(coords2)){
    return(geosphere::distm(coords)/1000)
  } else {
    return(geosphere::distm(coords, coords2)/1000)
  }
}

## scaled distance ----
#' @title Calculate a unit distance matrix
#' @family distm
#'
#' @details a distance matrix is fit and then divided by its maximum distance.
#'
#' @param coords a coordinate matrix with 2 columns and rows corresponding to
#' each location.
#' @param coords2 an optional coordinate matrix
#'
#' @param distm_FUN function used to calculate the distance matrix. This function
#' dictates the units of "max.dist"
#'
#' @return a distance matrix, scaled between 0 and 1, with a `max.dist`
#' attribute containing the maximum distance by which the initial values were
#' scaled.
#'
#' @examples
#' map.width = 3 # square map width
#' coords = expand.grid(x = 1:map.width, y = 1:map.width) # coordinate matrix
#' distm_scaled(coords) # calculate distance matrix
#'
#' @export
distm_scaled <- function(coords, coords2 = NULL, distm_FUN = "distm_km"){
  dist.f = match.fun(distm_FUN)
  D <- dist.f(coords, coords2)
  m = max(D)
  d = D/m
  attr(d, "max.dist") <- m
  return(d)
}

## Tapered-spherical covariance ----
#' @title Tapered-spherical distance-based covariance function
#'
#' @rdname covar_functions
#'
#' @param d a numeric vector or matrix of distances
#' @param theta distance beyond which covariances are forced to 0.
#' @param cor optional correlation parameter. If included, the covariance is
#' subtracted from \code{cor}.
#'
#' @details \code{covar_taper} calculates covariance v as follows:
#'
#' if \code{d <= theta}, then \code{v = ((1 - (d/theta))^2) * (1 + (d/(2 * theta)))}
#'
#' if \code{d > theta}, then \code{v = 0}
#'
#' @return a tapered-spherical transformation of d is returned.
#'
#' @examples
#'
#' # simulate dummy data
#' map.width = 5 # square map width
#' coords = expand.grid(x = 1:map.width, y = 1:map.width) # coordinate matrix
#'
#' # calculate distance
#' D = geosphere::distm(coords) # distance matrix
#'
#' # visualize covariance matrix
#' image(covar_taper(D, theta = .5*max(D)))
#'
#' # plot tapered covariance function
#' curve(covar_taper(x, theta = .5), from = 0, to = 1);abline(v = 0.5, lty = 2, col = "grey80")
#'
#' @export
covar_taper <- function(d, theta, cor = NULL){
  # conditional theta
  if (!missing(cor) && !is.null(cor)){
    theta <- exp(-theta)
  }
  # taper d, given theta
  x <- ifelse(test = d > theta,
              yes = 0,
              no = ((1 - (d/theta))^2) * (1 + (d/(2 * theta)))
  )
  # conditional return
  if (!missing(cor) && !is.null(cor)){
    return(cor - x) # taper.spherical.diff
  } else {
    return(x) # taper.spherical
  }
}

## Exponential covariance ----
#' @title Exponential distance-based covariance function
#'
#' @rdname covar_functions
#'
#' @param d a numeric vector or matrix of distances
#' @param range range parameter
#'
#' @details \code{covar_exp} calculates covariance v as follows:
#'
#' \code{v = exp(-d/range)}
#'
#' @return
#'
#' @examples
#'
#' # visualize covariance matrix
#' image(covar_exp(D, range = .2*max(D)))
#'
#' # plot exponential function with different ranges
#' curve(covar_exp(x, range = .2), from = 0, to = 1)
#' curve(covar_exp(x, range = .1), from = 0, to = 1, col = "blue", add = TRUE)
#' legend("topright", legend = c("range = 0.2", "range = 0.1"), col = c("black", "blue"), lty = 1)
#'
#' @export
covar_exp <- function(d, range){
  exp(-d/range)
}

## Exponential-power covariance ----
#' @title Exponential-power distance-based covariance function
#'
#' @rdname covar_functions
#'
#' @param d a numeric vector or matrix of distances
#' @param range range parameter
#' @param shape shape parameter
#'
#' @details \code{covar_exppow} calculates covariance v as follows:
#'
#' \code{v = exp(-(d/range)^2)}
#'
#' Note that \code{covar_exppow(..., shape = 1)} is equivalent to
#' \code{covar_exp()} but is needed as a separate function for use with \code{fitCor}.
#'
#' @return
#'
#' @examples
#'
#' # visualize Exponential covariance matrix
#' image(covar_exppow(D, range = .2*max(D), shape = 1))
#'
#' # visualize Exponential-power covariance matrix
#' image(covar_exppow(D, range = .2*max(D), shape = .5))
#'
#' # plot exponential power function with different shapes
#' curve(covar_exppow(x, range = .2, shape = 1), from = 0, to = 1)
#' curve(covar_exppow(x, range = .2, shape = .5), from = 0, to = 1, col = "blue", add = TRUE)
#' legend("topright", legend = c("shape = 1.0", "shape = 0.5"), col = c("black", "blue"), lty = 1)
#'
#' @export
covar_exppow <- function(d, range, shape){
  exp(-(d/range)^shape)
}

