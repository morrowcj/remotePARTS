#' Fit varcov matrix from distance matrix
#' @rdname fitV
#'
#' @param Dist n x n numeric distance matrix
#' @param spatialcor spatial correlation parameter(s)
#' @param method function with which to transform the distance matrix into the
#' varcov matrix
#'
#' @return n x n variance-covariance matrix
#'
#' @details \code{Dist} is transformed via \code{method}
#'
#' @examples
#'
#' @export
fitV <- function(Dist, spatialcor, method = "exponential") {

  if (method == "exponential") {
    return(exp(-Dist/spatialcor))
  } else if (method == "exponential-power") {
    return(exp(-(Dist/spatialcor[1])^spatialcor[2]))
  } else if (method == "taper-spherical") {
    return(taper_sphere(Dist, spatialcor))
  } else {
    stop(paste0("method '", method, "' not recognized."))
  }
}



#' Fit varcov matrix from distance matrix
#' @rdname fitV
#'
#' @details \code{fitV.switch()} is not exported and is the same as
#' \code{fitV()} but uses \code{switch()} instead of \code{if else} statements.
#' \code{fitV.switch()} will likely be removed in future implementations.
fitV.switch <- function(Dist, spatialcor, method = "exponential"){
  switch(method,
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
