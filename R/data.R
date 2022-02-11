# #' NDVI Remote sensing data for North America
# #'
# #' NDVI data for sites over time
# #'
# #' @format data frame with 348137 rows corresponding to sites and 37 columns:
# #' \describe{
# #'   \item{row}{grid row of the pixel}
# #'   \item{col}{grid column of the pixel}
# #'   \item{lng}{longitude of the pixel}
# #'   \item{lat}{latitude of the pixel}
# #'   \item{land}{dominant land class of the pixel}
# #'   \item{ndvi<t>}{ndvi value of the pixel during the year <t>}
# #' }
#
# "ndvi"

#' NDVI Remote sensing data for Alaska
#'
#' subset of the \code{ndvi} dataset containing only Alaska
#'
#' @format data frame with 31486 rows corresponding to sites and 38 columns:
#' \describe{
#'   \item{lng}{longitude of the pixel}
#'   \item{lat}{latitude of the pixel}
#'   \item{AR_coef}{pre-calculated AR REML coefficient}
#'   \item{CLS_coef}{pre-calculated CLS coefficient}
#'   \item{land}{dominant land class of the pixel}
#'   \item{land}{logical: is this land class rare?}
#'   \item{ndvi<t>}{ndvi value of the pixel during the year <t>}
#' }

"ndvi_AK"

#' NDVI Remote sensing data for Alaska
#'
#' subset of the \code{ndvi_AK} dataset containing 3000 random sites
#'
#' @format data frame with 3000 rows corresponding to sites and 38 columns:
#'
#' \describe{
#'   \item{lng}{longitude of the pixel}
#'   \item{lat}{latitude of the pixel}
#'   \item{AR_coef}{pre-calculated AR REML coefficient}
#'   \item{CLS_coef}{pre-calculated CLS coefficient}
#'   \item{land}{dominant land class of the pixel}
#'   \item{land}{logical: is this land class rare?}
#'   \item{ndvi<t>}{ndvi value of the pixel during the year <t>}
#' }

"ndvi_AK3000"

#' partitioned GLS results
#'
#' Example output from fitGLS_partition() fit to the \code{ndvi_AK} data set
#'
#' @format an S3 class "partGLS" object. See ?fitGLS_partition() for further
#' details

"partGLS_ndviAK"
