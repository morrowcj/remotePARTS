% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/optimize_GLS.R
\name{optim_GLS_func}
\alias{optim_GLS_func}
\alias{optimize_GLS}
\title{function to optimize GLS parameters r, a, and nugget}
\usage{
optim_GLS_func(
  par,
  y,
  modmat,
  D,
  verbose = FALSE,
  V.meth = "exponential-power"
)

optimize_GLS(
  formula,
  D,
  V.meth = "exponential-power",
  nugget = NA,
  spcor = NA,
  verbose = FALSE,
  contrasts = NULL,
  data,
  ...
)
}
\arguments{
\item{par}{parameters to optimize}

\item{y}{response vector}

\item{modmat}{predictor model matrix}

\item{D}{Distance matrix}

\item{verbose}{should the optimizer steps be printed to the console?}

\item{V.meth}{method passed to \code{fitV()} default: "exponential-power"}

\item{formula}{model formula}

\item{nugget}{NA: find maximum liklihood nugget}

\item{spcor}{NA: find maximum liklihood spatial correlation}

\item{contrasts}{possible contrasts object}

\item{data}{optional data to search for objects in formula}

\item{...}{additional arguments passed to \code{fitGLS2()}}
}
\value{
a remoteGLS object
}
\description{
function to optimize GLS parameters r, a, and nugget

Fit a GLS by estimating r, a, and nugget
}
\examples{
set.seed(916)

## load Alaska 3000 data
data("ndvi_AK3000")

## take a random subset of 100 pixels (to make example fast)
subsamp = sample.int(n = nrow(ndvi_AK3000), size = 100)

## subset the data: we now have 100 pixels, latitude, longitude, and land class
df = ndvi_AK3000[subsamp, c("lng", "lat", "land")] # subset the data

## simulate a response variable kappa
df$kappa = .5*df$lat + .2*df$lng + rnorm(100) #simulate response variable

## calculate distance matrix
D = geosphere::distm(df[, c("lng", "lat")])/1000 #distance in km

## fit the GLS, including ML spatial correlation and nugget
optimize_GLS(kappa ~ 0 + land, data = df, D = D, verbose = FALSE)
}
\seealso{
[fitGLS2()]
}
