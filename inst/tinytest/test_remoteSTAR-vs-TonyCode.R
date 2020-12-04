# This file is meant to test that the functions in this package produce the
# same output as the original code.

## Load Tony's Code
source("tony-og-code.R")
load("data-for-tests.rda", verbose = TRUE)
library(microbenchmark)

# Test Data Setup ----
load("../../data/ndvi_AK3000.rda", verbose = TRUE)


N = 100

set.seed(916)

subsamp = sample.int(n = nrow(ndvi_AK3000), size = N)
test.df = as.data.frame(ndvi_AK3000[subsamp, ])
test.grid = test.df[, c("row", "col")]
test.coord = test.df[, c("lng", "lat")]
test.land = test.df[, c("land")]

test.X = as.matrix(test.df[, grep(pattern = "ndvi", names(test.df))])
test.n = nrow(test.X); test.p = ncol(test.X)
test.t = scale(seq_len(ncol(test.X)))
test.D = geosphere::distm(test.coord)/1000

# CLS functions ----
## Speed Test
CLS.bench = microbenchmark(old.CLS <-  CLS.fit(test.X, test.t),
                           new.CLS <- cls_star(test.X, test.t),
                           times = 10L)

## Compare relevant output
expect_equivalent(old.CLS[, c("site", "mean", "c", "t", "p", "b", "MSE")],
new.CLS[, c("site", "mean", "Est", "t", "p", "x_t0.EST", "MSE")])

## Assign y
test.y = with(new.CLS, Est/mean)

# spatial correlation functions ----

spatcor.bench = microbenchmark(old.r <- spatialcor.fit(test.X, test.t, test.D,
                                                      fit.n.sample = 100),
                               new.r <- fit_spatialcor(test.X, test.t, dist = D,
                                                      location = test.coord,
                                                      fit.n = 100),
                               times = 10L)

expect_equivalent(old.r$spatialcor, new.r$spatialcor)
expect_equivalent(old.r$spatialcor.sigma, sigma(new.r$mod))

# V fitting functions ----
# V.bench = microbenchmark(old.V = V.fit(test.D, spatialcor = r))

