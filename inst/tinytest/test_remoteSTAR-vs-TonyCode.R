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
test.modmat = model.matrix(~ 0 + test.land)
test.n = nrow(test.X); test.p = ncol(test.X)
test.t = scale(seq_len(ncol(test.X)))
test.D = geosphere::distm(test.coord)/1000 # long time
test.form = "test.y ~ 0 + test.land"

# CLS functions ----
## benchmark
CLS.bench = microbenchmark(
  old.CLS = {old.CLS <-  CLS.fit(test.X, test.t)},
  new.CLS = {new.CLS <- cls_star(test.X, test.t)},
  times = 10L
)
## compare output
expect_equivalent(old.CLS[, c("site", "mean", "c", "t", "p", "b", "MSE")],
new.CLS[, c("site", "mean", "Est", "t", "p", "x_t0.EST", "MSE")])
## assign y
test.y = with(new.CLS, Est/mean)

# spatial correlation functions ----
## benchmark
spatcor.bench = microbenchmark(
  old.r = {old.r <- spatialcor.fit(test.X, test.t, test.D, fit.n.sample = N)},
  new.r = {new.r <- fit_spatialcor(test.X, test.t, dist = D,
                                   location = test.coord, fit.n = N)},
  times = 10L
)
## compare output
expect_equivalent(old.r$spatialcor, new.r$spatialcor)
expect_equivalent(old.r$spatialcor.sigma, sigma(new.r$mod))
## assign r
r = new.r$spatialcor


# V fitting functions ----
## benchmark
V.bench = microbenchmark(old.V = {old.V <-  V.fit(test.D, r)},
                         new.V = {new.V <- fitV(test.D, r)},
                         times = 10L)
## compare output
expect_equivalent(old.V, new.V)
## assign V
test.V = new.V


# nugget fitting functions ----
## benchmark
nug.bench = microbenchmark::microbenchmark(times = 10L,
  old.nug = {old.nug <- nugget.fit(test.form, test.df, test.V)},
  # new.nug_C = {new.nug_C <- optimize_nugget(X = model.matrix(~ 0 + test.land),
  #                              V = test.V, y = test.y)},
  new.nug = {new.nug_R <- fitNugget(X = test.modmat,
                                    V = test.V, y = test.y)}
) # C++ version wont' run for large data... why?
## compare output
expect_equivalent(old.nug, new.nug_R)
# expect_equivalent(old.nug, new.nug_C)

# GLS functions ----
## benchmark
GLS.bench = microbenchmark(times = 10L,
                           old.GLS = {old.GLS <- GLS.fit(test.form, data = test.df,  V = test.V,
                                              invcholV = invert_chol(test.V, new.nug_R))},
                           new.GLS = {new.GLS <- fitGLS(X = test.modmat,
                                             y = test.y, V = test.V, nugget = new.nug_R,
                                             X0 = cbind(rep(1, nrow(test.X))), save_xx = TRUE,
                                             threads = 1)})
## compare output
expect_equivalent(old.GLS$coef, new.GLS$betahat)
expect_equivalent(old.GLS$t, new.GLS$tstat)
expect_equivalent(old.GLS$SSE, new.GLS$SSE)
expect_equivalent(old.GLS$SSE0, new.GLS$SSE0)
expect_equivalent(old.GLS$xx, new.GLS$xx)
expect_equivalent(old.GLS$xx0, new.GLS$xx0)
expect_equivalent(old.GLS$F, new.GLS$Fstat)

if(FALSE){
  library(ggplot2)
  autoplot(CLS.bench)
  autoplot(V.bench)
  autoplot(nug.bench) # new is faster
  autoplot(GLS.bench) # new is marginally slower (N = 100, trheads = 4)
}
