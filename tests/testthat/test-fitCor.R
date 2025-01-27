# simulated dummy data
time.points = 30 # time series length
map.width = 8 # square map width
coords = expand.grid(x = 1:map.width, y = 1:map.width) # coordinate matrix

# spatiotemporal variables with autocorrelation (dput output)
X <- readRDS("fitCor_data_X.rds")
Z <- readRDS("fitCor_data_Z.rds")

# AR map to use
AR.map = fitAR_map(X, coords, formula = y ~ Z, X.list = list(Z = Z), resids.only = FALSE)

# using pre-defined covariance function
## exponential covariance
exp_fit <- fitCor(AR.map$residuals, coords, covar_FUN = "covar_exp", start = list(range = .1))

## exponential-power covariance
expow_fit <- fitCor(AR.map$residuals, coords, covar_FUN = "covar_exppow", start = list(range = .1, shape = .2))

# user-specified covariance function
custom_fit <- fitCor(AR.map$residuals, coords, covar_FUN = function(d, r){d^r}, start = list(r = .1))

# specify which pixels to use, for reproducibility
pixel_fit <- fitCor(AR.map$residuals, coords, index = 1:64) #all

test_that("fitCor spcor is right", {
  expect_equal(exp_fit$spcor, c(range = 1.1394), tolerance = 0.001)
  expect_equal(expow_fit$spcor, c(range = 383.628, shape = 0.1315), tolerance = 0.001)
  expect_equal(custom_fit$spcor, c(r = 0.3481), tolerance = 0.001)
  expect_equal(pixel_fit$spcor, c(r = 1.1394), tolerance = 0.001)
})
