# ---- fitAR ----
## Setup
# simulate dummy data
t = 1:30 # times series

# Random independent variable dput(round(Z, 2))
Z = c(-0.6, -0.59, -0.78, -1.9, 1.13, -0.44, 0.54, 0.01, -0.19, -0.26,
      1.94, -1.31, 0.66, -1.26, 0.11, 0.36, -0.13, -0.12, 0.58, 1.22,
      1.53, -1.77, 1.02, -0.44, 0.18, -1.17, -1.32, -1.83, -1.45, -1.28
)
x = .2*Z + (.05*t) # generate dependent effects
x[2:30] = x[2:30] + .2*x[1:29] # add autocorrelation

# fit the AR model, using Z as a covariate
AR = fitAR(x ~ Z)


# now using time as a covariate
AR_time <- fitAR(x ~ t)

# linear models for comparison
fmz <- lm(x ~ Z)
fmt <- lm(x ~ t)

# source variable from a dataframe
df = data.frame(y = x, t.scaled = t/30, Z = Z)
AR_both <- fitAR(y ~ t.scaled + Z, data = df)

## Tests

test_that("AR results are correct", {
  expect_equal(AR_both$coefficients, c(-.005, 1.75, 0.173), tolerance = .001, ignore_attr = TRUE)
})

test_that("AR output differs from lm output", {
  expect_false(any(fmz$coefficients == AR$coefficients))
  expect_false(any(fmt$coefficients == AR_time$coefficients))
})

test_that("AR methods work correctly", {
  expect_equal(AR$coefficients, coefficients(AR))
  expect_equal(AR$residuals, residuals(AR))
  expect_equal(AR$fitted.values, fitted(AR))
  expect_equal(AR$logLik, as.vector(logLik.remoteTS(AR)))
})

# ---- fitAR_map ----

# simulated dummy data
time.points = 9 # time series length
map.width = 5 # square map width
coords = expand.grid(x = 1:map.width, y = 1:map.width) # coordinate matrix

# spatiotemporal variables with autocorrelation (dput output)
X <- structure(c(0.16, 0.11, 0.13, 0.22, 0.32, 0.27, 0.3, 0.31, 0.29,
                 0.23, 0.26, 0.35, 0.26, 0.39, 0.39, 0.46, 0.41, 0.5, 0.46, 0.49,
                 0.46, 0.54, 0.56, 0.57, 0.65, 0.27, 0.18, 0.49, -0.1, 0.15, 0.72,
                 0.36, -0.08, 0.81, 0.44, 0.1, 0.14, 0.37, 0.37, -0.38, 0.32,
                 0.5, 0.34, 0.48, 0.4, 0.17, 0.28, 0.23, 0.35, 0.49, -0.43, 0.29,
                 0.67, 0.24, 0.75, 0.7, 0.05, 0.34, 0.56, 0.48, -0.1, 0.05, 0.68,
                 0.41, 0.2, 0.11, 0.57, 0.13, 0.32, -0.29, 0.18, -0.49, 0.79,
                 0.26, 0.67, 0.27, 0.53, 0.5, 0.6, 0.31, 0.68, 0.43, 0.37, 0.44,
                 0.67, -0.01, 0.48, -0.21, 0.21, 0.15, 0.38, 0.69, 0.34, 0.41,
                 0.19, 0.29, 0.34, 0.29, 0.5, 0.79, 0.2, -0.06, 0.22, 0.14, 0.36,
                 0.82, 0.84, 0.47, 0.28, 0.3, -0.04, 0.02, -0.59, 0.03, -0.23,
                 0.2, 1.29, -0.12, 0.92, 0.52, 0.56, 0.01, 0.16, 0.37, 0.47, 0.51,
                 0.24, -0.11, 0.38, 0.62, 1.05, 0.69, 0.79, 0.77, 0.56, -0.08,
                 0.35, -0.13, 0.98, 0.27, 0.19, 0.95, 0.06, 0.98, -0.16, 0.32,
                 -0.22, 0.51, 0.59, 0.59, 0.26, -0.04, 0.11, 0.32, 0.75, 0.67,
                 1.05, 0.78, 0.97, 0.58, 0.32, 0.93, 0.08, 0.55, 0.49, 0.26, 0.75,
                 0.04, 1.07, 0.51, 0.71, 0.22, 0.82, 1.03, 1.18, 0.19, 0.21, 0.23,
                 -0.51, 0.57, 0.98, 0.59, 0.44, 0.74, 0.82, 0.54, 1.53, 0.14,
                 0.28, 0.29, 0.82, 1.18, 0.02, 1.58, 0.46, 0.81, -0.17, 1.05,
                 0.55, 1.23, 0.75, -0.18, 0.37, 0.22, 0.36, 0.32, 0.64, 0.41,
                 0.63, 0.25, 0.67, 1.7, -0.01, -0.08, 0.41, 0.85, 1.13, 0.6, 0.67,
                 -0.13, 0.82, -0.24, 0.64, 0.73, 1.23), dim = c(25L, 9L))

Z <- structure(c(0.25, 0.3, 0.35, 0.4, 0.45, 0.45, 0.5, 0.55, 0.6,
                 0.65, 0.65, 0.7, 0.75, 0.8, 0.85, 0.85, 0.9, 0.95, 1, 1.05, 1.05,
                 1.1, 1.15, 1.2, 1.25, -0.94, -0.37, 1.76, 0.23, 0.35, 0.09, 2.93,
                 -0.08, 0.87, 1.46, -0.83, 0.44, 0.72, 0.44, -1.9, 1.51, 2.16,
                 0.46, 3.34, 0.6, 1.31, 0.61, 0.43, 0.47, 2.29, -1.64, 0.89, 0.92,
                 -1.67, 0.39, 1.38, 3.15, 0.64, 1.1, 2.26, -2.02, 1.65, -2.12,
                 2.15, -1.11, 1.23, 2.12, 0.32, 3.74, -1.22, 1.11, 0.14, 0.88,
                 1.72, 0.47, -1.72, -1.15, -0.05, -1.17, -0.82, 0.49, 2.13, 0.24,
                 0.35, 2.56, -2.87, 2.59, -3.35, 0.52, -1.04, 2.05, 3.46, -1.43,
                 4.38, 0.1, 0.04, -0.95, 2.48, 1.93, 0.66, -1.62, -1.69, 0.62,
                 -2.17, -0.65, 1.6, 1.03, 1.25, 0.08, 0.9, -2.75, 1.43, -5.22,
                 0.46, -1.48, 1.76, 5.48, -1.57, 5.17, -0.18, 0.67, -1.86, 1.2,
                 1.95, 2.55, -1.82, -3.44, -0.05, -1.78, -0.03, 2.14, 0.74, 1.44,
                 1.91, 0.14, -2.01, 2.21, -5.97, 1.64, 0.14, 2.09, 5.66, -2.33,
                 4.27, -1.34, 0.45, -2.31, 3.64, 0.71, 3.07, -1.52, -3.19, 0.13,
                 -2.64, 0.02, 2.34, 1.88, 1.28, 3.78, -0.27, -0.56, 3.82, -4.92,
                 0.27, 0.2, 2.13, 4.88, -1.07, 5.83, -0.69, 1.64, -1.31, 4.76,
                 1.09, 3.84, -3.03, -3.48, 1.43, -5.12, -0.21, 3.41, 0.74, 0.93,
                 3.93, 0.1, -0.42, 5.1, -5.74, -2.42, 0.09, 2.06, 4.81, 0.59,
                 4.87, -0.45, 1.61, -1.51, 4.38, 0.59, 3.69, -2.12, -3.68, 1.43,
                 -6.09, -1.7, 2.94, 0.4, 0.51, 2.37, 0.39, -0.72, 6.56, -5.11,
                 -2.05, 0.24, 1.03, 4.67, 1.53, 4.37, -0.11, 1.7, -2.98, 3.53,
                 1.82, 5.02), dim = c(25L, 9L))

# fit AR, showing all output
AR_map <- fitAR_map(X, coords, formula = y ~ t, resids.only = FALSE)

# with only residuals
AR_map_resids <- fitAR_map(X, coords, formula = y ~ t, resids.only = TRUE)

# fit AR with temporal and spatiotemporal predictors
AR_map2 <- fitAR_map(X, coords, formula = y ~ t + Z, X.list = list(t = 1:ncol(X),
                                                                   Z = Z), resids.only = FALSE)

## Tests
test_that("fitAR_map results are correct", {
  # values pulled from dput()
  expect_equal(AR_map$coefficients,
               structure(c(-0.07, 0.32, 0.11, 0.27, 0.32, 0.58, 0.22, 0.26,
                           0.4, 0.37, 0.21, 0.18, 0.28, 0.48, -0.12, 0.41, 0.37, 0.49, 0.24,
                           0.29, 0.16, 0.3, 0.34, 0.32, 0.44, 0.06, -0.03, 0.03, -0.02,
                           0.03, 0.02, 0.06, 0.03, 0.04, 0.03, 0.05, 0.17, -0.04, -0.03,
                           0.06, 0.05, 0.09, 0.01, 0.11, -0.01, 0.07, -0.06, 0.04, 0.05,
                           0.08),
                         dim = c(25L, 2L),
                         dimnames = list(NULL, c("(Intercept)", "t"))),
               tolerance = .1)
})

test_that("fitAR_map resids.only works properly", {
  expect_equal(AR_map$residuals, AR_map_resids$residuals)
  expect_null(AR_map_resids$coefficients)
})

test_that("fitAR_map methods work right", {
  expect_equal(AR_map$coefficients, coefficients(AR_map))
  expect_equal(AR_map$residuals, residuals(AR_map))
  expect_equal(AR_map$fitted.values, fitted(AR_map))
})

# ---- CLS ----

CLS <- fitCLS(x ~ t + Z, data = df)

test_that("fitCLS works correctly", {
  expect_equal(CLS$coefficients,
               c(`(Intercept)` = 0.122, x.0 = 1.865, t = -0.058, Z = -0.342),
               tolerance = 0.01)
})

