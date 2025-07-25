# Setup
## read data
data(ndvi_AK10000)
df = ndvi_AK10000[seq_len(200), ] # first 200 rows

## fit covariance matrix
V = covar_exp(distm_scaled(cbind(df$lng, df$lat)), range = .01)

## run GLS
(GLS = fitGLS(CLS_coef ~ 0 + land, data = df, V = V))

## with F-test calculations to compare with the NULL model
(GLS.F = fitGLS(CLS_coef ~ 0 + land, data = df, V = V, no.F = TRUE))

# Tests
test_that("fitGLS works correctly", {
  expect_equal(GLS$coefficients, c(landShrubland = 0.00165, landSavanna = 3e-05, landGrassland = 4e-05),
               tolerance = 0.001)
  expect_equal(GLS$covar_coef,
               structure(c(0, 0, 0, 0, 0, 0, 0, 0, 0), dim = c(3L, 3L),
                         dimnames = list(c("landShrubland", "landSavanna", "landGrassland"),
                                         c("landShrubland", "landSavanna", "landGrassland"))),
               tolerance = 0.01)
  expect_equal(GLS$pval_F, 0.00342, tolerance = 0.001)
  expect_equal(GLS$pval_t, c(landShrubland = 0, landSavanna = 0.942, landGrassland = 0.946),
               tolerance = 0.01)
})

test_that("fitGLS no.F works", {
  expect_null(GLS.F$pval_F)
})
