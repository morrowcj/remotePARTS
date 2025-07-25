test_that("distm_x works correctly", {
  # 3 lat-long pairs
  coords = data.frame(longitude = c(-89.377808, -89.727778, -89.852591),
                      latitude = c(  43.053619,  43.426389,  43.028171))

  # observed distances
  dists_observed = distm_km(coords)
  dists_ut = dists_observed[upper.tri(dists_observed)]

  # check against known distances (rounded)
  expect_equal(round(dists_ut, 4), c(50.2305, 38.7917, 45.3880))
  # check that the max distance attribute is correct
  expect_equal(round(attr(dists_observed, "max.dist"), 4), 50.2305)

  # scaled distances
  dists_scale = distm_scaled(coords)
  scaled_ut = dists_scale[upper.tri(dists_scale)]

  expect_equal(round(scaled_ut, 4), round(c(50.2305, 38.7917, 45.3880)/50.2305, 4))
  expect_equal(round(attr(dists_scale, "max.dist"), 4), 50.2305)
})

test_that("max distance works correctly", {
  coords = data.frame(longitude = c(-89.377808, -89.727778, -89.852591),
                      latitude = c(  43.053619,  43.426389,  43.028171))

  expect_equal(round(max_dist(coords), 4), 50.2305)
  expect_equal(max_dist(coords, "distm_scaled"), 1)
})

test_that("taper covar works correctly",{
  d = seq(0, 0.5, by = 0.1)

  expect_equal(round(covar_taper(d, 1), 2), c(1.00, 0.85, 0.70, 0.56, 0.43, 0.31))
  expect_equal(round(covar_taper(d, 0.3), 2), c(1.00, 0.52, 0.15, 0.00, 0.00, 0.00))
})

test_that("exp covar works correctly", {
  d = seq(0, 0.5, by = 0.1)

  expect_equal(round(covar_exp(d, .1), 2), c(1.00, 0.37, 0.14, 0.05, 0.02, 0.01))
  expect_equal(round(covar_exp(d, .2), 2),  c(1.00, 0.61, 0.37, 0.22, 0.14, 0.08))
})

test_that("exppow covar works correctly", {
  d = seq(0, 0.5, by = 0.1)

  expect_equal(round(covar_exppow(d, .1, 1), 2), round(covar_exp(d, .1), 2))
  expect_equal(round(covar_exppow(d, .1, .2), 2), c(1.00, 0.37, 0.32, 0.29, 0.27, 0.25))
})
