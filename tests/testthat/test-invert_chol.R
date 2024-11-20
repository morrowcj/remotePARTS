test_that("cholesky inversion is correct", {
  # positive definitive matrix
  M = matrix(c(2, -1, 0, -1, 2, -1, 0, -1, 2), nrow = 3, byrow = TRUE)

  # traditional solution
  chsky = chol(M)
  inv_chol = t(solve(chsky))

  # our invert_chol solution
  result = invert_chol(M)

  # expectation
  expect_equal(result[upper.tri(result, diag = TRUE)],
               inv_chol[upper.tri(inv_chol, diag = TRUE)])
})

test_that("invert_chol handles errors correctly", {
  # positive definitive matrix
  M = matrix(c(2, -1, 0, -1, 2, -1, 0, -1, 2), nrow = 3, byrow = TRUE)

  # non-matrix
  expect_error(invert_chol(as.data.frame(M)), "M is not of class 'matrix'")

  # integers
  Mint = matrix(as.integer(M), nrow = 3, byrow = TRUE)
  expect_error(invert_chol(Mint), "M is not of type 'double'")

  # non-positive definitive matrix
  sqr_sym = matrix(c(1, 2, 2, 1), nrow = 2, byrow = TRUE)
  expect_error(invert_chol(sqr_sym), "M is not positive definitive")
})
