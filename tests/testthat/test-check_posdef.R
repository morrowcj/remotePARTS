test_that("positive definitive assessed correctly", {

  # Square matrix, not symmetric, not positive definitive
  sqr_only = matrix(c(1, 2, 0, -3), nrow = 2, byrow = TRUE)
  expect_equal(check_posdef(sqr_only),
               c(TRUE, FALSE, FALSE),
               ignore_attr = TRUE)

  # Square and symmetric matrix, not pos def
  sqr_sym = matrix(c(1, 2, 2, 1), nrow = 2, byrow = TRUE)
  expect_equal(check_posdef(sqr_sym),
               c(TRUE, TRUE, FALSE),
               ignore_attr = TRUE)

  # positive definitive
  pos_def = matrix(c(2, -1, 0, -1, 2, -1, 0, -1, 2), nrow = 3, byrow = TRUE)
  expect_true(all(check_posdef(pos_def)))
})

test_that("correct structure of check_posdef output", {
  M = matrix(runif(9), nrow = 3)
  res = check_posdef(M)

  # correct dimensions
  expect_equal(length(res), 3)

  # column names are right
  expect_equal(names(res), c("sqr", "sym", "posdef"))

  # logical data type
  expect_true(is.logical(res) )
})

test_that("error handling of check_posdef", {
  M = as.data.frame(matrix(c(2, -1, 0, -1, 2, -1, 0, -1, 2), nrow = 3, byrow = TRUE))
  expect_error(check_posdef(M), "M must be a matrix")
})
