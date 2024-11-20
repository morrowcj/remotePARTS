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
