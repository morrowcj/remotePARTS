# Tests for invert_chol

## setup ----
M <- crossprod(matrix(data = 1:6, ncol = 2))
nug <- 0.25
M2 <- (1 - nug) * M + nug * diag(nrow(M))

## test helper function check_posdef() ----
expect_equivalent(check_posdef(matrix(1:9, 3)),
                  c(TRUE, FALSE, FALSE),
             info = "non positive definite input")
expect_equivalent(check_posdef(matrix(c(1,4,7,4,5,8,7,8,9), 3)),
                               c(TRUE,TRUE,FALSE),
                  info = "square and symmetric input")
expect_equivalent(check_posdef(M), c(TRUE, TRUE, TRUE),
                  info = "positive definite input")

## correct calculation ----
expect_equal(invert_chol(M),
             t(backsolve(chol(M), diag(nrow(M)))),
             info = "calc. w/no nugget")
# with nugget
expect_equal(invert_chol(M, nugget = nug),
             t(backsolve(chol(M2), diag(nrow(M2)))),
             info = "calc. w/nugget")

## error handling ----
expect_error(invert_chol(as.data.frame(M)),
             "M is not of class 'matrix'",
             info = "df input")
expect_error(invert_chol(matrix(1:9, 3),
                         "M is not of type 'double'"),
             info = "integer matrix input")
expect_error(invert_chol(crossprod(matrix(1:9, ncol = 3))),
             "M is not positive definite",
             info = "non positive definite input")


## invert_cholR
expect_true(FALSE, info = "need to depricate invert_cholR()")
