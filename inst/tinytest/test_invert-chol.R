# Tests for invert_chol

## setup ----
M <- crossprod(matrix(data = 1:6, ncol = 2))
nug <- 0.25
M2 <- (1 - nug) * M + nug * diag(nrow(M))

## correct calculation ----
# no nugget
expect_equal(invert_chol(M),
             t(backsolve(chol(M), diag(nrow(M)))))
expect_equal(invert_chol(as.data.frame(M)),
             t(backsolve(chol(M), diag(nrow(M)))))
# with nugget
expect_equal(invert_chol(M, nugget = nug),
             t(backsolve(chol(M2), diag(nrow(M2)))))

## error handling ----
expect_error(invert_chol(matrix(1:6, ncol = 2)),
             "not square")
expect_error(invert_chol(matrix(1:9, ncol = 3)),
             "not symmetric")
expect_error(invert_chol(crossprod(matrix(1:9, ncol = 3))),
             "not positive definite")
