## Test MatMult function

A <- matrix(rnorm(9), nrow = 3)
B <- matrix(rnorm(9), nrow = 3)

expect_equal(MatMult(A, B, 1), A %*% B)
