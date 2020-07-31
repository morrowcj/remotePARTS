# Test GLS_Chol.cpp

# n <- 40
#
# V.small <- V[1:n, 1:n]
# X.small <- Xmat[1:n, ]
# y.small <- rowSums(1.9 * V.small %*% X.small)
#
# ichol <- t(backsolve(chol(V.small), diag(n)))
#
# XX.small <- ichol %*% X.small
# yy.small <- ichol %*% y.small
#
# crossX <- crossprod(XX.small)
# crossY <- crossprod(XX.small, yy.small)
#
# b <- as.numeric(solve(crossX, crossY))
#
# save(V.small, X.small, y.small, ichol, XX.small, crossX, crossY, b,
#         file = "R/vignettes-and-examples/test-gls.rda")


load("R/vignettes-and-examples/test-gls.rda", verbose = TRUE)

Rcpp::sourceCpp("Cpp/GLS_Chol.cpp")

GLS_chol(X.small, V.small, y.small, 1)

solveone <- function(X, y, V){
  invcholV <- invert_cholR(V)
  xx <- invcholV %*% X
  yy <- invcholV %*% y
  coef <- as.numeric(solve(crossprod(xx), crossprod(xx,yy)))
}

(bench <- microbenchmark::microbenchmark(
  R = (beta.R <- solveone(X.small, y.small, V.small)),
  Cpp = (beta.C <- GLS_chol(X.small, V.small, y.small)),
  Cpp.5 = (beta.C5 <- GLS_chol(X.small, V.small, y.small, threads = 5)),
                               times = 1000L))
stopifnot(all.equal(beta.R, beta.C, beta.C5))
ggplot2::autoplot(bench)

(bigBench <- microbenchmark::microbenchmark(
  R = (beta.R <- solveone(Xmat, y, V)),
  Cpp = (beta.C <- GLS_chol(Xmat, V, y)),
  Cpp.5 = (beta.C5 <- GLS_chol(Xmat, V, y, threads = 5)), times = 10L
))
stopifnot(all.equal(beta.R, beta.C, beta.C5))
ggplot2::autoplot(bigBench)
