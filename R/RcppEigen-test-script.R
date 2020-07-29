# Rcpp Eigen Test Script

library(Rcpp);library(RcppEigen);library(inline)

ftrans <- cxxfunction(signature(AA = "matrix"), body = "
  using Eigen::Map;
  using Eigen::MatrixXi;
  // Map the integer matrix AA from R
  const Map<MatrixXi> A(as<Map<MatrixXi> >(AA));
  // evaluate and return the transpose of A
  const MatrixXi At(A.transpose());
  return wrap(At);
  ", plugin = "RcppEigen")

nr = 10000
nc = 10000

A <- matrix(1:(nr*nc), ncol = nc)

times.r <- numeric(10)
times.c <- numeric(10)
for(i in 1:10){
  times.r[i] <- system.time(At.R <- t(A))[3]
  times.c[i] <- system.time(At <- ftrans(A))[3]
}
cbind(R = summary(times.r), Cpp = summary(times.c))

stopifnot(all.equal(At.R, At))


# Test a leverage function
leverage <- function(X) {
  diag(X %*% solve(t(X) %*% X) %*% t(X))
}

leverage.cppEigen <- sourceCpp("R/CppSamp.cpp")
