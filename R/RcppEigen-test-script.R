# Rcpp Eigen Test Script

library(Rcpp);library(RcppEigen);library(inline)
library(ggplot2);library(microbenchmark)

# matrix transformation ----

sourceCpp("R/ftrans.cpp")
# sourceCpp("R/ftrans-alt.cpp")

nr = 1000
nc = 1000

A <- matrix(as.numeric(1:(nr*nc)), ncol = nc)

stopifnot(all.equal(ftrans_cpp(A), t(A)))

## Benchmark results

test.trans <- microbenchmark(trans.R = t(A),
                             trans.Cpp = ftrans_cpp(A),
                             transalt.Cpp = ftrans_cpp_alt(A),
                             times = 1000L)
test.trans # base R wins by a lot
autoplot(test.trans)


# Test a leverage function (matrix operations)----
leverage <- function(X) {
  diag(X %*% solve(t(X) %*% X) %*% t(X))
}

leverage.cppEigen <- sourceCpp("R/CppSamp.cpp")


data("attitude", package = "datasets")
model <- lm(rating ~ ., data = attitude)
X <- model.matrix(model)

stopifnot(all.equal(unname(leverage(X)), leverage_cpp(X)))

## Benchmark the results.
test <- microbenchmark(
  leverage_r         = leverage(X),
  leverage_cpp       = leverage_cpp(X),
  # leverage_cpp_short = leverage_cpp_short(X),
  times = 1000L
)

test # Rcpp wins by a mile
autoplot(test)
