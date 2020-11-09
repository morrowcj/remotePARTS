# Rcpp Eigen Test Script

## These commands are automatically set by [[Rcpp::plugins(openmp)]] in C++ file
# Sys.setenv("PKG_CXXFLAGS" = "-fopenmp")
# Sys.setenv("PKG_LIBS" = "-fopenmp")

library(Rcpp);library(RcppEigen);library(inline)
library(ggplot2);library(microbenchmark)

# matrix transformation ----

sourceCpp("code/Cpp/ftrans.cpp") # contains C++ transform functions

## Create test Matrix
nr = 1000; nc = 1000

A <- matrix(as.numeric(1:(nr*nc)), ncol = nc)

## compare functions numerically
stopifnot(all.equal(ftrans_cpp(A), t(A)))

## Benchmark results
test.trans <- microbenchmark(trans.R = t(A),
                             trans.Cpp = ftrans_cpp(A),
                             transalt.Cpp = ftrans_cpp_alt(A),
                             times = 1000L)
test.trans # base R wins by a lot
autoplot(test.trans)


# Test a leverage function (matrix operations)----

## R version
leverage <- function(X) {
  diag(X %*% solve(t(X) %*% X) %*% t(X))
}
## C++ version
sourceCpp("code/Cpp/CppSamp.cpp") # contains C++ leverage and coretest functions

## Load some test data
data("attitude", package = "datasets")
model <- lm(rating ~ ., data = attitude)
X <- model.matrix(model)

## compare numerically
stopifnot(all.equal(unname(leverage(X)), leverage_cpp(X)))

## Benchmark the results.
test <- microbenchmark(
  leverage_r         = leverage(X),
  leverage_cpp       = leverage_cpp(X),
  leverage_cpp_par   = leverage_cpp_par(X, 1), #one thread
  leverage_cpp_par5   = leverage_cpp_par(X, 5), #five threads
  # leverage_cpp_short = leverage_cpp_short(X),
  times = 1000L
)
test # Rcpp wins by a mile
autoplot(test)

## Test multiple cores with matrix multiplication of 500x500 matrices
test2 <- microbenchmark(one = corestest(1L),
                        three = corestest(3L),
                        five = corestest(5L),
                        seven = corestest(7L),
                        eight = corestest(8L),
                        times = 1000L)
test2 # Eigen is making use of multiple cores!!
autoplot(test2)

## (though the diff between 5 and 10 cores is negligible)
## for my machine 5 cores seems to be optimal for this type of operation
## It turns out that omp uses 'threads' instead of cores. However, my machine
## has 8 cores with 2 logical processors (threads) each. So maybe I'm mistaken.
## Perhaps it can only utilize each core, which is why it bottoms out at 7-8?
