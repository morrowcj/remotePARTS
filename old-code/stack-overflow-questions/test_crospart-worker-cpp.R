# setwd("code/stack-overflow-questions/")

# test_crosspart-worker-cpp.R

# setup ----
## source relevant file
  # normally build and load package instead of sourceCpp()
# if(!exists("crosspart_worker_cpp")){
  Rcpp::sourceCpp("crosspart-worker.cpp")
# }

## make reproducible
set.seed(75)

## set parameters
n = 500 # rows in original model matrix X
p = 32 # columns in original model matrix X
p0 = 1 # collumns in null model matrix X0

## calculate degrees of freedom
df2 = n - p - 1
df1 = n - p0 - 1 - df2

# Generate dummy data with arbitrary contents, but correct structure ----
## varcovar matrix
Vij <- matrix(abs(rnorm(2*n * 2*n)), nrow = 2*n)

## create list with only the components needed for GLS_worker_cpp()
make_example_list <- function(){
  L <-  list(
    xx = matrix(rnorm(n * p), nrow = n), # random matrix
    xx0 = matrix(rnorm(n * p0), nrow = n), # random null matrix
    # tInvCholV would normally this would be the inversion of a cholesky matrix
    tInvCholV = matrix(rnorm(n * n), nrow = n), # random matrix
    nugget = runif(n = 1, min = 0, max = 1) # random nugget
    )
  return(L)
}

## generate the lists
Li <- make_example_list()
Lj <- make_example_list()

# use the function ----
result <- crosspart_worker_cpp(xxi = Li$xx, xxj = Lj$xx,
                               xxi0 = Li$xx0, xxj0 = Lj$xx0,
                               tUinv_i = Li$tInvCholV,
                               tUinv_j = Lj$tInvCholV,
                               nug_i = Li$nugget,
                               nug_j = Lj$nugget,
                               Vij = Vij,
                               df1 = df1, df2 = df2)
