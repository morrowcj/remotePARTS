#include "function-declarations.h"
#include <omp.h>
// [[Rcpp::plugins(openmp)]]]


MatrixXd MatMult(const MapMatd& A, const MapMatd& B, int cores){
  if (cores > 1) {
    // Multicore functionality: works on windows 10 with src/Makevars.win
    omp_set_num_threads(cores);
  }
  return A * B;
}

/*** R
A <- matrix(rnorm(1000*1000), nrow = 1000, ncol = 1000)
B <- matrix(rnorm(1000*1000), nrow = 1000, ncol = 1000)
all.equal(MatMult(A, B, 1), (A %*% B))

# A <- matrix(rnorm(10000*10000), nrow = 10000, ncol = 10000)
# B <- matrix(rnorm(10000*10000), nrow = 10000, ncol = 10000)
# (bench <- microbenchmark::microbenchmark(one = MatMult(A, B, 1), two = MatMult(A, B, 2), five = MatMult(A, B, 5), times = 10L));autoplot(bench)

## Using multiple cores is NOT a linear relationship. There is communication
## overhead associated with windows multiprocessing. Only use if you have
## enough cores to leverage their benefit (>= 4).

# Unit: seconds
# expr      min       lq     mean   median       uq      max neval cld
# one  23.75706 24.48346 39.98477 41.29005 54.90878 55.01790    10  b
# two  54.75272 54.88268 55.22794 54.95965 55.93550 56.00776    10   c
# five 23.32917 24.09156 24.69750 24.49658 25.51383 26.23372    10 a
*/
