#ifndef FUNC_DECL
#define FUNC_DECL

#include "remoteSTAR_types.hpp"

/*
 * Helper Functions ----
 */

// t(A) %*% A
MatrixXd AtA(const MatrixXd& A);
// solve Ax + B
MatrixXd solve_cpp(const MatrixXd& A, const MatrixXd& B);
// solve Ax + I
MatrixXd solve_ident_cpp(const MatrixXd& A);

/*
 * Main Functions
 */

// calculate inverse cholesky matrix
MatrixXd invchol_cpp(const MapMatd& V, double nugget);
// fit GLS
List fitGLS_cpp(const MapMatd& X, const MapMatd& V, const MapMatd& y,
                const MapMatd& X0, double nugget, bool save_xx,
                const int threads);

/*
 * Test Functions ----
 */

// Function to test multicore multiplication

// function to test openmp matrix multiplication
//' Matrix multiplication, possibly with multiple cores
//'
//' @param A numeric matrix
//' @param B numeric matrix
//' @param cores integer number of cores to use
//'
// [[Rcpp::export]]
MatrixXd MatMult(const MapMatd& A, const MapMatd& B, int cores);

#endif // FUNC_DECL
