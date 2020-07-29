// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]
// #include <Eigen/Core>
#include <RcppEigen.h>
#include <omp.h>
#include <iostream>


using namespace Rcpp;
using Eigen::Map;
using Eigen::VectorXd;
using Eigen::MatrixXd;

typedef Map<MatrixXd> MapMatd;



// [[Rcpp::export]]
VectorXd leverage_cpp(const NumericMatrix& XX) {
  const MapMatd X(as<MapMatd>(XX));
  const MatrixXd H = X * (X.adjoint() * X).inverse() * X.adjoint();
  return H.diagonal();
}

// [[Rcpp::export]]
VectorXd leverage_cpp_par(const NumericMatrix& XX, int const threads = 5) {
  omp_set_num_threads(threads);
  Eigen::initParallel();
  const MapMatd X(as<MapMatd>(XX));
  const MatrixXd H = X * (X.adjoint() * X).inverse() * X.adjoint();
  return H.diagonal();
}

// [[Rcpp::export]]
int corestest(int const threads = 5){
  omp_set_num_threads(threads);
  Eigen::initParallel();
  int n = 500;
  MatrixXd A = MatrixXd::Ones(n,n);
  MatrixXd B = MatrixXd::Ones(n,n);
  MatrixXd C = MatrixXd::Ones(n,n);
  C.noalias() += A*B;
  return C.sum();
}
// omp actually only sets THREADS not cores.
