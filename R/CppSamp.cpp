// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

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

