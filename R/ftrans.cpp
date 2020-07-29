// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

using namespace Rcpp;
using Eigen::Map;
using Eigen::MatrixXd;
// using Rcpp::as;

typedef Map<MatrixXd> MapMatd;

// [[Rcpp::export]]
MatrixXd ftrans_cpp(const NumericMatrix& XX){
  // Map the numeric matrix XX from R
  const MapMatd A(as<MapMatd>(XX));
  // const Map<MatrixXi> A(as<Map<MatrixXi> >(AA));
  // const MatrixXd At(A.transpose());
  return A.transpose();
}
// evaluate and return the transpose of A



/* This next function is an alternative version of the matrix transformation
 * function.
 */

// [[Rcpp::export]]
MatrixXd ftrans_cpp_alt(Map<MatrixXd> X){
  return X.transpose();
}
// ftrans_cpp_alt(MatrixXd X) would make a copy of X. the Map<> is a pointer.

/* In the first function the '&' in the argument signifies a pointer but
 * 'NumericMatrix' is an Rcpp object while MatrixXd is an RcppEigen object.
 * The pointers are mapped differently (syntax)
 */
