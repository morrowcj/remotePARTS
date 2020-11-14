#include "function-declarations.hpp"

//' caculate t(A) \%*\% A
//'
//' @param A numeric matrix
//'
// [[Rcpp::export(.AtA)]]
MatrixXd AtA(const MatrixXd& A){
  int n(A.cols());
  return MatrixXd(n, n).setZero().selfadjointView<Lower>()
                       .rankUpdate(A.adjoint());
}

//' solve Ax = B for x
//'
//' @param A numeric matrix
//' @param B numeric matrix
//'
// [[Rcpp::export(.solve_cpp)]]
MatrixXd solve_cpp(const MatrixXd& A, const MatrixXd& B){
  return A.colPivHouseholderQr().solve(B);
}

//' solve Ax = I for x
//'
//' @param A numeric matrix
//'
// [[Rcpp::export(.solve_ident_cpp)]]
MatrixXd solve_ident_cpp(const MatrixXd& A){
  MatrixXd I = MatrixXd::Identity(A.rows(),A.cols());
  return A.colPivHouseholderQr().solve(I);
}
