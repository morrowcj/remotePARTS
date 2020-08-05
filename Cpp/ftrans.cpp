// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]
#include <RcppEigen.h>
#include <iostream>
#include <omp.h>

using namespace std;
using namespace Rcpp;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;
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

// [[Rcpp::export]]
MatrixXd test_solve(const MapMatd& X, const MapMatd& Y){
  // cout<<"X:\n"<<X<<endl;
  // cout<<"Y:\n"<<Y<<endl;

  // cout<<"Solution to Xb = Y:\n"<<X.colPivHouseholderQr().solve(Y)<<endl;

  cout<<"dim(X):\n"<<X.rows()<<", "<<X.cols()<<endl;
  cout<<"dim(Y):\n"<<Y.rows()<<", "<<Y.cols()<<endl;

  // MatrixXd a(1, 10);
  MatrixXd A = MatrixXd(10, 1).setOnes();
  cout<<"ones:\n"<<A<<endl;

  // MatrixXd A(3,3);
  // cout<<"A:\n"<<A<<endl;
  // MatrixXd B(3,3);
  // cout<<"B:\n"<<B<<endl;
  //
  // VectorXd V1(3);
  // cout<<"V1:\n"<<V1<<endl;
  // VectorXd V2(3);
  // cout<<"V2:\n"<<V1<<endl;
  //
  // VectorXd a(1);
  // cout<<"a:\n"<<a<<endl;
  // VectorXd b(1);
  // cout<<"b:\n"<<b<<endl;

  // cout<<"solution Ax = B:\n"<<A.colPivHouseholderQr().solve(B)<<endl;
  // cout<<"solution Ax = V1:\n"<<A.colPivHouseholderQr().solve(V1)<<endl;
  // cout<<"solution V1x = V2:\n"<<V1.colPivHouseholderQr().solve(V2)<<endl;
  // cout<<"solution ax = b:\n"<<a.colPivHouseholderQr().solve(b)<<endl;

  // return 0;
  return X.colPivHouseholderQr().solve(Y);
}

// [[Rcpp::export]]
MatrixXd solve_cpp(const MapMatd& X, const MapMatd Y){
  return test_solve(X, Y);
}

/***R
set.seed(798)

## 3x3 and 3x3
X = matrix(as.double(rnorm(9)), ncol = 3)
b <- matrix(rnorm(9), ncol = 3)
Y = X %*% b

test_solve(X, Y)

## 3x3 and 3x1
b <- matrix(rnorm(3), ncol = 1)
Y = X %*% b

test_solve(X, Y)

## 1x1 and 1x1
X = matrix(rnorm(1))
b = matrix(rnorm(1))
Y = X %*% b

test_solve(X, Y)

X <- matrix(as.double(1), nrow = 3, ncol = 1)
Y = X %*% b

test_solve(X, Y)

solve_cpp(X, Y)

*/

/* This next function is an alternative version of the matrix transformation
 * function.
 */

// [[Rcpp::export]]
MatrixXd ftrans_cpp_alt(const Map<MatrixXd> X){
  return X.transpose();
}
// ftrans_cpp_alt(MatrixXd X) would make a copy of X. the Map<> is a pointer.

/* In the first function the '&' in the argument signifies a pointer but
 * 'NumericMatrix' is an Rcpp object while MatrixXd is an RcppEigen object.
 * The pointers are mapped differently (syntax)
 */
