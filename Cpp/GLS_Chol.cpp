// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]
#include <iostream>
#include <RcppEigen.h>
// #include <omp.h>

using namespace std;
using namespace Rcpp;

// using Rcpp::as;
using Eigen::LLT;
using Eigen::Lower;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::Upper;
using Eigen::VectorXd;

typedef Map<MatrixXd> MapMatd;
typedef Map<MatrixXi> MapMati;
typedef Map<VectorXd> MapVecd;

// // [Rcpp::export]
// inline MatrixXd AtA(const MapMatd& A) {
//   int n(A.cols());
//   return MatrixXd(n,n).setZero().selfadjointView<Lower>()
//                       .rankUpdate(A.adjoint());
// }


// [[Rcpp::export]]
VectorXd GLS_chol(const NumericMatrix& XX,
                  const NumericMatrix& VV,
                  const NumericVector& yy,
                  const int threads = 1){

  // Threads for Parallel bits (not faster at 3000)
  // omp_set_num_threads(threads);
  // Eigen::initParallel();

  // Map the arguments to Eigen Objects
  const MapMatd X(as<MapMatd>(XX));
  // cout << "X:\n" << X <<endl;
  const MapMatd V(as<MapMatd>(VV));
  // cout << "V:\n" << V <<endl;
  const MapVecd y(as<MapVecd>(yy));
  // cout << "y:\n" << y <<endl;

  // Dimensions of X
  const int nX(X.rows()), pX(X.cols());
  // cout << "nX: " << nX <<" pX: " <<pX <<endl;

  // Dimensions of square matrix V
  const int nV(V.rows());
  // cout << "nV: " << nV << "\n" << endl;

  // chol decomp of V
  const LLT<MatrixXd> llt(V);
  // cout<<"chol:\n"<< MatrixXd(llt.matrixU())<<endl;

  // U inverse
  const MatrixXd UinvT(llt.matrixU().solve(MatrixXd::Identity(nV, nV)).adjoint());
  // cout << "U inverse:\n" << UinvT << endl;

  // Transpose of U inverse * X (and y)
  const MatrixXd XU(UinvT * X);
  // cout << "XU:\n" << XU <<endl;

  const MatrixXd yU(UinvT * y);
  // cout << "yU:\n" << yU <<endl;

  // crossprod(XU)
  const MatrixXd XtX(MatrixXd(XU.cols(), XU.cols()).setZero().selfadjointView<Lower>().
                       rankUpdate(XU.adjoint()));
  // cout << "XtX:\n" << XtX << endl;

  // crossprod(XU, yU)
  const MatrixXd yUy(XU.adjoint() * yU);
  // cout << "yUy:\n" << yUy << endl;
  //
  // Xtx * b = yUy i.e. solve(Xtx, yUy)
  const VectorXd betahat(XtX.colPivHouseholderQr().solve(yUy));
  // cout << "betahat:\n" << betahat << endl;

  //
  // return betahat;

  return betahat;
}
