// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]
#include <iostream>
#include <RcppEigen.h>
#include <math.h>
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

/* AtA(A)
 * inline code that computes t(A) %*% A
 */

// [[Rcpp::export]]
inline MatrixXd AtA(const MatrixXd& A) {
  int n(A.cols());
  return MatrixXd(n,n).setZero().selfadjointView<Lower>()
                      .rankUpdate(A.adjoint());
}

/***R
## Test the inline function AtA()
M <- matrix(as.double(1:6), ncol = 2)
crossprod(M) # native R
AtA(M) # cpp way
*/

/* cholCpp()
 * Example directly RcppEigen that computes chol(crossprod(A))
 * directly uses the AtA() function
 */

// [[Rcpp::export]]
List cholCpp(NumericMatrix& AA){
  const LLT<MatrixXd> llt(AtA(as<MapMatd>(AA)));

  return List::create(Named("L") = MatrixXd(llt.matrixL()),
                      Named("R") = MatrixXd(llt.matrixU()));
}

/***R
## Test the RcppEigen cholCpp() function
chol(crossprod(M)) # native R
cholCpp(M) # cpp
*/

/* tinvchol_cpp(V, nugget = 0)
 * computes the choleskey decomposition upper triangle (U) of V
 * and then returns the transpose: t(U)
 *
 * if a nugget is included... (nothing yet)
 */

// [[Rcpp::export]]
MatrixXd tinvchol_cpp(const MapMatd& V, int nugget = 0){
  // const MapMatd V(as<MapMatd>(VV)); //map V
  const int n = V.rows(); //dim of V

  // test
  // const MatrixXd ata(AtA(V));

  if (nugget == 0) {
    // no nugget case
    const LLT<MatrixXd> llt(V); //chol decomp of V: (U)
    const MatrixXd tUinv = llt.matrixU().solve(MatrixXd::Identity(n, n)).adjoint();
    return tUinv;
  } else {
    // with nugget
    const MatrixXd M = (1 - nugget) * V.array() + nugget * MatrixXd::Identity(n,n).array();
    const MatrixXd Mn = M.matrix();
    const LLT<MatrixXd> llt(Mn);
    const MatrixXd tUinv = llt.matrixU().solve(MatrixXd::Identity(n, n)).adjoint();
    return tUinv;
  }
}

/***R
## Test tinvchol_cpp()
V <- crossprod(M)

# no nugget (i.e. nugget = 0)
t(backsolve(chol(V), diag(nrow(V)))) # native R
tinvchol_cpp(V) # cpp

# with nugget
nug = .1
Vn = (1 - nug) * V + nug * diag(nrow(V))
t(backsolve(chol(Vn), diag(nrow(Vn)))) # R
tinvchol_cpp(Vn, nug) # cpp
*/

// [[Rcpp::export]]
List fitGLS_cpp(const MapMatd& X,
                const MapMatd& V,
                const MapMatd& y,
                const Nullable<NumericMatrix> X0_ = R_NilValue,
                const Nullable<NumericVector> y0_ = R_NilValue,
                int nugget = 0,
                const int threads = 1){

  const int nX = X.rows(), pX = X.cols(); // dimensions of X
  const MatrixXd tUinv = tinvchol_cpp(V); // transpose of chol(V) = t(inv(U))
  const MatrixXd xx = tUinv * X; // t(inv(U)) %*% X
  const MatrixXd yy = tUinv * y; // t(inv(U)) %*% y
  const MatrixXd varX = AtA(xx); // crossprod(xx)
  const MatrixXd XtY(xx.adjoint() * yy); // crossprod(xx, yy)


  // solve for betahat (using one specific solver - there are others)
  const VectorXd betahat(varX.colPivHouseholderQr().solve(XtY));

  // calculate some statistics
  int dft = nX - xx.cols();
  const VectorXd SSE = AtA(yy - xx * betahat); // SSE
  const VectorXd MSE = SSE/dft; // MSE
  const MatrixXd varXinv = varX.colPivHouseholderQr().solve(MatrixXd::Identity(varX.rows(), varX.cols()));
  const MatrixXd varcov = varXinv.array() * MSE.array()[0];
  // cout << "\ndiag(vcov):\n" << varcov.matrix().diagonal() << endl; // this works

  VectorXd se = varcov.matrix().diagonal();
  for (int i = 0; i < se.size(); i++){
    se[i] = std::sqrt(se[i]);
  }

  VectorXd tstat(betahat.size());
  for (int i = 0; i < betahat.size(); i++){
    tstat[i] = betahat[i] / se[i];
  }

  double logDetV = 0;
  for (int i = 0; i < tUinv.rows(); i++){
    logDetV += log(tUinv.diagonal()[i]);
  }
  logDetV *= -2;

  double logLik = -0.5 * (nX * log(2 * M_PI) + nX * log(dft * MSE.array()[0]/nX) + logDetV + nX);

  /*
   * Null Model
   */
  // if (y0_.isNotNull()){
  //   const MapVecd y0(as<MapVecd>(y0_));
  // } else {
  //
  // }

  const MatrixXd X0;
  if (X0_.isNotNull()){
    const MapMatd X0(as<MapMatd>(X0_));
  } else {
    const MatrixXd X0 = MatrixXd(nX,1).setOnes();
  }
  const MatrixXd xx0 = tUinv * X0;
  const MatrixXd varX0 = AtA(xx); // crossprod(xx0)
  const MatrixXd X0tY(xx0.adjoint() * yy); // crossprod(xx0, yy)
  // const VectorXd betahat0(varX0.colPivHouseholderQr().solve(X0tY)); // This line is breaking it
  // int df0 = betahat0.size();
  // const VectorXd SSE0 = AtA(yy - xx0 * betahat0); // SSE
  // double MSE0 = SSE0.array()[0]/(nX - xx.cols()); // MSE
  // double MSR = (SSE0.array()[0] - SSE.array()[0])/(xx.cols() - xx0.cols());
  // double logLik0 = -0.5 * (nX * log(2 * M_PI) + nX * log(dft * MSE0/nX) + logDetV + nX);

  // cout <<"\nsqrt(diag(vcov)):\n"<< se << endl;

  // se = vec1.sqrt();
  // const VectorXd se = sediag.sqrt();
  // const VectorXd tstat = betahat.array() / se.array()[0];

  // return a list of all needed values
  return List::create(Named("betahat") = betahat,
                      Named("VarX") = varX,
                      Named("SSE") = SSE,
                      Named("MSE") = MSE,
                      Named("varcov") = varcov.matrix(),
                      Named("SE") = se,
                      Named("tstat") = tstat,
                      Named("dft") = dft,
                      Named("logDetV") = logDetV,
                      Named("logLik") = logLik
                      // Named("betahat0") = betahat0,
                      // Named("SSE0") = SSE0,
                      // Named("MSE0") = MSE0,
                      // Named("MSR") = MSR,
                      // Named("df0") = df0
                      // Named("df0") = logLik0
                      );
}

/*** R
## test fitGLS_cpp()
load("../R/vignettes-and-examples/test-gls.rda", verbose = FALSE)
source("../R/fitGLS.R")

tmp.cpp <- fitGLS_cpp(X.small, V.small, y.small)

tmp.R <- fitGLS(X.small, V.small, y.small)

## compare results
check.equal <- function(string){
  all.equal(unlist(unname(tmp.R[[string]])), unlist(tmp.cpp[[string]]))
}

sapply(c("betahat","VarX", "SSE", "MSE", "varcov", "SE", "t.stat",
         "df.t", "logDetV", "logLik"),
       check.equal)

*/


