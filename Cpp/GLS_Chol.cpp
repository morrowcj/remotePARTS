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
using Eigen::VectorXi;

typedef Map<MatrixXd> MapMatd;
typedef Map<MatrixXi> MapMati;
typedef Map<VectorXd> MapVecd;

//==============================================================================
/* AtA(A) ----
 * inline code that computes t(A) %*% A
 */

// [[Rcpp::export]]
inline MatrixXd AtA(const MatrixXd& A) {
  int n(A.cols());
  return MatrixXd(n,n).setZero().selfadjointView<Lower>()
                      .rankUpdate(A.adjoint());
}

//==============================================================================
/* solve_cpp(A, B) ----
 * solves Ax = B for x.
 */

// [[Rcpp:export]]
inline MatrixXd solve_cpp(const MatrixXd& A, const MatrixXd& B){
  return A.colPivHouseholderQr().solve(B);
}

/***R
## Test the inline function AtA()
M <- matrix(as.double(1:6), ncol = 2)
crossprod(M) # native R
AtA(M) # cpp way
*/

//==============================================================================
/* tinvchol_cpp(V, nugget = 0)
 * computes the choleskey decomposition upper triangle (U) of V
 * and then returns the transpose: t(U)
 *
 * if a nugget is included... (nothing yet)
 */

// [[Rcpp::export]]
MatrixXd tinvchol_cpp(const MapMatd& V, int nugget = 0){
  const int n = V.rows(); //dim of V

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

//==============================================================================
/* fitGLS_cpp(X, V, y, X0, nugget = 0, threads = 1) ----
 *
 */

// [[Rcpp::export]]
List fitGLS_cpp(const MapMatd& X,
                const MapMatd& V,
                const MapMatd& y,
                const MapMatd& X0,
                int nugget = 0,
                const int threads = 1){

  const int nX = X.rows(), pX = X.cols(); // dimensions of X
  const MatrixXd tUinv = tinvchol_cpp(V, nugget); // transpose of chol(V) = t(inv(U))
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
  const MatrixXd varXinv = varX.colPivHouseholderQr().solve(
    MatrixXd::Identity(varX.rows(), varX.cols()));
  const MatrixXd varcov = varXinv.array() * MSE.array()[0];
  // cout << "\ndiag(vcov):\n" << varcov.matrix().diagonal() << endl; // this works

  // Vectorized operations
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

  double logLik = -0.5 * (nX * log(2 * M_PI) + nX *
                          log(dft * MSE.array()[0]/nX) + logDetV + nX);

  /* Null Model ----
   * The presence/absence of X0 should be handled entirely by R.
   * fitGLS_cpp() REQUIRES an X0 input.
   */

  const MatrixXd xx0 = tUinv * X0;
  const MatrixXd varX0 = AtA(xx0); // crossprod(xx0)
  const MatrixXd X0tY(xx0.adjoint() * yy); // crossprod(xx0, yy)
  const VectorXd betahat0(solve_cpp(varX0, X0tY));
  int df0 = betahat0.size();
  const VectorXd SSE0 = AtA(yy - xx0 * betahat0); // SSE
  double MSE0 = SSE0.array()[0]/(nX - xx.cols()); // MSE
  double MSR = (SSE0.array()[0] - SSE.array()[0])/(xx.cols() - xx0.cols());
  double logLik0 = -0.5 * (nX * log(2 * M_PI) + nX * log((nX - df0) * MSE0/nX) +
                           logDetV + nX);

  /* F test ----
   *
   */
  const double FF = (nX - xx.cols())/(xx.cols() - xx0.cols()) *
    (SSE0.array()[0] - SSE.array()[0]) / SSE.array()[0];
  VectorXi dfF(2);
  dfF(0) = xx.cols() - xx0.cols();
  dfF(1) = nX - xx.cols();

  // return a list of all needed values
  return List::create(Named("betahat") = betahat,
                      Named("VarX") = varX,
                      Named("SSE") = SSE,
                      Named("MSE") = MSE,
                      Named("varcov") = varcov.matrix(),
                      Named("SE") = se,
                      Named("tstat") = tstat,
                      Named("pval.t") = NA_REAL,
                      Named("dft") = dft,
                      Named("logDetV") = logDetV,
                      Named("logLik") = logLik,
                      Named("betahat0") = betahat0,
                      Named("SSE0") = SSE0,
                      Named("MSE0") = MSE0,
                      Named("MSR") = MSR,
                      Named("df0") = df0,
                      Named("logLik0") = logLik0,
                      Named("Fstat") = FF,
                      Named("pval.F") = NA_REAL,
                      Named("df.F") = dfF
                      );
}

/*** R
## test fitGLS_cpp()
load("../R/vignettes-and-examples/test-gls.rda", verbose = FALSE)
Xnull <- matrix(as.double(1), ncol = 1, nrow = nrow(X.small))
source("../R/fitGLS.R")

tmp.cpp <- fitGLS_cpp(X.small, V.small, y.small, Xnull)

tmp.R <- fitGLS(X.small, V.small, y.small)

## compare results
check.equal <- function(string){
  all.equal(unlist(unname(tmp.R[[string]])), unlist(tmp.cpp[[string]]))
}

sapply(c("betahat","VarX", "SSE", "MSE", "varcov", "SE", "t.stat",
         "df.t", "logDetV", "logLik", "betahat0", "SSE0", "MSE0", "MSR",
         "df0", "logLik0", "FF", "df.F"),
       check.equal)

*/


