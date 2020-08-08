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
# ## test fitGLS_cpp()
load("../R/vignettes-and-examples/test-gls.rda", verbose = FALSE)
Xnull <- matrix(as.double(1), ncol = 1, nrow = nrow(X.small))
# source("../R/fitGLS.R")
#
# tmp.cpp <- fitGLS_cpp(X.small, V.small, y.small, Xnull)
#
# tmp.R <- fitGLS(X.small, V.small, y.small)
#
# ## compare results
# check.equal <- function(string){
#   all.equal(unlist(unname(tmp.R[[string]])), unlist(tmp.cpp[[string]]))
# }
#
# sapply(c("betahat","VarX", "SSE", "MSE", "varcov", "SE", "t.stat",
#          "df.t", "logDetV", "logLik", "betahat0", "SSE0", "MSE0", "MSR",
#          "df0", "logLik0", "FF", "df.F"),
#        check.equal)

*/


//==============================================================================
/* LogLikGLS(XX, V, y, nugget)
 *
 */

// [[Rcpp::export]]
inline double LogLikGLS_cpp(const MapMatd& X,
                        const MapMatd& V,
                        const MapVecd& y,
                        int nugget = 0){

  int n = X.rows();

  const MatrixXd tUinv = tinvchol_cpp(V, nugget); // transpose of chol(V) = t(inv(U))
  const MatrixXd xx = tUinv * X; // t(inv(U)) %*% X
  const MatrixXd yy = tUinv * y; // t(inv(U)) %*% y
  const MatrixXd varX = AtA(xx); // crossprod(xx)
  const MatrixXd XtY(xx.adjoint() * yy); // crossprod(xx, yy)

  // solve for betahat (using one specific solver - there are others)
  const VectorXd bhat(varX.colPivHouseholderQr().solve(XtY));

  int df = n - xx.cols();
  const VectorXd SSE = AtA(yy - xx * bhat); // SSE
  const double MSE = SSE.array()[1]/df; // MSE

  // Log(Det(Vn))
  double logDetV = 0;
  for (int i = 0; i < tUinv.rows(); i++){
    logDetV += log(tUinv.diagonal()[i]);
  }
  logDetV *= -2;

  double logLik = -0.5*(n * log(2*M_PI) + n * log(df*MSE/n) + logDetV + n);

  return logLik;
}

/* Optimizer function
 * Currently uses the golden section search:
 *   https://chemicalstatistician.wordpress.com/2013/04/22/using-r-to-implement-the-golden-bisection-method/
 * However, the algorithm used by optimize() is faster:
 *   http://www.netlib.org/fmm/fmin.f
 * I should try and implement this algorithm instead
 */

// [[Rcpp::export]]
inline double nugOptim_cpp(const MapMatd& X,
                const MapMatd& V,
                const MapVecd& y,
                double lower = 0,
                double upper = 1,
                double tol = 0.00001){
  double mx; // declare x which maximise f(x)

  // golden ratio value
  double GR = 2/sqrt(5) + 1;

  // test bounds
  double x1 = upper - GR*(upper - lower);
  double x2 = lower + GR*(upper - lower);

  // current f() bounds
  double f1 = LogLikGLS_cpp(X, V, y, x1);
  double f2 = LogLikGLS_cpp(X, V, y, x2);

  int i = 0;
  while (abs(upper - lower ) > tol){
    i += 1;
    // cout << "current interval: [" << lower << ", " << upper << "]" << endl;
    if (f2 < f1) { // then the maximum is closer to x2
      // recycle new values according to the GR rule
      upper = x2;
      x2 = x1;
      f2 = f1;
      // calculate new values
      x1 = upper - GR*(upper - lower);
      f2 = LogLikGLS_cpp(X, V, y, x2);
    } else { // the maximum is closer to x1
      lower = x1;
      x1 = x2;
      f1 = f2;

      x2 = lower + GR*(upper - lower);
      f2 = LogLikGLS_cpp(X, V, y, x2);
    }
  }
  cout << "number of iterations: " << i << endl;

  mx = (upper + lower)/2;

  if(mx < tol){
    double f0 = LogLikGLS_cpp(X, V, y, 0);
    double fmx = LogLikGLS_cpp(X, V, y, mx);
    if (f0 > fmx){
      mx = 0;
    }
  }

  if(mx > 1 - tol){
    double f0 = LogLikGLS_cpp(X, V, y, 1);
    double fmx = LogLikGLS_cpp(X, V, y, mx);
    if (f0 > fmx){
      mx = 1;
    }
  }

  return mx;
}

/***R
# nugOptim_cpp(X.small, V.small, y.small, 0, 1, .00001) # much slower than using R's optimizer
*/
