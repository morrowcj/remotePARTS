// [[Rcpp::depends(RcppEigen)]]
// // [[Rcpp::interfaces(r, cpp)]]
#include <iostream>
#include <math.h>

// #include "remoteSTAR_types.h"
#include <RcppEigen.h>

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

// // [[Rcpp::plugins(openmp)]]
// #include <omp.h>

// using namespace Eigen;
// using Rcpp::as;

using namespace std;
using namespace Rcpp;

//==============================================================================
/* DistGeo() ----
 * calculate geographic distance matrix
 *
 * So far, this isn't working. I think it is becaus Rcpp::sourceCpp() doesn't
 * like that I'm using the external GeographicLib library.
 * https://stackoverflow.com/questions/50354302/rcpp-sourcecpp-undefined-symbol
 *
 * It may be best to use header files from another Rcpp package such as nngeo:
 * https://github.com/michaeldorman/nngeo/
 *
 * Check out this link for haversine great-circle distances instead:
 * https://www.btskinner.io/rworkshop/modules/hp_rcpp.html
 */

// MatrixXd DistGeo_cpp(const MapMatd& loc){
//   // Initialize locaiton objects
//   int n = loc.rows();
//   MatrixXd D = MatrixXd::Identity(n, n); // start with n x n Idnet. matrix
//
//   // Initialize 'geodesic' object
//   double a = 6378137, f = 1/298.257223563;
//   double lat_i, lon_i, lat_j, lon_j;
//   double azi_i, azi_j, sij;
//   struct geod_geodesic g;
//   geod_init(&g, a, f);
//
//   int i, j;
//   for(i = 1; i < n; ++i){
//     for(j = 0; j < i; ++j){
//       lat_i = loc(i, 0);
//       lon_i = loc(i, 1);
//       lat_j = loc(j, 0);
//       lon_j = loc(j, 1);
//       geod_inverse(&g, lat_i, lon_i, lat_j, lon_j, &sij, &azi_i, &azi_j);
//       D(i, j) = sij;
//       D(j, i) = sij;
//     }
//   }
//   return D;
// }

// MatrixXd DistGeo_cpp(const MapMatd& loc){
//   int n;
//   n = loc.rows();
//   MatrixXd D = MatrixXd::Identity(n, n); // start with n x n Idnet. matrix
//   int i, j;
//   i = 1, j = 2;
//   // Geodesic geod()
//   // Math::real lati, latj, loni, lonj;
//   double lati, latj, loni, lonj;
//   lati = loc(i, 1);
//   loni = loc(i, 2);
//   latj = loc(j, 1);
//   lonj = loc(j, 2);
//   Geodesic geod(Constants::WGS84_a(), Constants::WGS84_f());
//   double s12;
//   geod.Inverse(lati, loni, latj, lonj, s12);
//   Math::real distance = GeographicLib::Geodesic::Inverse(lati, loni, latj, lonj, s12);
//   D(i, j) = s12;
//   return D;
// }

//==============================================================================
/* AtA(A) ----
 * inline code that computes t(A) %*% A
 */

// [[Rcpp::export]]
inline Eigen::MatrixXd AtA(const MatrixXd& A) {
  int n(A.cols());
  return MatrixXd(n,n).setZero().selfadjointView<Lower>()
                      .rankUpdate(A.adjoint());
}

// // [[Rcpp::export]]
// inline Eigen::MatrixXd AAt(const MatrixXd& A) {
//   int m(A.rows());
//   return MatrixXd(m,m).setZero().selfadjointView<Lower>()
//                       .rankUpdate(A);
// }

//==============================================================================
/* solve_cpp(A, B) ----
 * solves Ax = B for x.
 */

//' solve Ax = B
//'
//' @param A numeric matrix
//' @param B numeric matrix
//'
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd solve_cpp(const MatrixXd& A, const MatrixXd& B){
  return A.colPivHouseholderQr().solve(B);
}

//' solve Ax = I
//'
//' @param A numeric matrix
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd solve_ident_cpp(const MatrixXd& A){
  MatrixXd I = MatrixXd::Identity(A.rows(),A.cols());
  return A.colPivHouseholderQr().solve(I);
}

/***R
# ## Test the inline function AtA()
# M <- matrix(as.double(1:6), ncol = 2)
# crossprod(M) # native R
# AtA(M) # cpp way
*/

//==============================================================================
/* tinvchol_cpp(V, nugget = 0)
 * computes the choleskey decomposition upper triangle (U) of V
 * and then returns the transpose: t(U)
 *
 * if a nugget is included... (nothing yet)
 */

//' Find the transposed inverse cholesky decomposition of V
//'
//' @param V numeric matrix
//' @param nugget numeric nugget to add to variance matrix
//'
//' @export
//'
//' @examples #TBA
// [[Rcpp::export]]
Eigen::MatrixXd tinvchol_cpp(const MapMatd& V, double nugget = 0.){

  double n = V.rows(); //dim of V

  if (abs(nugget) > 0) {
    // with nugget
    const MatrixXd M = (1 - nugget) * V.array() + nugget * MatrixXd::Identity(n,n).array();
    const MatrixXd Mn = M.matrix();
    const LLT<MatrixXd> llt(Mn);
    const MatrixXd tUinv = llt.matrixL().solve(MatrixXd::Identity(n, n));
    return tUinv;
  } else {
    // no nugget case
    const LLT<MatrixXd> llt(V); //chol decomp of V: (U)
    const MatrixXd tUinv = llt.matrixL().solve(MatrixXd::Identity(n, n));
    return tUinv;
  }
}

/***R
# ## Test tinvchol_cpp()
# V <- crossprod(M)
#
# # no nugget (i.e. nugget = 0)
# t(backsolve(chol(V), diag(nrow(V)))) # native R
# tinvchol_cpp(V) # cpp
#
# # with nugget
# nug = .1
# Vn = (1 - nug) * V + nug * diag(nrow(V))
# t(backsolve(chol(Vn), diag(nrow(Vn)))) # R
# tinvchol_cpp(Vn, nug) # cpp
*/

//==============================================================================
/* fitGLS_cpp(X, V, y, X0, nugget = 0, threads = 1) ----
 *
 */

//' Fit GLS to remote sensing data
//'
//' @details see `fitGLS()`
//'
//' @param X numeric matrix
//' @param V numeric matrix
//' @param y numeric vector
//' @param X0 numeric matrix
//' @param nugget numeric nugget to add to V
//' @param save_xx logical: should xx, xx0, and tInvCholV be returned? This
//' functionality is meant for use with the partitioned GLS whereby these
//' values are used to calculate cross-partition statistics.
//' @param threads integer indicating the number of threads to use. This current
//' version does not have multi-thread functionality so this argument does
//' nothing yet.
//'
//' @export
//' @examples #TBA
// [[Rcpp::export]]
List fitGLS_cpp(const MapMatd& X,
                const MapMatd& V,
                const MapMatd& y,
                const MapMatd& X0,
                double nugget = 0.,
                bool save_xx = false,
                const int threads = 1){

  int nX = X.rows(), pX = X.cols(); // dimensions of X
  MatrixXd tUinv = tinvchol_cpp(V, nugget); // transpose of chol(V) = t(inv(U))
  MatrixXd xx = tUinv * X; // t(inv(U)) %*% X
  MatrixXd yy = tUinv * y; // t(inv(U)) %*% y
  MatrixXd varX = AtA(xx); // crossprod(xx)
  MatrixXd XtY(xx.adjoint() * yy); // crossprod(xx, yy)

  // solve for betahat (using one specific solver - there are others)
  VectorXd betahat(varX.colPivHouseholderQr().solve(XtY));

  // calculate some statistics
  int dft = nX - xx.cols();
  VectorXd SSE = AtA(yy - xx * betahat); // SSE
  VectorXd MSE = SSE/dft; // MSE
  MatrixXd varXinv = varX.colPivHouseholderQr().solve(
    MatrixXd::Identity(varX.rows(), varX.cols()));
  MatrixXd varcov = varXinv.array() * MSE.array()[0];
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

  const MatrixXd varX0inv = varX0.colPivHouseholderQr().solve(
    MatrixXd::Identity(varX0.rows(), varX0.cols()));
  const MatrixXd varcov0 = MSE0 * varX0inv.array();

  VectorXd se0 = varcov0.matrix().diagonal();
  for (int i = 0; i < se0.size(); i++){
    se0[i] = std::sqrt(se0[i]);
  }

  double SSR = SSE0.array()[1] - SSE.array()[1];
  /* F test ----
   *
   */
  const double FF = (nX - xx.cols())/(xx.cols() - xx0.cols()) *
    (SSE0.array()[0] - SSE.array()[0]) / SSE.array()[0];
  VectorXi dfF(2);
  dfF(0) = xx.cols() - xx0.cols();
  dfF(1) = nX - xx.cols();

  // return a list of all needed values
  List res_list = List::create(Named("betahat") = betahat,
                               Named("SSE") = SSE,
                               Named("MSE") = MSE,
                               Named("SE") = se,
                               Named("dft") = dft,
                               Named("tstat") = tstat,
                               Named("pval.t") = NA_REAL,
                               Named("logLik") = logLik,
                               Named("betahat0") = betahat0,
                               Named("SSE0") = SSE0,
                               Named("MSE0") = MSE0,
                               Named("SE0") = se0,
                               Named("MSR") = MSR,
                               Named("df0") = df0,
                               Named("logLik0") = logLik0,
                               Named("df.F") = dfF,
                               Named("Fstat") = FF,
                               Named("pval.F") = NA_REAL);

  // handle the conditional return of large matrices
  if (save_xx){
    res_list.push_back(xx, "xx");
    res_list.push_back(xx0, "xx0");
    res_list.push_back(tUinv, "tInvCholV");
  } else{
    res_list.push_back(NA_REAL, "xx");
    res_list.push_back(NA_REAL, "xx0");
    res_list.push_back(NA_REAL, "tInvCholV");
  }

  return res_list;
}

/*** R
# # ## test fitGLS_cpp()
# load("..//R/vignettes-and-examples/test-gls.rda", verbose = FALSE)
# Xnull <- matrix(as.double(1), ncol = 1, nrow = nrow(X.small))
# source("../../R/fitGLS.R")
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

//' Caculate log-liklihood of GLS model
//'
//' @details this function is mostly meant to optimize the nugget for a paritular
//' set of data.
//'
//' Note: this function should be deprecated and simply added as functionality
//' to `fitGLS_cpp()`.
//'
//' @param nugget the nugget to add to V
//' @param X numeric matrix
//' @param V numeric matrix
//' @param y numeric vector
//'
//' @export
//' @examples #TBA
// [[Rcpp::export]]
inline double LogLikGLS_cpp(double nugget,
                        const MapMatd& X,
                        const MapMatd& V,
                        const MapMatd& y){

  const int nX = X.rows(); // dimensions of X
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
  double logDetV = 0;
  for (int i = 0; i < tUinv.rows(); i++){
    logDetV += log(tUinv.diagonal()[i]);
  }
  logDetV *= -2;

  double logLik = -0.5 * (nX * log(2 * M_PI) + nX *
                          log(dft * MSE.array()[0]/nX) + logDetV + nX);

  return logLik;
}

/* Optimizer function
 * Currently this function uses the same algorithm that optimize() derived from.
 * This code is a translation of the fortran fmin algorithm:
 * http://www.netlib.org/fmm/fmin.f
 */

//' Find the maximum likelihood estimate of the nugget
//'
//' @details this is the C++ version of `optimize()` which is specific to
//' finding the nugget value that maximizes the log-liklihood of `fitGLS_cpp()`
//'
//' Note: this function actually uses `LogLikGLS_cpp()` which should be swapped
//' for `fitGLS_cpp()` once the correct funcionality is added.
//'
//' @param X numeric matrix
//' @param V numeric matrix
//' @param y numeric vector
//' @param lower lower boundary for nugget search
//' @param upper upper boundary for nugget search
//' @param tol desired accuracy of nugget search
//' @param debug logical: debug mode?
//'
//' @export
//' @examples #TBA
// [[Rcpp::export]]
double optimizeNugget_cpp(const MapMatd& X, const MapMatd& V, const MapMatd& y,
                          double lower = 0, double upper = 1, double tol = .00001,
                          bool debug = false){

  // varible declaration
  double ax = lower;
  double bx = upper;
  double a,b,c,d,e,eps,xm,p,q,r,tol1,tol2,u,v,w;
  double fu,fv,fw,fx,x, f_min;

  // c is the squared inverse of the golden ratio
  c = 0.5*(3. - sqrt(5.));

  // eps is approximately sqrt of machine precision
  // eps = std::numeric_limits<double>::epsilon();
  // eps = sqrt(eps);
  // cout << "eps_std = "<<eps<<endl;

  eps = 1.;
  lab0:
    eps = eps/2.;
    tol1 = 1. + eps;
    if (tol1 > 1.) {goto lab0;}
    eps = sqrt(eps);
    // cout << "eps =" << eps<< endl;

  // initialization
  a = ax;
  b = bx;
  v = a + c*(b - a);
  w = v;
  x = v;
  e = 0.;
  fx = -LogLikGLS_cpp(x, X, V, y); //invert function to minimize
  fv = fx;
  fw = fx;

  int i = 0;
  // main loop start
  lab1:
    // cout << "Loop Start (lab1): iteration " <<i<<endl;
    if (debug) {
    i += 1;
      if (i >= 10) {
        cout << "breaking loop, too many iterations"<<endl;
        goto lab8;
      }
    }
    xm = 0.5*(a + b);
    tol1 = eps*abs(x) + tol/3.;
    tol2 = 2.*tol1;
    // check stoping criteria
    if (abs(x - xm) <= (tol2 - 0.5*(b - a))) {goto lab8;}
    // cout << "stop crit. not met: "<<abs(x-xm)<<" > "<<tol2-.5*(b-a)<<endl;
    // is golden section necessary?
    if (abs(e) <= tol1) {goto lab3;}
    // fit parabola
    r = (x - w)*(fx - fv);
    q = (x - v)*(fx - fw);
    p = (x - v)*q - (x - w)*r;
    q = 2.*(q - r);
    if (q > 0.) {p = -p;}
    q =  abs(q);
    r = e;
    e = d;
  lab2:
    // cout << "check parabola (lab2)" << endl;
    // is parabola acceptable?
    if (abs(p) >= abs(0.5*q*r)) {goto lab3;}
    if (p <= q*(a - x)) goto lab3;
    if (p >= q*(b - x)) goto lab3;
    // parabolic interpolation step
    d = p/q;
    u = x + d;
    // f must not be evaluated too close to ax or bx
    if ((u - a) < tol2) {d = copysign(tol1, xm - x);}
    if ((b - u) < tol2) {d = copysign(tol1, xm - x);}
    goto lab4;
  lab3:
    // cout << "golden section step (lab3)" <<endl;
    // golden section step
    if (x >= xm) {e = a - x;}
    if (x < xm) {e = b - x;}
    d = c*e;
  lab4:
    // cout << "check tolerance and update vars (lab4)" << endl;
    //f must not be evaluated too close to x
    if (abs(d) >= tol1) {u = x + d;}
    if (abs(d) < tol1) {u = x + copysign(tol1, d);}
    fu = -LogLikGLS_cpp(u, X, V, y);
    //update  a, b, v, w, and x
    if (fu > fx) {goto lab5;}
    if (u >= x) {a = x;}
    if (u < x) {b = x;}
    v = w;
    fv = fw;
    w = x;
    fw = fx;
    x = u;
    fx = fu;
    goto lab1;
  lab5:
    // cout << "conditional variable reset (lab5)" << endl;
    if (u < x) {a = u;}
    if (u >= x) {b = u;}
    if (fu <= fw) {goto lab6;}
    if (w == x) {goto lab6;}
    if (fu <= fv) {goto lab7;}
    if (v == x) {goto lab7;}
    if (v == w) {goto lab7;}
    goto lab1;
  lab6:
    // cout << "update function results (lab6)" << endl;
    v = w;
    fv = fw;
    w = u;
    fw = fu;
    goto lab1;
  lab7:
    // cout << "update function results alterante (lab7)" << endl;
    v = u;
    fv = fu;
    goto lab1;
  // end of main loop
  lab8:
    // cout << "return statement (lab8)" << endl;
    f_min = x;
    if (ax + tol >= f_min){
      if (fx <= -LogLikGLS_cpp(f_min, X, V, y)){
        return ax;
      }
    }
    return f_min;
}

/***R
# tol <- .00001
# system.time(tmp1 <- fitNugget(X.small, V.small, y.small, c(0, 1), tol))
# system.time(tmp2 <- optimizeNugget_cpp(X.small, V.small, y.small, 0, 1, tol))
# system.time(tmp3 <- fitNugget_Rcpp(X.small, V.small, y.small, c(0,1), tol))
#
# (vals <- c(tmp1, tmp2, tmp3))
# # xs <- seq(0, 1, length.out = 20)
# # fxs <- sapply(xs, function(x){LogLikGLS_cpp(x, X.small, V.small, y.small)})
# # plot(fxs ~ xs);abline(v = vals, col = c("red", "green", "blue"), lty = 1:3)

*/

/* GLS_worker_cpp()
 * This function will perform the GLS operations on each subset of X
 * For now, V_i needs to be provided with each X_i and y_i
 * because Vfit() is currently not implemented in C++. This is largely
 * due to the fact that since distGeo() needs to be done externally as well.
 */

//' Worker function 1 for paritioned GLS
//'
//' @details this function is the first of 2 (maybe 3) worker functions that,
//' together, perform the partitioned GLS analysis.
//'
//' This function is simply a wrapper for fitGLS_cpp() that finds the MLE nugget
//' and adds it to the output.
//'
//' NOTE: eventually, the worker functions will perform the analysis using
//' multiple cores but that has not yet been implemented.
//'
//' @param y numeric vector
//' @param X numeric matrix
//' @param V numeric matrix
//' @param X0 numeric matrix
//' @param save_xx logical: should xx, xx0, and tInvCholV be returned?
//'
//' @export
//' @examples #TBA
// [[Rcpp::export]]
List GLS_worker_cpp(const MapMatd& y,
                    const MapMatd& X,
                    const MapMatd& V,
                    const MapMatd& X0,
                    bool save_xx = false){

  // Estimate the nugget
  double nug = optimizeNugget_cpp(X, V, y);
  // run GLS
  List x_gls = fitGLS_cpp(X, V, y, X0, nug, save_xx, 1);
  // add the nugget to the output
  x_gls.push_back(nug, "nugget");

  return x_gls;
}

//' Worker function 2 for partitioned GLS
//'
//' @details this is the second worker function for the partitioned GLS analysis.
//'
//' NOTE: currently, there is no parallel functionality and the partitioned
//' form of the GLS is not implemented entirely in C++. Instead, the R function
//' fitGLS.partition_rcpp() weaves between R and C++ on a single core. While
//' this method is still much faster than the purely R implementation, migration
//' to entirely C++ will greatly improve speed further. This migration requires
//' calculating geographic distances with C++ which I've not yet written.
//'
//' Additionally, there seems to be a memory-related issue with this code. I've
//' successfully used this function when partitions have 100 or fewer rows (too
//' small). However, larger partitions cause a fatal error that causes a crash.
//'
//' @param xxi numeric matrix xx from  partition i
//' @param xxj numeric matrix xx from  partition j
//' @param xxi0 numeric matrix xx0 from  partition i
//' @param xxj0 numeric matrix xx0 from  partition j
//' @param tUinv_i numeric matrix tInvCholV from  partition i
//' @param tUinv_j numeric matrix tInvCholV from  partition j
//' @param Vsub numeric variance matrix for Xij (upper block)
//' @param df1 first degree of freedom
//' @param df2 second degree of freedom
//'
//' @export
//' @examples #TBA
// [[Rcpp::export]]
List crosspart_worker_cpp(const MapMatd& xxi,
                          const MapMatd& xxj,
                          const MapMatd& xxi0,
                          const MapMatd& xxj0,
                          const MapMatd& tUinv_i,
                          const MapMatd& tUinv_j,
                          const MapMatd& Vsub,
                          int df1,
                          int df2){
  // int N = Vij.cols(); // total rows
  // int np = N/2; // rows per partition
  int np = xxi.rows();

  // // scale the nuggets if nonzero
  // double nugget_i = nug_i == 0 ? 0 : (1 - nug_i) / nug_i;
  // double nugget_j = nug_j == 0 ? 0 : (1 - nug_j) / nug_j;
  // // combine the nuggets into a vector
  //    // equivalent to rep(c(nugget_i, nugget_j), each = np)
  // VectorXd nugget_vector(2*np);
  // for(int i = 0; i < np; i++){
  //   nugget_vector(i) = nug_i;
  // }
  // for(int j = np + 1; j < 2*np; j++){
  //   nugget_vector(j) = nug_j;
  // }
  // // // turn this into a diagonal matrix
  // // MatrixXd VDiag = nugget_vector.asDiagonal();
  // // then add Vij
  // VDiag = VDiag + Vij;

  // extract block matrix
  // MatrixXd Vsub = VDiag.block(1, np + 1, np, np);

  // Calculate some Statistics # this math is wrong.
  MatrixXd Rij = tUinv_i.adjoint() * Vsub * tUinv_j.adjoint();

  MatrixXd Hi = xxi * solve_ident_cpp(xxi.adjoint() * xxi) * xxi.adjoint();
  MatrixXd Hj = xxj * solve_ident_cpp(xxj.adjoint() * xxj) * xxj.adjoint();

  MatrixXd Hi0 = xxi0 * solve_ident_cpp(xxi0.adjoint() * xxi0) * xxi0.adjoint();
  MatrixXd Hj0 = xxj0 * solve_ident_cpp(xxj0.adjoint() * xxj0) * xxj0.adjoint();

  MatrixXd SiR = Hi - Hi0;
  MatrixXd SjR = Hj - Hj0;
  // Rcout << "SiR" <<SiR <<endl;

  MatrixXd npDiag = MatrixXd::Identity(np, np);
  MatrixXd SiE = npDiag - Hi;
  MatrixXd SjE = npDiag - Hj;

  // rSSR and rSSE were calculated incorrectly!
  VectorXd tmpRvec = (Rij * SjR * Rij.adjoint()).array();
  VectorXd SRvec = SiR.array();
  MatrixXd tmp_ijR = SRvec.adjoint() * tmpRvec;
  MatrixXd rSSRij = tmp_ijR.array()/df1;

  VectorXd tmpEvec = (Rij * SjE * Rij.adjoint()).array();
  VectorXd SEvec = SiE.array();
  MatrixXd tmp_ijE = SEvec.adjoint() * tmpEvec;
  MatrixXd rSSEij = tmp_ijE.array()/df2;

  // output
  List out_lst = List::create(Named("rSSRij") = rSSRij.matrix(),
                              Named("rSSEij") = rSSEij.matrix());

  return out_lst;
}
