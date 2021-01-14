#include "function-declarations.h"

//' Fit GLS to remote sensing data
//'
//' @details see `fitGLS()`
//'
//' @param X numeric matrix
//' @param V numeric matrix
//' @param y numeric vector
//' @param X0 numeric matrix
//' @param nugget numeric nugget to add to V
//' @param save_xx logical: should xx, xx0, and invcholV be returned? This
//' functionality is meant for use with the partitioned GLS whereby these
//' values are used to calculate cross-partition statistics.
//' @param threads integer indicating the number of threads to use. This current
//' version does not have multi-thread functionality so this argument does
//' nothing yet.
//'
//' @examples #TBA
// [[Rcpp::export(.fitGLS2_cpp)]]
void fitGLS2_cpp(List L,
                const MapMatd& X,
                const MapMatd& V,
                const MapMatd& y,
                const MapMatd& X0,
                double nugget,
                bool save_xx,
                const int threads){

  int nX = X.rows(), pX = X.cols(); // dimensions of X
  MatrixXd tUinv = invchol_cpp(V, nugget); // transpose of chol(V) = t(inv(U))
  MatrixXd xx = tUinv * X; // t(inv(U)) %*% X
  MatrixXd yy = tUinv * y; // t(inv(U)) %*% y
  MatrixXd varX = AtA(xx); // crossprod(xx)
  MatrixXd XtY(xx.adjoint() * yy); // crossprod(xx, yy)

  // solve for betahat (using one specific solver - there are others)
  VectorXd betahat(varX.colPivHouseholderQr().solve(XtY));

  // calculate some statistics
  int dft = nX - xx.cols();
  VectorXd SSE = AtA(yy - xx * betahat); // SSE
  // cout << "\nSSE:\n" << SSE << endl;
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

  MatrixXd xx0 = tUinv * X0;
  MatrixXd varX0 = AtA(xx0); // crossprod(xx0)
  MatrixXd X0tY(xx0.adjoint() * yy); // crossprod(xx0, yy)
  VectorXd betahat0(solve_cpp(varX0, X0tY));
  int df0 = betahat0.size();
  VectorXd SSE0 = AtA(yy - xx0 * betahat0); // SSE
  // cout << "\nSSE0:\n" << SSE0 << endl;
  double MSE0 = SSE0.array()[0]/(nX - xx.cols()); // MSE
  double MSR = (SSE0.array()[0] - SSE.array()[0])/(xx.cols() - xx0.cols());
  double logLik0 = -0.5 * (nX * log(2 * M_PI) + nX * log((nX - df0) * MSE0/nX) +
                           logDetV + nX);

  MatrixXd varX0inv = varX0.colPivHouseholderQr().solve(
    MatrixXd::Identity(varX0.rows(), varX0.cols()));
  MatrixXd varcov0 = MSE0 * varX0inv.array();

  VectorXd se0 = varcov0.matrix().diagonal();
  for (int i = 0; i < se0.size(); i++){
    se0[i] = std::sqrt(se0[i]);
  }

  VectorXd SSR = SSE0.array() - SSE.array();
  // cout << "\nSSR:\n" << SSR << endl;

  /* F test ----
   *
   */

  VectorXd FF = ((double)nX - (double)xx.cols())/ // force integers to doubles
    ((double)xx.cols() - (double)xx0.cols()) *
      SSR.array() / SSE.array();
  VectorXi dfF(2);
  dfF(0) = xx.cols() - xx0.cols();
  dfF(1) = nX - xx.cols();

  // replace elements in the list
    L["betahat"] = betahat;
    L["SSE"] = SSE;
    L["MSE"] = MSE;
    L["SE"] = se;
    L["dft"] = dft;
    L["logDetV"] = logDetV;
    L["tstat"] = tstat;
    L["logLik"] = logLik;
    L["betahat0"] = betahat0;
    L["SSE0"] = SSE0;
    L["MSE0"] = MSE0;
    L["SE0"] = se0;
    L["MSR"] = MSR;
    L["df0"] = df0;
    L["logLik0"] = logLik0;
    L["df.F"] = dfF;
    L["Fstat"] = FF;

  // handle the conditional return of large matrices
  if (save_xx){
    L["xx"] = xx;;
    L["xx0"] = xx0;
    L["invcholV"] = tUinv;
  }
  // return res_list;
}
