#include "function-declarations.h"

//' Fit GLS to remote sensing data
//'
//' @details see `fitGLS()`
//'
//' @param L empty remoteGLS list to be filled
//' @param X numeric matrix
//' @param V numeric matrix
//' @param y numeric vector
//' @param X0 numeric matrix
//' @param nugget numeric nugget to add to V
//' @param save_xx logical: should xx, xx0, and invcholV be returned? This
//' functionality is meant for use with the partitioned GLS whereby these
//' values are used to calculate cross-partition statistics.
//' @param LL_only logical: should only the log-liklihood be computed?
//' @param no_F logical: should calculations needed for F tests be skipped?
//' @param threads number of threads used by Eigen for matrix algebra
//'
//' @examples
// [[Rcpp::export(.fitGLS2_cpp)]]
void fitGLS2_cpp(List L,
                const MapMatd& X,
                const MapMatd& V,
                const MapMatd& y,
                const MapMatd& X0,
                double nugget,
                bool save_xx,
                bool LL_only,
                bool no_F,
                const int threads){

  int nX = X.rows(), pX = X.cols(); // dimensions of X
  MatrixXd tUinv = invchol_cpp(V, nugget); // transpose of chol(V) = t(inv(U))
  MatrixXd xx = tUinv * X; // t(inv(U)) %*% X
  MatrixXd yy = tUinv * y; // t(inv(U)) %*% y
  MatrixXd varX = AtA(xx); // crossprod(xx)
  MatrixXd XtY(xx.adjoint() * yy); // crossprod(xx, yy)

  // solve for betahat (using one specific solver - there are others)
  VectorXd betahat(varX.colPivHouseholderQr().solve(XtY));

  // calculate df
  int dft = nX - xx.cols();
  // calculate SS statistics
  VectorXd SSE = AtA(yy - xx * betahat); // SSE
  VectorXd MSE = SSE/dft; // MSE

  // Calculate Log-Liklihood ----
  double logDetV = 0;
  for (int i = 0; i < tUinv.rows(); i++){
    logDetV += log(tUinv.diagonal()[i]);
  }
  logDetV *= -2;
  double logLik = -0.5 * (nX * log(2 * M_PI) + nX *
                          log(dft * MSE.array()[0]/nX) + logDetV + nX);
  // add logLik to output list
  L["logLik"] = logLik;

  // exit the function if LL_only
  if (LL_only) {
    return;
  }

  // calculate varcov stats
  MatrixXd varXinv = varX.colPivHouseholderQr().solve(
    MatrixXd::Identity(varX.rows(), varX.cols()));
  MatrixXd varcov = varXinv.array() * MSE.array()[0];

  // caclulate standard errors
  VectorXd se = varcov.matrix().diagonal();
  for (int i = 0; i < se.size(); i++){
    se[i] = std::sqrt(se[i]);
  }
  // calculate t-statistic
  VectorXd tstat(betahat.size());
  for (int i = 0; i < betahat.size(); i++){
    tstat[i] = betahat[i] / se[i];
  }

  // Add all additional stats, up to t-test, into output list
  L["betahat"] = betahat;
  L["SSE"] = SSE;
  L["MSE"] = MSE;
  L["SE"] = se;
  L["dft"] = dft;
  L["logDetV"] = logDetV;
  L["tstat"] = tstat;
  // add xx and invcholV to output, if requested
  if (save_xx){
    L["xx"] = xx;
    L["invcholV"] = tUinv;
  }
  // break function before calculating F statistics, if asked
  if (no_F) {
    return;
  }

// **START F-test** ----
  // null model stats
  MatrixXd xx0 = tUinv * X0;
  MatrixXd varX0 = AtA(xx0); // crossprod(xx0)
  MatrixXd X0tY(xx0.adjoint() * yy); // crossprod(xx0, yy)
  VectorXd betahat0(solve_cpp(varX0, X0tY));
  // df
  int df0 = betahat0.size();
  // SS
  VectorXd SSE0 = AtA(yy - xx0 * betahat0); // SSE
  double MSE0 = SSE0.array()[0]/(nX - xx.cols()); // MSE
  double MSR = (SSE0.array()[0] - SSE.array()[0])/(xx.cols() - xx0.cols());
  // LL
  double logLik0 = -0.5 * (nX * log(2 * M_PI) + nX * log((nX - df0) * MSE0/nX) +
                           logDetV + nX);
  // covariance matrix
  MatrixXd varX0inv = varX0.colPivHouseholderQr().solve(
    MatrixXd::Identity(varX0.rows(), varX0.cols()));
  MatrixXd varcov0 = MSE0 * varX0inv.array();
  // SE
  VectorXd se0 = varcov0.matrix().diagonal();
  for (int i = 0; i < se0.size(); i++){
    se0[i] = std::sqrt(se0[i]);
  }

  VectorXd SSR = SSE0.array() - SSE.array();
  // cout << "\nSSR:\n" << SSR << endl;

  // Calculate F-statistic
  VectorXd FF = ((double)nX - (double)xx.cols())/ // force integers to doubles
    ((double)xx.cols() - (double)xx0.cols()) *
      SSR.array() / SSE.array();
  // F dfs
  VectorXi dfF(2);
  dfF(0) = xx.cols() - xx0.cols();
  dfF(1) = nX - xx.cols();

  // replace elements in the list

    L["betahat0"] = betahat0;
    L["SSE0"] = SSE0;
    L["MSE0"] = MSE0;
    L["SE0"] = se0;
    L["MSR"] = MSR;
    L["df0"] = df0;
    L["logLik0"] = logLik0;
    L["df.F"] = dfF;
    L["Fstat"] = FF;

  // add in xx0, if asked
  if (save_xx){
    L["xx0"] = xx0;
  }
  return;
}
