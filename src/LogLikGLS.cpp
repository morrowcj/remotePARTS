#include "function-declarations.h"

//' Caculate log-liklihood of GLS model
//'
//' @details this function is mostly meant to optimize the nugget for a paritular
//' set of data.
//'
//' Note: this function should be deprecated and simply added as functionality
//' to `.fitGLS_cpp()`.
//'
//' @param nugget the nugget to add to V
//' @param X numeric matrix
//' @param V numeric matrix
//' @param y numeric vector
//'
//' @examples #TBA
// [[Rcpp::export(.LogLikGLS_cpp)]]
double LogLikGLS_cpp(double nugget,
                     const MapMatd& X,
                     const MapMatd& V,
                     const MapMatd& y){

  int nX = X.rows(); // dimensions of X
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
  VectorXd MSE = SSE/dft; // MSE
  double logDetV = 0;
  for (int i = 0; i < tUinv.rows(); i++){
    logDetV += log(tUinv.diagonal()[i]);
  }
  logDetV *= -2;

  double logLik = -0.5 * (nX * log(2 * M_PI) + nX *
                          log(dft * MSE.array()[0]/nX) + logDetV + nX);

  return logLik;
}
