#include "function-declarations.hpp"

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

  const int nX = X.rows(); // dimensions of X
  const MatrixXd tUinv = invchol_cpp(V, nugget); // transpose of chol(V) = t(inv(U))
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
