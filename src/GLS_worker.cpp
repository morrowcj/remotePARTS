#include "function-declarations.h"

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
//' @param nug_l lower boundary for nugget optimization
//' @param nug_u upper boundary for nugget optimization
//' @param nug_tol tolerance of nugget optimization
//' @param save_xx logical: should xx, xx0, and invcholV be returned?
//'
//' @examples #TBA
// [[Rcpp::export(.GLS_worker_cpp)]]
List GLS_worker_cpp(const MapMatd& y,
                    const MapMatd& X,
                    const MapMatd& V,
                    const MapMatd& X0,
                    double nug_l,
                    double nug_u,
                    double nug_tol,
                    bool save_xx = false){

  // Estimate the nugget
  double nug = optimize_nugget_cpp(X, V, y, nug_l, nug_u, nug_tol, false);
  // run GLS
  List x_gls = fitGLS_cpp(X, V, y, X0, nug, save_xx, 1);
  // add the nugget to the output
  x_gls.push_back(nug, "nugget");

  return x_gls;
}
