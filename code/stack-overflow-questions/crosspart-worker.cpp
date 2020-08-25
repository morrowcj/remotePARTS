// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

using namespace Rcpp;

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;

typedef Map<MatrixXd> MapMatd;

// Helper function ----
//' solve Ax = B
//'
//' @param A numeric matrix
//' @param B numeric matrix
//'
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd solve_ident_cpp(const MatrixXd& A){
  MatrixXd I = MatrixXd::Identity(A.rows(),A.cols());
  return A.colPivHouseholderQr().solve(I);
}

// Main function ----
//' Worker function
//'
//' @param xxi n x p numeric matrix xx from  partition i
//' @param xxj n x p numeric matrix xx from  partition j
//' @param xxi0 n x p0 numeric matrix xx0 from  partition i
//' @param xxj0 n x p0 numeric matrix xx0 from  partition j
//' @param tUinv_i n x n numeric matrix tInvCholV from  partition i
//' @param tUinv_j n x n numeric matrix tInvCholV from  partition j
//' @param Vij n x n numeric variance matrix for Xij
//' @param df1 integer first degree of freedom
//' @param df2 integer second degree of freedom
//'
//' @export
//' @examples #TBA
// [[Rcpp::export]]
List crosspart_worker_cpp(const MapMatd& xxi, // n x p
                          const MapMatd& xxj, // n x p
                          const MapMatd& xxi0,
                          const MapMatd& xxj0,
                          const MapMatd& tUinv_i,
                          const MapMatd& tUinv_j,
                          double nug_i,
                          double nug_j,
                          const MapMatd& Vij,
                          int df1,
                          int df2){
  int N = Vij. cols(); // total rows
  int np = N/2; // rows per partition

  // scale the nuggets if nonzero
  double nugget_i = nug_i == 0 ? 0 : (1 - nug_i) / nug_i;
  double nugget_j = nug_j == 0 ? 0 : (1 - nug_j) / nug_j;
  // combine the nuggets into a vector
  // equivalent to rep(c(nugget_i, nugget_j), each = np)
  VectorXd nugget_vector(2*np);
  for(int i = 0; i < np; i++){
    nugget_vector(i) = nug_i;
  }
  for(int j = np + 1; j < 2*np; j++){
    nugget_vector(j) = nug_j;
  }
  // turn this into a diagonal matrix
  MatrixXd VDiag = nugget_vector.asDiagonal();
  // then add Vij
  VDiag = VDiag + Vij;
  // extract sub-block of varcovar matrix (only unique pairs)
  MatrixXd Vsub = VDiag.block(1, np + 1, np, np);

  // Calculate some Statistics
  MatrixXd Rij = tUinv_i.adjoint() * Vsub * tUinv_j.adjoint();

  MatrixXd Hi = xxi * solve_ident_cpp(xxi.adjoint() * xxi) * xxi.adjoint();
  MatrixXd Hj = xxj * solve_ident_cpp(xxj.adjoint() * xxj) * xxj.adjoint();

  MatrixXd Hi0 = xxi0 * solve_ident_cpp(xxi0.adjoint() * xxi0) * xxi0.adjoint();
  MatrixXd Hj0 = xxj0 * solve_ident_cpp(xxj0.adjoint() * xxj0) * xxj0.adjoint();

  MatrixXd SiR = Hi - Hi0;
  MatrixXd SjR = Hj - Hj0;

  MatrixXd npDiag = MatrixXd::Identity(np, np);
  MatrixXd SiE = npDiag - Hi;
  MatrixXd SjE = npDiag - Hj;

  MatrixXd tmp_ijR = SiR * (Rij * SjR * Rij.adjoint());
  MatrixXd rSSRij = tmp_ijR.array()/df1;

  MatrixXd tmp_ijE = SiE * (Rij * SjE * Rij.adjoint());
  MatrixXd rSSEij = tmp_ijE.array()/df2;

  // output
  List out_lst = List::create(Named("rSSRij") = rSSRij.matrix(),
                              Named("rSSEij") = rSSEij.matrix());

  return out_lst;
}
