#include "function-declarations.h"

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
//' @param invCholV_i numeric matrix invcholV from  partition i
//' @param invCholV_j numeric matrix invcholV from  partition j
//' @param Vsub numeric variance matrix for Xij (upper block)
//' @param nug_i nugget from partition i
//' @param nug_j nugget from partition j
//' @param df1 first degree of freedom
//' @param df2 second degree of freedom
//' @param Vcoef logical indicating if the coefficient covariance matrix
//' should be returned
//' @param ncores integer indicating nubmer of cores to use
//'
// [[Rcpp::export(.crosspart_worker_cpp)]]
List crosspart_worker_cpp(const MapMatd& xxi,
                          const MapMatd& xxj,
                          const MapMatd& xxi0,
                          const MapMatd& xxj0,
                          const MapMatd& invCholV_i,
                          const MapMatd& invCholV_j,
                          const MapMatd& Vsub,
                          double nug_i,
                          double nug_j,
                          int df1,
                          int df2,
                          bool Vcoef,
                          int ncores){
  Eigen::setNbThreads(ncores);

  // int N = Vsub.cols(); // total rows
  // int np = N/2; // rows per partition
  int np = xxi.rows();
  MatrixXd Vij;

  // // scale the nuggets if nonzero
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
  // // turn this into a diagonal matrix
  // // MatrixXd Vij = nugget_vector.asDiagonal();
  // // then add Vsub
    // // Note (18-Nov-2020): this needs to be fixed. The appropriate
    // // method is actually sqrt((1 - nug_i)*(1 - nug_j)) * Vij.array
  // Vij = Vij + Vsub;
  Vij = sqrt((1 - nug_i)*(1 - nug_j)) * Vsub.array();

  // Caclulate some satistics
  MatrixXd B = Vij * invCholV_j.adjoint(); //tcrossprod(Vij, invCholV_j)
  MatrixXd Rij = invCholV_i * B;

  MatrixXd Wi = solve_ident_cpp(xxi.adjoint() * xxi);
  MatrixXd Wj = solve_ident_cpp(xxj.adjoint() * xxj);

  MatrixXd Hi = xxi * Wi * xxi.adjoint();
  MatrixXd Hj = xxj * Wj * xxj.adjoint();

  MatrixXd Hi0 = xxi0 * solve_ident_cpp(xxi0.adjoint() * xxi0) * xxi0.adjoint();
  MatrixXd Hj0 = xxj0 * solve_ident_cpp(xxj0.adjoint() * xxj0) * xxj0.adjoint();

  MatrixXd SiR = Hi - Hi0;
  MatrixXd SjR = Hj - Hj0;
  // Rcout << "SiR" <<SiR <<endl;

  MatrixXd npDiag = MatrixXd::Identity(np, np);
  MatrixXd SiE = npDiag - Hi;
  MatrixXd SjE = npDiag - Hj;

  // Calculate rSSR ----
  // temporary R matrix
  MatrixXd R = Rij * SjR * Rij.adjoint();
  // turn R matrix into a column vector: not the .array() method does this wrong
  VectorXd Rvec(R.rows() * R.cols());
  int iter = 0;
  for(int c = 0; c < R.cols(); ++c){
    for(int r = 0; r < R.rows(); ++r){
      Rvec(iter) = R(r, c);
      ++iter;
    }
  }
  // turn SiR into a column vector
  VectorXd SiRvec(SiR.rows() * SiR.cols());
  iter = 0;
  for(int c = 0; c < SiR.cols(); ++c){
    for(int r = 0; r < SiR.rows(); ++r){
      SiRvec(iter) = SiR(r, c);
      ++iter;
    }
  }
  // calculate rSSRij
  MatrixXd rSSRij = (SiRvec.adjoint() * Rvec).array()/df1;

  // Calculate rSSE ----
  // temporary E matrix
  MatrixXd E = Rij * SjE * Rij.adjoint();
  // turn E into column vector
  VectorXd Evec(E.rows() * E.cols());
  iter = 0;
  for(int c = 0; c < E.cols(); ++c){
    for(int r = 0; r < E.rows(); ++r){
      Evec(iter) = E(r, c);
      ++iter;
    }
  }
  // turn SiE into a column vector
  VectorXd SiEvec(SiE.rows() * SiE.cols());
  iter = 0;
  for(int c = 0; c < SiE.cols(); ++c){
    for(int r = 0; r < SiE.rows(); ++r){
      SiEvec(iter) = SiE(r, c);
      ++iter;
    }
  }
  // calculate rSSEij
  MatrixXd rSSEij = (SiEvec.adjoint() * Evec).array()/df2;

  // calculate rcoef
  MatrixXd Vcoefij = Wi * (xxi.adjoint() * Rij * xxj) * Wj.adjoint();

  // -- Changes from Tony 06-June-2022 --
  MatrixXd sqrtdiag_i = pow(Wi.diagonal().array(), -0.5);
  MatrixXd sqrtdiag_j = pow(Wj.diagonal().array(), -0.5);
  MatrixXd rcoefij = sqrtdiag_i.asDiagonal() * Vcoefij * sqrtdiag_j.asDiagonal();
  // -- Older --
  // MatrixXd rcoefij = Vcoefij.array() * pow(Wi.array()*Wj.array(), -0.5);
  // MatrixXd rcoefij = Vcoefij.diagonal().array() *
  //   pow(Wi.diagonal().array() * Wj.diagonal().array(), -0.5);

  // output ----
  List out_lst = List::create(Named("Rij") = Rij,
                              Named("Hi") = Hi,
                              Named("Hj") = Hj,
                              Named("Hi0") = Hi0,
                              Named("Hj0") = Hj0,
                              Named("SiR") = SiR,
                              Named("SjR") = SjR,
                              Named("rcoefij") = rcoefij,
                              Named("rSSRij") = rSSRij.matrix(),
                              Named("rSSEij") = rSSEij.matrix());

  if(Vcoef){
    out_lst.push_back(Vcoefij, "Vcoefij");
  }

  return out_lst;
}
