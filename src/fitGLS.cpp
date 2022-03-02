#include "function-declarations.h"

//' Fit GLS to remote sensing data
//'
//' @param X numeric model matrix
//' @param V numeric covariance matrix
//' @param y numeric resposne vector
//' @param X0 numeric null model matrix
//' @param nugget numeric nugget
//' @param save_xx logical: should cross-partition stats be saved
//' @param save_invchol logical: should invCholV be returned?
//' @param LL_only logical: should only log-liklihood be computed?
//' @param no_F logical: should calculations F-test calculations be skipped?
//' @param optimize_nugget logical: should the ML nugget be obtained?
//' @param nug_l numeric lower value for nugget search
//' @param nug_u numeric upper value for nugget search
//' @param nug_tol numeric tolerance for nugget search
//' @param invCholV numeric inverse cholesky matrix
//' @param use_invCholV logical: should invCholV be used instead of V
//' @param ncores integer indicating cores to use
//'
//' @examples
//' #   data.file = system.file("extdata", "AK_ndvi_common-land.csv", package = "remotePARTS")
//' #
//' #   n = 1000
//' #
//' #   df = data.table::fread(data.file, nrows = n) # read first 1000 rows
//' #
//' #   ## format data
//' #   datalist = part_data(1, part_form = cls.coef ~ 0 + land, part_df = df,
//' #             part_mat = matrix(1:n, ncol = 1))
//' #
//' #   ## fit covariance matrix
//' #   V = covar_exp(distm_km(cbind(df$lng, df$lat)), range = .01)
//' #
//' #   # use V matrix
//' #   .fitGLS_cpp(X = datalist$X, V = V, y = datalist$y, X0 = datalist$X0, invCholV = diag(1),
//' #                nugget = 0.0, save_xx = FALSE, save_invchol = FALSE, LL_only = FALSE, no_F = FALSE,
//' #                use_invCholV = FALSE, optimize_nugget = FALSE, nug_l = 0.0, nug_u = 1.0, nug_tol = 1e-7)
//' #
//' #   # use inverse cholesky instead
//' #   .fitGLS_cpp(X = datalist$X, V = diag(1), y = datalist$y, X0 = datalist$X0, invCholV = invert_chol(V),
//' #                nugget = 0.0, save_xx = FALSE, save_invchol = FALSE, LL_only = FALSE, no_F = FALSE,
//' #                use_invCholV = TRUE, optimize_nugget = FALSE, nug_l = 0.0, nug_u = 1.0, nug_tol = 1e-7)
//' #
//' #   # optimize nugget
//' #   .fitGLS_cpp(X = datalist$X, V = V, y = datalist$y, X0 = datalist$X0, invCholV = diag(1),
//' #                nugget = 0.0, save_xx = FALSE, save_invchol = FALSE, LL_only = FALSE, no_F = FALSE,
//' #                use_invCholV = FALSE, optimize_nugget = TRUE, nug_l = 0.0, nug_u = 1.0, nug_tol = 1e-7)
//'
// [[Rcpp::export(.fitGLS_cpp)]]
List fitGLS_cpp(const MapMatd& X,
                 const MapMatd& V,
                 const MapMatd& y,
                 const MapMatd& X0,
                 double nugget,
                 bool save_xx,
                 bool save_invchol,
                 bool LL_only,
                 bool no_F,
                 bool optimize_nugget,
                 double nug_l,
                 double nug_u,
                 double nug_tol,
                 const MapMatd& invCholV,
                 bool use_invCholV,
                 int ncores){

  Eigen::setNbThreads(ncores);

  int nX = X.rows(), pX = X.cols(); // dimensions of X
  double nug;
  // Optimize nugget, if asked
  nug = (optimize_nugget) ? optimize_nugget_cpp(X, X0, V, y, nug_l, nug_u, nug_tol, invCholV, false, false, ncores) : nugget;

  // Inverse Cholesky Matrix
  MatrixXd iChol;
  iChol = (use_invCholV) ? invCholV : invchol_cpp(V, nug, ncores);
  // matrix multiplication
  MatrixXd xx = iChol * X; // t(inv(U)) %*% X
  MatrixXd yy = iChol * y; // t(inv(U)) %*% y
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
  for (int i = 0; i < iChol.rows(); i++){
    logDetV += log(iChol.diagonal()[i]);
  }
  logDetV *= -2;
  double logLik = -0.5 * (nX * log(2 * M_PI) + nX *
                          log(dft * MSE.array()[0]/nX) + logDetV + nX);
  // // add logLik to output list
  // L["logLik"] = logLik;

  // exit the function if LL_only
  if (LL_only) {
    List LL_list = List::create(Named("logLik") = logLik);
    return LL_list;
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

  // // Add all additional stats, up to t-test, into output list
  // L["coefficients"] = betahat;
  // L["SSE"] = SSE;
  // L["MSE"] = MSE;
  // L["SE"] = se;
  // L["df_t"] = dft;
  // L["logDetV"] = logDetV;
  // L["tstat"] = tstat;
  // // add xx and invcholV to output, if requested
  // if (save_xx){
  //   L["xx"] = xx;
  //   L["invcholV"] = iChol;
  // }
  // break function before calculating F statistics, if asked

  if (no_F) {
    List No_F = List::create(Named("coefficients") = betahat,
                             Named("covar_coef") = varcov,
                             Named("SSE") = SSE,
                             Named("MSE") = MSE,
                             Named("SE") = se,
                             Named("df_t") = dft,
                             Named("logDetV") = logDetV,
                             Named("tstat") = tstat,
                             Named("pval_t") = NA_REAL,
                             Named("logLik") = logLik,
                             Named("nugget") = nug);

    if (save_xx){
      No_F.push_back(xx, "xx");
    } else {
      No_F.push_back(NA_REAL, "xx");
    }

    if(save_invchol){
      No_F.push_back(iChol, "invcholV");
    } else {
      No_F.push_back(NA_REAL, "invcholV");
    }
    return No_F;
  }

  // **START F-test** ----
  // null model stats
  MatrixXd xx0 = iChol * X0;
  MatrixXd varX0 = AtA(xx0); // crossprod(xx0)
  MatrixXd X0tY(xx0.adjoint() * yy); // crossprod(xx0, yy)
  VectorXd betahat0(solve_cpp(varX0, X0tY));
  // df
  int df0 = betahat0.size();
  // SS
  VectorXd SSE0 = AtA(yy - xx0 * betahat0); // SSE
  double MSE0 = SSE0.array()[0]/(nX - xx.cols()); // MSE
  double MSR = (SSE0.array()[0] - SSE.array()[0])/(xx.cols() - xx0.cols());
  // logLik
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

  // L["coefficients0"] = betahat0;
  // L["SSE0"] = SSE0;
  // L["MSE0"] = MSE0;
  // L["SE0"] = se0;
  // L["MSR"] = MSR;
  // L["df0"] = df0;
  // L["logLik0"] = logLik0;
  // L["df_F"] = dfF;
  // L["Fstat"] = FF;
  //
  // // add in xx0, if asked
  // if (save_xx){
  //   L["xx0"] = xx0;
  // }

  List Full_out = List::create(Named("coefficients") = betahat,
                               // Named("covar_coef") = varcov,
                               Named("SSE") = SSE,
                               Named("MSE") = MSE,
                               Named("SE") = se,
                               Named("df_t") = dft,
                               Named("logDetV") = logDetV,
                               Named("tstat") = tstat,
                               Named("pval_t") = NA_REAL,
                               Named("logLik") = logLik,
                               Named("nugget") = nug,
                               Named("coefficients0") = betahat0,
                               Named("SSE0") = SSE0,
                               Named("MSE0") = MSE0,
                               Named("SE0") = se0,
                               Named("MSR") = MSR,
                               Named("df0") = df0,
                               Named("logLik0") = logLik0,
                               Named("df_F") = dfF,
                               Named("Fstat") = FF,
                               Named("pval_F") = NA_REAL);

  Full_out.push_back(varcov, "covar_coef"); // only way I could get this term to work

  if (save_xx){
    Full_out.push_back(xx, "xx");
    Full_out.push_back(xx0, "xx0");
  } else {
    Full_out.push_back(NA_REAL, "xx");
    Full_out.push_back(NA_REAL, "xx0");
  }

  if (save_invchol){
    Full_out.push_back(iChol, "invcholV");
  } else {
    Full_out.push_back(NA_REAL, "invcholV");
  }

  return Full_out;
}
