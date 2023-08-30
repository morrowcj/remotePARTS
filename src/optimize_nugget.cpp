#include "function-declarations.h"

//' Find the maximum likelihood estimate of the nugget
//'
//' @details this is the C++ version of `optimize()` which is specific to
//' finding the nugget value that maximizes the log-likelihood of `fitGLS_cpp()`
//' by minimizing the partial log likelihood (i.e., fitGLS_cpp(LL_only = TRUE)[["logLik"]] )
//'
//' This function is a translation from the forchan algorithm fmin into C++:
//' http://www.netlib.org/fmm/fmin.f
//'
//' @param X numeric model matrix
//' @param X0 numeric null model matrix (needed but not used)
//' @param V numeric covariance matrix
//' @param y numeric resposne vector
//' @param lower lower boundary for nugget search
//' @param upper upper boundary for nugget search
//' @param tol desired tolerance for nugget search
//' @param invchol numeric inverse cholesky matrix
//' @param use_invchol logical: should invchol be used instead of V?
//' @param debug logical: debug mode?
//' @param ncores integer indicating number of cores to use
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
//' #  .optimize_nugget_cpp(X = datalist$X, X0 = datalist$X0, V = V, y = datalist$y,
//' #                       lower = 0, upper = 1, tol = 1e-10, invchol = diag(1),
//' #                       use_invchol = FALSE, debug = TRUE)
//' #
//' #  .optimize_nugget_cpp(X = datalist$X, X0 = datalist$X0, V = diag(1), y = datalist$y,
//' #                       lower = 0, upper = 1, tol = 1e-10, invchol = invert_chol(V),
//' #                       use_invchol = TRUE, debug = TRUE)
//' #
//' #  warning("optimize_nugget_cpp CANNOT recycle the invchol matrix!! remove this functionality")
//'
// [[Rcpp::export(.optimize_nugget_cpp)]]
double optimize_nugget_cpp(const MapMatd& X, const MapMatd& X0, const MapMatd& V, const MapMatd& y,
                           double lower, double upper, double tol,
                           const MapMatd& invchol, bool use_invchol,
                           bool debug, int ncores){
  Eigen::setNbThreads(ncores);

  // varible declaration
  double ax = lower;
  double bx = upper;
  double a,b,c,d,e,eps,xm,p,q,r,tol1,tol2,u,v,w;
  double fv,fw, x, f_min;
  NumericVector fxV, fxV2, fuV;
  double fx, fx2, fu;
  List fxL, fuL, fxL2;

  // c is the squared inverse of the golden ratio
  c = 0.5*(3. - sqrt(5.));

  // eps is approximately sqrt of machine precision
  // eps = std::numeric_limits<double>::epsilon();
  // eps = sqrt(eps);
  // Rcout << "eps_std = "<<eps<<endl;

  eps = 1.;
  lab0:
    eps = eps/2.;
  tol1 = 1. + eps;
  if (tol1 > 1.) {goto lab0;}
  eps = sqrt(eps);
  // Rcout << "eps =" << eps<< endl;

  // initialization
  a = ax;
  b = bx;
  v = a + c*(b - a);
  w = v;
  x = v;
  e = 0.;
  // fx = -LogLikGLS_cpp(x, X, V, y, invchol, use_invchol); //invert function to minimize
  fxL = fitGLS_cpp(X, V, y, X0, x, false, false, true, false, false, 0., 0., 0., invchol, false, ncores);
  fxV = fxL["logLik"];
  fx = fxV[0] * -1;
  fv = fx;
  fw = fx;

  int i = 0;

  // main loop start
  lab1:
    // Rcout << "Loop Start (lab1): iteration " << i << endl;
    if (debug) {
      Rcout << "i = " << i << ", x = " << x << ", fx = " << fx << endl;
      i += 1;
      if (i >= 1000) {
      Rcout << "breaking loop, too many iterations"<<endl;
        goto lab8;
      }
    }
    xm = 0.5*(a + b);
    tol1 = eps*abs(x) + tol/3.;
    tol2 = 2.*tol1;
    // check stoping criteria
    if (abs(x - xm) <= (tol2 - 0.5*(b - a))) {goto lab8;}
    // Rcout << "stop crit. not met: "<<abs(x-xm)<<" > "<<tol2-.5*(b-a)<<endl;
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
    // Rcout << "check parabola (lab2)" << endl;
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
    // Rcout << "golden section step (lab3)" <<endl;
    // golden section step
    if (x >= xm) {e = a - x;}
    if (x < xm) {e = b - x;}
    d = c*e;

  lab4:
    // Rcout << "check tolerance and update vars (lab4)" << endl;
    //f must not be evaluated too close to x
    if (abs(d) >= tol1) {u = x + d;}
    if (abs(d) < tol1) {u = x + copysign(tol1, d);}
    // fu = -LogLikGLS_cpp(u, X, V, y, invchol, use_invchol);
    fuL = fitGLS_cpp(X, V, y, X0, u, false, false, true, false, false, 0., 0., 0., invchol, false, ncores);
    fuV = fuL["logLik"];
    fu = fuV[0] * -1;
    // Rcout << "u = " << u << " fu = " << fu << endl;
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
    // Rcout << "conditional variable reset (lab5)" << endl;
    if (u < x) {a = u;}
    if (u >= x) {b = u;}
    if (fu <= fw) {goto lab6;}
    if (w == x) {goto lab6;}
    if (fu <= fv) {goto lab7;}
    if (v == x) {goto lab7;}
    if (v == w) {goto lab7;}
    goto lab1;

  lab6:
    // Rcout << "update function results (lab6)" << endl;
    v = w;
    fv = fw;
    w = u;
    fw = fu;
    goto lab1;

  lab7:
    // Rcout << "update function results alternate (lab7)" << endl;
    v = u;
    fv = fu;
    goto lab1;

  // end of main loop
  lab8:
    // Rcout << "return statement (lab8)" << endl;
    // Rcout << "x = " << x << " fx = " << fx << "a = " << a << "ax = " << ax << endl;
    f_min = x;
    if (ax + tol >= f_min){
      // if (fx <= -LogLikGLS_cpp(f_min, X, V, y, invchol, use_invchol)){
      fxL2 = fitGLS_cpp(X, V, y, X0, f_min, false, false, true, false, false, 0., 0., 0., invchol, false, ncores);
      fxV2 = fxL2["logLik"];
      fx2 = fxV2[0] * -1;
      if (fx <= fx2){
        // Rcout << "returning starting value instead of f_min" << endl;
        // Rcout << "ax = " << ax << endl;
        return ax;
      }
    }
    // Rcout << "f_min = " << f_min << endl;
    return f_min;
}
