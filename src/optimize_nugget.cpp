#include "function-declarations.h"

//' Find the maximum likelihood estimate of the nugget
//'
//' @details this is the C++ version of `optimize()` which is specific to
//' finding the nugget value that maximizes the log-likelihood of `fitGLS_cpp()`
//'
//' This function is a translation from the forchan algorithm fmin into C++:
//' http://www.netlib.org/fmm/fmin.f
//'
//' Note: this function actually uses `LogLikGLS_cpp()` which should be swapped
//' for `fitGLS_cpp()` once the correct functionality is added to the latter.
//'
//' @param X numeric matrix
//' @param V numeric matrix
//' @param y numeric vector
//' @param lower lower boundary for nugget search
//' @param upper upper boundary for nugget search
//' @param tol desired accuracy of nugget search
//' @param debug logical: debug mode?
//'
//' @examples #TBA
// [[Rcpp::export(.optimize_nugget_cpp)]]
double optimize_nugget_cpp(const MapMatd& X, const MapMatd& V, const MapMatd& y,
                          double lower, double upper,
                          double tol, bool debug){

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
      cout << "x = " << x << " fx = " << fx << endl;
      i += 1;
      if (i >= 100) {
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
                // cout << "returning starting value instead of f_min" << endl;
                return ax;
              }
            }
            return f_min;
}
