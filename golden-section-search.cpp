  // [[Rcpp::export]]
  inline double nugOptim_cpp(const MapMatd& X,
                  const MapMatd& V,
                  const MapVecd& y,
                  double lower = 0,
                  double upper = 1,
                  double tol = 0.00001){
    double mx; // declare x which maximise f(x)

    // golden ratio value
    double GR = 2/sqrt(5) + 1;

    // test bounds
    double x1 = upper - GR*(upper - lower);
    double x2 = lower + GR*(upper - lower);

    // current f() bounds
    double f1 = LogLikGLS_cpp(x1, X, V, y);
    double f2 = LogLikGLS_cpp(x2, X, V, y);

    int i = 0;
    while (abs(upper - lower ) > tol){
      i += 1;
      // cout << "current interval: [" << lower << ", " << upper << "]" << endl;
      if (f2 < f1) { // then the maximum is closer to x2
        // recycle new values according to the GR rule
        upper = x2;
        x2 = x1;
        f2 = f1;
        // calculate new values
        x1 = upper - GR*(upper - lower);
        f2 = LogLikGLS_cpp(x2, X, V, y);
      } else { // the maximum is closer to x1
        lower = x1;
        x1 = x2;
        f1 = f2;

        x2 = lower + GR*(upper - lower);
        f2 = LogLikGLS_cpp(x2, X, V, y);
      }
    }
    cout << "number of iterations: " << i << endl;

    mx = (upper + lower)/2;

    if(mx < tol){
      double f0 = LogLikGLS_cpp(0, X, V, y);
      double fmx = LogLikGLS_cpp(mx, X, V, y);
      if (f0 > fmx){
        mx = 0;
      }
    }

    if(mx > 1 - tol){
      double f0 = LogLikGLS_cpp(1, X, V, y);
      double fmx = LogLikGLS_cpp(mx, X, V, y);
      if (f0 > fmx){
        mx = 1;
      }
    }

    return mx;
  }
