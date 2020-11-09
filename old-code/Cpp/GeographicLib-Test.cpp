// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]
#include <iostream>
#include <RcppEigen.h>
#include <GeographicLib/Geodesic.hpp>
// #include <GeographicLib/Math.hpp>
// #include <GeographicLib/Constants.hpp>

using namespace GeographicLib;
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
double GDist_cpp(double lat1, double lon1, double lat2, double lon2){
    double s12;
    Geodesic geod(Constants::WGS84_a(), Constants::WGS84_f());
    geod.Inverse(lat1, lon1, lat2, lon2, s12);
    return s12;
}

// /*** R
// lat1 = 40.6; lon1 = -73.8 # JFK Airport
// lat2 = 51.6; lon2 = -0.5 # LHR Airport
// GDist_cpp(lat1, lon1, lat2, lon2)
// */

int main(){
  double out = GDist_cpp(40.6, -73.8, 51.6, -.5);
  cout <<"distance: " << out << endl;
  return 0;
}
