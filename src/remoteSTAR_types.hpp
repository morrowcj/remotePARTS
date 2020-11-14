// [[Rcpp::depends(RcppEigen)]]
#pragma once

#include <RcppEigen.h>
#include <math.h>
#include <iostream>

using Eigen::LLT;
using Eigen::Lower;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::Upper;
using Eigen::VectorXd;
using Eigen::VectorXi;

typedef Map<MatrixXd> MapMatd;
typedef Map<MatrixXi> MapMati;
typedef Map<VectorXd> MapVecd;

using namespace std;
using namespace Rcpp;
