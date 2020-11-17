#ifndef FUNC_DECL
#define FUNC_DECL

#include "remoteSTAR_types.h"

/*
 * Helper Functions ----
 */

// t(A) %*% A
MatrixXd AtA(const MatrixXd& A);
// solve Ax + B
MatrixXd solve_cpp(const MatrixXd& A, const MatrixXd& B);
// solve Ax + I
MatrixXd solve_ident_cpp(const MatrixXd& A);

/*
 * Main Functions
 */

// calculate inverse cholesky matrix
MatrixXd invchol_cpp(const MapMatd& V, double nugget);
// fit GLS
List fitGLS_cpp(const MapMatd& X, const MapMatd& V, const MapMatd& y,
                const MapMatd& X0, double nugget, bool save_xx,
                const int threads);
// calculate log-likelihood only
double LogLikGLS_cpp(double nugget, const MapMatd& X, const MapMatd& V,
                     const MapMatd& y);
// find the maximum likelihood nugget
double optimize_nugget_cpp(const MapMatd& X, const MapMatd& V, const MapMatd& y,
                          double lower, double upper, double tol,
                          bool debug);
// worker function to execute fitGLS()
List GLS_worker_cpp(const MapMatd& y, const MapMatd& X, const MapMatd& V,
                    const MapMatd& X0, double nug_l, double nug_u,
                    double nug_tol, bool save_xx);
// worker function to compare partitions from GLS_worker_cpp
List crosspart_worker_cpp(const MapMatd& xxi, const MapMatd& xxj,
                          const MapMatd& xxi0, const MapMatd& xxj0,
                          const MapMatd& tUinv_i, const MapMatd& tUinv_j,
                          const MapMatd& Vsub, int df1, int df2);

#endif // FUNC_DECL
