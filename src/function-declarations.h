#ifndef FUNC_DECL
#define FUNC_DECL

#include "remotePARTS_types.h"

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
MatrixXd invchol_cpp(const MapMatd& V, double nugget, int ncores);

double optimize_nugget_cpp(const MapMatd& X, const MapMatd& X0, const MapMatd& V,
                           const MapMatd& y, double lower, double upper, double tol,
                           const MapMatd& invCholV, bool use_invCholV, bool debug,
                           int ncores);

// fit GLS
List fitGLS_cpp(const MapMatd& X, const MapMatd& V, const MapMatd& y,
                const MapMatd& X0, double nugget, bool save_xx, bool save_invchol,
                bool LL_only, bool no_F, bool optimize_nugget, double nug_l,
                double nug_u, double nug_tol, const MapMatd& invCholV,
                bool use_invCholV, int ncores);

#endif // FUNC_DECL

#ifdef _OPENMP
#include <omp.h>
#endif
