#pragma once
#include <Eigen/Dense>

namespace sctoolbox {

// Gauss-Jacobi quadrature data for SC routines. BETA holds the turning
// angles of the finite singularities; columns with beta(j) <= -1 (or NaN)
// are left as zero, matching +sctool/scqdata.m. Returns an
// nqpts x 2*(n+1) matrix [qnode, qwght].
Eigen::MatrixXd scqdata(const Eigen::VectorXd& beta, int nqpts);

}  // namespace sctoolbox
