#pragma once
#include <Eigen/Dense>

namespace sctoolbox {

// Perturbed Cholesky decomposition (Dennis & Schnabel Algorithm A5.5.2).
// maxoffl=0 means H is assumed positive definite (the only case exercised
// by real call sites in nehook/nemodel; the maxoffl!=0 branch of the
// original MATLAB references an undefined variable and is dead code).
void nechdcmp(const Eigen::MatrixXd& H, double maxoffl, Eigen::MatrixXd& L, double& maxadd);

}  // namespace sctoolbox
