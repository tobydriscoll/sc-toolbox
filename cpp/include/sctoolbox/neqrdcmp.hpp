#pragma once
#include <Eigen/Dense>

namespace sctoolbox {

// QR decomposition encoded as Householder rotations (Dennis & Schnabel
// Algorithm A3.2.1). Square matrices only. M is overwritten in place with
// the encoded factorization, matching +sctool/neqrdcmp.m's [M,M1,M2,sing].
void neqrdcmp(Eigen::MatrixXd& M, Eigen::VectorXd& M1, Eigen::VectorXd& M2, int& sing);

}  // namespace sctoolbox
