#pragma once
#include <Eigen/Dense>
#include <vector>

namespace sctoolbox {

// Port of @stripmap/private/stquadh.m. Recursively subdivides the
// (assumed horizontal) integration interval until the "alpha rule" is
// satisfied for singularities lying between the endpoints, then
// delegates to stquad for each safe sub-interval.
Eigen::VectorXcd stquadh(const Eigen::VectorXcd& z1, const Eigen::VectorXcd& z2, const std::vector<int>& sing1,
                         const Eigen::VectorXcd& z, const Eigen::VectorXd& beta, const Eigen::MatrixXd& qdat);

}  // namespace sctoolbox
