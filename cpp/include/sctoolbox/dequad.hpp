#pragma once
#include <Eigen/Dense>
#include <vector>

namespace sctoolbox {

// Port of @extermap/private/dequad.m. z holds the finite prevertices
// (origin excluded -- it is handled internally via an implicit beta=-2
// term). sing1 follows MATLAB's 1-indexed convention directly:
// sing1[k] == 0 means z1(k) is not a singularity, otherwise sing1[k] is
// the 1-indexed position of z1(k) within z.
Eigen::VectorXcd dequad(const Eigen::VectorXcd& z1, const Eigen::VectorXcd& z2, const std::vector<int>& sing1,
                        const Eigen::VectorXcd& z, const Eigen::VectorXd& beta, const Eigen::MatrixXd& qdat);

}  // namespace sctoolbox
