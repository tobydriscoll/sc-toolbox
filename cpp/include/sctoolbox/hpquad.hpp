#pragma once
#include <Eigen/Dense>
#include <vector>

namespace sctoolbox {

// Port of @hplmap/private/hpquad.m's single-sided (6-argument) form.
// sing1 follows MATLAB's 1-indexed convention directly: sing1[k] == 0
// means z1(k) is not a singularity, otherwise sing1[k] is the 1-indexed
// position of z1(k) within z (the finite prevertices only).
Eigen::VectorXcd hpquad(const Eigen::VectorXcd& z1, const Eigen::VectorXcd& z2, const std::vector<int>& sing1,
                        const Eigen::VectorXcd& z, const Eigen::VectorXd& beta, const Eigen::MatrixXd& qdat);

}  // namespace sctoolbox
