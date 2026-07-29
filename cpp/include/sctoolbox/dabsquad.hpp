#pragma once
#include <Eigen/Dense>
#include <vector>

namespace sctoolbox {

// Port of +sctool/dabsquad.m. sing1 follows MATLAB's 1-indexed convention
// directly: sing1[k] == 0 means z1(k) is not a singularity, otherwise
// sing1[k] is the 1-indexed position of z1(k) within z. Returns a real
// vector (the absolute-value side-length integral is real by construction).
Eigen::VectorXd dabsquad(const Eigen::VectorXcd& z1, const Eigen::VectorXcd& z2, const std::vector<int>& sing1,
                         const Eigen::VectorXcd& z, const Eigen::VectorXd& beta, const Eigen::MatrixXd& qdat);

}  // namespace sctoolbox
