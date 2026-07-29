#pragma once
#include <Eigen/Dense>
#include <vector>

namespace sctoolbox {

// Port of @crdiskmap/private/crimap0.m.
Eigen::VectorXcd crimap0(const Eigen::VectorXcd& wp, const Eigen::VectorXcd& z, const Eigen::VectorXd& beta,
                         const Eigen::Vector2cd& aff, const Eigen::MatrixXd& qdat,
                         const std::vector<double>& options = {});

}  // namespace sctoolbox
