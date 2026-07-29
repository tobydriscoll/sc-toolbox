#pragma once
#include <Eigen/Dense>
#include <vector>

#include "sctoolbox/dinvmap.hpp"

namespace sctoolbox {

// Port of @stripmap/private/stinvmap.m.
InvMapResult stinvmap(const Eigen::VectorXcd& wp, const Eigen::VectorXcd& w, const Eigen::VectorXd& beta,
                       const Eigen::VectorXcd& z, std::complex<double> c, const Eigen::MatrixXd& qdat,
                       const Eigen::VectorXcd& z0 = {}, const std::vector<double>& options = {});

}  // namespace sctoolbox
