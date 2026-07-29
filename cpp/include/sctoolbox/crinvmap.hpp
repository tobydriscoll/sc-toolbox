#pragma once
#include <Eigen/Dense>
#include <vector>

#include "sctoolbox/qgraph.hpp"

namespace sctoolbox {

// Port of @crdiskmap/private/crinvmap.m.
Eigen::VectorXcd crinvmap(const Eigen::VectorXcd& wp, const Eigen::VectorXcd& w, const Eigen::VectorXd& beta,
                          const Eigen::VectorXd& cr, const Eigen::MatrixXcd& aff, const Eigen::VectorXcd& wcfix,
                          const QGraph& Q, const Eigen::MatrixXd& qdat, const std::vector<double>& options = {});

}  // namespace sctoolbox
