#pragma once
#include <Eigen/Dense>

#include "sctoolbox/qgraph.hpp"

namespace sctoolbox {

// Port of @crdiskmap/private/crmap.m. wcfix/aff as in crderiv.
Eigen::VectorXcd crmap(const Eigen::VectorXcd& zp, const Eigen::VectorXcd& w, const Eigen::VectorXd& beta,
                       const Eigen::VectorXd& cr, const Eigen::MatrixXcd& aff, const Eigen::VectorXcd& wcfix,
                       const QGraph& Q, const Eigen::MatrixXd& qdat);

}  // namespace sctoolbox
