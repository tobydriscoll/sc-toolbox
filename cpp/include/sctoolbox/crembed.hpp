#pragma once
#include <Eigen/Dense>

#include "sctoolbox/qgraph.hpp"

namespace sctoolbox {

// Port of @crdiskmap/private/crembed.m. qnum is 0-indexed.
Eigen::VectorXcd crembed(const Eigen::VectorXd& cr, const QGraph& Q, int qnum);

}  // namespace sctoolbox
