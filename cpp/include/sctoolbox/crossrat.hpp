#pragma once
#include <Eigen/Dense>

#include "sctoolbox/qgraph.hpp"

namespace sctoolbox {

// Port of @crdiskmap/private/crossrat.m: the n-3 target crossratios implied
// by the polygon w and its quadrilateral graph Q.
Eigen::VectorXcd crossrat(const Eigen::VectorXcd& w, const QGraph& Q);

}  // namespace sctoolbox
