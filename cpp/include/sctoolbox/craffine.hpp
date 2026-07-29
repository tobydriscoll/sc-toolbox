#pragma once
#include <Eigen/Dense>

#include "sctoolbox/qgraph.hpp"

namespace sctoolbox {

// Port of @crdiskmap/private/craffine.m, restricted to the case where w is
// fully known (no NaN entries to be deduced) -- the only form needed by
// crparam, which always has the complete target polygon available.
// Returns the (n-3) x 2 affine table.
Eigen::MatrixXcd craffine(const Eigen::VectorXcd& w, const Eigen::VectorXd& beta, const Eigen::VectorXd& cr,
                          const QGraph& Q, double tol = 1e-8);

}  // namespace sctoolbox
