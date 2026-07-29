#pragma once
#include <Eigen/Dense>

#include "sctoolbox/qgraph.hpp"

namespace sctoolbox {

// Port of @crdiskmap/private/crderiv.m. wcfix is the 5-element vector
// [quadnum(1-indexed, as in MATLAB), mt1, mt2, mt3, mt4] returned by
// crfixwc; aff is the n3 x 2 affine table (one row per quadrilateral).
Eigen::VectorXcd crderiv(const Eigen::VectorXcd& zp, const Eigen::VectorXd& beta, const Eigen::VectorXd& cr,
                         const Eigen::MatrixXcd& aff, const Eigen::VectorXcd& wcfix, const QGraph& Q);

}  // namespace sctoolbox
