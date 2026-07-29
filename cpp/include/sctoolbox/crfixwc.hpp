#pragma once
#include <Eigen/Dense>

#include "sctoolbox/qgraph.hpp"

namespace sctoolbox {

// Port of @crdiskmap/private/crfixwc.m (the wc-given form only; the
// interactive ginput-based prompt for omitted wc is not applicable here).
// Returns the 5-element wcfix vector [quadnum(1-indexed as double), mt1..mt4].
Eigen::VectorXcd crfixwc(const Eigen::VectorXcd& w, const Eigen::VectorXd& beta, const Eigen::VectorXd& cr,
                         const Eigen::MatrixXcd& aff, const QGraph& Q, std::complex<double> wc);

}  // namespace sctoolbox
