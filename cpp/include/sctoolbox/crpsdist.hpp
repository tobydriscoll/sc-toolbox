#pragma once
#include <Eigen/Dense>
#include <array>

namespace sctoolbox {

// Port of @crdiskmap/private/crpsdist.m: distance from each point in pts to
// the line segment described by segment's two endpoints.
Eigen::VectorXd crpsdist(const std::array<std::complex<double>, 2>& segment, const Eigen::VectorXcd& pts);

}  // namespace sctoolbox
