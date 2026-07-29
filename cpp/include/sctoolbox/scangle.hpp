#pragma once
#include <Eigen/Dense>

namespace sctoolbox {

// Turning angles of a polygon with vertices W (port of +sctool/scangle.m).
// A vertex's angle is NaN if it or either neighbor is at infinity.
Eigen::VectorXd scangle(const Eigen::VectorXcd& w);

}  // namespace sctoolbox
