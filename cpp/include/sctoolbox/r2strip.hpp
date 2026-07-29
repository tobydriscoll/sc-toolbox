#pragma once
#include <Eigen/Dense>

namespace sctoolbox {

struct R2Strip {
    Eigen::VectorXcd yp, yprime;
};

// Port of @rectmap/private/r2strip.m. Maps from the rectangle (with
// prevertices z, only the corners of which matter) to the strip
// 0 <= Im z <= 1, via log(sn(z|m))/pi with m = exp(-2*pi*L).
R2Strip r2strip(const Eigen::VectorXcd& zp, const Eigen::VectorXcd& z, double L);

}  // namespace sctoolbox
