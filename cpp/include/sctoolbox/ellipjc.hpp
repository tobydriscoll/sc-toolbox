#pragma once
#include <Eigen/Dense>

namespace sctoolbox {

struct EllipJC {
    Eigen::VectorXcd sn, cn, dn;
};

// Port of @rectmap/private/ellipjc.m. u may be a vector; L is a scalar
// parameter (m = exp(-2*pi*L)). If mIsM is true, L is interpreted directly
// as the recursive call's "m" argument (mirrors MATLAB's 3-argument
// recursive-call convention, where the high/low rectangle-half folding is
// skipped).
EllipJC ellipjc(const Eigen::VectorXcd& u, double L, bool mIsM = false);

}  // namespace sctoolbox
