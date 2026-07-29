#pragma once
#include <Eigen/Dense>
#include <array>
#include <limits>

namespace sctoolbox {

struct ScaddvtxResult {
    Eigen::VectorXcd w;
    Eigen::VectorXd beta;
};

// Insert a new vertex immediately after w[pos] (0-indexed). Port of
// +sctool/scaddvtx.m. `window` bounds the new vertex when it must be
// extended away from an infinite neighbor (re-min, re-max, im-min, im-max).
ScaddvtxResult scaddvtx(const Eigen::VectorXcd& w, const Eigen::VectorXd& beta, int pos,
                         std::array<double, 4> window = {-std::numeric_limits<double>::infinity(),
                                                           std::numeric_limits<double>::infinity(),
                                                           -std::numeric_limits<double>::infinity(),
                                                           std::numeric_limits<double>::infinity()});

}  // namespace sctoolbox
