#pragma once
#include <Eigen/Dense>
#include <array>

namespace sctoolbox {

struct StParamResult {
    Eigen::VectorXcd z;       // length N, original vertex order, with -Inf/+Inf at the ends
    std::complex<double> c;
    Eigen::MatrixXd qdat;     // N+1 column-pair size, original vertex order
};

// Port of @stripmap/private/stparam.m, restricted to the z0=[] (automatic
// initial guess) calling form, matching dparam/hpparam/deparam's
// restriction. `ends` is the 1-indexed 2-vector of vertex positions
// mapping to the strip's left/right ends, as in the .m source.
StParamResult stparam(const Eigen::VectorXcd& w, const Eigen::VectorXd& beta, std::array<int, 2> ends,
                      double tol = 1e-8, int method = 2);

}  // namespace sctoolbox
