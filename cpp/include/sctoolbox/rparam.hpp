#pragma once
#include <Eigen/Dense>
#include <array>

namespace sctoolbox {

struct RParamResult {
    Eigen::VectorXcd z;  // length n, rectangle prevertices, original vertex order
    std::complex<double> c;
    double L;
    Eigen::MatrixXd qdat;  // n+1 column-pair size, original vertex order
};

// Port of @rectmap/private/rparam.m, restricted to the z0=[] (automatic
// initial guess) calling form, matching dparam/hpparam/deparam/stparam's
// restriction. `cnr` is the 1-indexed 4-vector of corner positions, as in
// the .m source (first two entries describe a long side, CCW order).
RParamResult rparam(const Eigen::VectorXcd& w, const Eigen::VectorXd& beta, std::array<int, 4> cnr,
                    double tol = 1e-8, int method = 2);

}  // namespace sctoolbox
