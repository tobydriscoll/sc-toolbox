#pragma once
#include <Eigen/Dense>

namespace sctoolbox {

struct HpParamResult {
    Eigen::VectorXcd z;
    std::complex<double> c;
    Eigen::MatrixXd qdat;
};

// Port of @hplmap/private/hpparam.m, restricted to the z0=[] (automatic
// initial guess) calling form, matching dparam's restriction.
HpParamResult hpparam(const Eigen::VectorXcd& w, const Eigen::VectorXd& beta, double tol = 1e-8, int method = 2);

}  // namespace sctoolbox
