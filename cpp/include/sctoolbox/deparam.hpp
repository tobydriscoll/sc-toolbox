#pragma once
#include <Eigen/Dense>

namespace sctoolbox {

struct DeParamResult {
    Eigen::VectorXcd z;
    std::complex<double> c;
    Eigen::MatrixXd qdat;
};

// Port of @extermap/private/deparam.m, restricted to the z0=[] (automatic
// initial guess) calling form, matching dparam/hpparam's restriction.
DeParamResult deparam(const Eigen::VectorXcd& w, const Eigen::VectorXd& beta, double tol = 1e-8, int method = 2);

}  // namespace sctoolbox
