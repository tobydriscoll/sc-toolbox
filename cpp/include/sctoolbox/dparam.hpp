#pragma once
#include <Eigen/Dense>

namespace sctoolbox {

struct DParamResult {
    Eigen::VectorXcd z;
    std::complex<double> c;
    Eigen::MatrixXd qdat;
};

// Port of @diskmap/private/dparam.m, restricted to the z0=[] (automatic
// initial guess) calling form -- the only form exercised by golden tests
// and by dinvmap/etc. callers in this codebase. `tol` and `method` mirror
// scmapopt's Tolerance/SolverMethod (method: 1=line search, 2=trust
// region/hook step, matching nesolve's details(2) convention).
DParamResult dparam(const Eigen::VectorXcd& w, const Eigen::VectorXd& beta, double tol = 1e-8, int method = 2);

}  // namespace sctoolbox
