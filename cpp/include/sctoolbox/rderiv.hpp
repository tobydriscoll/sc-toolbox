#pragma once
#include <Eigen/Dense>
#include <complex>

namespace sctoolbox {

// Port of @rectmap/private/rderiv.m. If zs is empty, it is computed
// internally as r2strip(z,z,L) snapped to the strip edges (as MATLAB does
// when called with 5 args).
Eigen::VectorXcd rderiv(const Eigen::VectorXcd& zp, const Eigen::VectorXcd& z, const Eigen::VectorXd& beta,
                        std::complex<double> c, double L, const Eigen::VectorXcd& zs = {});

}  // namespace sctoolbox
