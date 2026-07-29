#pragma once
#include <Eigen/Dense>
#include <complex>

namespace sctoolbox {

// Port of @diskmap/private/dderiv.m. f'(zp) = c * prod_i (1 - zp/z(i))^beta(i).
Eigen::VectorXcd dderiv(const Eigen::VectorXcd& zp, const Eigen::VectorXcd& z, const Eigen::VectorXd& beta,
                        std::complex<double> c = 1.0);

// Overload also returning the second derivative (dderiv.m's nargout>1 branch),
// used by dinvmap's Newton iteration.
Eigen::VectorXcd dderiv(const Eigen::VectorXcd& zp, const Eigen::VectorXcd& z, const Eigen::VectorXd& beta,
                        std::complex<double> c, Eigen::VectorXcd& d2f);

}  // namespace sctoolbox
