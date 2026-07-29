#pragma once
#include <Eigen/Dense>
#include <complex>

namespace sctoolbox {

// Port of @hplmap/private/hpderiv.m. f'(zp) = c * prod_i (zp - z(i))^beta(i),
// over the finite prevertices only (infinite entries of z are stripped).
Eigen::VectorXcd hpderiv(const Eigen::VectorXcd& zp, const Eigen::VectorXcd& z, const Eigen::VectorXd& beta,
                         std::complex<double> c = 1.0);

}  // namespace sctoolbox
