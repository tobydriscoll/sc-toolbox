#pragma once
#include <Eigen/Dense>
#include <complex>

namespace sctoolbox {

// Port of @extermap/private/dederiv.m. f'(zp) = c * prod_i (1 - zp/z(i))^beta(i) * zp^(-2),
// where the zp^(-2) factor accounts for the implicit singularity at the origin.
Eigen::VectorXcd dederiv(const Eigen::VectorXcd& zp, const Eigen::VectorXcd& z, const Eigen::VectorXd& beta,
                         std::complex<double> c = 1.0);

}  // namespace sctoolbox
