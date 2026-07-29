#pragma once
#include <Eigen/Dense>
#include <complex>

namespace sctoolbox {

// Port of @rectmap/private/rmap.m.
Eigen::VectorXcd rmap(const Eigen::VectorXcd& zp, const Eigen::VectorXcd& w, const Eigen::VectorXd& beta,
                      const Eigen::VectorXcd& z, std::complex<double> c, double L, const Eigen::MatrixXd& qdat);

}  // namespace sctoolbox
