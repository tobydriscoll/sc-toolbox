#pragma once
#include <Eigen/Dense>

namespace sctoolbox {

// Port of @hplmap/private/hpmap.m. Evaluates the SC half-plane map at points zp.
Eigen::VectorXcd hpmap(const Eigen::VectorXcd& zp, const Eigen::VectorXcd& w, const Eigen::VectorXd& beta,
                       const Eigen::VectorXcd& z, std::complex<double> c, const Eigen::MatrixXd& qdat);

}  // namespace sctoolbox
