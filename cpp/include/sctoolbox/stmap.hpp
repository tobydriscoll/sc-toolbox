#pragma once
#include <Eigen/Dense>

namespace sctoolbox {

// Port of @stripmap/private/stmap.m. Evaluates the SC strip map at points zp.
Eigen::VectorXcd stmap(const Eigen::VectorXcd& zp, const Eigen::VectorXcd& w, const Eigen::VectorXd& beta,
                       const Eigen::VectorXcd& z, std::complex<double> c, const Eigen::MatrixXd& qdat);

}  // namespace sctoolbox
