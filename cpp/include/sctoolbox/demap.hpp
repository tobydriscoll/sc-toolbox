#pragma once
#include <Eigen/Dense>

namespace sctoolbox {

// Port of @extermap/private/demap.m. Evaluates the SC exterior map at points zp.
Eigen::VectorXcd demap(const Eigen::VectorXcd& zp, const Eigen::VectorXcd& w, const Eigen::VectorXd& beta,
                       const Eigen::VectorXcd& z, std::complex<double> c, const Eigen::MatrixXd& qdat);

}  // namespace sctoolbox
