#pragma once
#include <Eigen/Dense>

namespace sctoolbox {

// Port of @diskmap/private/dmap.m. Evaluates the SC disk map at points zp.
Eigen::VectorXcd dmap(const Eigen::VectorXcd& zp, const Eigen::VectorXcd& w, const Eigen::VectorXd& beta,
                      const Eigen::VectorXcd& z, std::complex<double> c, const Eigen::MatrixXd& qdat);

}  // namespace sctoolbox
