#pragma once
#include <Eigen/Dense>

namespace sctoolbox {

// Port of @crdiskmap/private/crmap0.m. aff is the 2-element complex affine
// table row (scale, translation) for this embedding.
Eigen::VectorXcd crmap0(const Eigen::VectorXcd& zp, const Eigen::VectorXcd& z, const Eigen::VectorXd& beta,
                        const Eigen::Vector2cd& aff, const Eigen::MatrixXd& qdat);

}  // namespace sctoolbox
