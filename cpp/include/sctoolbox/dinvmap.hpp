#pragma once
#include <Eigen/Dense>
#include <vector>

namespace sctoolbox {

struct InvMapResult {
    Eigen::VectorXcd zp;
    std::vector<int> flag;  // indices (0-based) where Newton failed to converge
};

// Port of @diskmap/private/dinvmap.m.
InvMapResult dinvmap(const Eigen::VectorXcd& wp, const Eigen::VectorXcd& w, const Eigen::VectorXd& beta,
                      const Eigen::VectorXcd& z, std::complex<double> c, const Eigen::MatrixXd& qdat,
                      const Eigen::VectorXcd& z0 = {}, const std::vector<double>& options = {});

}  // namespace sctoolbox
