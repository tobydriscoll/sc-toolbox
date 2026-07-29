#pragma once
#include <Eigen/Dense>
#include <vector>

namespace sctoolbox {

struct RCorners {
    Eigen::VectorXcd w;
    Eigen::VectorXd beta;
    Eigen::VectorXcd z;
    std::vector<int> corners;  // 0-indexed, into the renumbered w/beta/z above
};

// Port of @rectmap/private/rcorners.m.
RCorners rcorners(const Eigen::VectorXcd& w, const Eigen::VectorXd& beta, const Eigen::VectorXcd& z);

}  // namespace sctoolbox
