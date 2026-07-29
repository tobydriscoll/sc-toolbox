#pragma once
#include <Eigen/Dense>
#include <string>
#include <vector>

namespace sctoolbox {

struct ScfixResult {
    Eigen::VectorXcd w;
    Eigen::VectorXd beta;
    std::vector<int> aux;
};

// Port of scfix.m. `type` is one of "hp", "d", "de", "st", "r". `aux`
// (1-indexed vertex positions, matching MATLAB's convention) holds the
// strip-end pair for "st" or the four corner positions for "r"; required
// for those two types since the interactive (mouse-click) fallback used by
// the original MATLAB for unsupplied strip ends has no headless equivalent
// here.
ScfixResult scfix(const std::string& type, Eigen::VectorXcd w, Eigen::VectorXd beta,
                   std::vector<int> aux = {});

}  // namespace sctoolbox
