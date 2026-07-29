#pragma once
#include <Eigen/Dense>

namespace sctoolbox {

// The "quadrilateral graph" built by @crdiskmap/private/crqgraph.m from a
// polygon triangulation. All fields are 0-indexed (MATLAB's edge/qlvert/
// qledge hold 1-indexed vertex/edge ids; adjacent is a logical matrix).
struct QGraph {
    Eigen::MatrixXi edge;      // 2 x (2n-3), 0-indexed vertex ids; diagonals (interior edges) in cols 0..n3-1
    Eigen::MatrixXi qlvert;    // 4 x n3, 0-indexed vertex ids
    Eigen::MatrixXi qledge;    // 4 x n3, 0-indexed edge ids (range 0..2n-4)
    Eigen::MatrixXi adjacent;  // n3 x n3, 0/1
};

}  // namespace sctoolbox
