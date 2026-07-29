#pragma once
#include <Eigen/Dense>

namespace sctoolbox {

// Mirrors @crdiskmap/private/crtriang.m's [edge,triedge,edgetri] triple
// directly: all fields hold 1-indexed vertex/edge/triangle ids, with 0
// as the "unset" sentinel, exactly as MATLAB's zero-initialized arrays.
// edge is 2 x (2n-3) (interior diagonals in columns 1..n-3, 1-indexed);
// triedge is 3 x (n-2); edgetri is 2 x (2n-3) (edgetri(2,k)==0 means edge
// k is a boundary edge, member of only one triangle).
struct CrTriangulation {
    Eigen::MatrixXi edge;
    Eigen::MatrixXi triedge;
    Eigen::MatrixXi edgetri;
};

// Port of @crdiskmap/private/crtriang.m: triangulates the simple polygon w.
CrTriangulation crtriang(const Eigen::VectorXcd& w);

}  // namespace sctoolbox
