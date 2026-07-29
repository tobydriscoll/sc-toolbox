#pragma once
#include "sctoolbox/crtriang.hpp"

namespace sctoolbox {

// Port of @crdiskmap/private/crcdt.m: constrained Delaunay triangulation
// (edge-flipping) of a polygon given an initial triangulation.
CrTriangulation crcdt(const Eigen::VectorXcd& w, CrTriangulation t);

}  // namespace sctoolbox
