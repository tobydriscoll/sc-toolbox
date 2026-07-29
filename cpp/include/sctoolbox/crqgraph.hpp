#pragma once
#include "sctoolbox/crtriang.hpp"
#include "sctoolbox/qgraph.hpp"

namespace sctoolbox {

// Port of @crdiskmap/private/crqgraph.m's 4-argument form (triangulation
// already computed). Converts the 1-indexed CrTriangulation into the
// 0-indexed QGraph used by crembed/crspread/crgather/crderiv/crmap/crinvmap.
QGraph crqgraph(const Eigen::VectorXcd& w, const CrTriangulation& tri);

}  // namespace sctoolbox
