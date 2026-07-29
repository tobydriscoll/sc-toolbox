#pragma once
#include <vector>

#include "sctoolbox/crtriang.hpp"

namespace sctoolbox {

struct CrSplitResult {
    Eigen::VectorXcd w;
    std::vector<bool> orig;  // true at entries corresponding to an original vertex
    CrTriangulation tri;     // final constrained Delaunay triangulation
};

// Port of @crdiskmap/private/crsplit.m. Phase 1 (chopping very sharp
// corners) and phase 2's split-detection (geodesic-distance check for
// narrow channels) are fully implemented. If phase 2 actually finds a
// channel needing a split, the dynamic re-triangulation surgery that
// follows (@crdiskmap/private/crsplit.m lines ~103-143: inserting new
// edges/triangles into the CDT in place) is not implemented -- no existing
// golden test exercises it (both crdiskmap test polygons need zero splits,
// confirmed directly against MATLAB), and guessing at untested mesh-surgery
// logic is worse than failing loudly. Throws std::runtime_error in that case.
CrSplitResult crsplit(const Eigen::VectorXcd& w);

}  // namespace sctoolbox
