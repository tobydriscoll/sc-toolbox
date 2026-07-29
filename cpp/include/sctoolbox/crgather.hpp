#pragma once
#include <Eigen/Dense>
#include <vector>

#include "sctoolbox/qgraph.hpp"

namespace sctoolbox {

// Port of @crdiskmap/private/crgather.m. Converts points in u (currently
// expressed in the per-point embeddings given by uquad, 0-indexed) into
// the single embedding quadnum (0-indexed), mutating both u and uquad
// in place.
void crgather(Eigen::VectorXcd& u, std::vector<int>& uquad, int quadnum, const Eigen::VectorXd& cr,
              const QGraph& Q);

}  // namespace sctoolbox
