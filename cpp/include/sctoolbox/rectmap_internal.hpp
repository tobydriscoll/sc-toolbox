#pragma once
#include <Eigen/Dense>
#include <utility>

namespace sctoolbox {

// 0-indexed positions i such that imag(z[i]) != imag(z[(i+1)%n]), matching
// MATLAB's `find(diff(imag(z([1:n 1]))))` (1-indexed) minus 1. Used by both
// rderiv and rmap to locate where to splice in the strip's +/-Inf ends.
std::pair<int, int> rectStripEnds(const Eigen::VectorXcd& z);

// Builds the (n+3)-column-pair qdat used by stripmap.evaluate/deriv from the
// original (n+1)-column-pair qdat returned by scqdata(beta,...), splicing in
// scqdata's built-in "neutral" (beta=0) column at the two inserted +/-Inf
// slots in the augmented n+2-prevertex array. The result needs n+2+1 = n+3
// column-pairs (not n+2) because stmap/stquad index qdat assuming the same
// "(prevertex count)+1" neutral-column convention scqdata itself uses, where
// "prevertex count" here is the augmented count n+2. Matches MATLAB's
// `idx = [1:ends(1) n+1 ends(1)+1:ends(2) n+1 ends(2)+1:n n+1]` in
// rmap.m/rderiv.m, whose trailing "n+1" filler is not actually redundant.
Eigen::MatrixXd rectAugQdat(const Eigen::MatrixXd& qdat, int n, int e0, int e1);

}  // namespace sctoolbox
