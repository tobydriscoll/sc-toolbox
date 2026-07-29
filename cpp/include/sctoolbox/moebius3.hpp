#pragma once
#include <Eigen/Dense>
#include <array>

namespace sctoolbox {

// Port of the @moebius/moebius.m (Z,W) constructor's "all finite" branch
// only: the Moebius transformation mapping the 3-vector z to w, as
// coefficients [c1,c2,c3,c4] for (c1 + c2*u) / (c3 + c4*u). crspread/
// crgather only ever call moebius with finite points (the rectangle-frame
// embedding corners), so the Inf-handling branches of moebius.m are not
// needed here.
std::array<std::complex<double>, 4> moebius3(const Eigen::Vector3cd& z, const Eigen::Vector3cd& w);

}  // namespace sctoolbox
