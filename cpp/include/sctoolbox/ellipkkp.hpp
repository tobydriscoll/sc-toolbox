#pragma once

namespace sctoolbox {

// Port of @rectmap/private/ellipkkp.m. Returns {K, Kp} for parameter
// m = exp(-2*pi*L), 0 < L < inf, via the AGM method.
struct EllipKKp {
    double K;
    double Kp;
};
EllipKKp ellipkkp(double L);

}  // namespace sctoolbox
