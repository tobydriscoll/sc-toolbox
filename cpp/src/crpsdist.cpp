#include "sctoolbox/crpsdist.hpp"

#include <cmath>

namespace sctoolbox {

Eigen::VectorXd crpsdist(const std::array<std::complex<double>, 2>& segment, const Eigen::VectorXcd& pts) {
    const int n = static_cast<int>(pts.size());
    Eigen::VectorXd d(n);
    if (n == 0) return d;

    const std::complex<double> diffSeg = segment[1] - segment[0];
    const std::complex<double> rot = diffSeg / std::abs(diffSeg);  // complex sign(diffSeg)
    const double xmax = std::abs(diffSeg);

    for (int i = 0; i < n; ++i) {
        const std::complex<double> p = (pts(i) - segment[0]) / rot;
        const double re = p.real();
        if (re >= 0.0 && re <= xmax) {
            d(i) = std::abs(p.imag());
        } else if (re < 0.0) {
            d(i) = std::abs(p);
        } else {
            d(i) = std::abs(p - xmax);
        }
    }
    return d;
}

}  // namespace sctoolbox
