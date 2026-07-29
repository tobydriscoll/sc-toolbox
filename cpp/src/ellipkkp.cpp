#include "sctoolbox/ellipkkp.hpp"

#include <cmath>
#include <limits>

namespace sctoolbox {

namespace {
double agmK(double a0, double b0) {
    double mm = 1.0;
    int i1 = 0;
    while (mm > std::numeric_limits<double>::epsilon()) {
        const double a1 = (a0 + b0) / 2.0;
        const double b1 = std::sqrt(a0 * b0);
        const double c1 = (a0 - b0) / 2.0;
        ++i1;
        mm = std::pow(2.0, i1) * c1 * c1;
        a0 = a1;
        b0 = b1;
    }
    return M_PI / (2.0 * a0);
}
}  // namespace

EllipKKp ellipkkp(double L) {
    EllipKKp r;
    if (L > 10.0) {
        r.K = M_PI / 2.0;
        r.Kp = M_PI * L + std::log(4.0);
        return r;
    }

    const double m = std::exp(-2.0 * M_PI * L);
    r.K = agmK(1.0, std::sqrt(1.0 - m));
    if (m == 1.0) r.K = std::numeric_limits<double>::infinity();

    r.Kp = agmK(1.0, std::sqrt(m));
    if (m == 0.0) r.Kp = std::numeric_limits<double>::infinity();

    return r;
}

}  // namespace sctoolbox
