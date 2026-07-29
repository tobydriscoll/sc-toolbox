#include "sctoolbox/r2strip.hpp"

#include <cmath>

#include "sctoolbox/ellipjc.hpp"

namespace sctoolbox {

R2Strip r2strip(const Eigen::VectorXcd& zp, const Eigen::VectorXcd& /*z*/, double L) {
    const int n = static_cast<int>(zp.size());
    const EllipJC e = ellipjc(zp, L);

    R2Strip r;
    r.yp.resize(n);
    r.yprime.resize(n);
    for (int i = 0; i < n; ++i) {
        std::complex<double> sn = e.sn(i);
        sn = std::complex<double>(sn.real(), std::max(sn.imag(), 0.0));

        std::complex<double> yp = std::log(sn) / M_PI;
        yp = std::complex<double>(yp.real(), std::max(0.0, yp.imag()));
        yp = std::complex<double>(yp.real(), std::min(1.0, yp.imag()));
        r.yp(i) = yp;

        r.yprime(i) = e.cn(i) * e.dn(i) / sn / M_PI;
    }
    return r;
}

}  // namespace sctoolbox
