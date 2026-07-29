#include "sctoolbox/rderiv.hpp"

#include <cmath>
#include <limits>

#include "sctoolbox/r2strip.hpp"
#include "sctoolbox/rectmap_internal.hpp"
#include "sctoolbox/stderiv.hpp"

namespace sctoolbox {

Eigen::VectorXcd rderiv(const Eigen::VectorXcd& zp, const Eigen::VectorXcd& z, const Eigen::VectorXd& beta,
                        std::complex<double> c, double L, const Eigen::VectorXcd& zsIn) {
    const int n = static_cast<int>(z.size());

    Eigen::VectorXcd zs = zsIn;
    if (zs.size() == 0) {
        zs = r2strip(z, z, L).yp;
        for (int i = 0; i < n; ++i) zs(i) = std::complex<double>(zs(i).real(), std::round(zs(i).imag()));
    }

    const R2Strip strip = r2strip(zp, z, L);
    const Eigen::VectorXcd& F = strip.yp;
    const Eigen::VectorXcd& dF = strip.yprime;

    const auto [e0, e1] = rectStripEnds(z);
    const double inf = std::numeric_limits<double>::infinity();

    Eigen::VectorXcd zsAug(n + 2);
    Eigen::VectorXd bsAug(n + 2);
    int k = 0;
    for (int i = 0; i <= e0; ++i) {
        zsAug(k) = zs(i);
        bsAug(k) = beta(i);
        ++k;
    }
    zsAug(k) = std::complex<double>(inf, 0.0);
    bsAug(k) = 0.0;
    ++k;
    for (int i = e0 + 1; i <= e1; ++i) {
        zsAug(k) = zs(i);
        bsAug(k) = beta(i);
        ++k;
    }
    zsAug(k) = std::complex<double>(-inf, 0.0);
    bsAug(k) = 0.0;
    ++k;
    for (int i = e1 + 1; i < n; ++i) {
        zsAug(k) = zs(i);
        bsAug(k) = beta(i);
        ++k;
    }

    const Eigen::VectorXcd dG = stderiv(F, zsAug, bsAug);

    const int p = static_cast<int>(zp.size());
    Eigen::VectorXcd fprime(p);
    for (int i = 0; i < p; ++i) fprime(i) = c * dF(i) * dG(i);
    return fprime;
}

}  // namespace sctoolbox
