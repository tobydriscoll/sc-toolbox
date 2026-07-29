#include "sctoolbox/rmap.hpp"

#include <cmath>
#include <limits>

#include "sctoolbox/r2strip.hpp"
#include "sctoolbox/rcorners.hpp"
#include "sctoolbox/rectmap_internal.hpp"
#include "sctoolbox/stmap.hpp"

namespace sctoolbox {

Eigen::VectorXcd rmap(const Eigen::VectorXcd& zpIn, const Eigen::VectorXcd& wIn, const Eigen::VectorXd& betaIn,
                      const Eigen::VectorXcd& zIn, std::complex<double> c, double L, const Eigen::MatrixXd& qdat) {
    if (zpIn.size() == 0) return {};

    const int n = static_cast<int>(wIn.size());
    const RCorners rc = rcorners(wIn, betaIn, zIn);
    const Eigen::VectorXcd& w = rc.w;
    const Eigen::VectorXd& beta = rc.beta;
    const Eigen::VectorXcd& z = rc.z;

    const double Kp = z.array().imag().maxCoeff();

    Eigen::VectorXcd zs = r2strip(z, z, L).yp;
    for (int i = 0; i < n; ++i) zs(i) = std::complex<double>(zs(i).real(), std::round(zs(i).imag()));

    const auto [e0, e1] = rectStripEnds(zs);
    const double inf = std::numeric_limits<double>::infinity();
    const std::complex<double> nanc(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN());

    Eigen::VectorXcd zsAug(n + 2), wsAug(n + 2);
    Eigen::VectorXd bsAug(n + 2);
    int k = 0;
    for (int i = 0; i <= e0; ++i) {
        zsAug(k) = zs(i);
        wsAug(k) = w(i);
        bsAug(k) = beta(i);
        ++k;
    }
    zsAug(k) = std::complex<double>(inf, 0.0);
    wsAug(k) = nanc;
    bsAug(k) = 0.0;
    ++k;
    for (int i = e0 + 1; i <= e1; ++i) {
        zsAug(k) = zs(i);
        wsAug(k) = w(i);
        bsAug(k) = beta(i);
        ++k;
    }
    zsAug(k) = std::complex<double>(-inf, 0.0);
    wsAug(k) = nanc;
    bsAug(k) = 0.0;
    ++k;
    for (int i = e1 + 1; i < n; ++i) {
        zsAug(k) = zs(i);
        wsAug(k) = w(i);
        bsAug(k) = beta(i);
        ++k;
    }

    const Eigen::MatrixXd qdatAug = rectAugQdat(qdat, n, e0, e1);

    const int p = static_cast<int>(zpIn.size());
    Eigen::VectorXcd zp = zpIn;
    for (int i = 0; i < p; ++i) {
        if (std::abs(zp(i)) < 2.0 * std::numeric_limits<double>::epsilon()) {
            zp(i) += 100.0 * std::numeric_limits<double>::epsilon();
        }
        if (std::abs(zp(i) - std::complex<double>(0.0, Kp)) < 2.0 * std::numeric_limits<double>::epsilon()) {
            zp(i) -= std::complex<double>(0.0, 100.0 * std::numeric_limits<double>::epsilon() * Kp);
        }
    }

    const Eigen::VectorXcd yp = r2strip(zp, z, L).yp;
    return stmap(yp, wsAug, bsAug, zsAug, c, qdatAug);
}

}  // namespace sctoolbox
