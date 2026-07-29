#include "sctoolbox/hplmap.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#include "sctoolbox/hpderiv.hpp"
#include "sctoolbox/hpinvmap.hpp"
#include "sctoolbox/hpmap.hpp"
#include "sctoolbox/hpparam.hpp"
#include "sctoolbox/hpquad.hpp"
#include "sctoolbox/scfix.hpp"

namespace sctoolbox {

namespace {
Polygon scfixedPolygon(const Polygon& poly) {
    const ScfixResult fixed = scfix("hp", poly.vertex(), poly.angle().array() - 1.0);
    Eigen::VectorXd alpha = fixed.beta.array() + 1.0;
    return Polygon(fixed.w, alpha);
}
bool isInf(std::complex<double> v) { return !std::isfinite(v.real()) || !std::isfinite(v.imag()); }
}  // namespace

HplMap::HplMap(Polygon poly, double tol) : poly_(scfixedPolygon(poly)), tol_(tol) {
    const Eigen::VectorXcd w = poly_.vertex();
    const Eigen::VectorXd beta = poly_.angle().array() - 1.0;
    const HpParamResult r = hpparam(w, beta, tol_);
    z_ = r.z;
    c_ = r.c;
    qdat_ = r.qdat;
}

Eigen::VectorXcd HplMap::eval(const Eigen::VectorXcd& zp) const {
    const Eigen::VectorXcd w = poly_.vertex();
    const Eigen::VectorXd beta = poly_.angle().array() - 1.0;
    const int p = static_cast<int>(zp.size());
    Eigen::VectorXcd wp = Eigen::VectorXcd::Constant(
        p, std::complex<double>(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()));

    std::vector<int> idx;
    for (int i = 0; i < p; ++i)
        if (zp(i).imag() > -std::numeric_limits<double>::epsilon()) idx.push_back(i);
    if (idx.empty()) return wp;

    Eigen::VectorXcd zpActive(idx.size());
    for (size_t i = 0; i < idx.size(); ++i) zpActive(i) = zp(idx[i]);
    const Eigen::VectorXcd wpActive = hpmap(zpActive, w, beta, z_, c_, qdat_);
    for (size_t i = 0; i < idx.size(); ++i) wp(idx[i]) = wpActive(i);
    return wp;
}

InvMapResult HplMap::evalinv(const Eigen::VectorXcd& wp, const Eigen::VectorXcd& z0) const {
    const Eigen::VectorXcd w = poly_.vertex();
    const Eigen::VectorXd beta = poly_.angle().array() - 1.0;
    return hpinvmap(wp, w, beta, z_, c_, qdat_, z0, {0.0, accuracy()});
}

Eigen::VectorXcd HplMap::evaldiff(const Eigen::VectorXcd& zp) const {
    const Eigen::VectorXd beta = poly_.angle().array() - 1.0;
    return hpderiv(zp, z_, beta, c_);
}

double HplMap::accuracy() const {
    if (acc_ >= 0.0) return acc_;

    const Eigen::VectorXcd w = poly_.vertex();
    const Eigen::VectorXd beta = poly_.angle().array() - 1.0;
    const int n = static_cast<int>(w.size());

    std::vector<int> idxAll;  // 1-indexed positions in 1..n-1 where w is finite
    for (int p = 1; p <= n - 1; ++p)
        if (!isInf(w(p - 1))) idxAll.push_back(p);

    const int m = static_cast<int>(idxAll.size()) - 1;
    Eigen::VectorXcd z1(m), z2(m), mid(m);
    std::vector<int> sing1(m), sing2(m);
    for (int k = 0; k < m; ++k) {
        const int a = idxAll[k], b = idxAll[k + 1];
        const std::complex<double> za = z_(a - 1), zb = z_(b - 1);
        z1(k) = za;
        z2(k) = zb;
        sing1[k] = a;
        sing2[k] = b;
        mid(k) = (za + zb) / 2.0 + std::complex<double>(0.0, std::abs(za - zb) / 2.0);
    }

    const Eigen::VectorXcd zFinite = z_.head(n - 1);
    const Eigen::VectorXd betaFinite = beta.head(n - 1);
    const Eigen::VectorXcd I = hpquad(z1, mid, sing1, zFinite, betaFinite, qdat_) -
                                hpquad(z2, mid, sing2, zFinite, betaFinite, qdat_);

    double maxAbs = 0.0;
    for (int k = 0; k < m; ++k) {
        const int a = idxAll[k], b = idxAll[k + 1];
        maxAbs = std::max(maxAbs, std::abs(c_ * I(k) - (w(b - 1) - w(a - 1))));
    }
    acc_ = maxAbs;
    return acc_;
}

}  // namespace sctoolbox
