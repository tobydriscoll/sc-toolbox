#include "sctoolbox/diskmap.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#include "sctoolbox/dderiv.hpp"
#include "sctoolbox/dmap.hpp"
#include "sctoolbox/dparam.hpp"
#include "sctoolbox/dquad.hpp"
#include "sctoolbox/scfix.hpp"

namespace sctoolbox {

namespace {
Polygon scfixedPolygon(const Polygon& poly) {
    const ScfixResult fixed = scfix("d", poly.vertex(), poly.angle().array() - 1.0);
    Eigen::VectorXd alpha = fixed.beta.array() + 1.0;
    return Polygon(fixed.w, alpha);
}
}  // namespace

DiskMap::DiskMap(Polygon poly, double tol) : poly_(scfixedPolygon(poly)), tol_(tol) {
    const Eigen::VectorXcd w = poly_.vertex();
    const Eigen::VectorXd beta = poly_.angle().array() - 1.0;
    const DParamResult r = dparam(w, beta, tol_);
    z_ = r.z;
    c_ = r.c;
    qdat_ = r.qdat;
}

Eigen::VectorXcd DiskMap::eval(const Eigen::VectorXcd& zp) const {
    const Eigen::VectorXcd w = poly_.vertex();
    const Eigen::VectorXd beta = poly_.angle().array() - 1.0;
    const int p = static_cast<int>(zp.size());
    Eigen::VectorXcd wp = Eigen::VectorXcd::Constant(
        p, std::complex<double>(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()));

    std::vector<int> idx;
    for (int i = 0; i < p; ++i)
        if (std::abs(zp(i)) <= 1.0 + std::numeric_limits<double>::epsilon()) idx.push_back(i);
    if (idx.empty()) return wp;

    Eigen::VectorXcd zpActive(idx.size());
    for (size_t i = 0; i < idx.size(); ++i) zpActive(i) = zp(idx[i]);
    const Eigen::VectorXcd wpActive = dmap(zpActive, w, beta, z_, c_, qdat_);
    for (size_t i = 0; i < idx.size(); ++i) wp(idx[i]) = wpActive(i);
    return wp;
}

InvMapResult DiskMap::evalinv(const Eigen::VectorXcd& wp, const Eigen::VectorXcd& z0) const {
    const Eigen::VectorXcd w = poly_.vertex();
    const Eigen::VectorXd beta = poly_.angle().array() - 1.0;
    return dinvmap(wp, w, beta, z_, c_, qdat_, z0, {0.0, accuracy()});
}

Eigen::VectorXcd DiskMap::evaldiff(const Eigen::VectorXcd& zp) const {
    const Eigen::VectorXd beta = poly_.angle().array() - 1.0;
    return dderiv(zp, z_, beta, c_);
}

double DiskMap::accuracy() const {
    if (acc_ >= 0.0) return acc_;

    const Eigen::VectorXcd w = poly_.vertex();
    const Eigen::VectorXd beta = poly_.angle().array() - 1.0;
    const int n = static_cast<int>(w.size());

    std::vector<int> finite;
    for (int i = 0; i < n; ++i)
        if (!std::isinf(w(i).real()) && !std::isinf(w(i).imag())) finite.push_back(i);

    const int m = static_cast<int>(finite.size());
    Eigen::VectorXcd z1(m), z2(m), wf1(m), wf2(m);
    std::vector<int> sing1(m), sing2(m);
    for (int k = 0; k < m; ++k) {
        const int a = finite[k];
        const int b = finite[(k + 1) % m];
        z1(k) = z_(a);
        z2(k) = z_(b);
        sing1[k] = a + 1;
        sing2[k] = b + 1;
        wf1(k) = w(a);
        wf2(k) = w(b);
    }
    const Eigen::VectorXcd mid = Eigen::VectorXcd::Zero(m);

    const Eigen::VectorXcd I = dquad(z1, mid, sing1, z_, beta, qdat_) - dquad(z2, mid, sing2, z_, beta, qdat_);

    double maxAbs = 0.0;
    for (int k = 0; k < m; ++k) maxAbs = std::max(maxAbs, std::abs(c_ * I(k) - (wf2(k) - wf1(k))));
    acc_ = maxAbs;
    return acc_;
}

}  // namespace sctoolbox
