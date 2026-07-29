#include "sctoolbox/extermap.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#include "sctoolbox/dederiv.hpp"
#include "sctoolbox/deinvmap.hpp"
#include "sctoolbox/demap.hpp"
#include "sctoolbox/dequad.hpp"
#include "sctoolbox/deparam.hpp"
#include "sctoolbox/scfix.hpp"

namespace sctoolbox {

namespace {
bool isInf(std::complex<double> v) { return !std::isfinite(v.real()) || !std::isfinite(v.imag()); }

// w=flipud(vertex(poly)); beta=flipud(1-angle(poly))
std::pair<Eigen::VectorXcd, Eigen::VectorXd> deConvention(const Polygon& poly) {
    return {poly.vertex().reverse(), 1.0 - poly.angle().reverse().array()};
}

Polygon scfixedExterPolygon(const Polygon& poly) {
    auto [w, beta] = deConvention(poly);
    const ScfixResult fixed = scfix("de", w, beta);
    const Eigen::VectorXcd polyW = fixed.w.reverse();
    const Eigen::VectorXd polyAlpha = 1.0 - fixed.beta.reverse().array();
    return Polygon(polyW, polyAlpha);
}
}  // namespace

ExterMap::ExterMap(Polygon poly, double tol) : poly_(scfixedExterPolygon(poly)), tol_(tol) {
    auto [w, beta] = deConvention(poly_);
    const DeParamResult r = deparam(w, beta, tol_);
    z_ = r.z;
    c_ = r.c;
    qdat_ = r.qdat;
}

Eigen::VectorXcd ExterMap::eval(const Eigen::VectorXcd& zp) const {
    auto [w, beta] = deConvention(poly_);
    const int p = static_cast<int>(zp.size());
    Eigen::VectorXcd wp = Eigen::VectorXcd::Constant(
        p, std::complex<double>(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()));

    std::vector<int> idx;
    for (int i = 0; i < p; ++i)
        if (std::abs(zp(i)) <= 1.0 + std::numeric_limits<double>::epsilon()) idx.push_back(i);
    if (idx.empty()) return wp;

    Eigen::VectorXcd zpActive(idx.size());
    for (size_t i = 0; i < idx.size(); ++i) zpActive(i) = zp(idx[i]);
    const Eigen::VectorXcd wpActive = demap(zpActive, w, beta, z_, c_, qdat_);
    for (size_t i = 0; i < idx.size(); ++i) wp(idx[i]) = wpActive(i);
    return wp;
}

InvMapResult ExterMap::evalinv(const Eigen::VectorXcd& wp, const Eigen::VectorXcd& z0) const {
    auto [w, beta] = deConvention(poly_);
    return deinvmap(wp, w, beta, z_, c_, qdat_, z0, {0.0, accuracy()});
}

Eigen::VectorXcd ExterMap::evaldiff(const Eigen::VectorXcd& zp) const {
    auto [w, beta] = deConvention(poly_);
    (void)w;
    return dederiv(zp, z_, beta, c_);
}

double ExterMap::accuracy() const {
    if (acc_ >= 0.0) return acc_;

    auto [w, beta] = deConvention(poly_);
    const int n = static_cast<int>(w.size());

    std::vector<int> idxAll;  // 1-indexed finite positions
    for (int p = 1; p <= n; ++p)
        if (!isInf(w(p - 1))) idxAll.push_back(p);
    const int m = static_cast<int>(idxAll.size());

    Eigen::VectorXcd z1(m), z2(m), mid(m);
    std::vector<int> sing1(m), sing2(m);
    for (int k = 0; k < m; ++k) {
        const int a = idxAll[k], b = idxAll[(k + 1) % m];
        const std::complex<double> za = z_(a - 1), zb = z_(b - 1);
        z1(k) = za;
        z2(k) = zb;
        sing1[k] = a;
        sing2[k] = b;
        const double dtheta = std::fmod(std::arg(zb / za) + 2.0 * M_PI, 2.0 * M_PI);
        mid(k) = za * std::exp(std::complex<double>(0.0, dtheta / 2.0));
    }

    const Eigen::VectorXcd I = dequad(z1, mid, sing1, z_, beta, qdat_) - dequad(z2, mid, sing2, z_, beta, qdat_);

    double maxAbs = 0.0;
    for (int k = 0; k < m; ++k) {
        const int a = idxAll[k], b = idxAll[(k + 1) % m];
        maxAbs = std::max(maxAbs, std::abs(c_ * I(k) - (w(b - 1) - w(a - 1))));
    }
    acc_ = maxAbs;
    return acc_;
}

}  // namespace sctoolbox
