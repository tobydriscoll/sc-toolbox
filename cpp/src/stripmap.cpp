#include "sctoolbox/stripmap.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#include "sctoolbox/scfix.hpp"
#include "sctoolbox/stderiv.hpp"
#include "sctoolbox/stinvmap.hpp"
#include "sctoolbox/stmap.hpp"
#include "sctoolbox/stparam.hpp"
#include "sctoolbox/stquad.hpp"
#include "sctoolbox/stquadh.hpp"

namespace sctoolbox {

namespace {
bool isInf(std::complex<double> v) { return !std::isfinite(v.real()) || !std::isfinite(v.imag()); }

Polygon scfixedPolygon(const Polygon& poly, std::array<int, 2>& ends) {
    std::vector<int> aux{ends[0], ends[1]};
    const ScfixResult fixed = scfix("st", poly.vertex(), poly.angle().array() - 1.0, aux);
    ends[0] = fixed.aux[0];
    ends[1] = fixed.aux[1];
    Eigen::VectorXd alpha = fixed.beta.array() + 1.0;
    return Polygon(fixed.w, alpha);
}
}  // namespace

StripMap::StripMap(Polygon poly, std::array<int, 2> ends, double tol)
    : poly_(scfixedPolygon(poly, ends)), tol_(tol) {
    const Eigen::VectorXcd w = poly_.vertex();
    const Eigen::VectorXd beta = poly_.angle().array() - 1.0;
    const StParamResult r = stparam(w, beta, ends, tol_);
    z_ = r.z;
    c_ = r.c;
    qdat_ = r.qdat;
}

Eigen::VectorXcd StripMap::eval(const Eigen::VectorXcd& zp) const {
    const Eigen::VectorXcd w = poly_.vertex();
    const Eigen::VectorXd beta = poly_.angle().array() - 1.0;
    const int p = static_cast<int>(zp.size());
    Eigen::VectorXcd wp = Eigen::VectorXcd::Constant(
        p, std::complex<double>(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()));

    std::vector<int> idx;
    for (int i = 0; i < p; ++i) {
        const double im = zp(i).imag();
        if (im > -std::numeric_limits<double>::epsilon() && im < 1.0 + std::numeric_limits<double>::epsilon())
            idx.push_back(i);
    }
    if (idx.empty()) return wp;

    Eigen::VectorXcd zpActive(idx.size());
    for (size_t i = 0; i < idx.size(); ++i) zpActive(i) = zp(idx[i]);
    const Eigen::VectorXcd wpActive = stmap(zpActive, w, beta, z_, c_, qdat_);
    for (size_t i = 0; i < idx.size(); ++i) wp(idx[i]) = wpActive(i);
    return wp;
}

InvMapResult StripMap::evalinv(const Eigen::VectorXcd& wp, const Eigen::VectorXcd& z0) const {
    const Eigen::VectorXcd w = poly_.vertex();
    const Eigen::VectorXd beta = poly_.angle().array() - 1.0;
    return stinvmap(wp, w, beta, z_, c_, qdat_, z0, {0.0, accuracy()});
}

Eigen::VectorXcd StripMap::evaldiff(const Eigen::VectorXcd& zp) const {
    const Eigen::VectorXd beta = poly_.angle().array() - 1.0;
    return stderiv(zp, z_, beta, c_);
}

double StripMap::accuracy() const {
    if (acc_ >= 0.0) return acc_;

    const Eigen::VectorXcd w = poly_.vertex();
    const Eigen::VectorXd beta = poly_.angle().array() - 1.0;
    const int n = static_cast<int>(w.size());

    int end1 = -1, end2 = -1;  // 0-indexed
    for (int i = 0; i < n; ++i) {
        if (isInf(z_(i)) && z_(i).real() < 0) end1 = i;
        if (isInf(z_(i)) && z_(i).real() > 0) end2 = i;
    }

    std::vector<bool> bot(n, false), top(n, false);
    for (int i = 0; i < n; ++i) {
        if (isInf(z_(i))) continue;
        if (z_(i).imag() == 0.0) bot[i] = true;
        else top[i] = true;
    }
    bot[end1] = bot[end2] = false;
    top[end1] = top[end2] = false;

    std::vector<int> idxbot, idxtop;  // 0-indexed
    for (int i = 0; i < n; ++i)
        if (bot[i] && !isInf(w(i))) idxbot.push_back(i);
    for (int i = 0; i < n; ++i)
        if (top[i] && !isInf(w(i))) idxtop.push_back(i);

    std::vector<std::pair<int, int>> pairs;
    for (size_t k = 0; k + 1 < idxbot.size(); ++k) pairs.emplace_back(idxbot[k], idxbot[k + 1]);
    for (size_t k = 0; k + 1 < idxtop.size(); ++k) pairs.emplace_back(idxtop[k], idxtop[k + 1]);

    // Cross-strip check: vertex right after end1, paired with the idxtop
    // entry whose real part is closest.
    const int afterEnd1 = (end1 + 1) % n;
    int bestK = 0;
    double bestD = std::numeric_limits<double>::infinity();
    for (size_t k = 0; k < idxtop.size(); ++k) {
        const double d = std::abs(z_(idxtop[k]).real() - z_(afterEnd1).real());
        if (d < bestD) {
            bestD = d;
            bestK = static_cast<int>(k);
        }
    }
    pairs.emplace_back(afterEnd1, idxtop[bestK]);

    const int m = static_cast<int>(pairs.size());
    Eigen::VectorXcd I(m);
    for (int p = 0; p < m; ++p) {
        const int a = pairs[p].first, b = pairs[p].second;
        const std::complex<double> zl = z_(a), zr = z_(b);
        const bool s2 = (b - a == 1);
        if (s2) {
            Eigen::VectorXcd zlv(1), zrv(1), midv(1);
            std::vector<int> singL{a + 1}, singR{b + 1};
            zlv(0) = zl;
            zrv(0) = zr;
            midv(0) = (zl + zr) / 2.0;
            const Eigen::VectorXcd I1 = stquadh(zlv, midv, singL, z_, beta, qdat_);
            const Eigen::VectorXcd I2 = stquadh(zrv, midv, singR, z_, beta, qdat_);
            I(p) = I1(0) - I2(0);
        } else {
            Eigen::VectorXcd zlv(1), zrv(1), mid1v(1), mid2v(1);
            std::vector<int> singL{a + 1}, singR{b + 1}, zeroSing{0};
            zlv(0) = zl;
            zrv(0) = zr;
            mid1v(0) = std::complex<double>(zl.real(), 0.5);
            mid2v(0) = std::complex<double>(zr.real(), 0.5);
            const Eigen::VectorXcd I1 = stquad(zlv, mid1v, singL, z_, beta, qdat_);
            const Eigen::VectorXcd I2 = stquadh(mid1v, mid2v, zeroSing, z_, beta, qdat_);
            const Eigen::VectorXcd I3 = stquad(zrv, mid2v, singR, z_, beta, qdat_);
            I(p) = I1(0) + I2(0) - I3(0);
        }
    }

    double maxAbs = 0.0;
    for (int p = 0; p < m; ++p) {
        const int a = pairs[p].first, b = pairs[p].second;
        maxAbs = std::max(maxAbs, std::abs(c_ * I(p) - (w(b) - w(a))));
    }
    acc_ = maxAbs;
    return acc_;
}

}  // namespace sctoolbox
