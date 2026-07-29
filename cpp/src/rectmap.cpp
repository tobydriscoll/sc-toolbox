#include "sctoolbox/rectmap.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#include "sctoolbox/isinpoly.hpp"
#include "sctoolbox/r2strip.hpp"
#include "sctoolbox/rcorners.hpp"
#include "sctoolbox/rderiv.hpp"
#include "sctoolbox/rectmap_internal.hpp"
#include "sctoolbox/rinvmap.hpp"
#include "sctoolbox/rmap.hpp"
#include "sctoolbox/rparam.hpp"
#include "sctoolbox/scfix.hpp"
#include "sctoolbox/stquad.hpp"

namespace sctoolbox {

namespace {
bool isInf(std::complex<double> v) { return !std::isfinite(v.real()) || !std::isfinite(v.imag()); }

Polygon scfixedPolygon(const Polygon& poly, std::array<int, 4>& cnr) {
    std::vector<int> aux{cnr[0], cnr[1], cnr[2], cnr[3]};
    const ScfixResult fixed = scfix("r", poly.vertex(), poly.angle().array() - 1.0, aux);
    for (int i = 0; i < 4; ++i) cnr[i] = fixed.aux[i];
    Eigen::VectorXd alpha = fixed.beta.array() + 1.0;
    return Polygon(fixed.w, alpha);
}
}  // namespace

RectMap::RectMap(Polygon poly, std::array<int, 4> cnr, double tol) : poly_(scfixedPolygon(poly, cnr)), tol_(tol) {
    const Eigen::VectorXcd w = poly_.vertex();
    const Eigen::VectorXd beta = poly_.angle().array() - 1.0;
    const RParamResult r = rparam(w, beta, cnr, tol_);
    z_ = r.z;
    c_ = r.c;
    L_ = r.L;
    qdat_ = r.qdat;
}

std::array<int, 4> RectMap::corners() const {
    const int n = static_cast<int>(z_.size());
    const double K = z_.real().maxCoeff();
    const double Kp = z_.array().imag().maxCoeff();
    const std::array<std::complex<double>, 4> rect = {std::complex<double>(K, 0.0), std::complex<double>(K, Kp),
                                                       std::complex<double>(-K, Kp), std::complex<double>(-K, 0.0)};
    std::array<int, 4> corner;
    for (int c = 0; c < 4; ++c) {
        double best = std::numeric_limits<double>::infinity();
        int bestI = 0;
        for (int i = 0; i < n; ++i) {
            const double d = std::abs(z_(i) - rect[c]);
            if (d < best) {
                best = d;
                bestI = i;
            }
        }
        corner[c] = bestI;
    }
    return corner;
}

Eigen::VectorXcd RectMap::eval(const Eigen::VectorXcd& zp) const {
    const Eigen::VectorXcd w = poly_.vertex();
    const Eigen::VectorXd beta = poly_.angle().array() - 1.0;
    const std::array<int, 4> cnr = corners();
    Eigen::VectorXcd zr(4);
    for (int i = 0; i < 4; ++i) zr(i) = z_(cnr[i]);

    const int p = static_cast<int>(zp.size());
    Eigen::VectorXcd wp = Eigen::VectorXcd::Constant(
        p, std::complex<double>(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()));

    const Eigen::VectorXd inPoly = isinpoly(zp, zr, tol_);
    std::vector<int> idx;
    for (int i = 0; i < p; ++i)
        if (inPoly(i) != 0.0) idx.push_back(i);
    if (idx.empty()) return wp;

    Eigen::VectorXcd zpActive(idx.size());
    for (size_t i = 0; i < idx.size(); ++i) zpActive(i) = zp(idx[i]);
    const Eigen::VectorXcd wpActive = rmap(zpActive, w, beta, z_, c_, L_, qdat_);
    for (size_t i = 0; i < idx.size(); ++i) wp(idx[i]) = wpActive(i);
    return wp;
}

InvMapResult RectMap::evalinv(const Eigen::VectorXcd& wp, const Eigen::VectorXcd& z0) const {
    const Eigen::VectorXcd w = poly_.vertex();
    const Eigen::VectorXd beta = poly_.angle().array() - 1.0;
    return rinvmap(wp, w, beta, z_, c_, L_, qdat_, z0, {0.0, accuracy()});
}

Eigen::VectorXcd RectMap::evaldiff(const Eigen::VectorXcd& zp) const {
    const Eigen::VectorXd beta = poly_.angle().array() - 1.0;
    return rderiv(zp, z_, beta, c_, L_);
}

double RectMap::accuracy() const {
    if (acc_ >= 0.0) return acc_;

    const Eigen::VectorXcd wOrig = poly_.vertex();
    const Eigen::VectorXd betaOrig = poly_.angle().array() - 1.0;
    const RCorners rc = rcorners(wOrig, betaOrig, z_);
    const Eigen::VectorXcd& w = rc.w;
    const Eigen::VectorXd& beta = rc.beta;
    const Eigen::VectorXcd& z = rc.z;
    const int n = static_cast<int>(w.size());
    const int corner3 = rc.corners[2];

    Eigen::VectorXcd zs = r2strip(z, z, L_).yp;
    for (int i = 0; i < n; ++i) zs(i) = std::complex<double>(zs(i).real(), std::round(zs(i).imag()));

    std::vector<int> idxbot, idxtop;
    for (int i = 0; i < corner3; ++i)
        if (!isInf(w(i))) idxbot.push_back(i);
    for (int i = corner3; i < n; ++i)
        if (!isInf(w(i))) idxtop.push_back(i);

    std::vector<std::pair<int, int>> pairs;
    std::vector<std::complex<double>> mids;
    for (size_t k = 0; k + 1 < idxbot.size(); ++k) {
        const int a = idxbot[k], b = idxbot[k + 1];
        pairs.emplace_back(a, b);
        mids.emplace_back(((zs(a) + zs(b)) / 2.0).real(), 0.5);
    }
    for (size_t k = 0; k + 1 < idxtop.size(); ++k) {
        const int a = idxtop[k], b = idxtop[k + 1];
        pairs.emplace_back(a, b);
        mids.emplace_back(((zs(a) + zs(b)) / 2.0).real(), 0.5);
    }

    int bestK = 0;
    double bestD = std::numeric_limits<double>::infinity();
    for (size_t k = 0; k < idxtop.size(); ++k) {
        const double d = std::abs(zs(idxtop[k]).real() - zs(0).real());
        if (d < bestD) {
            bestD = d;
            bestK = static_cast<int>(k);
        }
    }
    pairs.emplace_back(0, idxtop[bestK]);
    mids.push_back((zs(0) + zs(idxtop[bestK])) / 2.0);

    const auto [e0, e1] = rectStripEnds(zs);
    const double inf = std::numeric_limits<double>::infinity();
    Eigen::VectorXcd zq(n + 2);
    Eigen::VectorXd bq(n + 2);
    {
        int k = 0;
        for (int i = 0; i <= e0; ++i) {
            zq(k) = zs(i);
            bq(k) = beta(i);
            ++k;
        }
        zq(k) = std::complex<double>(inf, 0.0);
        bq(k) = 0.0;
        ++k;
        for (int i = e0 + 1; i <= e1; ++i) {
            zq(k) = zs(i);
            bq(k) = beta(i);
            ++k;
        }
        zq(k) = std::complex<double>(-inf, 0.0);
        bq(k) = 0.0;
        ++k;
        for (int i = e1 + 1; i < n; ++i) {
            zq(k) = zs(i);
            bq(k) = beta(i);
            ++k;
        }
    }
    const Eigen::MatrixXd qdataAug = rectAugQdat(qdat_, n, e0, e1);

    const int m = static_cast<int>(pairs.size());
    Eigen::VectorXcd I(m);
    for (int p = 0; p < m; ++p) {
        const int a = pairs[p].first, b = pairs[p].second;
        const int aShift = a + (a > e0 ? 1 : 0) + (a > e1 ? 1 : 0);
        const int bShift = b + (b > e0 ? 1 : 0) + (b > e1 ? 1 : 0);
        Eigen::VectorXcd zlv(1), zrv(1), midv(1);
        std::vector<int> singL{aShift + 1}, singR{bShift + 1};
        zlv(0) = zs(a);
        zrv(0) = zs(b);
        midv(0) = mids[p];
        const Eigen::VectorXcd I1 = stquad(zlv, midv, singL, zq, bq, qdataAug);
        const Eigen::VectorXcd I2 = stquad(zrv, midv, singR, zq, bq, qdataAug);
        I(p) = I1(0) - I2(0);
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
