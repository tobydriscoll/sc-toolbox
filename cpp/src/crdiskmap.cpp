#include "sctoolbox/crdiskmap.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "sctoolbox/crderiv.hpp"
#include "sctoolbox/crembed.hpp"
#include "sctoolbox/crfixwc.hpp"
#include "sctoolbox/crinvmap.hpp"
#include "sctoolbox/crmap.hpp"
#include "sctoolbox/crossrat.hpp"
#include "sctoolbox/crparam.hpp"
#include "sctoolbox/crquad.hpp"
#include "sctoolbox/isinpoly.hpp"
#include "sctoolbox/scfix.hpp"

namespace sctoolbox {

CrDiskMap::CrDiskMap(Polygon poly, double tol) : poly_(poly), tol_(tol) {
    const ScfixResult fixed = scfix("d", poly.vertex(), poly.angle().array() - 1.0);
    const CrParamResult r = crparam(fixed.w, fixed.beta, tol_);
    poly_ = Polygon(r.w, r.beta.array() + 1.0);
    cr_ = r.cr;
    aff_ = r.aff;
    Q_ = r.Q;
    qdat_ = r.qdat;

    const std::complex<double> wc = (r.w(Q_.qlvert(0, 0)) + r.w(Q_.qlvert(1, 0)) + r.w(Q_.qlvert(2, 0))) / 3.0;
    wcfix_ = crfixwc(r.w, r.beta, r.cr, r.aff, Q_, wc);
}

Eigen::VectorXcd CrDiskMap::eval(const Eigen::VectorXcd& zp) const {
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
    const Eigen::VectorXcd wpActive = crmap(zpActive, w, beta, cr_, aff_, wcfix_, Q_, qdat_);
    for (size_t i = 0; i < idx.size(); ++i) wp(idx[i]) = wpActive(i);
    return wp;
}

Eigen::VectorXcd CrDiskMap::evalinv(const Eigen::VectorXcd& wp) const {
    const Eigen::VectorXcd w = poly_.vertex();
    const Eigen::VectorXd beta = poly_.angle().array() - 1.0;
    const int p = static_cast<int>(wp.size());
    Eigen::VectorXcd zp = Eigen::VectorXcd::Constant(
        p, std::complex<double>(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()));

    const Eigen::VectorXd inPoly = isinpoly(wp, w);
    std::vector<int> idx;
    for (int i = 0; i < p; ++i)
        if (inPoly(i) != 0.0) idx.push_back(i);
    if (idx.empty()) return zp;

    Eigen::VectorXcd wpActive(idx.size());
    for (size_t i = 0; i < idx.size(); ++i) wpActive(i) = wp(idx[i]);
    const Eigen::VectorXcd zpActive = crinvmap(wpActive, w, beta, cr_, aff_, wcfix_, Q_, qdat_, {0.0, accuracy()});
    for (size_t i = 0; i < idx.size(); ++i) zp(idx[i]) = zpActive(i);
    return zp;
}

Eigen::VectorXcd CrDiskMap::evaldiff(const Eigen::VectorXcd& zp) const {
    const Eigen::VectorXd beta = poly_.angle().array() - 1.0;
    return crderiv(zp, beta, cr_, aff_, wcfix_, Q_);
}

double CrDiskMap::accuracy() const {
    if (acc_ >= 0.0) return acc_;

    const Eigen::VectorXcd w = poly_.vertex();
    const Eigen::VectorXd beta = poly_.angle().array() - 1.0;
    const int n = static_cast<int>(w.size());
    const int n3 = n - 3;

    const Eigen::VectorXcd crtarget = crossrat(w, Q_);
    Eigen::VectorXcd crimage(n3);

    for (int k = 0; k < n3; ++k) {
        const Eigen::VectorXcd prever = crembed(cr_, Q_, k);
        const Eigen::Vector4i idx = Q_.qlvert.col(k);
        Eigen::VectorXcd z4(4);
        std::vector<int> sing4(4);
        for (int i = 0; i < 4; ++i) {
            z4(i) = prever(idx(i));
            sing4[i] = idx(i) + 1;
        }
        const Eigen::VectorXcd wq = -crquad(z4, sing4, prever, beta, qdat_);
        crimage(k) = ((wq(1) - wq(0)) * (wq(3) - wq(2))) / ((wq(2) - wq(1)) * (wq(0) - wq(3)));
    }

    double maxAbs = 0.0;
    for (int k = 0; k < n3; ++k) maxAbs = std::max(maxAbs, std::abs(crimage(k) - crtarget(k)));
    acc_ = maxAbs;
    return acc_;
}

}  // namespace sctoolbox
