#include "sctoolbox/rinvmap.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#include "sctoolbox/findz0.hpp"
#include "sctoolbox/ode45.hpp"
#include "sctoolbox/r2strip.hpp"
#include "sctoolbox/rcorners.hpp"
#include "sctoolbox/rderiv.hpp"
#include "sctoolbox/rmap.hpp"
#include "sctoolbox/scinvopt.hpp"

namespace sctoolbox {

namespace {
Eigen::VectorXcd gather(const Eigen::VectorXcd& v, const std::vector<bool>& mask) {
    std::vector<std::complex<double>> out;
    for (int i = 0; i < v.size(); ++i)
        if (!mask[i]) out.push_back(v(i));
    return Eigen::Map<const Eigen::VectorXcd>(out.data(), static_cast<int>(out.size()));
}

void scatter(Eigen::VectorXcd& v, const std::vector<bool>& mask, const Eigen::VectorXcd& vals) {
    int k = 0;
    for (int i = 0; i < v.size(); ++i)
        if (!mask[i]) v(i) = vals(k++);
}

struct Rect {
    double xmin, xmax, ymin, ymax;
};

Eigen::VectorXcd rectproject(const Eigen::VectorXcd& zp, const Rect& rect) {
    Eigen::VectorXcd out(zp.size());
    for (int i = 0; i < zp.size(); ++i) {
        const double re = std::max(std::min(zp(i).real(), rect.xmax), rect.xmin);
        const double im = std::max(std::min(zp(i).imag(), rect.ymax), rect.ymin);
        out(i) = std::complex<double>(re, im);
    }
    return out;
}
}  // namespace

InvMapResult rinvmap(const Eigen::VectorXcd& wpIn, const Eigen::VectorXcd& wIn, const Eigen::VectorXd& betaIn,
                      const Eigen::VectorXcd& zIn, std::complex<double> c, double L, const Eigen::MatrixXd& qdat,
                      const Eigen::VectorXcd& z0In, const std::vector<double>& options) {
    const int n = static_cast<int>(wIn.size());
    const RCorners rc = rcorners(wIn, betaIn, zIn);
    const Eigen::VectorXcd& w = rc.w;
    const Eigen::VectorXd& beta = rc.beta;
    const Eigen::VectorXcd& z = rc.z;

    Rect rect;
    rect.xmin = std::numeric_limits<double>::infinity();
    rect.xmax = -std::numeric_limits<double>::infinity();
    rect.ymin = std::numeric_limits<double>::infinity();
    rect.ymax = -std::numeric_limits<double>::infinity();
    for (int idx : rc.corners) {
        rect.xmin = std::min(rect.xmin, z(idx).real());
        rect.xmax = std::max(rect.xmax, z(idx).real());
        rect.ymin = std::min(rect.ymin, z(idx).imag());
        rect.ymax = std::max(rect.ymax, z(idx).imag());
    }

    Eigen::VectorXcd zs = r2strip(z, z, L).yp;
    for (int i = 0; i < n; ++i) zs(i) = std::complex<double>(zs(i).real(), std::round(zs(i).imag()));

    const int lenwpTotal = static_cast<int>(wpIn.size());
    Eigen::VectorXcd zp = Eigen::VectorXcd::Zero(lenwpTotal);

    const InvOpt opt = scinvopt(options);

    std::vector<bool> done(lenwpTotal, false);
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < lenwpTotal; ++i) {
            if (!done[i] && std::abs(wpIn(i) - w(j)) < 3.0 * std::numeric_limits<double>::epsilon()) {
                zp(i) = z(j);
                done[i] = true;
            }
        }
    }
    int lenwp = 0;
    for (bool d : done) if (!d) ++lenwp;
    if (lenwp == 0) return InvMapResult{zp, {}};

    Eigen::VectorXcd zn;

    if (opt.ode) {
        Eigen::VectorXcd z0, w0;
        if (z0In.size() == 0) {
            const Eigen::VectorXcd wpActive = gather(wpIn, done);
            auto mapfun = [&](const Eigen::VectorXcd& zq) { return rmap(zq, w, beta, z, c, L, qdat); };
            FindZ0Result r = findz0("r", wpActive, mapfun, w, beta, z, c, qdat);
            z0 = r.z0;
            w0 = r.w0;
        } else if (z0In.size() == 1) {
            Eigen::VectorXcd allz0 = Eigen::VectorXcd::Constant(lenwpTotal, z0In(0));
            Eigen::VectorXcd allw0 = rmap(allz0, w, beta, z, c, L, qdat);
            z0 = gather(allz0, done);
            w0 = gather(allw0, done);
        } else {
            Eigen::VectorXcd allw0 = rmap(z0In, w, beta, z, c, L, qdat);
            z0 = gather(z0In, done);
            w0 = gather(allw0, done);
        }

        const double odetol = std::max(opt.tol, opt.newton ? 1e-3 : 0.0);
        const Eigen::VectorXcd scale = gather(wpIn, done) - w0;

        Eigen::VectorXd y0(2 * lenwp);
        y0.head(lenwp) = z0.real();
        y0.tail(lenwp) = z0.imag();

        OdeFun odefun = [&](double, const Eigen::VectorXd& y) {
            const Eigen::VectorXcd zq = y.head(lenwp) + std::complex<double>(0, 1) * y.tail(lenwp);
            const Eigen::VectorXcd f = scale.cwiseQuotient(rderiv(zq, z, beta, c, L, zs));
            Eigen::VectorXd zdot(2 * lenwp);
            zdot.head(lenwp) = f.real();
            zdot.tail(lenwp) = f.imag();
            return zdot;
        };

        const Eigen::VectorXd y1 = ode45(odefun, 0.0, 1.0, y0, odetol, odetol);
        Eigen::VectorXcd zpActive = y1.head(lenwp) + std::complex<double>(0, 1) * y1.tail(lenwp);
        zpActive = rectproject(zpActive, rect);
        scatter(zp, done, zpActive);
    }

    if (opt.newton) {
        zn = zp;

        int k = 0;
        Eigen::VectorXcd Factive;
        while (true) {
            int remaining = 0;
            for (bool d : done) if (!d) ++remaining;
            if (remaining == 0 || k >= opt.maxiter) break;

            const Eigen::VectorXcd znActive = gather(zn, done);
            const Eigen::VectorXcd wpActive = gather(wpIn, done);
            Factive = wpActive - rmap(znActive, w, beta, z, c, L, qdat);
            const Eigen::VectorXcd dF = rderiv(znActive, z, beta, c, L, zs);
            Eigen::VectorXcd znNew = znActive + Factive.cwiseQuotient(dF);
            znNew = rectproject(znNew, rect);
            scatter(zn, done, znNew);

            int idx = 0;
            for (int i = 0; i < lenwpTotal; ++i) {
                if (!done[i]) {
                    if (std::abs(Factive(idx)) < opt.tol) done[i] = true;
                    ++idx;
                }
            }
            ++k;
        }
        zp = zn;
    }

    std::vector<int> flag;
    for (int i = 0; i < lenwpTotal; ++i)
        if (!done[i]) flag.push_back(i);
    return InvMapResult{zp, flag};
}

}  // namespace sctoolbox
