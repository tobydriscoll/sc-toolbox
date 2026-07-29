#include "sctoolbox/hpinvmap.hpp"

#include <cmath>
#include <limits>

#include "sctoolbox/findz0.hpp"
#include "sctoolbox/hpderiv.hpp"
#include "sctoolbox/hpmap.hpp"
#include "sctoolbox/ode45.hpp"
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
}  // namespace

InvMapResult hpinvmap(const Eigen::VectorXcd& wp, const Eigen::VectorXcd& w, const Eigen::VectorXd& beta,
                       const Eigen::VectorXcd& z, std::complex<double> c, const Eigen::MatrixXd& qdat,
                       const Eigen::VectorXcd& z0In, const std::vector<double>& options) {
    const int n = static_cast<int>(w.size());
    const int lenwpTotal = static_cast<int>(wp.size());
    Eigen::VectorXcd zp = Eigen::VectorXcd::Zero(lenwpTotal);

    const InvOpt opt = scinvopt(options);

    std::vector<bool> done(lenwpTotal, false);
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < lenwpTotal; ++i) {
            if (!done[i] && std::abs(wp(i) - w(j)) < 3.0 * std::numeric_limits<double>::epsilon()) {
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
            const Eigen::VectorXcd wpActive = gather(wp, done);
            auto mapfun = [&](const Eigen::VectorXcd& zq) { return hpmap(zq, w, beta, z, c, qdat); };
            FindZ0Result r = findz0("hp", wpActive, mapfun, w, beta, z, c, qdat);
            z0 = r.z0;
            w0 = r.w0;
        } else if (z0In.size() == 1) {
            Eigen::VectorXcd allz0 = Eigen::VectorXcd::Constant(lenwpTotal, z0In(0));
            Eigen::VectorXcd allw0 = hpmap(allz0, w, beta, z, c, qdat);
            z0 = gather(allz0, done);
            w0 = gather(allw0, done);
        } else {
            Eigen::VectorXcd allw0 = hpmap(z0In, w, beta, z, c, qdat);
            z0 = gather(z0In, done);
            w0 = gather(allw0, done);
        }

        const double odetol = std::max(opt.tol, opt.newton ? 1e-3 : 0.0);
        const Eigen::VectorXcd scale = gather(wp, done) - w0;

        Eigen::VectorXd y0(2 * lenwp);
        y0.head(lenwp) = z0.real();
        y0.tail(lenwp) = z0.imag();

        OdeFun odefun = [&](double, const Eigen::VectorXd& y) {
            Eigen::VectorXd yim = y.tail(lenwp).cwiseMax(0.0);
            const Eigen::VectorXcd zq = y.head(lenwp) + std::complex<double>(0, 1) * yim;
            const Eigen::VectorXcd f = scale.cwiseQuotient(hpderiv(zq, z, beta, c));
            Eigen::VectorXd zdot(2 * lenwp);
            zdot.head(lenwp) = f.real();
            zdot.tail(lenwp) = f.imag();
            return zdot;
        };

        const Eigen::VectorXd y1 = ode45(odefun, 0.0, 1.0, y0, odetol, odetol);
        Eigen::VectorXcd zpActive = y1.head(lenwp) + std::complex<double>(0, 1) * y1.tail(lenwp);
        for (int i = 0; i < lenwp; ++i)
            if (zpActive(i).imag() < 0.0) zpActive(i) = std::complex<double>(zpActive(i).real(), 0.0);
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
            const Eigen::VectorXcd wpActive = gather(wp, done);
            Factive = wpActive - hpmap(znActive, w, beta, z, c, qdat);
            const Eigen::VectorXcd dF = hpderiv(znActive, z, beta, c);
            const Eigen::VectorXcd znNewRaw = znActive + Factive.cwiseQuotient(dF);
            Eigen::VectorXcd znNew(znNewRaw.size());
            for (int i = 0; i < znNewRaw.size(); ++i)
                znNew(i) = std::complex<double>(znNewRaw(i).real(), std::max(0.0, znNewRaw(i).imag()));
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
