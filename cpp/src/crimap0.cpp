#include "sctoolbox/crimap0.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#include "sctoolbox/crmap0.hpp"
#include "sctoolbox/dderiv.hpp"
#include "sctoolbox/ode45.hpp"
#include "sctoolbox/scinvopt.hpp"

namespace sctoolbox {

Eigen::VectorXcd crimap0(const Eigen::VectorXcd& wp, const Eigen::VectorXcd& z, const Eigen::VectorXd& beta,
                         const Eigen::Vector2cd& aff, const Eigen::MatrixXd& qdat,
                         const std::vector<double>& options) {
    const int lenwp = static_cast<int>(wp.size());
    Eigen::VectorXcd zp = Eigen::VectorXcd::Zero(lenwp);

    const InvOpt opt = scinvopt(options);

    std::vector<int> maskIdx;
    for (int i = 0; i < beta.size(); ++i)
        if (std::abs(beta(i)) > std::numeric_limits<double>::epsilon()) maskIdx.push_back(i);
    const int n2 = static_cast<int>(maskIdx.size());
    Eigen::VectorXcd z2(n2);
    Eigen::VectorXd beta2(n2);
    for (int i = 0; i < n2; ++i) {
        z2(i) = z(maskIdx[i]);
        beta2(i) = beta(maskIdx[i]);
    }

    const Eigen::VectorXcd w0 = Eigen::VectorXcd::Constant(lenwp, aff(1));
    std::vector<bool> done(lenwp, false);

    if (opt.ode) {
        const double odetol = std::max(opt.tol, opt.newton ? 1e-2 : 0.0);
        const Eigen::VectorXcd scale = wp - w0;

        const Eigen::VectorXd y0 = Eigen::VectorXd::Zero(2 * lenwp);
        OdeFun odefun = [&](double, const Eigen::VectorXd& y) {
            const Eigen::VectorXcd zq = y.head(lenwp) + std::complex<double>(0, 1) * y.tail(lenwp);
            const Eigen::VectorXcd f = scale.cwiseQuotient(dderiv(zq, z2, beta2, aff(0)));
            Eigen::VectorXd zdot(2 * lenwp);
            zdot.head(lenwp) = f.real();
            zdot.tail(lenwp) = f.imag();
            return zdot;
        };

        const Eigen::VectorXd y1 = ode45(odefun, 0.0, 1.0, y0, odetol, odetol);
        zp = y1.head(lenwp) + std::complex<double>(0, 1) * y1.tail(lenwp);
        for (int i = 0; i < lenwp; ++i)
            if (std::abs(zp(i)) > 1.0) zp(i) /= std::abs(zp(i));
    }

    if (opt.newton) {
        Eigen::VectorXcd zn = opt.ode ? zp : Eigen::VectorXcd::Zero(lenwp);

        int k = 0;
        Eigen::VectorXcd F;
        while (true) {
            int remaining = 0;
            for (bool d : done)
                if (!d) ++remaining;
            if (remaining == 0 || k >= 16) break;  // crimap0.m hardcodes maxiter=16

            std::vector<int> active;
            for (int i = 0; i < lenwp; ++i)
                if (!done[i]) active.push_back(i);
            const int mm = static_cast<int>(active.size());

            Eigen::VectorXcd znActive(mm), wpActive(mm);
            for (int j = 0; j < mm; ++j) {
                znActive(j) = zn(active[j]);
                wpActive(j) = wp(active[j]);
            }

            F = wpActive - crmap0(znActive, z, beta, aff, qdat);

            Eigen::VectorXcd dF(mm);
            for (int j = 0; j < mm; ++j) {
                std::complex<double> logsum(0.0, 0.0);
                for (int t = 0; t < n2; ++t) logsum += beta2(t) * std::log(1.0 - znActive(j) / z2(t));
                dF(j) = aff(0) * std::exp(logsum);
            }

            for (int j = 0; j < mm; ++j) zn(active[j]) = znActive(j) + F(j) / dF(j);
            for (int j = 0; j < mm; ++j)
                if (std::abs(F(j)) < opt.tol) done[active[j]] = true;
            ++k;
        }
        zp = zn;
    }

    return zp;
}

}  // namespace sctoolbox
