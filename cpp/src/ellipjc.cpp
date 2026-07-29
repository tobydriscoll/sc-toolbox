#include "sctoolbox/ellipjc.hpp"

#include <cmath>

#include "sctoolbox/ellipkkp.hpp"

namespace sctoolbox {

EllipJC ellipjc(const Eigen::VectorXcd& uIn, double L, bool mIsM) {
    const int n = static_cast<int>(uIn.size());
    Eigen::VectorXcd u = uIn;
    std::vector<bool> high(n, false);
    double m;

    if (!mIsM) {
        const EllipKKp kk = ellipkkp(L);
        for (int i = 0; i < n; ++i) {
            if (u(i).imag() > kk.Kp / 2.0) {
                high[i] = true;
                u(i) = std::complex<double>(0.0, kk.Kp) - u(i);
            }
        }
        m = std::exp(-2.0 * M_PI * L);
    } else {
        m = L;
    }

    EllipJC r;
    r.sn.resize(n);
    r.cn.resize(n);
    r.dn.resize(n);

    const double eps = std::numeric_limits<double>::epsilon();
    if (m < 4.0 * eps) {
        for (int i = 0; i < n; ++i) {
            const std::complex<double> sinu = std::sin(u(i));
            const std::complex<double> cosu = std::cos(u(i));
            r.sn(i) = sinu + m / 4.0 * (sinu * cosu - u(i)) * cosu;
            r.cn(i) = cosu + m / 4.0 * (-sinu * cosu + u(i)) * sinu;
            r.dn(i) = 1.0 + m / 4.0 * (cosu * cosu - sinu * sinu - 1.0);
        }
    } else {
        double kappa;
        if (m > 1e-3) {
            kappa = (1.0 - std::sqrt(1.0 - m)) / (1.0 + std::sqrt(1.0 - m));
        } else {
            // polyval([132,42,14,5,2,1,0], m/4)
            const double x = m / 4.0;
            kappa = (((((132.0 * x + 42.0) * x + 14.0) * x + 5.0) * x + 2.0) * x + 1.0) * x + 0.0;
        }
        const double mu = kappa * kappa;
        Eigen::VectorXcd v(n);
        for (int i = 0; i < n; ++i) v(i) = u(i) / (1.0 + kappa);

        const EllipJC r1 = ellipjc(v, mu, true);
        for (int i = 0; i < n; ++i) {
            const std::complex<double> denom = 1.0 + kappa * r1.sn(i) * r1.sn(i);
            r.sn(i) = (1.0 + kappa) * r1.sn(i) / denom;
            r.cn(i) = r1.cn(i) * r1.dn(i) / denom;
            r.dn(i) = (1.0 - kappa * r1.sn(i) * r1.sn(i)) / denom;
        }
    }

    const std::complex<double> I(0.0, 1.0);
    for (int i = 0; i < n; ++i) {
        if (!high[i]) continue;
        const std::complex<double> snh = r.sn(i), cnh = r.cn(i), dnh = r.dn(i);
        const double sqm = std::sqrt(m);
        r.sn(i) = -1.0 / (sqm * snh);
        r.cn(i) = I * dnh / (sqm * snh);
        r.dn(i) = I * cnh / snh;
    }

    return r;
}

}  // namespace sctoolbox
