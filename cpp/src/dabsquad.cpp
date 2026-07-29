#include "sctoolbox/dabsquad.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace sctoolbox {

Eigen::VectorXd dabsquad(const Eigen::VectorXcd& z1, const Eigen::VectorXcd& z2, const std::vector<int>& sing1in,
                         const Eigen::VectorXcd& z, const Eigen::VectorXd& beta, const Eigen::MatrixXd& qdat) {
    const int nqpts = static_cast<int>(qdat.rows());
    const int n = static_cast<int>(z.size());
    const int m = static_cast<int>(z1.size());

    Eigen::VectorXd argz(n);
    for (int i = 0; i < n; ++i) argz(i) = std::arg(z(i));

    std::vector<int> sing1 = sing1in;
    if (sing1.empty()) sing1.assign(m, 0);

    Eigen::VectorXd I = Eigen::VectorXd::Zero(m);

    for (int k = 0; k < m; ++k) {
        if (z1(k) == z2(k)) continue;

        double arga = std::arg(z1(k));
        double argb = std::arg(z2(k));
        const double ang21 = std::arg(z2(k) / z1(k));
        if ((argb - arga) * ang21 < 0.0) {
            const double s = (ang21 > 0.0) ? 1.0 : (ang21 < 0.0 ? -1.0 : 0.0);
            argb += 2.0 * M_PI * s;
        }

        const std::complex<double> za = z1(k);
        const std::complex<double> zb = z2(k);
        const int sng = sing1[k];  // 1-indexed; 0 = none
        const int sng0 = sng - 1;

        double mindist = std::numeric_limits<double>::infinity();
        for (int i = 0; i < n; ++i) {
            if (i == sng0) continue;
            mindist = std::min(mindist, std::abs(z(i) - za));
        }
        double dist = std::min(1.0, 2.0 * mindist / std::abs(zb - za));
        double argr = arga + dist * (argb - arga);

        const int ind = (sng0 >= 0) ? sng0 : n;

        Eigen::VectorXd nd(nqpts), wt(nqpts);
        for (int q = 0; q < nqpts; ++q) {
            nd(q) = ((argr - arga) * qdat(q, ind) + argr + arga) / 2.0;
            wt(q) = (std::abs(argr - arga) / 2.0) * qdat(q, ind + n + 1);
        }

        Eigen::MatrixXd terms(n, nqpts);
        bool coincident = false;
        for (int i = 0; i < n; ++i) {
            for (int q = 0; q < nqpts; ++q) {
                double theta = std::fmod(nd(q) - argz(i) + 2.0 * M_PI, 2.0 * M_PI);
                if (theta > M_PI) theta = 2.0 * M_PI - theta;
                const double t = 2.0 * std::sin(theta / 2.0);
                terms(i, q) = t;
                if (t == 0.0) coincident = true;
            }
        }

        if (coincident) {
            I(k) = 0.0;
            continue;
        }

        if (sng0 >= 0) {
            for (int q = 0; q < nqpts; ++q) terms(sng0, q) /= std::abs(nd(q) - arga);
            for (int q = 0; q < nqpts; ++q) wt(q) *= std::pow(std::abs(argr - arga) / 2.0, beta(sng0));
        }

        double sum0 = 0.0;
        for (int q = 0; q < nqpts; ++q) {
            double logsum = 0.0;
            for (int i = 0; i < n; ++i) logsum += beta(i) * std::log(terms(i, q));
            sum0 += std::exp(logsum) * wt(q);
        }
        I(k) = sum0;

        while (dist < 1.0) {
            const double argl = argr;
            const std::complex<double> zl = std::exp(std::complex<double>(0.0, argl));
            double mindist2 = std::numeric_limits<double>::infinity();
            for (int i = 0; i < n; ++i) mindist2 = std::min(mindist2, std::abs(z(i) - zl));
            dist = std::min(1.0, 2.0 * mindist2 / std::abs(zl - zb));
            argr = argl + dist * (argb - argl);

            Eigen::VectorXd nd2(nqpts), wt2(nqpts);
            for (int q = 0; q < nqpts; ++q) {
                nd2(q) = ((argr - argl) * qdat(q, n) + argr + argl) / 2.0;
                wt2(q) = (std::abs(argr - argl) / 2.0) * qdat(q, 2 * n + 1);
            }

            double seg = 0.0;
            for (int q = 0; q < nqpts; ++q) {
                double logsum = 0.0;
                for (int i = 0; i < n; ++i) {
                    double theta = std::fmod(nd2(q) - argz(i) + 2.0 * M_PI, 2.0 * M_PI);
                    if (theta > M_PI) theta = 2.0 * M_PI - theta;
                    logsum += beta(i) * std::log(2.0 * std::sin(theta / 2.0));
                }
                seg += std::exp(logsum) * wt2(q);
            }
            I(k) += seg;
        }
    }
    return I;
}

}  // namespace sctoolbox
