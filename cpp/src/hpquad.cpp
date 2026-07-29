#include "sctoolbox/hpquad.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>

namespace sctoolbox {

namespace {
// Evaluates exp(sum_i beta(i)*log(terms(i,q))) . wt(q) (dot product over
// quadrature points q), where terms(i,q) = nd(q) - z(i).
std::complex<double> integrateBlock(const Eigen::VectorXcd& nd, const Eigen::VectorXcd& wt,
                                     const Eigen::VectorXcd& z, const Eigen::VectorXd& beta, int sng0) {
    const int n = static_cast<int>(z.size());
    const int nqpts = static_cast<int>(nd.size());

    std::complex<double> total(0.0, 0.0);
    for (int q = 0; q < nqpts; ++q) {
        std::complex<double> logsum(0.0, 0.0);
        for (int i = 0; i < n; ++i) {
            std::complex<double> term = nd(q) - z(i);
            if (i == sng0) term /= std::abs(term);
            logsum += beta(i) * std::log(term);
        }
        total += std::exp(logsum) * wt(q);
    }
    return total;
}
}  // namespace

Eigen::VectorXcd hpquad(const Eigen::VectorXcd& z1, const Eigen::VectorXcd& z2, const std::vector<int>& sing1,
                        const Eigen::VectorXcd& z, const Eigen::VectorXd& beta, const Eigen::MatrixXd& qdat) {
    const int n = static_cast<int>(z.size());
    const int m = static_cast<int>(z1.size());

    Eigen::VectorXcd I = Eigen::VectorXcd::Zero(m);

    for (int k = 0; k < m; ++k) {
        if (z1(k) == z2(k)) continue;
        const std::complex<double> za = z1(k);
        const std::complex<double> zb = z2(k);
        const int sng = sing1.empty() ? 0 : sing1[k];  // MATLAB 1-indexed; 0 = none
        const int sng0 = sng - 1;                      // 0-indexed singularity row, or -1

        double mindist = std::numeric_limits<double>::infinity();
        for (int i = 0; i < n; ++i) {
            if (i == sng0) continue;
            mindist = std::min(mindist, std::abs(z(i) - za));
        }
        double dist = std::min(1.0, 2.0 * mindist / std::abs(zb - za));
        std::complex<double> zr = za + dist * (zb - za);

        const int ind = (sng0 >= 0) ? sng0 : n;  // 0-indexed qdat column
        const int nqpts = static_cast<int>(qdat.rows());

        Eigen::VectorXcd nd(nqpts);
        Eigen::VectorXcd wt(nqpts);
        for (int q = 0; q < nqpts; ++q) {
            nd(q) = ((zr - za) * qdat(q, ind) + zr + za) / 2.0;
            wt(q) = (zr - za) / 2.0 * qdat(q, ind + n + 1);
        }

        bool coincident = false;
        for (int q = 0; q < nqpts && !coincident; ++q)
            for (int i = 0; i < n; ++i)
                if (nd(q) - z(i) == std::complex<double>(0.0, 0.0)) coincident = true;

        if (coincident) {
            I(k) = 0.0;
            continue;
        }

        if (sng0 >= 0) {
            for (int q = 0; q < nqpts; ++q) wt(q) *= std::pow(std::abs(zr - za) / 2.0, beta(sng0));
        }
        I(k) = integrateBlock(nd, wt, z, beta, sng0);

        while (dist < 1.0) {
            const std::complex<double> zl = zr;
            double mindist2 = std::numeric_limits<double>::infinity();
            for (int i = 0; i < n; ++i) mindist2 = std::min(mindist2, std::abs(z(i) - zl));
            dist = std::min(1.0, 2.0 * mindist2 / std::abs(zl - zb));
            zr = zl + dist * (zb - zl);

            Eigen::VectorXcd nd2(nqpts), wt2(nqpts);
            for (int q = 0; q < nqpts; ++q) {
                nd2(q) = ((zr - zl) * qdat(q, n) + zr + zl) / 2.0;
                wt2(q) = (zr - zl) / 2.0 * qdat(q, 2 * n + 1);
            }
            I(k) += integrateBlock(nd2, wt2, z, beta, -1);
        }
    }
    return I;
}

}  // namespace sctoolbox
