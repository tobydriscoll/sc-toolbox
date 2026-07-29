#include "sctoolbox/dmap.hpp"

#include <cmath>
#include <limits>

#include "sctoolbox/dquad.hpp"

namespace sctoolbox {

Eigen::VectorXcd dmap(const Eigen::VectorXcd& zp, const Eigen::VectorXcd& w, const Eigen::VectorXd& beta,
                      const Eigen::VectorXcd& z, std::complex<double> c, const Eigen::MatrixXd& qdat) {
    const int n = static_cast<int>(z.size());
    const double tol = std::pow(10.0, -static_cast<double>(qdat.rows()));
    const int p = static_cast<int>(zp.size());

    Eigen::VectorXcd wp = Eigen::VectorXcd::Zero(p);

    std::vector<int> sing(p);
    std::vector<double> dist(p);
    for (int k = 0; k < p; ++k) {
        double best = std::numeric_limits<double>::infinity();
        int bestI = 0;
        for (int i = 0; i < n; ++i) {
            const double d = std::abs(zp(k) - z(i));
            if (d < best) {
                best = d;
                bestI = i;
            }
        }
        dist[k] = best;
        sing[k] = bestI;  // 0-indexed
    }

    std::vector<bool> vertex(p), bad(p);
    for (int k = 0; k < p; ++k) vertex[k] = dist[k] < tol;
    for (int k = 0; k < p; ++k) {
        if (vertex[k]) wp(k) = w(sing[k]);
    }

    std::vector<bool> atinf(n);
    for (int i = 0; i < n; ++i) atinf[i] = !std::isfinite(w(i).real()) || !std::isfinite(w(i).imag());

    bool anyBad = false;
    for (int k = 0; k < p; ++k) {
        bad[k] = atinf[sing[k]] && !vertex[k];
        if (bad[k]) anyBad = true;
    }

    std::complex<double> wc(0.0, 0.0);
    if (anyBad) {
        const bool wnm1Finite = std::isfinite(w(n - 2).real()) && std::isfinite(w(n - 2).imag());
        if (wnm1Finite) {
            Eigen::VectorXcd z1(1), z2(1);
            z1(0) = z(n - 2);
            z2(0) = 0.0;
            std::vector<int> sng{n - 1};  // 1-indexed position n-1 -> 0-indexed n-2
            const Eigen::VectorXcd I = dquad(z1, z2, sng, z, beta, qdat);
            wc = w(n - 2) + c * I(0);
        } else {
            Eigen::VectorXcd z1(1), z2(1);
            z1(0) = z(n - 1);
            z2(0) = 0.0;
            std::vector<int> sng{n};  // 1-indexed position n -> 0-indexed n-1
            const Eigen::VectorXcd I = dquad(z1, z2, sng, z, beta, qdat);
            wc = w(n - 1) + c * I(0);
        }
    }

    std::vector<int> normalIdx, badIdx;
    for (int k = 0; k < p; ++k) {
        if (!bad[k] && !vertex[k]) normalIdx.push_back(k);
        if (bad[k]) badIdx.push_back(k);
    }

    if (!normalIdx.empty()) {
        const int m = static_cast<int>(normalIdx.size());
        Eigen::VectorXcd zs(m), zpv(m);
        std::vector<int> sng(m);
        for (int j = 0; j < m; ++j) {
            const int k = normalIdx[j];
            zs(j) = z(sing[k]);
            zpv(j) = zp(k);
            sng[j] = sing[k] + 1;  // back to 1-indexed
        }
        const Eigen::VectorXcd I = dquad(zs, zpv, sng, z, beta, qdat);
        for (int j = 0; j < m; ++j) {
            const int k = normalIdx[j];
            wp(k) = w(sing[k]) + c * I(j);
        }
    }

    if (!badIdx.empty()) {
        const int m = static_cast<int>(badIdx.size());
        Eigen::VectorXcd zpv(m), zeroV(m);
        std::vector<int> sng(m, 0);
        for (int j = 0; j < m; ++j) {
            const int k = badIdx[j];
            zpv(j) = zp(k);
            zeroV(j) = 0.0;
        }
        const Eigen::VectorXcd I = dquad(zpv, zeroV, sng, z, beta, qdat);
        for (int j = 0; j < m; ++j) {
            const int k = badIdx[j];
            wp(k) = wc - c * I(j);
        }
    }

    return wp;
}

}  // namespace sctoolbox
