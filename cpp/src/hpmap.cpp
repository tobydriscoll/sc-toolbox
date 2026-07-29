#include "sctoolbox/hpmap.hpp"

#include <cmath>
#include <limits>

#include "sctoolbox/hpquad.hpp"

namespace sctoolbox {

namespace {
bool isInf(std::complex<double> v) { return !std::isfinite(v.real()) || !std::isfinite(v.imag()); }
}  // namespace

Eigen::VectorXcd hpmap(const Eigen::VectorXcd& zp, const Eigen::VectorXcd& w, const Eigen::VectorXd& beta,
                       const Eigen::VectorXcd& z, std::complex<double> c, const Eigen::MatrixXd& qdat) {
    const int n = static_cast<int>(w.size());
    const double tol = std::pow(10.0, -static_cast<double>(qdat.rows()));
    const int p = static_cast<int>(zp.size());

    const Eigen::VectorXcd zFinite = z.head(n - 1);
    const Eigen::VectorXd betaFinite = beta.head(n - 1);

    Eigen::VectorXcd wp = Eigen::VectorXcd::Zero(p);

    std::vector<int> sing(p);  // 1-indexed
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
        sing[k] = bestI + 1;
    }

    std::vector<bool> vertex(p, false);
    for (int k = 0; k < p; ++k) {
        vertex[k] = dist[k] < tol;
        if (vertex[k]) wp(k) = w(sing[k] - 1);
    }
    for (int k = 0; k < p; ++k) {
        if (isInf(zp(k))) {
            wp(k) = w(n - 1);
            vertex[k] = true;
        }
    }

    std::vector<bool> atinf(n);
    for (int i = 0; i < n; ++i) atinf[i] = isInf(w(i));

    std::vector<bool> bad(p, false);
    bool anyBad = false;
    for (int k = 0; k < p; ++k) {
        bad[k] = atinf[sing[k] - 1] && !vertex[k];
        if (bad[k]) anyBad = true;
    }

    std::vector<std::complex<double>> mid(p);
    if (anyBad) {
        for (int k = 0; k < p; ++k) {
            if (!bad[k]) continue;
            const double direcn = (zp(k) - z(sing[k] - 1)).real();
            const double sgn = (direcn > 0.0) ? 1.0 : ((direcn < 0.0) ? -1.0 : 0.0);
            const int adj = static_cast<int>(sgn) + (direcn == 0.0 ? 1 : 0);
            sing[k] += adj;
            mid[k] = (z(sing[k] - 1) + zp(k)) / 2.0;
        }
    }

    Eigen::VectorXcd zs(p), ws(p);
    for (int k = 0; k < p; ++k) {
        zs(k) = z(sing[k] - 1);
        ws(k) = w(sing[k] - 1);
    }

    std::vector<int> normalIdx, badIdx;
    for (int k = 0; k < p; ++k) {
        if (!bad[k] && !vertex[k]) normalIdx.push_back(k);
        if (bad[k]) badIdx.push_back(k);
    }

    if (!normalIdx.empty()) {
        const int m = static_cast<int>(normalIdx.size());
        Eigen::VectorXcd zsv(m), zpv(m);
        std::vector<int> sngv(m);
        for (int j = 0; j < m; ++j) {
            const int k = normalIdx[j];
            zsv(j) = zs(k);
            zpv(j) = zp(k);
            sngv[j] = sing[k];
        }
        const Eigen::VectorXcd I = hpquad(zsv, zpv, sngv, zFinite, betaFinite, qdat);
        for (int j = 0; j < m; ++j) {
            const int k = normalIdx[j];
            wp(k) = ws(k) + c * I(j);
        }
    }

    if (!badIdx.empty()) {
        const int m = static_cast<int>(badIdx.size());
        Eigen::VectorXcd zsv(m), midv(m), zpv(m);
        std::vector<int> sngv(m), zeroV(m, 0);
        for (int j = 0; j < m; ++j) {
            const int k = badIdx[j];
            zsv(j) = zs(k);
            midv(j) = mid[k];
            zpv(j) = zp(k);
            sngv[j] = sing[k];
        }
        const Eigen::VectorXcd I1 = hpquad(zsv, midv, sngv, zFinite, betaFinite, qdat);
        const Eigen::VectorXcd I2 = hpquad(zpv, midv, zeroV, zFinite, betaFinite, qdat);
        for (int j = 0; j < m; ++j) {
            const int k = badIdx[j];
            wp(k) = ws(k) + c * (I1(j) - I2(j));
        }
    }

    return wp;
}

}  // namespace sctoolbox
