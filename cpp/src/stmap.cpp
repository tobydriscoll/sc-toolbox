#include "sctoolbox/stmap.hpp"

#include <cmath>
#include <limits>

#include "sctoolbox/stquad.hpp"
#include "sctoolbox/stquadh.hpp"

namespace sctoolbox {

namespace {
bool isInf(std::complex<double> v) { return !std::isfinite(v.real()) || !std::isfinite(v.imag()); }
}  // namespace

Eigen::VectorXcd stmap(const Eigen::VectorXcd& zp, const Eigen::VectorXcd& w, const Eigen::VectorXd& beta,
                       const Eigen::VectorXcd& z, std::complex<double> c, const Eigen::MatrixXd& qdat) {
    const int n = static_cast<int>(w.size());
    const double tol = std::pow(10.0, -static_cast<double>(qdat.rows()));
    const int p = static_cast<int>(zp.size());

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

    int negInfIdx = -1, posInfIdx = -1;
    for (int i = 0; i < n; ++i) {
        if (z(i) == std::complex<double>(-std::numeric_limits<double>::infinity(), 0.0)) negInfIdx = i;
        if (z(i) == std::complex<double>(std::numeric_limits<double>::infinity(), 0.0)) posInfIdx = i;
    }

    std::vector<bool> leftend(p, false), rightend(p, false);
    for (int k = 0; k < p; ++k) {
        if (isInf(zp(k)) && zp(k).real() < 0) {
            leftend[k] = true;
            wp(k) = w(negInfIdx);
        }
        if (isInf(zp(k)) && zp(k).real() > 0) {
            rightend[k] = true;
            wp(k) = w(posInfIdx);
        }
        vertex[k] = vertex[k] || leftend[k] || rightend[k];
    }

    std::vector<bool> atinf(n);
    for (int i = 0; i < n; ++i) atinf[i] = isInf(w(i));

    std::vector<bool> bad(p, false);
    bool anyBad = false;
    for (int k = 0; k < p; ++k) {
        bad[k] = atinf[sing[k] - 1] && !vertex[k];
        if (bad[k]) anyBad = true;
    }

    std::vector<std::complex<double>> mid1(p), mid2(p);
    if (anyBad) {
        Eigen::VectorXcd zf = z;
        for (int i = 0; i < n; ++i)
            if (atinf[i]) zf(i) = std::complex<double>(std::numeric_limits<double>::infinity(), 0.0);

        for (int k = 0; k < p; ++k) {
            if (!bad[k]) continue;
            double best = std::numeric_limits<double>::infinity();
            int bestI = 0;
            for (int i = 0; i < n; ++i) {
                const double d = std::abs(zp(k) - zf(i));
                if (d < best) {
                    best = d;
                    bestI = i;
                }
            }
            sing[k] = bestI + 1;
            mid1[k] = std::complex<double>(z(bestI).real(), 0.5);
            mid2[k] = std::complex<double>(zp(k).real(), 0.5);
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
        const Eigen::VectorXcd I = stquad(zsv, zpv, sngv, z, beta, qdat);
        for (int j = 0; j < m; ++j) {
            const int k = normalIdx[j];
            wp(k) = ws(k) + c * I(j);
        }
    }

    if (!badIdx.empty()) {
        const int m = static_cast<int>(badIdx.size());
        Eigen::VectorXcd zsv(m), mid1v(m), mid2v(m), zpv(m);
        std::vector<int> sngv(m), zeroV(m, 0);
        for (int j = 0; j < m; ++j) {
            const int k = badIdx[j];
            zsv(j) = zs(k);
            mid1v(j) = mid1[k];
            mid2v(j) = mid2[k];
            zpv(j) = zp(k);
            sngv[j] = sing[k];
        }
        const Eigen::VectorXcd I1 = stquad(zsv, mid1v, sngv, z, beta, qdat);
        const Eigen::VectorXcd I2 = stquadh(mid1v, mid2v, zeroV, z, beta, qdat);
        const Eigen::VectorXcd I3 = stquad(zpv, mid2v, zeroV, z, beta, qdat);
        for (int j = 0; j < m; ++j) {
            const int k = badIdx[j];
            wp(k) = ws(k) + c * (I1(j) + I2(j) - I3(j));
        }
    }

    return wp;
}

}  // namespace sctoolbox
