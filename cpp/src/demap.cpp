#include "sctoolbox/demap.hpp"

#include <cmath>
#include <limits>

#include "sctoolbox/dequad.hpp"

namespace sctoolbox {

Eigen::VectorXcd demap(const Eigen::VectorXcd& zp, const Eigen::VectorXcd& w, const Eigen::VectorXd& beta,
                       const Eigen::VectorXcd& z, std::complex<double> c, const Eigen::MatrixXd& qdat) {
    const int n = static_cast<int>(w.size());
    const double tol = std::pow(10.0, -static_cast<double>(qdat.rows()));
    const int p = static_cast<int>(zp.size());

    Eigen::VectorXcd wp = Eigen::VectorXcd::Zero(p);

    std::vector<int> sing(p);  // 0-indexed
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
        sing[k] = bestI;
    }

    std::vector<bool> vertex(p, false);
    for (int k = 0; k < p; ++k) {
        vertex[k] = dist[k] < tol;
        if (vertex[k]) wp(k) = w(sing[k]);
    }
    for (int k = 0; k < p; ++k) {
        if (std::abs(zp(k)) < tol) {
            wp(k) = std::complex<double>(std::numeric_limits<double>::infinity(), 0.0);
            vertex[k] = true;
        }
    }

    for (int k = 0; k < p; ++k) {
        if (!vertex[k]) wp(k) = w(sing[k]);
    }

    for (int k = 0; k < p; ++k) {
        if (vertex[k]) continue;
        const double abszpk = std::abs(zp(k));
        std::complex<double> zold = z(sing[k]);
        int sng = sing[k] + 1;  // 1-indexed, only valid for the first sub-segment
        bool done = false;
        while (!done) {
            const double d = std::min(1.0, 2.0 * abszpk / std::abs(zp(k) - zold));
            const std::complex<double> znew = zold + d * (zp(k) - zold);

            Eigen::VectorXcd z1(1), z2(1);
            z1(0) = zold;
            z2(0) = znew;
            std::vector<int> sngv{sng};
            const Eigen::VectorXcd I = dequad(z1, z2, sngv, z, beta, qdat);
            wp(k) += c * I(0);

            done = !(d < 1.0);
            zold = znew;
            sng = 0;
        }
    }

    return wp;
}

}  // namespace sctoolbox
