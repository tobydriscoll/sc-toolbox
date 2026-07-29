#include "sctoolbox/crquad.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>

namespace sctoolbox {

Eigen::VectorXcd crquad(const Eigen::VectorXcd& z1, const std::vector<int>& sing1, const Eigen::VectorXcd& z,
                        const Eigen::VectorXd& beta, const Eigen::MatrixXd& qdat) {
    const int n = static_cast<int>(z.size());
    const int m = static_cast<int>(z1.size());
    const int nqpts = static_cast<int>(qdat.rows());
    const double eps_ = std::numeric_limits<double>::epsilon();

    std::vector<bool> ignore(n);
    for (int i = 0; i < n; ++i) ignore[i] = std::abs(beta(i)) < eps_;

    std::vector<int> keep;
    for (int i = 0; i < n; ++i)
        if (!ignore[i]) keep.push_back(i);
    const int nkeep = static_cast<int>(keep.size());

    Eigen::VectorXcd I = Eigen::VectorXcd::Zero(m);

    for (int k = 0; k < m; ++k) {
        if (std::abs(z1(k)) <= eps_) continue;

        const int sng = sing1.empty() ? 0 : sing1[k];  // MATLAB 1-indexed; 0 = none
        const int sng0 = sng - 1;

        double mindist = std::numeric_limits<double>::infinity();
        for (int i = 0; i < n; ++i) {
            if (ignore[i] || i == sng0) continue;
            mindist = std::min(mindist, std::abs(z(i) - z1(k)));
        }
        const int panels = std::max(1, static_cast<int>(std::ceil(-std::log(mindist) / std::log(2.0))));

        const int qcol0 = (sng0 >= 0) ? sng0 : n;

        int sngK = -1;  // position of sng0 within `keep`, or -1 if not normalized
        if (sng0 >= 0 && !ignore[sng0]) {
            int cnt = 0;
            for (int i = 0; i < sng0; ++i)
                if (!ignore[i]) ++cnt;
            sngK = cnt;
        }

        std::complex<double> za = z1(k);
        std::complex<double> zb(0.0, 0.0);
        int qcol = qcol0;
        double h = 0.0;

        for (int j = 1; j <= panels; ++j) {
            if (j == 1) {
                h = std::pow(2.0, 1 - panels);
            } else {
                h += std::pow(2.0, j - panels - 1);
                za = zb;
                qcol = n;
            }
            zb = z1(k) * (1.0 - h);

            Eigen::VectorXcd nd(nqpts), wt(nqpts);
            for (int q = 0; q < nqpts; ++q) {
                nd(q) = ((zb - za) * qdat(q, qcol) + zb + za) / 2.0;
                wt(q) = (zb - za) / 2.0 * qdat(q, qcol + n + 1);
            }

            bool coincident = false;
            for (int q = 0; q + 1 < nqpts; ++q)
                if (nd(q) == nd(q + 1)) {
                    coincident = true;
                    break;
                }
            for (int q = 0; q < nqpts && !coincident; ++q)
                for (int row = 0; row < nkeep; ++row)
                    if (1.0 - nd(q) / z(keep[row]) == std::complex<double>(0.0, 0.0)) coincident = true;

            if (coincident) {
                I(k) = 0.0;
                continue;
            }

            double wtScale = 1.0;
            if (sngK >= 0 && qcol < n) wtScale = std::pow(std::abs(zb - za) / 2.0, beta(keep[sngK]));

            std::complex<double> contribution(0.0, 0.0);
            for (int q = 0; q < nqpts; ++q) {
                std::complex<double> logsum(0.0, 0.0);
                for (int row = 0; row < nkeep; ++row) {
                    std::complex<double> term = 1.0 - nd(q) / z(keep[row]);
                    if (sngK >= 0 && qcol < n && row == sngK) term /= std::abs(term);
                    logsum += beta(keep[row]) * std::log(term);
                }
                contribution += std::exp(logsum) * wt(q);
            }
            I(k) += contribution * wtScale;
        }
    }
    return I;
}

}  // namespace sctoolbox
