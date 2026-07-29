#include "sctoolbox/stquadh.hpp"

#include <cmath>
#include <limits>

#include "sctoolbox/stquad.hpp"

namespace sctoolbox {

Eigen::VectorXcd stquadh(const Eigen::VectorXcd& z1, const Eigen::VectorXcd& z2, const std::vector<int>& sing1,
                         const Eigen::VectorXcd& z, const Eigen::VectorXd& beta, const Eigen::MatrixXd& qdat) {
    const int n = static_cast<int>(z.size());
    const int m = static_cast<int>(z1.size());
    const double alf = 0.75;

    Eigen::VectorXcd I = Eigen::VectorXcd::Zero(m);

    for (int k = 0; k < m; ++k) {
        if (z1(k) == z2(k)) continue;

        const std::complex<double> za = z1(k);
        const std::complex<double> zb = z2(k);
        const int sng = sing1.empty() ? 0 : sing1[k];  // 1-indexed; 0 = none
        const int sng0 = sng - 1;

        const double d = zb.real() - za.real();
        const double sgnD = (d > 0.0) ? 1.0 : ((d < 0.0) ? -1.0 : 0.0);

        double Lmin = std::numeric_limits<double>::infinity();
        for (int i = 0; i < n; ++i) {
            const bool zinf = !std::isfinite(z(i).real()) || !std::isfinite(z(i).imag());
            const double dx = (z(i).real() - za.real()) * sgnD;
            const double dy = std::abs(z(i).imag() - za.imag());
            const bool toright = (dx > 0.0) && !zinf;
            bool active = (dx > dy / alf) && toright;
            if (sng0 >= 0 && i == sng0) active = false;

            if (active) {
                const double x = dx;
                const double y = dy;
                const double Lcand = (x - std::sqrt((alf * x) * (alf * x) - (1 - alf * alf) * y * y)) / (1 - alf * alf);
                Lmin = std::min(Lmin, Lcand);
            } else if (toright) {
                Lmin = std::min(Lmin, dy / alf);
            }
        }

        if (Lmin < std::abs(d)) {
            const std::complex<double> zmid = za + Lmin * sgnD;
            Eigen::VectorXcd za1(1), zmid1(1);
            za1(0) = za;
            zmid1(0) = zmid;
            std::vector<int> sngv{sng};
            const Eigen::VectorXcd I1 = stquad(za1, zmid1, sngv, z, beta, qdat);

            Eigen::VectorXcd zmid2(1), zb1(1);
            zmid2(0) = zmid;
            zb1(0) = zb;
            std::vector<int> noSing{0};
            const Eigen::VectorXcd I2 = stquadh(zmid2, zb1, noSing, z, beta, qdat);

            I(k) = I1(0) + I2(0);
        } else {
            Eigen::VectorXcd za1(1), zb1(1);
            za1(0) = za;
            zb1(0) = zb;
            std::vector<int> sngv{sng};
            const Eigen::VectorXcd I1 = stquad(za1, zb1, sngv, z, beta, qdat);
            I(k) = I1(0);
        }
    }
    return I;
}

}  // namespace sctoolbox
