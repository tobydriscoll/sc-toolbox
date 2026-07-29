#include "sctoolbox/stderiv.hpp"

#include <cmath>
#include <complex>
#include <vector>

namespace sctoolbox {

Eigen::VectorXcd stderiv(const Eigen::VectorXcd& zp, const Eigen::VectorXcd& z, const Eigen::VectorXd& beta,
                         std::complex<double> c, int j) {
    const int nFull = static_cast<int>(z.size());
    std::vector<int> ends;  // 0-indexed positions of the infinite entries
    for (int i = 0; i < nFull; ++i)
        if (!std::isfinite(z(i).real()) || !std::isfinite(z(i).imag())) ends.push_back(i);

    double theta = beta(ends[1]) - beta(ends[0]);
    if (z(ends[0]).real() < 0.0) theta = -theta;

    int jAdj = j;
    if (j > 0) {
        if (j > ends[0] + 1) jAdj -= 1;
        if (j > ends[1] + 1) jAdj -= 1;
    }

    const int n = nFull - 2;
    Eigen::VectorXcd zred(n);
    Eigen::VectorXd bred(n);
    int idx = 0;
    for (int i = 0; i < nFull; ++i) {
        if (i == ends[0] || i == ends[1]) continue;
        zred(idx) = z(i);
        bred(idx) = beta(i);
        ++idx;
    }

    const int npts = static_cast<int>(zp.size());
    Eigen::VectorXcd fprime(npts);
    const double log2 = 0.69314718055994531;
    const std::complex<double> iUnit(0.0, 1.0);

    for (int p = 0; p < npts; ++p) {
        std::complex<double> sumterm(0.0, 0.0);
        for (int i = 0; i < n; ++i) {
            std::complex<double> term = -M_PI / 2.0 * (zp(p) - zred(i));
            const bool lower = (zred(i).imag() == 0.0);
            if (lower) term = -term;
            const double rt = term.real();
            std::complex<double> tval;
            if (std::abs(rt) > 40.0) {
                const double s = (rt > 0.0) - (rt < 0.0);
                tval = s * (term - iUnit * (M_PI / 2.0)) - log2;
            } else {
                tval = std::log(-iUnit * std::sinh(term));
            }
            if (jAdj > 0 && (i + 1) == jAdj) tval -= std::log(std::abs(zp(p) - zred(i)));
            sumterm += tval * bred(i);
        }
        fprime(p) = c * std::exp(M_PI / 2.0 * theta * zp(p) + sumterm);
    }
    return fprime;
}

}  // namespace sctoolbox
