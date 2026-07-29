#include "sctoolbox/hpderiv.hpp"

#include <vector>

namespace sctoolbox {

Eigen::VectorXcd hpderiv(const Eigen::VectorXcd& zp, const Eigen::VectorXcd& z, const Eigen::VectorXd& beta,
                         std::complex<double> c) {
    const int npts = static_cast<int>(zp.size());
    const int n = static_cast<int>(z.size());

    std::vector<std::complex<double>> zf;
    std::vector<double> bf;
    for (int i = 0; i < n; ++i) {
        if (std::isfinite(z(i).real()) && std::isfinite(z(i).imag())) {
            zf.push_back(z(i));
            bf.push_back(beta(i));
        }
    }
    const int m = static_cast<int>(zf.size());

    Eigen::VectorXcd fprime(npts);
    for (int p = 0; p < npts; ++p) {
        std::complex<double> logsum(0.0, 0.0);
        for (int i = 0; i < m; ++i) {
            std::complex<double> term = zp(p) - zf[i];
            logsum += bf[i] * std::log(term);
        }
        fprime(p) = c * std::exp(logsum);
    }
    return fprime;
}

}  // namespace sctoolbox
