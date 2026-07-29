#include "sctoolbox/dderiv.hpp"

namespace sctoolbox {

Eigen::VectorXcd dderiv(const Eigen::VectorXcd& zp, const Eigen::VectorXcd& z, const Eigen::VectorXd& beta,
                        std::complex<double> c) {
    const int npts = static_cast<int>(zp.size());
    const int n = static_cast<int>(z.size());
    Eigen::VectorXcd fprime(npts);

    for (int p = 0; p < npts; ++p) {
        std::complex<double> logsum(0.0, 0.0);
        for (int i = 0; i < n; ++i) {
            std::complex<double> term = 1.0 - zp(p) / z(i);
            logsum += beta(i) * std::log(term);
        }
        fprime(p) = c * std::exp(logsum);
    }
    return fprime;
}

Eigen::VectorXcd dderiv(const Eigen::VectorXcd& zp, const Eigen::VectorXcd& z, const Eigen::VectorXd& beta,
                        std::complex<double> c, Eigen::VectorXcd& d2f) {
    const int npts = static_cast<int>(zp.size());
    const int n = static_cast<int>(z.size());
    Eigen::VectorXcd fprime(npts);
    d2f = Eigen::VectorXcd::Zero(npts);

    for (int p = 0; p < npts; ++p) {
        std::complex<double> logsum(0.0, 0.0);
        for (int i = 0; i < n; ++i) {
            std::complex<double> term = 1.0 - zp(p) / z(i);
            logsum += beta(i) * std::log(term);
            d2f(p) -= (beta(i) / z(i)) / term;
        }
        fprime(p) = c * std::exp(logsum);
    }
    d2f = d2f.cwiseProduct(fprime);
    return fprime;
}

}  // namespace sctoolbox
