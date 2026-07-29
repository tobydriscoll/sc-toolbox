#include "sctoolbox/scangle.hpp"

#include <cmath>
#include <complex>
#include <limits>

namespace sctoolbox {

Eigen::VectorXd scangle(const Eigen::VectorXcd& w) {
    const int n = static_cast<int>(w.size());
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double pi = M_PI;

    Eigen::VectorXd beta = Eigen::VectorXd::Constant(n, nan);

    auto isInf = [](const std::complex<double>& z) {
        return std::isinf(z.real()) || std::isinf(z.imag());
    };

    for (int i = 0; i < n; ++i) {
        const int prev = (i - 1 + n) % n;
        const int next = (i + 1) % n;
        if (isInf(w(i)) || isInf(w(prev)) || isInf(w(next))) continue;

        const std::complex<double> dw = w(i) - w(prev);
        const std::complex<double> dwshift = w(next) - w(i);
        double b = std::arg(dw * std::conj(dwshift)) / pi;
        if (std::abs(b + 1.0) < 1e-12) b = 1.0;
        beta(i) = b;
    }

    return beta;
}

}  // namespace sctoolbox
