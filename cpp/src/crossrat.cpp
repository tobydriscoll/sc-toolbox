#include "sctoolbox/crossrat.hpp"

namespace sctoolbox {

Eigen::VectorXcd crossrat(const Eigen::VectorXcd& w, const QGraph& Q) {
    const int n3 = static_cast<int>(Q.qlvert.cols());
    Eigen::VectorXcd cr(n3);
    for (int k = 0; k < n3; ++k) {
        const std::complex<double> w1 = w(Q.qlvert(0, k));
        const std::complex<double> w2 = w(Q.qlvert(1, k));
        const std::complex<double> w3 = w(Q.qlvert(2, k));
        const std::complex<double> w4 = w(Q.qlvert(3, k));
        cr(k) = ((w2 - w1) * (w4 - w3)) / ((w3 - w2) * (w1 - w4));
    }
    return cr;
}

}  // namespace sctoolbox
