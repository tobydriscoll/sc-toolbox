#include "sctoolbox/moebius3.hpp"

namespace sctoolbox {

std::array<std::complex<double>, 4> moebius3(const Eigen::Vector3cd& z, const Eigen::Vector3cd& w) {
    const std::complex<double> t1 = -(z(1) - z(0)) * (w(2) - w(1));
    const std::complex<double> t2 = -(z(2) - z(1)) * (w(1) - w(0));
    const std::complex<double> A1 = w(2) * z(0) * t2 - w(0) * z(2) * t1;
    const std::complex<double> A2 = w(0) * t1 - w(2) * t2;
    const std::complex<double> A3 = z(0) * t2 - z(2) * t1;
    const std::complex<double> A4 = t1 - t2;
    return {A1, A2, A3, A4};
}

}  // namespace sctoolbox
