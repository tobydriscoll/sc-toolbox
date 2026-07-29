#include "sctoolbox/nefdjac.hpp"

#include <algorithm>
#include <cmath>

namespace sctoolbox {

Eigen::MatrixXd nefdjac(const Fvec& fvec, const Eigen::VectorXd& fc, Eigen::VectorXd xc,
                         const Eigen::VectorXd& sx, const Eigen::VectorXd& details, int& nofun) {
    const int n = static_cast<int>(fc.size());
    const double sqrteta = std::sqrt(details(12));
    Eigen::MatrixXd J(n, n);

    for (int j = 0; j < n; ++j) {
        const double xcj = xc(j);
        const double sign = (xcj > 0.0 ? 1.0 : (xcj < 0.0 ? -1.0 : 0.0)) + (xcj == 0.0 ? 1.0 : 0.0);
        const double stepsizej0 = sqrteta * std::max(std::abs(xcj), 1.0 / sx(j)) * sign;
        const double tempj = xcj;
        xc(j) = xcj + stepsizej0;
        const double stepsizej = xc(j) - tempj;
        const Eigen::VectorXd fj = fvec(xc);
        ++nofun;
        J.col(j) = (fj - fc) / stepsizej;
        xc(j) = tempj;
    }
    return J;
}

}  // namespace sctoolbox
