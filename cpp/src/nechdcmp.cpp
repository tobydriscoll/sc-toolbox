#include "sctoolbox/nechdcmp.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace sctoolbox {

void nechdcmp(const Eigen::MatrixXd& H, double maxoffl, Eigen::MatrixXd& L, double& maxadd) {
    const int n = static_cast<int>(H.rows());
    const double eps = std::numeric_limits<double>::epsilon();
    L = Eigen::MatrixXd::Zero(n, n);

    const double minl = std::pow(eps, 0.25) * maxoffl;

    double minl2 = 0.0;
    if (maxoffl == 0.0) {
        maxoffl = std::sqrt(H.diagonal().maxCoeff());
        minl2 = std::sqrt(eps) * maxoffl;
    }

    maxadd = 0.0;

    for (int j = 0; j < n; ++j) {
        if (j == 0) {
            L(j, j) = H(j, j);
        } else {
            L(j, j) = H(j, j) - L.row(j).head(j) * L.row(j).head(j).transpose();
        }
        double minljj = 0.0;
        for (int i = j + 1; i < n; ++i) {
            if (j == 0) {
                L(i, j) = H(j, i);
            } else {
                L(i, j) = H(j, i) - L.row(i).head(j) * L.row(j).head(j).transpose();
            }
            minljj = std::max(std::abs(L(i, j)), minljj);
        }
        minljj = std::max(minljj / maxoffl, minl);
        if (L(j, j) > minljj * minljj) {
            L(j, j) = std::sqrt(L(j, j));
        } else {
            if (minljj < minl2) {
                minljj = minl2;
            }
            maxadd = std::max(maxadd, minljj * minljj - L(j, j));
            L(j, j) = minljj;
        }
        for (int i = j + 1; i < n; ++i) {
            L(i, j) = L(i, j) / L(j, j);
        }
    }
}

}  // namespace sctoolbox
