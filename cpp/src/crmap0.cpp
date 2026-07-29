#include "sctoolbox/crmap0.hpp"

#include <limits>
#include <vector>

#include "sctoolbox/crquad.hpp"

namespace sctoolbox {

Eigen::VectorXcd crmap0(const Eigen::VectorXcd& zp, const Eigen::VectorXcd& z, const Eigen::VectorXd& beta,
                        const Eigen::Vector2cd& aff, const Eigen::MatrixXd& qdat) {
    const int n = static_cast<int>(z.size());
    const int np = static_cast<int>(zp.size());
    const double eps = std::numeric_limits<double>::epsilon();

    std::vector<int> sing(np, 0);  // 1-indexed; 0 = none
    for (int k = 0; k < np; ++k) {
        for (int i = 0; i < n; ++i) {
            if (std::abs(z(i) - zp(k)) < 10.0 * eps) {
                sing[k] = i + 1;
                break;
            }
        }
    }

    const Eigen::VectorXcd I = crquad(zp, sing, z, beta, qdat);
    Eigen::VectorXcd wp(np);
    for (int k = 0; k < np; ++k) wp(k) = -aff(0) * I(k) + aff(1);
    return wp;
}

}  // namespace sctoolbox
