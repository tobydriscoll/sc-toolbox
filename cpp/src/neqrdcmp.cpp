#include "sctoolbox/neqrdcmp.hpp"

#include <cmath>

namespace sctoolbox {

void neqrdcmp(Eigen::MatrixXd& M, Eigen::VectorXd& M1, Eigen::VectorXd& M2, int& sing) {
    const int n = static_cast<int>(M.rows());
    M1 = Eigen::VectorXd::Zero(n);
    M2 = Eigen::VectorXd::Zero(n);
    sing = 0;

    for (int k = 0; k < n - 1; ++k) {
        // MATLAB: eta = max(M(k:n,k)) -- max of values, not abs, ported verbatim.
        double eta = M.col(k).segment(k, n - k).maxCoeff();
        if (eta == 0.0) {
            M1(k) = 0.0;
            M2(k) = 0.0;
            sing = 1;
        } else {
            M.col(k).segment(k, n - k) /= eta;
            const double mkk = M(k, k);
            const double sign = (mkk > 0.0 ? 1.0 : (mkk < 0.0 ? -1.0 : 0.0)) + (mkk == 0.0 ? 1.0 : 0.0);
            const double sigma = sign * M.col(k).segment(k, n - k).norm();
            M(k, k) = M(k, k) + sigma;
            M1(k) = sigma * M(k, k);
            M2(k) = -eta * sigma;
            Eigen::VectorXd colk = M.col(k).segment(k, n - k);
            Eigen::RowVectorXd tau = (colk.transpose() * M.block(k, k + 1, n - k, n - k - 1)) / M1(k);
            M.block(k, k + 1, n - k, n - k - 1) -= colk * tau;
        }
    }
    M2(n - 1) = M(n - 1, n - 1);
}

}  // namespace sctoolbox
