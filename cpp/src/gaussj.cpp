#include "sctoolbox/gaussj.hpp"

#include <Eigen/Eigenvalues>
#include <cmath>

namespace sctoolbox {

void gaussj(int n, double alf, double bet, Eigen::VectorXd& z, Eigen::VectorXd& w) {
    const double apb = alf + bet;

    Eigen::VectorXd a(n);
    Eigen::VectorXd b(n > 1 ? n - 1 : 0);

    // a(1)/b(1) use closed forms that avoid a 0/0 when alf=bet=0 (n=1
    // case), unlike the general recurrence used for index >= 2.
    a(0) = (bet - alf) / (apb + 2.0);
    if (n > 1) {
        b(0) = std::sqrt(4.0 * (1.0 + alf) * (1.0 + bet) /
                          ((apb + 3.0) * (apb + 2.0) * (apb + 2.0)));
    }
    for (int k = 2; k <= n; ++k) {
        const double N = static_cast<double>(k);
        a(k - 1) = apb * (bet - alf) / ((apb + 2.0 * N) * (apb + 2.0 * N - 2.0));
    }
    for (int k = 2; k <= n - 1; ++k) {
        const double N = static_cast<double>(k);
        const double apb2N = apb + 2.0 * N;
        b(k - 1) = std::sqrt(4.0 * N * (N + alf) * (N + bet) * (N + apb) /
                              ((apb2N * apb2N - 1.0) * apb2N * apb2N));
    }

    Eigen::VectorXd nodes(n);
    Eigen::VectorXd firstRow(n);  // V(1,:) in MATLAB: first component of each eigenvector
    if (n > 1) {
        Eigen::MatrixXd T = Eigen::MatrixXd::Zero(n, n);
        for (int i = 0; i < n; ++i) T(i, i) = a(i);
        for (int i = 0; i < n - 1; ++i) {
            T(i, i + 1) = b(i);
            T(i + 1, i) = b(i);
        }
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(T);
        nodes = es.eigenvalues();              // ascending, like MATLAB's sort(diag(D))
        firstRow = es.eigenvectors().row(0);
    } else {
        nodes(0) = a(0);
        firstRow(0) = 1.0;
    }

    const double c = std::pow(2.0, apb + 1.0) * std::tgamma(alf + 1.0) *
                      std::tgamma(bet + 1.0) / std::tgamma(apb + 2.0);

    z = nodes;
    w = c * firstRow.array().square();
}

}  // namespace sctoolbox
