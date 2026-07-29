#include "sctoolbox/crderiv.hpp"

#include <cmath>
#include <limits>
#include <set>
#include <vector>

#include "sctoolbox/crembed.hpp"
#include "sctoolbox/crspread.hpp"
#include "sctoolbox/dderiv.hpp"

namespace sctoolbox {

Eigen::VectorXcd crderiv(const Eigen::VectorXcd& zp, const Eigen::VectorXd& beta, const Eigen::VectorXd& cr,
                         const Eigen::MatrixXcd& aff, const Eigen::VectorXcd& wcfix, const QGraph& Q) {
    const int m = static_cast<int>(zp.size());
    Eigen::VectorXcd fp = Eigen::VectorXcd::Zero(m);

    const int quadnum = static_cast<int>(std::lround(wcfix(0).real())) - 1;  // wcfix(1) is 1-indexed
    const std::complex<double> mt1 = wcfix(1), mt2 = wcfix(2), mt3 = wcfix(3), mt4 = wcfix(4);

    Eigen::VectorXcd zl(m), d0(m);
    for (int k = 0; k < m; ++k) {
        const std::complex<double> denom = mt3 * zp(k) + mt4;
        zl(k) = (mt1 * zp(k) + mt2) / denom;
        d0(k) = (mt1 * mt4 - mt2 * mt3) / (denom * denom);
    }

    const CrSpreadResult spread = crspread(zl, quadnum, cr, Q);

    std::vector<int> idx(m);
    for (int k = 0; k < m; ++k) {
        double best = std::numeric_limits<double>::infinity();
        int bestQ = 0;
        for (int q = 0; q < static_cast<int>(cr.size()); ++q) {
            const double a = std::abs(spread.ul(q, k));
            if (a < best) {
                best = a;
                bestQ = q;
            }
        }
        idx[k] = bestQ;
    }

    std::set<int> uniqueQ(idx.begin(), idx.end());
    for (int q : uniqueQ) {
        const Eigen::VectorXcd z = crembed(cr, Q, q);

        std::vector<int> maskIdx;
        for (int k = 0; k < m; ++k)
            if (idx[k] == q) maskIdx.push_back(k);
        const int mm = static_cast<int>(maskIdx.size());

        Eigen::VectorXcd zlv(mm);
        for (int j = 0; j < mm; ++j) zlv(j) = spread.ul(q, maskIdx[j]);

        const Eigen::VectorXcd dfp = dderiv(zlv, z, beta);

        for (int j = 0; j < mm; ++j) {
            const int k = maskIdx[j];
            fp(k) = aff(q, 0) * d0(k) * spread.dl(q, k) * dfp(j);
        }
    }

    return fp;
}

}  // namespace sctoolbox
