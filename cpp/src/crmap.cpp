#include "sctoolbox/crmap.hpp"

#include <cmath>
#include <limits>
#include <set>
#include <vector>

#include "sctoolbox/crembed.hpp"
#include "sctoolbox/crmap0.hpp"
#include "sctoolbox/crspread.hpp"

namespace sctoolbox {

Eigen::VectorXcd crmap(const Eigen::VectorXcd& zp, const Eigen::VectorXcd& /*w*/, const Eigen::VectorXd& beta,
                       const Eigen::VectorXd& cr, const Eigen::MatrixXcd& aff, const Eigen::VectorXcd& wcfix,
                       const QGraph& Q, const Eigen::MatrixXd& qdat) {
    const int m = static_cast<int>(zp.size());
    Eigen::VectorXcd wp = Eigen::VectorXcd::Zero(m);

    const int quadnum = static_cast<int>(std::lround(wcfix(0).real())) - 1;
    const std::complex<double> mt1 = wcfix(1), mt2 = wcfix(2), mt3 = wcfix(3), mt4 = wcfix(4);

    Eigen::VectorXcd zl(m);
    for (int k = 0; k < m; ++k) zl(k) = (mt1 * zp(k) + mt2) / (mt3 * zp(k) + mt4);

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

        const Eigen::Vector2cd affq(aff(q, 0), aff(q, 1));
        const Eigen::VectorXcd wpv = crmap0(zlv, z, beta, affq, qdat);

        for (int j = 0; j < mm; ++j) wp(maskIdx[j]) = wpv(j);
    }

    return wp;
}

}  // namespace sctoolbox
