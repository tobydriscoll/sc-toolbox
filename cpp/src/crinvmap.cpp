#include "sctoolbox/crinvmap.hpp"

#include <cmath>

#include "sctoolbox/crembed.hpp"
#include "sctoolbox/crgather.hpp"
#include "sctoolbox/crimap0.hpp"
#include "sctoolbox/isinpoly.hpp"

namespace sctoolbox {

Eigen::VectorXcd crinvmap(const Eigen::VectorXcd& wp, const Eigen::VectorXcd& w, const Eigen::VectorXd& beta,
                          const Eigen::VectorXd& cr, const Eigen::MatrixXcd& aff, const Eigen::VectorXcd& wcfix,
                          const QGraph& Q, const Eigen::MatrixXd& qdat, const std::vector<double>& options) {
    const int n3 = static_cast<int>(cr.size());
    const int lenwp = static_cast<int>(wp.size());
    const double tol = std::pow(10.0, -static_cast<double>(qdat.rows()));

    Eigen::VectorXcd zp = Eigen::VectorXcd::Zero(lenwp);
    std::vector<int> quadnum(lenwp, -1);  // 0-indexed; -1 = unassigned

    for (int q = 0; q < n3; ++q) {
        std::vector<int> idx;
        for (int i = 0; i < lenwp; ++i)
            if (quadnum[i] < 0) idx.push_back(i);
        if (idx.empty()) break;

        Eigen::VectorXcd wpIdx(idx.size());
        for (size_t j = 0; j < idx.size(); ++j) wpIdx(static_cast<int>(j)) = wp(idx[j]);

        Eigen::VectorXcd quadVerts(4);
        for (int i = 0; i < 4; ++i) quadVerts(i) = w(Q.qlvert(i, q));

        const Eigen::VectorXd mask = isinpoly(wpIdx, quadVerts, tol);

        std::vector<int> activeIdx;
        for (int j = 0; j < static_cast<int>(idx.size()); ++j)
            if (mask(j) != 0.0) activeIdx.push_back(idx[j]);

        if (!activeIdx.empty()) {
            const Eigen::VectorXcd z = crembed(cr, Q, q);

            Eigen::VectorXcd wpActive(activeIdx.size());
            for (size_t j = 0; j < activeIdx.size(); ++j) wpActive(static_cast<int>(j)) = wp(activeIdx[j]);

            const Eigen::Vector2cd affq(aff(q, 0), aff(q, 1));
            const Eigen::VectorXcd zpActive = crimap0(wpActive, z, beta, affq, qdat, options);

            for (size_t j = 0; j < activeIdx.size(); ++j) {
                zp(activeIdx[j]) = zpActive(static_cast<int>(j));
                quadnum[activeIdx[j]] = q;
            }
        }

        bool allAssigned = true;
        for (int qn : quadnum)
            if (qn < 0) allAssigned = false;
        if (allAssigned) break;
    }

    const int targetQuad = static_cast<int>(std::lround(wcfix(0).real())) - 1;
    crgather(zp, quadnum, targetQuad, cr, Q);

    const std::complex<double> mt1 = wcfix(1), mt2 = wcfix(2), mt3 = wcfix(3), mt4 = wcfix(4);
    Eigen::VectorXcd result(lenwp);
    for (int i = 0; i < lenwp; ++i) result(i) = (-mt4 * zp(i) + mt2) / (mt3 * zp(i) - mt1);
    return result;
}

}  // namespace sctoolbox
