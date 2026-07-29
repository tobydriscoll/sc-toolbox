#include "sctoolbox/crfixwc.hpp"

#include "sctoolbox/crembed.hpp"
#include "sctoolbox/crimap0.hpp"
#include "sctoolbox/isinpoly.hpp"
#include "sctoolbox/scqdata.hpp"

namespace sctoolbox {

Eigen::VectorXcd crfixwc(const Eigen::VectorXcd& w, const Eigen::VectorXd& beta, const Eigen::VectorXd& cr,
                         const Eigen::MatrixXcd& aff, const QGraph& Q, std::complex<double> wc) {
    const int n = static_cast<int>(cr.size()) + 3;
    const int n3 = n - 3;

    int quadnum = 0;
    for (int q = 0; q < n3; ++q) {
        Eigen::VectorXcd quadVerts(4);
        for (int i = 0; i < 4; ++i) quadVerts(i) = w(Q.qlvert(i, q));
        Eigen::VectorXcd wcv(1);
        wcv(0) = wc;
        const Eigen::VectorXd r = isinpoly(wcv, quadVerts, 1e-4);
        if (r(0) != 0.0) {
            quadnum = q;
            break;
        }
    }

    const Eigen::VectorXcd z = crembed(cr, Q, quadnum);
    const Eigen::Vector2cd affq(aff(quadnum, 0), aff(quadnum, 1));
    const Eigen::MatrixXd qdat = scqdata(beta, 8);
    Eigen::VectorXcd wcv(1);
    wcv(0) = wc;
    const Eigen::VectorXcd zcv = crimap0(wcv, z, beta, affq, qdat);
    const std::complex<double> zc = zcv(0);

    const std::complex<double> zn = z(n - 1);
    const std::complex<double> ratio = (1.0 - std::conj(zc) * zn) / (zn - zc);
    const std::complex<double> a = ratio / std::abs(ratio);

    Eigen::VectorXcd wcfix(5);
    wcfix(0) = static_cast<double>(quadnum + 1);
    wcfix(1) = 1.0;
    wcfix(2) = zc * a;
    wcfix(3) = std::conj(zc);
    wcfix(4) = a;
    return wcfix;
}

}  // namespace sctoolbox
