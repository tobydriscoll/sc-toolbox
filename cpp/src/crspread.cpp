#include "sctoolbox/crspread.hpp"

#include <cmath>
#include <vector>

#include "sctoolbox/moebius3.hpp"

namespace sctoolbox {

namespace {
Eigen::Vector4cd quadCorners(double r) {
    const double f1 = std::sqrt(1.0 / (r + 1.0));
    const double f2 = std::sqrt(r / (r + 1.0));
    Eigen::Vector4cd c;
    c << std::complex<double>(-f1, -f2), std::complex<double>(-f1, f2), std::complex<double>(f1, f2),
        std::complex<double>(f1, -f2);
    return c;
}
}  // namespace

CrSpreadResult crspread(const Eigen::VectorXcd& uIn, int quadnum, const Eigen::VectorXd& cr, const QGraph& Q) {
    const int n3 = static_cast<int>(cr.size());
    const int m = static_cast<int>(uIn.size());

    Eigen::MatrixXcd ul = Eigen::MatrixXcd::Zero(n3, m);
    Eigen::MatrixXcd dl = Eigen::MatrixXcd::Zero(n3, m);
    ul.row(quadnum) = uIn.transpose();
    dl.row(quadnum) = Eigen::RowVectorXcd::Ones(m);

    Eigen::MatrixXcd zr = Eigen::MatrixXcd::Zero(4, n3);
    zr.col(quadnum) = quadCorners(cr(quadnum));

    std::vector<bool> done(n3, false);
    done[quadnum] = true;
    std::vector<bool> todo(n3, false);
    for (int i = 0; i < n3; ++i) todo[i] = Q.adjacent(i, quadnum) != 0;

    while (true) {
        bool anyNotDone = false;
        for (bool d : done)
            if (!d) anyNotDone = true;
        if (!anyNotDone) break;

        int q = -1;
        for (int i = 0; i < n3; ++i)
            if (todo[i]) {
                q = i;
                break;
            }
        if (q < 0) break;

        zr.col(q) = quadCorners(cr(q));

        int qn = -1;
        for (int i = 0; i < n3; ++i)
            if (done[i] && Q.adjacent(i, q) != 0) {
                qn = i;
                break;
            }

        std::vector<int> i1, i2;
        for (int a = 0; a < 4; ++a)
            for (int b = 0; b < 4; ++b)
                if (Q.qlvert(a, q) == Q.qlvert(b, qn)) {
                    i1.push_back(a);
                    i2.push_back(b);
                }

        Eigen::Vector3cd src, dst;
        for (int t = 0; t < 3; ++t) {
            src(t) = zr(i2[t], qn);
            dst(t) = zr(i1[t], q);
        }
        const auto mt = moebius3(src, dst);

        for (int p = 0; p < m; ++p) {
            const std::complex<double> uqn = ul(qn, p);
            const std::complex<double> den = mt[3] * uqn + mt[2];
            ul(q, p) = (mt[1] * uqn + mt[0]) / den;
            dl(q, p) = dl(qn, p) * (mt[1] * mt[2] - mt[0] * mt[3]) / (den * den);
        }

        done[q] = true;
        todo[q] = false;
        for (int i = 0; i < n3; ++i)
            if (Q.adjacent(i, q) != 0 && !done[i]) todo[i] = true;
    }

    return {ul, dl};
}

}  // namespace sctoolbox
