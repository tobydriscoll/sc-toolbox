#include "sctoolbox/crgather.hpp"

#include <cmath>
#include <limits>

#include "sctoolbox/moebius3.hpp"

namespace sctoolbox {

namespace {
void crgatherRec(Eigen::VectorXcd& u, std::vector<int>& uquad, int quadnum, const Eigen::VectorXd& cr,
                  const QGraph& Q, Eigen::MatrixXcd& zr) {
    const int n3 = static_cast<int>(cr.size());
    const double r = cr(quadnum);
    const double f1 = std::sqrt(1.0 / (r + 1.0));
    const double f2 = std::sqrt(r / (r + 1.0));
    zr(0, quadnum) = {-f1, -f2};
    zr(1, quadnum) = {-f1, f2};
    zr(2, quadnum) = {f1, f2};
    zr(3, quadnum) = {f1, -f2};

    std::vector<int> nbrs;
    for (int q = 0; q < n3; ++q)
        if (Q.adjacent(quadnum, q) != 0 && std::isnan(zr(0, q).real())) nbrs.push_back(q);

    for (int q : nbrs) {
        crgatherRec(u, uquad, q, cr, Q, zr);

        std::vector<int> i1, i2;
        for (int a = 0; a < 4; ++a)
            for (int b = 0; b < 4; ++b)
                if (Q.qlvert(a, quadnum) == Q.qlvert(b, q)) {
                    i1.push_back(a);
                    i2.push_back(b);
                }

        Eigen::Vector3cd src, dst;
        for (int t = 0; t < 3; ++t) {
            src(t) = zr(i2[t], q);
            dst(t) = zr(i1[t], quadnum);
        }
        const auto mt = moebius3(src, dst);

        for (int p = 0; p < static_cast<int>(uquad.size()); ++p) {
            if (uquad[p] == q) {
                const std::complex<double> up = u(p);
                u(p) = (mt[1] * up + mt[0]) / (mt[3] * up + mt[2]);
                uquad[p] = quadnum;
            }
        }
    }
}
}  // namespace

void crgather(Eigen::VectorXcd& u, std::vector<int>& uquad, int quadnum, const Eigen::VectorXd& cr,
              const QGraph& Q) {
    const int n3 = static_cast<int>(cr.size());
    const std::complex<double> nanc(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN());
    Eigen::MatrixXcd zr = Eigen::MatrixXcd::Constant(4, n3, nanc);
    crgatherRec(u, uquad, quadnum, cr, Q, zr);
}

}  // namespace sctoolbox
