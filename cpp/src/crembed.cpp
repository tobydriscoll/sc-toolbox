#include "sctoolbox/crembed.hpp"

#include <cmath>
#include <limits>
#include <vector>

namespace sctoolbox {

Eigen::VectorXcd crembed(const Eigen::VectorXd& cr, const QGraph& Q, int qnum) {
    const int n3 = static_cast<int>(cr.size());
    const int n = n3 + 3;
    const double eps = std::numeric_limits<double>::epsilon();

    Eigen::VectorXcd z = Eigen::VectorXcd::Zero(n);

    const double r = cr(qnum);
    const double f1 = std::sqrt(1.0 / (r + 1.0));
    const double f2 = std::sqrt(r / (r + 1.0));
    const std::complex<double> corner[4] = {
        {-f1, -f2},
        {-f1, f2},
        {f1, f2},
        {f1, -f2},
    };
    const Eigen::Vector4i idx0 = Q.qlvert.col(qnum);
    for (int i = 0; i < 4; ++i) z(idx0(i)) = corner[i];

    std::vector<bool> vtxdone(n, false);
    for (int i = 0; i < 4; ++i) vtxdone[idx0(i)] = true;

    const int numEdges = 2 * n - 3;
    std::vector<bool> edgedone(numEdges, false);
    edgedone[qnum] = true;
    for (int e = n3; e < numEdges; ++e) edgedone[e] = true;  // boundary edges

    std::vector<bool> edgetodo(numEdges, false);
    const Eigen::Vector4i qle0 = Q.qledge.col(qnum);
    for (int i = 0; i < 4; ++i) edgetodo[qle0(i)] = !edgedone[qle0(i)];

    while (true) {
        int e = -1;
        for (int i = 0; i < numEdges; ++i)
            if (edgetodo[i]) {
                e = i;
                break;
            }
        if (e < 0) break;

        Eigen::Vector4i idx = Q.qlvert.col(e);
        if (!vtxdone[idx(1)]) {
            const Eigen::Vector4i tmp = idx;
            idx << tmp(2), tmp(3), tmp(0), tmp(1);
        }

        if (std::abs(z(idx(1)) - z(idx(0))) < 5.0 * eps) {
            z(idx(3)) = z(idx(1));
        } else {
            const std::complex<double> fac = -cr(e) * (z(idx(2)) - z(idx(1))) / (z(idx(1)) - z(idx(0)));
            z(idx(3)) = (z(idx(2)) + fac * z(idx(0))) / (1.0 + fac);
            z(idx(3)) = z(idx(3)) / std::abs(z(idx(3)));
        }

        vtxdone[idx(3)] = true;
        edgedone[e] = true;
        edgetodo[e] = false;
        const Eigen::Vector4i qle = Q.qledge.col(e);
        for (int i = 0; i < 4; ++i) edgetodo[qle(i)] = !edgedone[qle(i)];
    }

    return z;
}

}  // namespace sctoolbox
