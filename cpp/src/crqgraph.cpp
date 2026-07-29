#include "sctoolbox/crqgraph.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

#include "sctoolbox/scangle.hpp"

namespace sctoolbox {

QGraph crqgraph(const Eigen::VectorXcd& w, const CrTriangulation& tri) {
    const int N = static_cast<int>(tri.edge.cols());
    const int n3 = (N + 3) / 2 - 3;

    Eigen::MatrixXi qlvert = Eigen::MatrixXi::Zero(4, n3);
    Eigen::MatrixXi qledge = Eigen::MatrixXi::Zero(4, n3);
    Eigen::MatrixXi T = Eigen::MatrixXi::Zero(n3, n3);

    for (int e = 1; e <= n3; ++e) {
        const int t1tri = tri.edgetri(0, e - 1);
        const int t2tri = tri.edgetri(1, e - 1);
        std::array<int, 3> t1, t2;
        for (int i = 0; i < 3; ++i) {
            t1[i] = tri.triedge(i, t1tri - 1);
            t2[i] = tri.triedge(i, t2tri - 1);
        }

        int i1pos = -1, i2pos = -1;
        for (int i = 0; i < 3; ++i) {
            if (t1[i] == e) i1pos = i;
            if (t2[i] == e) i2pos = i;
        }
        std::array<int, 3> t1r, t2r;
        for (int k = 0; k < 3; ++k) t1r[k] = t1[(i1pos + k) % 3];
        for (int k = 0; k < 3; ++k) t2r[k] = t2[(i2pos + k) % 3];

        Eigen::VectorXcd mids1(3), mids2(3);
        for (int k = 0; k < 3; ++k) {
            const std::complex<double> a = w(tri.edge(0, t1r[k] - 1) - 1);
            const std::complex<double> b = w(tri.edge(1, t1r[k] - 1) - 1);
            mids1(k) = (a + b) / 2.0;
            const std::complex<double> c = w(tri.edge(0, t2r[k] - 1) - 1);
            const std::complex<double> d = w(tri.edge(1, t2r[k] - 1) - 1);
            mids2(k) = (c + d) / 2.0;
        }
        if (scangle(mids1).sum() > 0.0) std::swap(t1r[1], t1r[2]);
        if (scangle(mids2).sum() > 0.0) std::swap(t2r[1], t2r[2]);

        qledge(0, e - 1) = t1r[1];
        qledge(1, e - 1) = t1r[2];
        qledge(2, e - 1) = t2r[1];
        qledge(3, e - 1) = t2r[2];

        for (int k = 0; k < 4; ++k) {
            const int val = qledge(k, e - 1);
            if (val <= n3) {
                T(val - 1, e - 1) = 1;
                T(e - 1, val - 1) = 1;
            }
        }

        qlvert(0, e - 1) = tri.edge(0, e - 1);
        qlvert(2, e - 1) = tri.edge(1, e - 1);
        const int q1 = qlvert(0, e - 1), q3 = qlvert(2, e - 1);

        std::vector<int> e1vals, e2vals;
        for (int k = 0; k < 3; ++k) {
            e1vals.push_back(tri.edge(0, t1r[k] - 1));
            e1vals.push_back(tri.edge(1, t1r[k] - 1));
            e2vals.push_back(tri.edge(0, t2r[k] - 1));
            e2vals.push_back(tri.edge(1, t2r[k] - 1));
        }
        std::sort(e1vals.begin(), e1vals.end());
        std::sort(e2vals.begin(), e2vals.end());
        const std::array<int, 3> e1{e1vals[0], e1vals[2], e1vals[4]};
        const std::array<int, 3> e2{e2vals[0], e2vals[2], e2vals[4]};

        for (int v : e1)
            if (v != q1 && v != q3) qlvert(1, e - 1) = v;
        for (int v : e2)
            if (v != q1 && v != q3) qlvert(3, e - 1) = v;

        Eigen::VectorXcd qpts(4);
        for (int k = 0; k < 4; ++k) qpts(k) = w(qlvert(k, e - 1) - 1);
        if (std::abs(scangle(qpts).sum() - 2.0) > 1e-8) {
            const int tmp = qlvert(1, e - 1);
            qlvert(1, e - 1) = qlvert(3, e - 1);
            qlvert(3, e - 1) = tmp;
        }
    }

    QGraph Q;
    Q.edge = tri.edge.array() - 1;
    Q.qlvert = qlvert.array() - 1;
    Q.qledge = qledge.array() - 1;
    Q.adjacent = T;
    return Q;
}

}  // namespace sctoolbox
