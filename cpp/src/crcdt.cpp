#include "sctoolbox/crcdt.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace sctoolbox {

CrTriangulation crcdt(const Eigen::VectorXcd& w, CrTriangulation tri) {
    Eigen::MatrixXi& edge = tri.edge;
    Eigen::MatrixXi& triedge = tri.triedge;
    Eigen::MatrixXi& edgetri = tri.edgetri;
    const int numedge = static_cast<int>(edge.cols());

    std::vector<bool> interior(numedge);
    for (int i = 0; i < numedge; ++i) interior[i] = (edgetri(1, i) != 0);

    std::vector<bool> done(numedge, true);
    for (int i = 0; i < numedge; ++i)
        if (interior[i]) done[i] = false;

    int quadvtx[5];  // 1-indexed slots 1..4 used

    while (true) {
        int e = -1;
        for (int i = 0; i < numedge; ++i)
            if (!done[i]) {
                e = i + 1;
                break;
            }
        if (e < 0) break;

        const int t1tri = edgetri(0, e - 1);
        const int t2tri = edgetri(1, e - 1);
        int t1[3], t2[3];
        for (int i = 0; i < 3; ++i) {
            t1[i] = triedge(i, t1tri - 1);
            t2[i] = triedge(i, t2tri - 1);
        }

        quadvtx[1] = edge(0, e - 1);
        quadvtx[3] = edge(1, e - 1);

        std::vector<int> e1vals, e2vals;
        for (int i = 0; i < 3; ++i) {
            e1vals.push_back(edge(0, t1[i] - 1));
            e1vals.push_back(edge(1, t1[i] - 1));
            e2vals.push_back(edge(0, t2[i] - 1));
            e2vals.push_back(edge(1, t2[i] - 1));
        }
        std::sort(e1vals.begin(), e1vals.end());
        std::sort(e2vals.begin(), e2vals.end());
        const std::vector<int> e1{e1vals[0], e1vals[2], e1vals[4]};
        const std::vector<int> e2{e2vals[0], e2vals[2], e2vals[4]};

        for (int v : e1)
            if (v != quadvtx[1] && v != quadvtx[3]) quadvtx[2] = v;
        for (int v : e2)
            if (v != quadvtx[1] && v != quadvtx[3]) quadvtx[4] = v;

        const std::complex<double> edgew = w(quadvtx[3] - 1) - w(quadvtx[1] - 1);
        double alpha[4];
        alpha[0] = std::arg(edgew / (w(quadvtx[2] - 1) - w(quadvtx[1] - 1)));
        alpha[1] = std::arg((w(quadvtx[4] - 1) - w(quadvtx[1] - 1)) / edgew);
        alpha[2] = std::arg(-edgew / (w(quadvtx[4] - 1) - w(quadvtx[3] - 1)));
        alpha[3] = std::arg(-(w(quadvtx[2] - 1) - w(quadvtx[3] - 1)) / edgew);
        double sumAlpha = 0.0;
        for (int i = 0; i < 4; ++i) sumAlpha += std::abs(alpha[i]);

        if (sumAlpha < M_PI) {
            edge(0, e - 1) = quadvtx[2];
            edge(1, e - 1) = quadvtx[4];

            int i1pos = -1, i2pos = -1;
            for (int i = 0; i < 3; ++i) {
                if (t1[i] == e) i1pos = i;
                if (t2[i] == e) i2pos = i;
            }
            int i1 = (i1pos + 1) % 3;
            int i2 = (i2pos + 1) % 3;

            bool shareEndpoint = false;
            for (int p = 0; p < 2; ++p)
                for (int q = 0; q < 2; ++q)
                    if (edge(p, t1[i1] - 1) == edge(q, t2[i2] - 1)) shareEndpoint = true;
            if (shareEndpoint) i2 = (i2pos + 2) % 3;

            triedge(i1, t1tri - 1) = t2[i2];
            const int edgeT2I2 = t2[i2];
            for (int r = 0; r < 2; ++r)
                if (edgetri(r, edgeT2I2 - 1) == t2tri) edgetri(r, edgeT2I2 - 1) = t1tri;

            triedge(i2, t2tri - 1) = t1[i1];
            const int edgeT1I1 = t1[i1];
            for (int r = 0; r < 2; ++r)
                if (edgetri(r, edgeT1I1 - 1) == t1tri) edgetri(r, edgeT1I1 - 1) = t2tri;

            for (int i = 0; i < 3; ++i) {
                const int ee = t1[i];
                done[ee - 1] = !interior[ee - 1];
            }
            for (int i = 0; i < 3; ++i) {
                const int ee = t2[i];
                done[ee - 1] = !interior[ee - 1];
            }
        }
        done[e - 1] = true;
    }

    return tri;
}

}  // namespace sctoolbox
